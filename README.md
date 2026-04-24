# NextPass Ephemeris

Compute satellite/target pass trajectories from pre-computed az/el ephemeris files. Determines rise/set times, checks antenna rate limits, proposes feasible tracking windows with slew checking, and generates antenna control commands.

## Files

- `NextPass.py` — Core module. Parses ephemeris files, finds passes, checks rate limits, generates commands.
- `plot_pass.py` — Plotting script. Visualizes pass trajectory, rates, violations, and trackable windows.
- `test_nextpass.py` — Unit tests covering pass detection, rate/keyhole violations, and the 0/360 azimuth wrap.

## Ephemeris Format

Input is a text file with pre-computed topocentric az/el positions in sexagesimal (DD:MM:SS) format. The target name is read from the `NAME =` header line.

```
FORMAT = EPHEMERIS
VELDEF = VRAD-TOP
COORDMODE = azel
HEAD = date utc az el
NAME = MyTarget
2026-03-20 17:50:00.000   15:00:00.0   -20:00:00.0
2026-03-20 17:51:00.000   17:00:00.0   -16:00:00.0
2026-03-20 17:52:00.000   19:00:00.0   -12:00:00.0
...
```

The module also supports RA/DEC ephemeris files in J2000 coordinates. Set `COORDMODE = J2000` (or `radec`) and provide RA (HH:MM:SS) and DEC (DD:MM:SS):

```
FORMAT = EPHEMERIS
VELDEF = VRAD-TOP
COORDMODE = J2000
HEAD = date utc ra dec dra ddec
NAME = LRO
2026-Feb-27 23:00 07:50:26.72 +24:17:12.7 -1404.6334 -1545.7100
2026-Feb-27 23:05 07:50:18.38 +24:15:36.9 -1563.5493 -736.2550
...
```

RA/DEC coordinates are automatically converted to AzEl using the observer location (defined at the top of `NextPass.py`). Both date formats (`2026-03-20` and `2026-Feb-27`) are supported.

The ephemeris can span beyond the rise and set times (include negative elevations). The module will find the actual horizon crossings and rate-allowable trajectories.

## Dependencies

`NextPass.py` uses `astropy` for the observer location and for RA/DEC to AzEl coordinate conversion. The plotting script additionally needs `matplotlib` and `numpy`.

```
pip install astropy matplotlib numpy
```

## Command Line Usage

The target name is taken from the ephemeris `NAME =` header; it is not a CLI argument.

```bash
python NextPass.py ephemeris.txt

# Custom minimum elevation and search range
python NextPass.py ephemeris.txt --min-el 5.0 --range 24.0

# Specify start time (otherwise uses the ephemeris start time)
python NextPass.py ephemeris.txt --start-time "2026-03-20 17:55:00"
python NextPass.py ephemeris.txt --start-time now
```

Example output:

```
Now: 2026-03-19 20:34:49.198247+00:00
Time until  (hh:mm:ss): 21:22:10.801607
Time/Target of aquisition 2026-03-20 17:57:00+00:00 32.0 8.0
Pass at midpoint 2026-03-20 18:03:00+00:00 140.0 52.0
Time/Target at set  2026-03-20 18:09:00+00:00 244.0 8.0
Pass duration (s) 720.0
Az wrap: CWwrap
Rate limits OK

Commands:
Slew(Location("AzEl", 32.0, 8.0))
WaitFor('17:57:00')
Track("MyTarget", None, 720)
```

## Plotting

```bash
python plot_pass.py ephemeris.txt
```

Produces a 5-panel figure with elevation vs time, azimuth vs time, az/el polar plot, elevation rate vs time, and azimuth rate vs time. Rate limit violations are shaded in red across the time-series panels. The plotting script also prints the same pass summary and commands as the CLI above.

Additional options:

```bash
python plot_pass.py ephemeris.txt --min-el 5.0 --range 24.0 \
                                  --start-time now \
                                  --max-az-rate 1.0 --max-el-rate 0.5
```

## Using from Another Script

### Basic Usage

```python
from datetime import datetime, timezone
from NextPass import GetNextPass

now = datetime.now(timezone.utc)
rtn = GetNextPass('ephemeris.txt', minEl=5.0, now=now)

print(rtn.start_time(), rtn.start_az(), rtn.start_el())
print(rtn.end_time(), rtn.end_az(), rtn.end_el())
print(rtn.pass_duration())
print(rtn.which_wrap())
```

`GetNextPass` signature: `GetNextPass(ephem_file, minEl, now=None, trange_hrs=12.0, max_az_rate=MAX_AZ_RATE, max_el_rate=MAX_EL_RATE)`. The `now` argument accepts a `datetime`, the string `'now'`, the string `'ephem'` (use the ephemeris start time — also the default when `now=None`), or a datetime string like `'2026-03-20 17:55:00'`.

### Getting Commands as Strings for use in AstrID

`generate_commands()` returns a list of command strings. `print_commands()` prints them and returns the list.

```python
rtn = GetNextPass('ephemeris.txt', 5.0, now)

# Get without printing
commands = rtn.generate_commands()

# Print and get
commands = rtn.print_commands()
```

Each trackable window produces three commands which can be executed by AstrID to track the allowable trajectory of the given source:

```python
['Slew(Location("AzEl", 32.0, 8.0))',
 "WaitFor('17:57:00')",
 'Track("MyTarget", None, 720)']
```

The target name used in the `Track(...)` command is taken from the ephemeris `NAME =` header.

### Executing Commands Against Your Control System

For more granularity, an AstrID script can iterate the trackable windows directly rather than parsing strings:

```python
from NextPass import GetNextPass

rtn = GetNextPass('ephemeris.txt', 5.0, now)

for ws, we in rtn.trackable_windows:
    az, el = rtn._orbit.SatPositionAt(ws)
    dur = (we - ws).total_seconds()

    loc = Location("AzEl", az, el)
    Slew(loc)
    WaitFor(ws.strftime('%H:%M:%S'))
    Track(loc, None, dur)
```

### Checking for Rate Limit Violations

Violations are accounted for in the trajectory suggestions but you can query them directly as below.

```python
rtn = GetNextPass('ephemeris.txt', 5.0, now)

if rtn.has_violations():
    print(rtn.rate_warning)

    for v in rtn.violations:
        print(f"Violation: {v['start']} - {v['end']} ({v['reason']})")

    for ws, we in rtn.trackable_windows:
        dur = (we - ws).total_seconds()
        print(f"Safe window: {ws} - {we} ({dur:.0f}s)")
```

Each violation dict includes `start`, `end`, `peak_az_rate`, `peak_el_rate`, `peak_elev`, and a human-readable `reason` string.

### Accessing the Sampled Trajectory

The trajectory is pre-computed and attached to the result for plotting or further analysis:

```python
rtn = GetNextPass('ephemeris.txt', 5.0, now)

traj = rtn.trajectory
# traj['secs']     - seconds from rise
# traj['az']       - azimuth in degrees (unwrapped across 0/360; may go below 0 or above 360)
# traj['el']       - elevation in degrees
# traj['az_rate']  - azimuth rate in deg/sec
# traj['el_rate']  - elevation rate in deg/sec
# traj['times']    - datetime objects
```

The azimuth array is unwrapped so it stays continuous across the 0/360 boundary — convenient for plotting and differencing. Take `az % 360` if you need wrapped values.

## Rate Limit Checking

The module checks antenna rate limits (default: 0.3 deg/s elevation, 0.6 deg/s azimuth) plus an elevation keyhole (default 80°, above which the telescope cannot observe). A configurable safety buffer (`RATE_LIMIT_BUFFER = 0.9` by default) is applied to the rate limits to leave headroom for acceleration and braking.

Three scenarios are handled:

1. **No violations** — full pass is one trackable window.
2. **Violations with feasible windows** — the pass is split into trackable segments, with slew feasibility checked between them. If the antenna can't slew to a window's start position in time, the start is delayed to the earliest reachable point, or the window is dropped entirely.
3. **Entire pass exceeds limits** — no commands are generated.

Custom rate limits can be passed to `GetNextPass`:

```python
rtn = GetNextPass('ephemeris.txt', 5.0, now,
                  max_az_rate=1.0, max_el_rate=0.5)
```

### Azimuth 0/360 Wrap Handling

Rate and slew calculations use a wrap-aware angular difference so that a motion from e.g. 359° to 1° is treated as a 2° change rather than a 358° change. This prevents phantom rate violations at the 0/360 seam and allows the slew-feasibility check to take the short way around.

## Testing

Unit tests live in `test_nextpass.py`. Each test synthesizes a small in-memory ephemeris file engineered to trigger (or not trigger) one specific code path, so the suite is self-contained and does not need external data.

### Running the tests

From the directory containing `NextPass.py` and `test_nextpass.py`:

```bash
# Run all tests, verbose
python -m unittest test_nextpass.py -v

# Run all tests, quiet
python -m unittest test_nextpass.py

# Run one test class
python -m unittest test_nextpass.TestAzimuthWrap -v

# Run a single test method
python -m unittest test_nextpass.TestAzRateLimit.test_fast_azimuth_motion_flagged -v

# Run the file directly (uses the __main__ block)
python test_nextpass.py
```

If you add more test files later, `python -m unittest discover` will find and run them all.

### What the tests cover

| Test class | Scenario |
|---|---|
| `TestAzimuthWrap` | A slow pass that crosses the 0/360 boundary produces no phantom rate violation and a single trackable window. Also covers the `_az_diff` helper directly. |
| `TestAzRateLimit` | Azimuth motion above the effective limit (0.6 × buffer) is flagged with `az_rate` in the reason. |
| `TestElRateLimit` | Elevation motion above the effective limit (0.3 × buffer) is flagged with `el_rate` in the reason. |
| `TestElKeyhole` | Elevation above 80° is flagged as a keyhole violation even when rates are small. |
| `TestElMinHorizon` | Rise and set times are trimmed to the `minEl` crossings, not the ephemeris extremes. |
| `TestNormalPass` | A gentle pass well within all limits yields exactly one trackable window spanning rise to set, with no violations. |
| `TestMultipleTrackableWindows` | A mid-pass violation splits the pass into two or more non-overlapping trackable windows, and `generate_commands()` emits three commands per window. |

## Telescope Location

The observer location is set at the top of `NextPass.py` using astropy's site registry:

```python
from astropy.coordinates import EarthLocation
site = EarthLocation.of_site("Green Bank Telescope")
```

To use a different site, replace that line with another registered site name, or construct an `EarthLocation` directly from lat/lon/elevation:

```python
site = EarthLocation(lat=38.050417 * u.deg,
                     lon=-81.104889 * u.deg,
                     height=517.0 * u.m)
```

This is used for RA/DEC to AzEl conversion when processing J2000 ephemeris files. It does not affect AzEl ephemeris files since those already contain pre-computed topocentric positions.