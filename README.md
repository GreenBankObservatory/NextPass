# NextPass Ephemeris

Compute satellite/target pass trajectories from pre-computed az/el ephemeris files. Determines rise/set times, checks antenna rate limits, proposes feasible tracking windows with slew checking, and generates antenna control commands.

## Files

- `NextPass_ephem.py` — Core module. Parses ephemeris files, finds passes, checks rate limits, generates commands.
- `plot_pass.py` — Plotting script. Visualizes pass trajectory, rates, violations, and trackable windows.

## Ephemeris Format

Input is a text file with pre-computed topocentric az/el positions in sexagesimal (DD:MM:SS) format:

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

The module also supports RA/DEC ephemeris files in J2000 coordinates. Set `COORDMODE = J2000` and provide RA (HH:MM:SS) and DEC (DD:MM:SS):
 
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
 
RA/DEC coordinates are automatically converted to AzEl using the observer location (defined at the top of `NextPass_ephem.py`). Both date formats (`2026-03-20` and `2026-Feb-27`) are supported.

The ephemeris can span beyond the rise and set times (include negative elevations). The module will find the actual horizon crossings and rate allowable trajectories.

## Dependencies

`NextPass_ephem.py` requires `astropy` for RA/DEC to AzEl coordinate conversion (only loaded when processing J2000 ephemeris files). The plotting script additionally needs `matplotlib` and `numpy`.

```
pip install matplotlib numpy
```

## Command Line Usage

```bash
python NextPass_ephem.py ephemeris.txt TargetName

# Custom minimum elevation and uptime range
python NextPass_ephem.py ephemeris.txt TargetName --min-el 5.0 --range 24.0
```

Example output:

```
Now: 2026-03-19 20:34:49.198247+00:00
Time until rise: 21:22:10.801607
Time/Target of rise 2026-03-20 17:57:00+00:00 32.0 8.0
Pass at midpoint 2026-03-20 18:03:00+00:00 140.0 52.0
Time/Target at set  2026-03-20 18:09:00+00:00 244.0 8.0
Pass duration 720.0
Az wrap: CWwrap

Commands:
Slew(Location("AzEl", 32.0, 8.0))
WaitFor('17:57:00')
Track(Location("AzEl", 32.0, 8.0), None, 301)
Slew(Location("AzEl", 181.7, 45.3))
WaitFor('18:04:16')
Track(Location("AzEl", 181.7, 45.3), None, 284)
```

## Plotting

```bash
python plot_pass.py ephemeris.txt TargetName
```

Produces a 5-panel figure (`pass_plot.png`) with elevation vs time, azimuth vs time, az/el polar plot, elevation rate vs time, and azimuth rate vs time. Rate limit violations are shaded in red across all panels. This also prints all output from CL usage above.

## Using from Another Script

### Basic Usage

```python
from datetime import datetime, timezone
from NextPass_ephem import GetNextPass

now = datetime.now(timezone.utc)
rtn = GetNextPass('ephemeris.txt', 'TargetName', 5.0, now, 24.0)

print(rtn.start_time(), rtn.start_az(), rtn.start_el())
print(rtn.end_time(), rtn.end_az(), rtn.end_el())
print(rtn.pass_duration())
print(rtn.which_wrap())
```

### Getting Commands as Strings for use in AstrID

`generate_commands()` returns a list of command strings. `print_commands()` prints them and returns the list.

```python
rtn = GetNextPass('ephemeris.txt', 'TargetName', 5.0, now, 24.0)

# Get without printing
commands = rtn.generate_commands()

# Print and get
commands = rtn.print_commands()
```

Each trackable window produces three commands which can be executed by AstrID to track allowable trajectory of given source:

```python
['Slew(Location("AzEl", 32.0, 8.0))',
 "WaitFor('17:57:00')",
 'Track(Location("AzEl", 32.0, 8.0), None, 301)']
```

### Executing Commands Against Your Control System

For more granularity, AstrID script could iterate the trackable windows directly rather than parsing strings:

```python
from NextPass_ephem import GetNextPass

rtn = GetNextPass('ephemeris.txt', 'TargetName', 5.0, now, 24.0)

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
rtn = GetNextPass('ephemeris.txt', 'TargetName', 5.0, now, 24.0)

if rtn.has_violations():
    print(rtn.rate_warning)

    for v in rtn.violations:
        print(f"Violation: {v['start']} - {v['end']} ({v['reason']})")

    for ws, we in rtn.trackable_windows:
        dur = (we - ws).total_seconds()
        print(f"Safe window: {ws} - {we} ({dur:.0f}s)")
```

### Accessing the Sampled Trajectory

The trajectory is pre-computed and attached to the result for plotting or further analysis:

```python
rtn = GetNextPass('ephemeris.txt', 'TargetName', 5.0, now, 24.0)

traj = rtn.trajectory
# traj['secs']     - seconds from rise
# traj['az']       - azimuth in degrees
# traj['el']       - elevation in degrees
# traj['az_rate']  - azimuth rate in deg/sec
# traj['el_rate']  - elevation rate in deg/sec
# traj['times']    - datetime objects
```

## Rate Limit Checking

The module checks antenna rate limits (default: 0.3 deg/s elevation, 0.6 deg/s azimuth) and handles three scenarios:

1. **No violations** — full pass is one trackable window.
2. **Violations with feasible windows** — the pass is split into trackable segments, with slew feasibility checked between them. If the antenna can't slew to a window's start position in time, the start is delayed to the earliest reachable point, or the window is dropped entirely.
3. **Entire pass exceeds limits** — no commands are generated.

Custom rate limits can be passed to `GetNextPass`:

```python
rtn = GetNextPass('ephemeris.txt', 'TargetName', 5.0, now, 24.0,
                  max_az_rate=1.0, max_el_rate=0.5)
```

## Telescope Location
 
The observer location is defined as module-level constants at the top of `NextPass_ephem.py`:
 
```python
OBSERVER_LAT = 38.050417
OBSERVER_LON = -81.104889
OBSERVER_ELEV = 517.0
```
 
This is used for RA/DEC to AzEl conversion when processing J2000 ephemeris files. It does not affect AzEl ephemeris files since those already contain pre-computed topocentric positions.