# -*- coding: utf-8 -*-
from __future__ import print_function
from datetime import datetime, timedelta, tzinfo
from astropy.coordinates import EarthLocation
import warnings
from astropy.utils.iers import conf as iers_conf
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time
import astropy.units as u


# UTC tzinfo singleton (Python 2.7 has no datetime.timezone)
class _UTC(tzinfo):
    ZERO = timedelta(0)

    def utcoffset(self, dt):
        return self.ZERO

    def tzname(self, dt):
        return "UTC"

    def dst(self, dt):
        return self.ZERO

UTC = _UTC()


# Default antenna rate limits (deg/sec)
MAX_EL_RATE = 0.3
MAX_AZ_RATE = 0.6
RATE_LIMIT_BUFFER = 0.9

# Observer location
site = EarthLocation.of_site("Green Bank Telescope")

def _az_diff(a1, a0):
    """Shortest signed difference a1 - a0 in degrees, wrapped to (-180, 180]."""
    return (a1 - a0 + 180.0) % 360.0 - 180.0

def _unwrap_az(az_list):
    """Make an azimuth sequence continuous by removing 360-degree jumps."""
    if not az_list:
        return az_list
    out = [az_list[0]]
    for a in az_list[1:]:
        out.append(out[-1] + _az_diff(a, out[-1] % 360.0))
    return out


class NextPass:
    def __init__(self, rise_t, rise_az, rise_el, set_t, set_az, set_el):
        self.rise_t = rise_t
        self.rise_az = rise_az
        self.rise_el = rise_el
        self.set_t = set_t
        self.set_az = set_az
        self.set_el = set_el
        self.violations = []
        self.trackable_windows = []
        self.rate_warning = None

    def start_time(self):
        return self.rise_t
    def start_time_datetime(self):
        return self.rise_t

    def start_az(self):
        return float(self.rise_az)
    def start_el(self):
        return float(self.rise_el)

    def end_time(self):
        return self.set_t
    def end_time_datetime(self):
        return self.set_t

    def end_az(self):
        return float(self.set_az)
    def end_el(self):
        return float(self.set_el)

    def pass_duration(self, now=None):
        """Returns the pass duration in seconds"""
        delta = self.end_time_datetime() - self.start_time_datetime()
        return delta.total_seconds()

    def time_until_rise(self, when=None):
        if when is None:
            when = datetime.now(UTC)
        return self.rise_t - when

    def setmidpt(self, t, azz, ell):
        self.midp_t = t
        self.midp_az = azz
        self.midp_el = ell

    def midpoint_time(self):
        return self.midp_t
    def midpoint_az(self):
        return self.midp_az
    def midpoint_el(self):
        return self.midp_el

    def has_violations(self):
        return len(self.violations) > 0

    def generate_commands(self):
        """Generate Slew/WaitFor/Track command sequence for trackable windows.

        Returns a list of command strings.
        """
        commands = []
        name = getattr(self._orbit, 'name', None) or 'target'
        for ws, we in self.trackable_windows:
            az, el = self._orbit.SatPositionAt(ws)
            dur = (we - ws).total_seconds()
            commands.append('Slew(Location("AzEl", {0:.1f}, {1:.1f}))'.format(az, el))
            commands.append("WaitFor('{0}')".format(ws.strftime('%H:%M:%S')))
            commands.append(
                'Track("{0}", None, {1:.0f})'.format(name, dur)
            )
        return commands

    def print_commands(self):
        """Print and return the command sequence."""
        commands = self.generate_commands()
        for cmd in commands:
            print(cmd)
        return commands

    def quadrant(self, x):
        x = (x + 360.0) % 360.0
        if x >= 0 and x < 90:
            return 0
        elif x >= 90 and x < 180:
            return 1
        elif x >= 180 and x < 270:
            return 2
        else:
            return 3

    def which_wrap(self):
        pts = [self.rise_az, self.midp_az, self.set_az]
        is_NW = False
        is_NE = False
        is_SW = False
        is_SE = False
        for pt in pts:
            q = self.quadrant(pt)
            if q == 0:
                is_NE = True
            elif q == 1:
                is_SE = True
            elif q == 2:
                is_SW = True
            else:
                is_NW = True
        all_N = not is_SW and not is_SE
        all_S = not is_NW and not is_NE
        if all_N or all_S:
            return "Auto"
        if is_SW:
            return "CWwrap"
        if is_SE:
            return "CCWwrap"
        return "Auto"


def sexagesimal_to_decimal(s):
    """Convert DD:MM:SS.s string to decimal degrees."""
    s = s.strip()
    sign = -1 if s.startswith('-') else 1
    s = s.lstrip('+-')
    parts = s.split(':')
    d = float(parts[0])
    m = float(parts[1]) if len(parts) > 1 else 0.0
    sec = float(parts[2]) if len(parts) > 2 else 0.0
    return sign * (d + m / 60.0 + sec / 3600.0)


def hms_to_degrees(s):
    """Convert HH:MM:SS.s string to degrees (RA: hours -> degrees)."""
    return sexagesimal_to_decimal(s) * 15.0


def _parse_datetime(date_str, time_str):
    """Try multiple datetime formats."""
    combined = date_str + ' ' + time_str
    for fmt in ('%Y-%m-%d %H:%M:%S.%f', '%Y-%m-%d %H:%M:%S',
                '%Y-%b-%d %H:%M:%S.%f', '%Y-%b-%d %H:%M:%S',
                '%Y-%b-%d %H:%M'):
        try:
            return datetime.strptime(combined, fmt)
        except ValueError:
            continue
    return None


def _radec_to_azel(times_utc, ra_degs, dec_degs, lat, lon, elev_m):
    """Convert lists of RA/DEC (degrees, J2000) to AzEl for an observer.

    Returns lists of (az_deg, el_deg).
    """
    warnings.filterwarnings('ignore', module='astropy')
    iers_conf.auto_download = False

    observer = EarthLocation(lat=lat * u.deg, lon=lon * u.deg, height=elev_m * u.m)

    az_list = []
    el_list = []
    for dt, ra, dec in zip(times_utc, ra_degs, dec_degs):
        t = Time(dt, scale='utc')
        coord = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame='icrs')
        altaz = coord.transform_to(AltAz(obstime=t, location=observer))
        az_list.append(altaz.az.deg)
        el_list.append(altaz.alt.deg)

    return az_list, el_list

def _gappt_to_azel(times_utc, ra_degs, dec_degs, lat, lon, elev_m):
    """Convert lists of RA/DEC (degrees, GAPPT/apparent) to AzEl for an observer.

    Returns lists of (az_deg, el_deg).
    """
    from astropy.coordinates import GCRS
    warnings.filterwarnings('ignore', module='astropy')
    iers_conf.auto_download = False

    observer = EarthLocation(lat=lat * u.deg, lon=lon * u.deg, height=elev_m * u.m)

    az_list = []
    el_list = []
    for dt, ra, dec in zip(times_utc, ra_degs, dec_degs):
        t = Time(dt, scale='utc')
        coord = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame=GCRS(obstime=t))
        altaz = coord.transform_to(AltAz(obstime=t, location=observer))
        az_list.append(altaz.az.deg)
        el_list.append(altaz.alt.deg)

    return az_list, el_list

def parse_ephemeris(filename, observer_lat=site.lat.deg, observer_lon=site.lon.deg,
                    observer_elev=site.height.value):
    """Parse an ephemeris text file (azel or J2000 RA/DEC).

    Reads COORDMODE from the header to determine format.
    For J2000, converts RA/DEC to AzEl using the observer location.

    Returns a tuple of (points, ephem_name) where points is a list of
    (datetime_utc, az_deg, el_deg) tuples and ephem_name is the NAME
    from the header (or None if not present).
    """
    coordmode = 'azel'
    ephem_name = None
    header_lines = []
    data_lines = []

    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if '=' in line:
                key, val = line.split('=', 1)
                key = key.strip().upper()
                val_raw = val.strip()
                if key == 'COORDMODE':
                    coordmode = val_raw.lower()
                elif key == 'NAME':
                    ephem_name = val_raw
                header_lines.append(line)
                continue
            data_lines.append(line)

    points = []

    if coordmode == 'azel':
        for line in data_lines:
            parts = line.split()
            if len(parts) < 4:
                continue
            dt = _parse_datetime(parts[0], parts[1])
            if dt is None:
                continue
            dt = dt.replace(tzinfo=UTC)
            az = sexagesimal_to_decimal(parts[2])
            el = sexagesimal_to_decimal(parts[3])
            points.append((dt, az, el))

    elif coordmode in ('j2000', 'radec', 'gappt'):
        times_utc = []
        ra_degs = []
        dec_degs = []
        for line in data_lines:
            parts = line.split()
            if len(parts) < 4:
                continue
            dt = _parse_datetime(parts[0], parts[1])
            if dt is None:
                continue
            dt = dt.replace(tzinfo=UTC)
            ra = hms_to_degrees(parts[2])
            dec = sexagesimal_to_decimal(parts[3])
            times_utc.append(dt)
            ra_degs.append(ra)
            dec_degs.append(dec)

        if times_utc:
            if coordmode == 'gappt':
                az_list, el_list = _gappt_to_azel(
                    times_utc, ra_degs, dec_degs,
                    observer_lat, observer_lon, observer_elev
                )
            else:
                az_list, el_list = _radec_to_azel(
                    times_utc, ra_degs, dec_degs,
                    observer_lat, observer_lon, observer_elev
                )
            for dt, az, el in zip(times_utc, az_list, el_list):
                points.append((dt, az, el))
    else:
        raise Exception("Unknown COORDMODE: {0} (expected azel, j2000, radec, or gappt)".format(coordmode))

    return points, ephem_name


def interp_at(points, t):
    """Linearly interpolate az/el at time t from the ephemeris points."""
    if t <= points[0][0]:
        return points[0][1], points[0][2]
    if t >= points[-1][0]:
        return points[-1][1], points[-1][2]
    for i in range(len(points) - 1):
        t0, az0, el0 = points[i]
        t1, az1, el1 = points[i + 1]
        if t0 <= t <= t1:
            span = (t1 - t0).total_seconds()
            if span == 0:
                return az0, el0
            frac = (t - t0).total_seconds() / span
            # Wrap-aware az interpolation: take the short way around the circle
            az_delta = _az_diff(az1, az0)
            az = (az0 + frac * az_delta) % 360.0
            el = el0 + frac * (el1 - el0)
            return az, el
    return points[-1][1], points[-1][2]


class OrbitPasses:
    """Load a pre-computed az/el ephemeris and provide position lookups."""

    def __init__(self, ephem_file, name=None):
        points, ephem_name = parse_ephemeris(ephem_file)
        # Prefer name from ephem file header; fall back to caller-supplied name
        self.name = ephem_name or name
        self.points = points
        if not self.points:
            raise Exception("No data points found in {0}".format(ephem_file))

    def SatPositionAt(self, t):
        """Return (az_deg, el_deg) at datetime t via interpolation."""
        return interp_at(self.points, t)

    def VelocityAt(self, t, dt_sec=1.0):
        """Return (az_rate, el_rate) in deg/sec at datetime t via central difference."""
        delta = timedelta(seconds=dt_sec)
        az0, el0 = interp_at(self.points, t - delta)
        az1, el1 = interp_at(self.points, t + delta)
        az_rate = _az_diff(az1, az0) / (2.0 * dt_sec)    # <-- changed
        el_rate = (el1 - el0) / (2.0 * dt_sec)
        return az_rate, el_rate

    def sample_trajectory(self, start_t, end_t, n_points=500):
        """Sample positions and velocities across a time window.

        Returns a dict with lists: 'times', 'secs', 'az', 'el', 'az_rate', 'el_rate'
        where 'secs' is seconds from start_t.
        """
        duration = (end_t - start_t).total_seconds()
        n_points = min(max(int(duration), 100), n_points)

        times = []
        secs = []
        azs = []
        els = []
        az_rates = []
        el_rates = []
        for i in range(n_points + 1):
            t = start_t + timedelta(seconds=i * duration / n_points)
            az, el = self.SatPositionAt(t)
            az_rate, el_rate = self.VelocityAt(t)
            times.append(t)
            secs.append((t - start_t).total_seconds())
            azs.append(az)
            els.append(el)
            az_rates.append(az_rate)
            el_rates.append(el_rate)

        azs = _unwrap_az(azs)

        return {
            'times': times,
            'secs': secs,
            'az': azs,
            'el': els,
            'az_rate': az_rates,
            'el_rate': el_rates,
        }

    def check_rate_limits(self, start_t, end_t, max_az_rate=MAX_AZ_RATE,
                          max_el_rate=MAX_EL_RATE, max_el_keyhole=80.0,
                          sample_sec=1.0, min_window_sec=10.0):
        """Check if antenna rate limits or keyhole are exceeded during a time window.

        A violation occurs when the az or el rate exceeds its limit, or when
        the elevation exceeds the keyhole (telescope cannot observe there).

        Returns:
            violations: list of dicts with 'start', 'end', 'peak_az_rate',
                        'peak_el_rate', 'peak_elev', 'reason' for each violation window.
            trackable_windows: list of (start, end) tuples where rates and elevation
                               are within limits -- the safe windows for tracking.
        """
        # add buffer in case antena motion isnt full captured in the motion, ex. braking or accelerating
        max_az_rate = max_az_rate * RATE_LIMIT_BUFFER
        max_el_rate = max_el_rate * RATE_LIMIT_BUFFER

        duration = (end_t - start_t).total_seconds()
        n_samples = max(int(duration / sample_sec), 2)
        step = duration / n_samples

        # Sample rates and elevation across the pass
        exceeded = []
        for i in range(n_samples + 1):
            t = start_t + timedelta(seconds=i * step)
            az, el = self.SatPositionAt(t)
            az_rate, el_rate = self.VelocityAt(t)
            over = (abs(az_rate) > max_az_rate
                    or abs(el_rate) > max_el_rate
                    or el > max_el_keyhole)
            exceeded.append((t, over, az_rate, el_rate, el))

        # Group contiguous violations into raw windows
        raw_violations = []
        in_violation = False
        v_start = None
        peak_az = 0.0
        peak_el = 0.0
        peak_elev = 0.0

        for t, over, az_r, el_r, el in exceeded:
            if over and not in_violation:
                in_violation = True
                v_start = t
                peak_az = abs(az_r)
                peak_el = abs(el_r)
                peak_elev = el
            elif over and in_violation:
                peak_az = max(peak_az, abs(az_r))
                peak_el = max(peak_el, abs(el_r))
                peak_elev = max(peak_elev, el)
            elif not over and in_violation:
                in_violation = False
                raw_violations.append({
                    'start': v_start, 'end': t,
                    'peak_az_rate': peak_az, 'peak_el_rate': peak_el,
                    'peak_elev': peak_elev,
                })

        if in_violation:
            raw_violations.append({
                'start': v_start, 'end': exceeded[-1][0],
                'peak_az_rate': peak_az, 'peak_el_rate': peak_el,
                'peak_elev': peak_elev,
            })

        # Merge violations separated by gaps shorter than min_window_sec
        if len(raw_violations) <= 1:
            merged = raw_violations
        else:
            merged = [raw_violations[0].copy()]
            for v in raw_violations[1:]:
                gap = (v['start'] - merged[-1]['end']).total_seconds()
                if gap < min_window_sec:
                    merged[-1]['end'] = v['end']
                    merged[-1]['peak_az_rate'] = max(merged[-1]['peak_az_rate'],
                                                     v['peak_az_rate'])
                    merged[-1]['peak_el_rate'] = max(merged[-1]['peak_el_rate'],
                                                     v['peak_el_rate'])
                    merged[-1]['peak_elev'] = max(merged[-1]['peak_elev'],
                                                  v['peak_elev'])
                else:
                    merged.append(v.copy())

        # Add reason strings
        violations = []
        for v in merged:
            reasons = []
            if v['peak_az_rate'] > max_az_rate:
                reasons.append("az_rate={0:.4f} > {1}".format(v['peak_az_rate'], max_az_rate))
            if v['peak_el_rate'] > max_el_rate:
                reasons.append("el_rate={0:.4f} > {1}".format(v['peak_el_rate'], max_el_rate))
            if v['peak_elev'] > max_el_keyhole:
                reasons.append("elev={0:.1f} > keyhole {1}".format(v['peak_elev'], max_el_keyhole))
            v['reason'] = ', '.join(reasons)
            violations.append(v)

        # Build trackable windows from gaps between merged violations
        trackable_windows = []
        cursor = start_t
        for v in violations:
            if (v['start'] - cursor).total_seconds() >= min_window_sec:
                trackable_windows.append((cursor, v['start']))
            cursor = v['end']
        if (end_t - cursor).total_seconds() >= min_window_sec:
            trackable_windows.append((cursor, end_t))

        # Check slew feasibility between windows
        feasible_windows = []
        for i, (ws, we) in enumerate(trackable_windows):
            if i == 0:
                feasible_windows.append((ws, we))
                continue

            prev_end = feasible_windows[-1][1] if feasible_windows else ws
            prev_az, prev_el = self.SatPositionAt(prev_end)
            tgt_az, tgt_el = self.SatPositionAt(ws)
            gap_sec = (ws - prev_end).total_seconds()

            if gap_sec > 0:
                az_slew_rate = abs(_az_diff(tgt_az, prev_az)) / gap_sec   # <-- changed
                el_slew_rate = abs(tgt_el - prev_el) / gap_sec
            else:
                az_slew_rate = float('inf')
                el_slew_rate = float('inf')

            if az_slew_rate <= max_az_rate and el_slew_rate <= max_el_rate:
                feasible_windows.append((ws, we))
            else:
                new_start = None
                window_dur = (we - ws).total_seconds()
                for offset in range(int(sample_sec), int(window_dur), int(max(sample_sec, 1))):
                    candidate = ws + timedelta(seconds=offset)
                    cand_az, cand_el = self.SatPositionAt(candidate)
                    slew_time = (candidate - prev_end).total_seconds()
                    if slew_time <= 0:
                        continue
                    az_sr = abs(_az_diff(cand_az, prev_az)) / slew_time   
                    el_sr = abs(cand_el - prev_el) / slew_time
                    if az_sr <= max_az_rate and el_sr <= max_el_rate:
                        new_start = candidate
                        break

                if new_start and (we - new_start).total_seconds() >= min_window_sec:
                    feasible_windows.append((new_start, we))

        return violations, feasible_windows

def GetNextPass(ephem_file, minEl, now=None, max_az_rate=MAX_AZ_RATE, max_el_rate=MAX_EL_RATE):
    """Find the next pass and check rate limits.

    Always returns a NextPass object. If rate limits are exceeded,
    the violations and trackable_windows are stored on the object
    and a warning is printed. Check rtn.has_violations() to see
    if there were issues.
    """
    orbit = OrbitPasses(ephem_file)
    points = orbit.points

    # Resolve start time
    if now is None or (isinstance(now, str) and now.strip().lower() == 'ephem'):
        now = points[0][0]
    elif isinstance(now, str) and now.strip().lower() == 'now':
        now = datetime.now(UTC)
    elif isinstance(now, str):
        parsed = _parse_datetime(*now.split(' ', 1))
        if parsed is None:
            raise ValueError("Could not parse start time: '{0}'".format(now))
        now = parsed.replace(tzinfo=UTC)

    orbit = OrbitPasses(ephem_file)
    points = orbit.points

    # Find rise
    rise = None
    rise_idx = None
    for i, (t, az, el) in enumerate(points):
        if t < now:
            continue
        if el >= minEl:
            if i == 0 or points[i - 1][2] < minEl:
                rise = (t, az, el)
                rise_idx = i
                break
            else:
                rise = (t, az, el)
                rise_idx = i
                break

    if rise is None:
        raise Exception("No rise event found in ephemeris data")

    # Find set
    sett = None
    for i in range(rise_idx + 1, len(points)):
        t, az, el = points[i]
        if el < minEl:
            sett = points[i - 1]
            break

    if sett is None:
        sett = points[-1]

    risetime = rise[0]
    settime = sett[0]

    rtn = NextPass(risetime, rise[1], rise[2], settime, sett[1], sett[2])
    rtn._orbit = orbit

    # Find midpoint
    mid_offset = timedelta(seconds=rtn.pass_duration() / 2.0)
    t2 = rtn.start_time() + mid_offset
    azz, ell = orbit.SatPositionAt(t2)
    rtn.setmidpt(t2, azz, ell)

    # Check rate limits
    violations, trackable_windows = orbit.check_rate_limits(
        risetime, settime, max_az_rate, max_el_rate
    )

    rtn.violations = violations
    if violations:
        rtn.trackable_windows = trackable_windows

        msg = "WARNING: Antenna rate limits exceeded during pass for {0}\n".format(orbit.name)
        msg += "  Pass: {0} - {1} UTC\n".format(risetime.strftime('%H:%M:%S'), settime.strftime('%H:%M:%S'))
        msg += "  Limits: az_rate={0} deg/s, el_rate={1} deg/s\n".format(max_az_rate, max_el_rate)
        msg += "\n  Violations ({0}):\n".format(len(violations))
        for i, v in enumerate(violations, 1):
            msg += "    {0}. {1} - {2} ".format(i, v['start'].strftime('%H:%M:%S'), v['end'].strftime('%H:%M:%S'))
            msg += "({0})\n".format(v['reason'])
        msg += "\n  Proposed trackable windows ({0}):\n".format(len(trackable_windows))
        if trackable_windows:
            for i, (ws, we) in enumerate(trackable_windows, 1):
                dur = (we - ws).total_seconds()
                s_az, s_el = orbit.SatPositionAt(ws)
                e_az, e_el = orbit.SatPositionAt(we)
                msg += "    {0}. TRACK {1} - {2} ".format(i, ws.strftime('%H:%M:%S'), we.strftime('%H:%M:%S'))
                msg += "({0:.0f}s) ".format(dur)
                msg += "Az={0:.1f}->{1:.1f} El={2:.1f}->{3:.1f}\n".format(s_az, e_az, s_el, e_el)
        else:
            msg += "    ** No trackable windows -- pass exceeds rate limits "
            msg += "or antenna cannot slew to reachable positions **\n"

        rtn.rate_warning = msg
        print(msg)
    else:
        rtn.trackable_windows = [(risetime, settime)]

    # Attach sampled trajectory for plotting
    rtn.trajectory = orbit.sample_trajectory(risetime, settime)

    return rtn


def print_summary(rtn):
    """Print pass summary info and commands."""
    print("Now:", datetime.now(UTC))
    print("Time until  (hh:mm:ss):", rtn.time_until_rise(None))
    print("Time/Target of aquisition", rtn.start_time(), rtn.start_az(), rtn.start_el())
    print("Pass at midpoint", rtn.midpoint_time(), rtn.midpoint_az(), rtn.midpoint_el())
    print("Time/Target at set ", rtn.end_time(), rtn.end_az(), rtn.end_el())
    print("Pass duration (s)", rtn.pass_duration())
    print("Az wrap:", rtn.which_wrap())
    if rtn.has_violations():
        print("Rate limits EXCEEDED -- see warnings above")
    else:
        print("Rate limits OK")
    print()
    print("Commands:")
    rtn.print_commands()


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Compute satellite pass from ephemeris')
    parser.add_argument('ephem_file', help='Ephemeris text file')
    #parser.add_argument('name', nargs='?', default='target', help='Target name')
    parser.add_argument('--min-el', type=float, default=10.0, help='Minimum elevation (deg)')
    parser.add_argument("--start-time", default=None,
                        help="Start time (YYYY-MM-DD HH:MM:SS UTC), 'now', "
                             "or omit to use ephemeris start time")
    args = parser.parse_args()

    rtn = GetNextPass(args.ephem_file, args.min_el, args.start_time)
    print_summary(rtn)