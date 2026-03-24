from datetime import datetime, timedelta, timezone


# Default antenna rate limits (deg/sec)
MAX_EL_RATE = 0.3
MAX_AZ_RATE = 0.6

# Observer location
OBSERVER_LAT = 38.050417
OBSERVER_LON = -81.104889
OBSERVER_ELEV = 517.0


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
            when = datetime.now(timezone.utc)
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
        for ws, we in self.trackable_windows:
            az, el = self._orbit.SatPositionAt(ws)
            dur = (we - ws).total_seconds()
            commands.append(f'Slew(Location("AzEl", {az:.1f}, {el:.1f}))')
            commands.append(f"WaitFor('{ws.strftime('%H:%M:%S')}')")
            commands.append(
                f'Track(Location("AzEl", {az:.1f}, {el:.1f}), None, {dur:.0f})'
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
    parts = s.strip().split(':')
    d = float(parts[0])
    m = float(parts[1]) if len(parts) > 1 else 0.0
    sec = float(parts[2]) if len(parts) > 2 else 0.0
    sign = -1 if d < 0 else 1
    return sign * (abs(d) + m / 60.0 + sec / 3600.0)


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
    import warnings
    warnings.filterwarnings('ignore', module='astropy')
    from astropy.utils.iers import conf as iers_conf
    iers_conf.auto_download = False
    from astropy.coordinates import SkyCoord, EarthLocation, AltAz
    from astropy.time import Time
    import astropy.units as u

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


def parse_ephemeris(filename, observer_lat=OBSERVER_LAT, observer_lon=OBSERVER_LON,
                    observer_elev=OBSERVER_ELEV):
    """Parse an ephemeris text file (azel or J2000 RA/DEC).

    Reads COORDMODE from the header to determine format.
    For J2000, converts RA/DEC to AzEl using the observer location.

    Returns a list of (datetime_utc, az_deg, el_deg) tuples.
    """
    coordmode = 'azel'
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
                val = val.strip().lower()
                if key == 'COORDMODE':
                    coordmode = val
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
            dt = dt.replace(tzinfo=timezone.utc)
            az = sexagesimal_to_decimal(parts[2])
            el = sexagesimal_to_decimal(parts[3])
            points.append((dt, az, el))

    elif coordmode in ('j2000', 'radec'):
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
            dt = dt.replace(tzinfo=timezone.utc)
            ra = hms_to_degrees(parts[2])
            dec = sexagesimal_to_decimal(parts[3])
            times_utc.append(dt)
            ra_degs.append(ra)
            dec_degs.append(dec)

        if times_utc:
            az_list, el_list = _radec_to_azel(
                times_utc, ra_degs, dec_degs,
                observer_lat, observer_lon, observer_elev
            )
            for dt, az, el in zip(times_utc, az_list, el_list):
                points.append((dt, az, el))
    else:
        raise Exception(f"Unknown COORDMODE: {coordmode}")

    return points


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
            az = az0 + frac * (az1 - az0)
            el = el0 + frac * (el1 - el0)
            return az, el
    return points[-1][1], points[-1][2]


class OrbitPasses:
    """Load a pre-computed az/el ephemeris and provide position lookups."""

    def __init__(self, ephem_file, name=None):
        self.name = name
        self.points = parse_ephemeris(ephem_file)
        if not self.points:
            raise Exception(f"No data points found in {ephem_file}")

    def SatPositionAt(self, t):
        """Return (az_deg, el_deg) at datetime t via interpolation."""
        return interp_at(self.points, t)

    def VelocityAt(self, t, dt_sec=1.0):
        """Return (az_rate, el_rate) in deg/sec at datetime t via central difference."""
        delta = timedelta(seconds=dt_sec)
        az0, el0 = interp_at(self.points, t - delta)
        az1, el1 = interp_at(self.points, t + delta)
        az_rate = (az1 - az0) / (2.0 * dt_sec)
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

        return {
            'times': times,
            'secs': secs,
            'az': azs,
            'el': els,
            'az_rate': az_rates,
            'el_rate': el_rates,
        }

    def check_rate_limits(self, start_t, end_t, max_az_rate=MAX_AZ_RATE,
                          max_el_rate=MAX_EL_RATE, sample_sec=1.0,
                          min_window_sec=10.0):
        """Check if antenna rate limits are exceeded during a time window.

        Returns:
            violations: list of dicts with 'start', 'end', 'peak_az_rate',
                        'peak_el_rate', 'reason' for each violation window.
            trackable_windows: list of (start, end) tuples where rates are
                               within limits — the safe windows for tracking.
        """
        duration = (end_t - start_t).total_seconds()
        n_samples = max(int(duration / sample_sec), 2)
        step = duration / n_samples

        # Sample rates across the pass
        exceeded = []
        for i in range(n_samples + 1):
            t = start_t + timedelta(seconds=i * step)
            az_rate, el_rate = self.VelocityAt(t)
            over = abs(az_rate) > max_az_rate or abs(el_rate) > max_el_rate
            exceeded.append((t, over, az_rate, el_rate))

        # Group contiguous violations into raw windows
        raw_violations = []
        in_violation = False
        v_start = None
        peak_az = 0.0
        peak_el = 0.0

        for t, over, az_r, el_r in exceeded:
            if over and not in_violation:
                in_violation = True
                v_start = t
                peak_az = abs(az_r)
                peak_el = abs(el_r)
            elif over and in_violation:
                peak_az = max(peak_az, abs(az_r))
                peak_el = max(peak_el, abs(el_r))
            elif not over and in_violation:
                in_violation = False
                raw_violations.append({
                    'start': v_start, 'end': t,
                    'peak_az_rate': peak_az, 'peak_el_rate': peak_el,
                })

        if in_violation:
            raw_violations.append({
                'start': v_start, 'end': exceeded[-1][0],
                'peak_az_rate': peak_az, 'peak_el_rate': peak_el,
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
                else:
                    merged.append(v.copy())

        # Add reason strings
        violations = []
        for v in merged:
            reasons = []
            if v['peak_az_rate'] > max_az_rate:
                reasons.append(f"az_rate={v['peak_az_rate']:.4f} > {max_az_rate}")
            if v['peak_el_rate'] > max_el_rate:
                reasons.append(f"el_rate={v['peak_el_rate']:.4f} > {max_el_rate}")
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
                az_slew_rate = abs(tgt_az - prev_az) / gap_sec
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
                    az_sr = abs(cand_az - prev_az) / slew_time
                    el_sr = abs(cand_el - prev_el) / slew_time
                    if az_sr <= max_az_rate and el_sr <= max_el_rate:
                        new_start = candidate
                        break

                if new_start and (we - new_start).total_seconds() >= min_window_sec:
                    feasible_windows.append((new_start, we))

        return violations, feasible_windows


def GetNextPass(ephem_file, name, minEl, now=None, trange_hrs=12.0,
                max_az_rate=MAX_AZ_RATE, max_el_rate=MAX_EL_RATE):
    """Find the next pass and check rate limits.

    Always returns a NextPass object. If rate limits are exceeded,
    the violations and trackable_windows are stored on the object
    and a warning is printed. Check rtn.has_violations() to see
    if there were issues.
    """
    if now is None:
        now = datetime.now(timezone.utc)

    orbit = OrbitPasses(ephem_file, name)
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

        msg = f"WARNING: Antenna rate limits exceeded during pass for {name}\n"
        msg += f"  Pass: {risetime.strftime('%H:%M:%S')} - {settime.strftime('%H:%M:%S')} UTC\n"
        msg += f"  Limits: az_rate={max_az_rate} deg/s, el_rate={max_el_rate} deg/s\n"
        msg += f"\n  Violations ({len(violations)}):\n"
        for i, v in enumerate(violations, 1):
            msg += f"    {i}. {v['start'].strftime('%H:%M:%S')} - {v['end'].strftime('%H:%M:%S')} "
            msg += f"({v['reason']})\n"
        msg += f"\n  Proposed trackable windows ({len(trackable_windows)}):\n"
        if trackable_windows:
            for i, (ws, we) in enumerate(trackable_windows, 1):
                dur = (we - ws).total_seconds()
                s_az, s_el = orbit.SatPositionAt(ws)
                e_az, e_el = orbit.SatPositionAt(we)
                msg += f"    {i}. TRACK {ws.strftime('%H:%M:%S')} - {we.strftime('%H:%M:%S')} "
                msg += f"({dur:.0f}s) "
                msg += f"Az={s_az:.1f}->{e_az:.1f} El={s_el:.1f}->{e_el:.1f}\n"
        else:
            msg += "    ** No trackable windows — pass exceeds rate limits "
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
    print("Now:", datetime.now(timezone.utc))
    print("Time until rise:", rtn.time_until_rise(None))
    print("Time/Target of rise", rtn.start_time(), rtn.start_az(), rtn.start_el())
    print("Pass at midpoint", rtn.midpoint_time(), rtn.midpoint_az(), rtn.midpoint_el())
    print("Time/Target at set ", rtn.end_time(), rtn.end_az(), rtn.end_el())
    print("Pass duration", rtn.pass_duration())
    print("Az wrap:", rtn.which_wrap())
    if rtn.has_violations():
        print("Rate limits EXCEEDED — see warnings above")
    else:
        print("Rate limits OK")
    print()
    print("Commands:")
    rtn.print_commands()


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Compute satellite pass from ephemeris')
    parser.add_argument('ephem_file', help='Ephemeris text file')
    parser.add_argument('name', nargs='?', default='target', help='Target name')
    parser.add_argument('--min-el', type=float, default=10.0, help='Minimum elevation (deg)')
    parser.add_argument('--range', type=float, default=24.0, dest='trange_hrs',
                        help='Search range (hours)')
    args = parser.parse_args()

    rtn = GetNextPass(args.ephem_file, args.name, args.min_el, None, args.trange_hrs)
    print_summary(rtn)