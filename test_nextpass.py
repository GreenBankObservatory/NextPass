# -*- coding: utf-8 -*-
"""Unit tests for NextPass.py. Run: python -m unittest test_nextpass -v"""
from __future__ import print_function
import os, shutil, tempfile, unittest
from datetime import datetime, timedelta
from NextPass import GetNextPass, MAX_AZ_RATE, MAX_EL_RATE, RATE_LIMIT_BUFFER, UTC, _az_diff

EFF_AZ = MAX_AZ_RATE * RATE_LIMIT_BUFFER
EFF_EL = MAX_EL_RATE * RATE_LIMIT_BUFFER
T0 = datetime(2026, 1, 1, 0, 0, 0, tzinfo=UTC)


def _dms(x):
    s = '-' if x < 0 else ''
    a = abs(x); d = int(a); mf = (a - d) * 60.0; m = int(mf); sec = (mf - m) * 60.0
    return "{0}{1:02d}:{2:02d}:{3:05.2f}".format(s, d, m, sec)


def write_ephem(path, name, samples, step=1.0):
    with open(path, 'w') as f:
        f.write("FORMAT = EPHEMERIS\nCOORDMODE = azel\nNAME = {0}\nHEAD = date utc az el\n".format(name))
        for i, (az, el) in enumerate(samples):
            t = T0 + timedelta(seconds=i * step)
            f.write("{0} {1} {2} {3}\n".format(
                t.strftime('%Y-%m-%d'), t.strftime('%H:%M:%S.%f')[:-3],
                _dms(az % 360.0), _dms(el)))


class Base(unittest.TestCase):
    def setUp(self): self.tmp = tempfile.mkdtemp(prefix="nextpass_")
    def tearDown(self): shutil.rmtree(self.tmp, ignore_errors=True)
    def path(self, name): return os.path.join(self.tmp, name + ".ephem")


class TestNextPass(Base):

    def test_az_diff_helper(self):
        # _az_diff returns signed short-way difference across the 0/360 seam
        self.assertAlmostEqual(_az_diff(1.0, 359.0), 2.0)
        self.assertAlmostEqual(_az_diff(359.0, 1.0), -2.0)
        self.assertAlmostEqual(_az_diff(10.0, 5.0), 5.0)

    def test_wrap_crossing_zero(self):
        # Slow pass crossing 0/360 must not trigger a phantom rate violation
        p = self.path("wrap")
        write_ephem(p, "wrap", [(355.0 + 0.03 * i, 15.0 + 0.05 * i) for i in range(300)])
        r = GetNextPass(p, 10.0)
        self.assertFalse(r.has_violations())
        self.assertEqual(len(r.trackable_windows), 1)
        self.assertLess(max(abs(x) for x in r.trajectory['az_rate']), EFF_AZ)

    def test_az_rate_violation(self):
        # Az motion far above limit (1.0 deg/s vs 0.54) must be flagged
        p = self.path("az")
        write_ephem(p, "az", [(100.0 + 1.0 * i, 30.0) for i in range(200)])
        r = GetNextPass(p, 10.0)
        self.assertTrue(r.has_violations())
        self.assertIn("az_rate", " ".join(v['reason'] for v in r.violations))

    def test_el_rate_violation(self):
        # El motion far above limit (0.5 deg/s vs 0.27) must be flagged
        p = self.path("el")
        write_ephem(p, "el", [(90.0, min(15.0 + 0.5 * i, 70.0)) for i in range(100)])
        r = GetNextPass(p, 10.0)
        self.assertTrue(r.has_violations())
        self.assertIn("el_rate", " ".join(v['reason'] for v in r.violations))

    def test_keyhole_violation(self):
        # Elevation > 80 deg must be flagged as keyhole even when rates are tiny
        p = self.path("kh")
        write_ephem(p, "kh", [(180.0, min(70.0 + 0.05 * i, 85.0)) for i in range(400)])
        r = GetNextPass(p, 10.0)
        self.assertTrue(r.has_violations())
        self.assertIn("keyhole", " ".join(v['reason'] for v in r.violations))

    def test_min_horizon(self):
        # Rise/set must be trimmed to minEl crossings, not the ephemeris extremes
        p = self.path("hz"); n = 200
        samples = [(90.0 + 0.05 * i,
                    -5.0 + (i / (n / 2.0)) * 35.0 if i <= n / 2
                    else 30.0 - ((i - n / 2.0) / (n / 2.0)) * 35.0) for i in range(n)]
        write_ephem(p, "hz", samples)
        r = GetNextPass(p, 10.0)
        self.assertGreaterEqual(r.start_el(), 10.0)
        self.assertGreaterEqual(r.end_el(), 10.0)
        self.assertLess(r.pass_duration(), (n - 1) * 1.0)

    def test_normal_pass(self):
        # Gentle pass within all limits yields exactly one full-span trackable window
        p = self.path("ok")
        samples = [(45.0 + 0.1 * i, 20.0 + 30.0 * (1 - (2 * i / 600.0 - 1) ** 2)) for i in range(600)]
        write_ephem(p, "ok", samples)
        r = GetNextPass(p, 10.0)
        self.assertFalse(r.has_violations())
        self.assertEqual(len(r.trackable_windows), 1)
        self.assertEqual(r.trackable_windows[0], (r.start_time(), r.end_time()))

    def test_multiple_windows(self):
        # Mid-pass violation splits into >=2 non-overlapping windows with 3 commands each
        p = self.path("split"); samples = []; az = 30.0
        for i in range(500):
            az += 1.5 if 200 <= i < 260 else 0.05
            samples.append((az, 30.0))
        write_ephem(p, "split", samples)
        r = GetNextPass(p, 10.0)
        self.assertTrue(r.has_violations())
        self.assertGreaterEqual(len(r.trackable_windows), 2)
        for a, b in zip(r.trackable_windows, r.trackable_windows[1:]):
            self.assertLessEqual(a[1], b[0])
        self.assertEqual(len(r.generate_commands()), 3 * len(r.trackable_windows))


if __name__ == "__main__":
    unittest.main(verbosity=2)