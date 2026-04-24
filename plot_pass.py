# -*- coding: utf-8 -*-
from __future__ import print_function
import sys
import warnings
warnings.filterwarnings('ignore', category=DeprecationWarning, module='matplotlib')
warnings.filterwarnings('ignore', message='.*dubious year.*')
try:
    from numpy import VisibleDeprecationWarning
    warnings.filterwarnings('ignore', category=VisibleDeprecationWarning)
except ImportError:
    pass
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
from NextPass import GetNextPass, MAX_AZ_RATE, MAX_EL_RATE, UTC
def plot_pass(ephem_file, min_el=10.0, start=None, max_az_rate=MAX_AZ_RATE, max_el_rate=MAX_EL_RATE,
              output=None):
    rtn = GetNextPass(ephem_file, min_el, start, max_az_rate, max_el_rate)
    print("Rise:  {0} Az={1:.1f} El={2:.1f}".format(rtn.start_time(), rtn.start_az(), rtn.start_el()))
    print("Mid:   {0} Az={1:.1f} El={2:.1f}".format(rtn.midpoint_time(), rtn.midpoint_az(), rtn.midpoint_el()))
    print("Set:   {0} Az={1:.1f} El={2:.1f}".format(rtn.end_time(), rtn.end_az(), rtn.end_el()))
    print("Duration: {0:.0f}s  Wrap: {1}".format(rtn.pass_duration(), rtn.which_wrap()))
    print()
    print("Commands:")
    rtn.print_commands()
    # Unpack pre-computed trajectory
    traj = rtn.trajectory
    secs = traj['secs']
    azs = traj['az']
    els = traj['el']
    az_rates = traj['az_rate']
    el_rates = traj['el_rate']
    fig = plt.figure(figsize=(15, 9))
    fig.suptitle("{0} -- {1} UTC".format(rtn._orbit.name, rtn.start_time().strftime('%Y-%m-%d %H:%M:%S')))
    # Helper to shade violation regions on an axis
    t0 = rtn.start_time()
    def shade_violations(ax):
        for v in rtn.violations:
            s0 = (v['start'] - t0).total_seconds()
            s1 = (v['end'] - t0).total_seconds()
            ax.axvspan(s0, s1, alpha=0.2, color='red', label='Rate exceeded')
        handles, labels = ax.get_legend_handles_labels()
        seen = set()
        unique = [(h, l) for h, l in zip(handles, labels) if l not in seen and not seen.add(l)]
        if unique:
            ax.legend(*zip(*unique))
    # Row 1: El vs Time
    ax1 = fig.add_subplot(2, 3, 1)
    ax1.plot(secs, els)
    ax1.axhline(y=min_el, color='r', linestyle='--', label='Min El ({0})'.format(min_el))
    ax1.set_xlabel('Seconds from rise')
    ax1.set_ylabel('Elevation (deg)')
    ax1.set_title('Elevation vs Time')
    shade_violations(ax1)
    ax1.grid(True)
    # Row 1: Az vs Time
    ax2 = fig.add_subplot(2, 3, 2)
    ax2.plot(secs, azs)
    ax2.set_xlabel('Seconds from rise')
    ax2.set_ylabel('Azimuth (deg)')
    ax2.set_title('Azimuth vs Time')
    shade_violations(ax2)
    ax2.grid(True)
    # Row 1: Polar plot (az/el)
    ax_polar = fig.add_subplot(2, 3, 3, projection='polar')
    az_rad = np.radians(azs)
    ax_polar.plot(az_rad, els)
    ax_polar.plot(az_rad[0], els[0], 'go', markersize=8, label='Rise')
    ax_polar.plot(az_rad[-1], els[-1], 'rs', markersize=8, label='Set')
    ax_polar.set_theta_zero_location('N')
    ax_polar.set_theta_direction(-1)
    ax_polar.set_title('\nAz/El Polar')
    ax_polar.legend(loc='lower right')
    # Row 2: El rate vs Time
    ax3 = fig.add_subplot(2, 3, 4)
    ax3.plot(secs, el_rates)
    ax3.axhline(y=0, color='k', linestyle='-', linewidth=0.5)
    ax3.axhline(y=MAX_EL_RATE, color='r', linestyle='--', label='Max ({0} deg/s)'.format(MAX_EL_RATE))
    ax3.axhline(y=-MAX_EL_RATE, color='r', linestyle='--')
    ax3.set_xlabel('Seconds from rise')
    ax3.set_ylabel('El Rate (deg/s)')
    ax3.set_title('Elevation Rate vs Time')
    shade_violations(ax3)
    ax3.grid(True)
    # Row 2: Az rate vs Time
    ax4 = fig.add_subplot(2, 3, 5)
    ax4.plot(secs, az_rates)
    ax4.axhline(y=0, color='k', linestyle='-', linewidth=0.5)
    ax4.axhline(y=MAX_AZ_RATE, color='r', linestyle='--', label='Max ({0} deg/s)'.format(MAX_AZ_RATE))
    ax4.axhline(y=-MAX_AZ_RATE, color='r', linestyle='--')
    ax4.set_xlabel('Seconds from rise')
    ax4.set_ylabel('Az Rate (deg/s)')
    ax4.set_title('Azimuth Rate vs Time')
    shade_violations(ax4)
    ax4.grid(True)
    plt.tight_layout()
    plt.show()
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Plot satellite pass from ephemeris')
    parser.add_argument('ephem_file', help='Ephemeris text file')
    #parser.add_argument('name', nargs='?', default='target', help='Target name')
    parser.add_argument('--min-el', type=float, default=10.0, help='Minimum elevation (deg)')
    parser.add_argument('--start-time', default=None,
                        help="'now', 'ephem', or YYYY-MM-DD HH:MM:SS (default: ephem start)")
    parser.add_argument('--max-az-rate', type=float, default=MAX_AZ_RATE,
                        help='Max azimuth rate deg/s (default: {0})'.format(MAX_AZ_RATE))
    parser.add_argument('--max-el-rate', type=float, default=MAX_EL_RATE,
                        help='Max elevation rate deg/s (default: {0})'.format(MAX_EL_RATE))
    parser.add_argument('-o', '--output', default=None,
                        help='Save plot to file instead of displaying')
    args = parser.parse_args()

    plot_pass(args.ephem_file, args.min_el, args.start_time, args.max_az_rate, args.max_el_rate, args.output)