import sys
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime, timezone
from NextPass import GetNextPass, MAX_AZ_RATE, MAX_EL_RATE

def plot_pass(ephem_file, target_name, min_el=5.0, trange_hrs=24.0):
    now = datetime.now(timezone.utc)
    rtn = GetNextPass(ephem_file, target_name, min_el, now, trange_hrs)

    print(f"Rise:  {rtn.start_time()} Az={rtn.start_az():.1f} El={rtn.start_el():.1f}")
    print(f"Mid:   {rtn.midpoint_time()} Az={rtn.midpoint_az():.1f} El={rtn.midpoint_el():.1f}")
    print(f"Set:   {rtn.end_time()} Az={rtn.end_az():.1f} El={rtn.end_el():.1f}")
    print(f"Duration: {rtn.pass_duration():.0f}s  Wrap: {rtn.which_wrap()}")
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
    fig.suptitle(f"{target_name} — {rtn.start_time().strftime('%Y-%m-%d %H:%M:%S')} UTC")

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
    ax1.axhline(y=min_el, color='r', linestyle='--', label=f'Min El ({min_el}°)')
    ax1.set_xlabel('Seconds from rise')
    ax1.set_ylabel('Elevation (°)')
    ax1.set_title('Elevation vs Time')
    shade_violations(ax1)
    ax1.grid(True)

    # Row 1: Az vs Time
    ax2 = fig.add_subplot(2, 3, 2)
    ax2.plot(secs, azs)
    ax2.set_xlabel('Seconds from rise')
    ax2.set_ylabel('Azimuth (°)')
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
    ax_polar.set_title('Az/El Polar', pad=15)
    ax_polar.legend(loc='lower right')

    # Row 2: El rate vs Time
    ax3 = fig.add_subplot(2, 3, 4)
    ax3.plot(secs, el_rates)
    ax3.axhline(y=0, color='k', linestyle='-', linewidth=0.5)
    ax3.axhline(y=MAX_EL_RATE, color='r', linestyle='--', label=f'Max ({MAX_EL_RATE}°/s)')
    ax3.axhline(y=-MAX_EL_RATE, color='r', linestyle='--')
    ax3.set_xlabel('Seconds from rise')
    ax3.set_ylabel('El Rate (°/s)')
    ax3.set_title('Elevation Rate vs Time')
    shade_violations(ax3)
    ax3.grid(True)

    # Row 2: Az rate vs Time
    ax4 = fig.add_subplot(2, 3, 5)
    ax4.plot(secs, az_rates)
    ax4.axhline(y=0, color='k', linestyle='-', linewidth=0.5)
    ax4.axhline(y=MAX_AZ_RATE, color='r', linestyle='--', label=f'Max ({MAX_AZ_RATE}°/s)')
    ax4.axhline(y=-MAX_AZ_RATE, color='r', linestyle='--')
    ax4.set_xlabel('Seconds from rise')
    ax4.set_ylabel('Az Rate (°/s)')
    ax4.set_title('Azimuth Rate vs Time')
    shade_violations(ax4)
    ax4.grid(True)

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    ephem_file = sys.argv[1] if len(sys.argv) > 1 else 'test_ephem.txt'
    target_name = sys.argv[2] if len(sys.argv) > 2 else 'target'
    plot_pass(ephem_file, target_name)
