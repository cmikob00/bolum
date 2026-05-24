'''
bolum plotting file
'''
from bolum.constants import *

# plotting function — simplified to single luminosity
def plot_outputs(basename, header, times, taus, lums, alts, frag_times, frag_alts,
                 velocities, accels, qs, T_stags, p_stags, T_surfs, init_strength,
                 E_rad_totals, peak_power, peak_rad_i, E_rad_tot, E_event):

    owtname = basename + f'/bolide_plots.pdf'
    pdf_pages = PdfPages(owtname)

    # 1. luminosity vs time (linear/linear)
    fig1 = plt.figure(figsize=(9.0, 6.5))
    xmin = -0.1
    xmax = 1.05 * times[-1]
    ymin = 0.8 * lums[0]
    ymax = 2. * peak_power
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.plot(times, lums, color='blue', label='Total Luminosity')
    plt.xscale('linear')
    plt.yscale('linear')
    plt.xlabel('Time (seconds)')
    plt.ylabel('Luminosity (Watts)')
    plt.gcf().text(.20, .89, header, fontsize=10, color='black')
    if frag_times:
        plt.axvline(frag_times[0], linestyle='--', color='green', alpha=0.6, label=f"Frag at {frag_alts[0]/1000.:.1f} km and {frag_times[0]:.2f} sec")
    plt.gcf().text(.14, .74, f"Total Rad E = {E_rad_tot / jkt:.3f} kt", fontsize=10, color='black')
    plt.gcf().text(.14, .71, f"Total Rad E = {E_rad_tot:.3e} J", fontsize=10, color='black')
    plt.gcf().text(.14, .68, f"Approx Event E = {E_event / jkt:.3f} kt", fontsize=10, color='black')
    plt.gcf().text(.14, .65, f"Approx Event E = {E_event:.3e} J", fontsize=10, color='black')
    plt.axhline(peak_power, linestyle='--', color='red', label=f"Peak = {peak_power:.3e} W")
    plt.title(f'Bolide Light Curve: Luminosity vs Time\n')
    plt.grid(True)
    plt.legend(loc='upper left')
    plt.gcf().text(.10, .95, 'UNCLASSIFIED', fontsize=10, color='green')
    plt.gcf().text(.80, .05, 'UNCLASSIFIED', fontsize=10, color='green')
    pdf_pages.savefig(fig1)

    # 2. luminosity in Watts vs time in log/log scale
    fig2 = plt.figure(figsize=(9.0, 6.5))
    xmin = 0.003
    xmax = 1.05 * times[-1]
    ymin = 0.8 * lums[0]
    ymax = 2. * peak_power
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.plot(times, lums, color='blue', label='Total Luminosity')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Time (seconds)')
    plt.ylabel('Luminosity (Watts)')
    plt.gcf().text(.20, .89, header, fontsize=10, color='black')
    if frag_times:
         plt.axvline(frag_times[0], linestyle='--', color='green', alpha=0.6, label=f"Frag at {frag_alts[0]/1000.:.1f} km and {frag_times[0]:.2f} sec")
    plt.gcf().text(.14, .74, f"Total Rad E = {E_rad_tot / jkt:.3f} kt", fontsize=10, color='black')
    plt.gcf().text(.14, .71, f"Total Rad E = {E_rad_tot:.3e} J", fontsize=10, color='black')
    plt.gcf().text(.14, .68, f"Approx Event E = {E_event / jkt:.3f} kt", fontsize=10, color='black')
    plt.gcf().text(.14, .65, f"Approx Event E = {E_event:.3e} J", fontsize=10, color='black')
    plt.axhline(peak_power, linestyle='--', color='red', label=f"Peak Power = {peak_power:.3e} W")
    plt.title(f'Bolide Light Curve: Luminosity vs Time\n')
    plt.grid(True)
    plt.legend(loc='upper left')
    plt.gcf().text(.10, .95, 'UNCLASSIFIED', fontsize=10, color='green')
    plt.gcf().text(.80, .05, 'UNCLASSIFIED', fontsize=10, color='green')
    pdf_pages.savefig(fig2)

    # 3. radiant intensity in W/sr vs time in linear/linear scale
    fig3 = plt.figure(figsize=(9.0, 6.5))
    # convert luminosity in watts to W/sr
    radiant_intensities = lums / (4. * pi)
    xmin = -0.1
    xmax = 1.05 * times[-1]
    ymin = 0.8 * radiant_intensities[0]
    ymax = 2. * peak_rad_i
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.plot(times, radiant_intensities, color='blue', label='Total Rad I')
    plt.xscale('linear')
    plt.yscale('linear')
    plt.xlabel('Time (seconds)')
    plt.ylabel('Radiant Intensity (W/sr)')
    plt.gcf().text(.20, .89, header, fontsize=10, color='black')
    if frag_times:
         plt.axvline(frag_times[0], linestyle='--', color='green', alpha=0.6, label=f"Frag at {frag_alts[0]/1000.:.1f} km and {frag_times[0]:.2f} sec")
    plt.gcf().text(.14, .74, f"Total Rad E = {E_rad_tot / jkt:.3f} kt", fontsize=10, color='black')
    plt.gcf().text(.14, .71, f"Total Rad E = {E_rad_tot:.3e} J", fontsize=10, color='black')
    plt.gcf().text(.14, .68, f"Approx Event E = {E_event / jkt:.3f} kt", fontsize=10, color='black')
    plt.gcf().text(.14, .65, f"Approx Event E = {E_event:.3e} J", fontsize=10, color='black')
    plt.axhline(peak_rad_i, linestyle='--', color='red', label=f"Peak Rad I = {peak_rad_i:.3e} W/sr")
    plt.title(f'Bolide Light Curve: Radiant Intensity vs Time\n')
    plt.grid(True)
    plt.legend(loc='upper left')
    plt.gcf().text(.10, .95, 'UNCLASSIFIED', fontsize=10, color='green')
    plt.gcf().text(.80, .05, 'UNCLASSIFIED', fontsize=10, color='green')
    pdf_pages.savefig(fig3)

    # 4. radiant intensity in W/sr vs time in log/log scale
    fig4 = plt.figure(figsize=(9.0, 6.5))
    # convert luminosity in watts to W/sr
    radiant_intensities = lums / (4. * pi)
    xmin = 0.003
    xmax = 1.05 * times[-1]
    ymin = 0.8 * radiant_intensities[0]
    ymax = 2. * peak_rad_i
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.plot(times, radiant_intensities, color='blue', label='Total Rad I')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Time (seconds)')
    plt.ylabel('Radiant Intensity (W/sr)')
    plt.gcf().text(.20, .89, header, fontsize=10, color='black')
    if frag_times:
        plt.axvline(frag_times[0], linestyle='--', color='green', alpha=0.6, label=f"Frag at {frag_alts[0]/1000.:.1f} km and {frag_times[0]:.2f} sec")
    plt.gcf().text(.14, .74, f"Total Rad E = {E_rad_tot / jkt:.3f} kt", fontsize=10, color='black')
    plt.gcf().text(.14, .71, f"Total Rad E = {E_rad_tot:.3e} J", fontsize=10, color='black')
    plt.gcf().text(.14, .68, f"Approx Event E = {E_event / jkt:.3f} kt", fontsize=10, color='black')
    plt.gcf().text(.14, .65, f"Approx Event E = {E_event:.3e} J", fontsize=10, color='black')
    plt.axhline(peak_rad_i, linestyle='--', color='red', label=f"Peak Rad I = {peak_rad_i:.3e} W/sr")
    plt.title(f'Bolide Light Curve: Radiant Intensity vs Time\n')
    plt.grid(True)
    plt.legend(loc='upper left')
    plt.gcf().text(.10, .95, 'UNCLASSIFIED', fontsize=10, color='green')
    plt.gcf().text(.80, .05, 'UNCLASSIFIED', fontsize=10, color='green')
    pdf_pages.savefig(fig4)

    # luminosity vs altitude
    fig5 = plt.figure(figsize=(9.0, 6.5))
    ymin = 0.0
    ymax = 100.
    xmin = 0.8 * lums[0]
    xmax = 2. * peak_power
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.plot(lums, alts / 1000., color='blue', label='Luminosity vs Altitude')
    plt.xscale('linear')
    plt.xlabel('Luminosity (Watts)')
    plt.ylabel('Altitude (km)')
    plt.gcf().text(.20, .89, header, fontsize=10, color='black')
    if frag_alts:
        plt.axhline(frag_alts[0] / 1000., linestyle='--', color='green', alpha=0.6, label=f"Frag at {frag_alts[0]/1000.:.1f} km and {frag_times[0]:.2f} sec")
    plt.gcf().text(.65, .74, f"Total Rad E = {E_rad_tot / jkt:.3f} kt", fontsize=10, color='black')
    plt.gcf().text(.65, .71, f"Total Rad E = {E_rad_tot:.3e} J", fontsize=10, color='black')
    plt.gcf().text(.65, .68, f"Approx Event E = {E_event / jkt:.3f} kt", fontsize=10, color='black')
    plt.gcf().text(.65, .65, f"Approx Event E = {E_event:.3e} J", fontsize=10, color='black')
    plt.axvline(peak_power, linestyle='--', color='red', label=f'Peak Power = {peak_power:.3e} W')
    plt.title(f'Bolide Light Curve: Altitude vs Luminosity\n')
    plt.grid(True)
    plt.legend(loc='upper right')
    plt.gcf().text(.10, .95, 'UNCLASSIFIED', fontsize=10, color='green')
    plt.gcf().text(.80, .05, 'UNCLASSIFIED', fontsize=10, color='green')
    pdf_pages.savefig(fig5)

    # radiant intensity in W/sr vs altitude
    fig6 = plt.figure(figsize=(9.0, 6.5))
    # convert luminosity in watts to W/sr
    radiant_intensities = lums / (4. * pi)
    ymin = 0.0
    ymax = 100.
    xmin = 0.8 * radiant_intensities[0]
    xmax = 2. * peak_rad_i
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.plot(radiant_intensities, alts / 1000., color='blue', label='Rad I vs Altitude')
    plt.xscale('linear')
    plt.xlabel('Radiant Intensity (W/sr)')
    plt.ylabel('Altitude (km)')
    plt.gcf().text(.20, .89, header, fontsize=10, color='black')
    if frag_alts:
        plt.axhline(frag_alts[0] / 1000., linestyle='--', color='green', label=f"Frag at {frag_alts[0]/1000.:.1f} km and {frag_times[0]:.2f} sec")
    plt.gcf().text(.65, .74, f"Total Rad E = {E_rad_tot / jkt:.3f} kt", fontsize=10, color='black')
    plt.gcf().text(.65, .71, f"Total Rad E = {E_rad_tot:.3e} J", fontsize=10, color='black')
    plt.gcf().text(.65, .68, f"Approx Event E = {E_event / jkt:.3f} kt", fontsize=10, color='black')
    plt.gcf().text(.65, .65, f"Approx Event E = {E_event:.3e} J", fontsize=10, color='black')
    plt.axvline(peak_rad_i, linestyle='--', color='red', label=f'Peak Rad I = {peak_rad_i:.3e} W/sr')
    plt.title(f'Bolide Light Curve: Altitude vs Radiant Intensity\n')
    plt.grid(True)
    plt.legend(loc='upper right')
    plt.gcf().text(.10, .95, 'UNCLASSIFIED', fontsize=10, color='green')
    plt.gcf().text(.80, .05, 'UNCLASSIFIED', fontsize=10, color='green')
    pdf_pages.savefig(fig6)

    # altitude vs time
    fig7 = plt.figure(figsize=(9.0, 6.5))
    plt.plot(times, alts / 1000.)
    plt.xscale('log')
    plt.xlabel('Time (seconds))')
    plt.ylabel('Altitude (km)')
    plt.gcf().text(.20, .89, header, fontsize=10, color='black')
    if frag_times:
        plt.axvline(frag_times[0], linestyle='--', color='green', alpha=0.6)
        plt.gcf().text(.14, .86, f"Frag at {frag_alts[0]/1000.:.1f} km", fontsize=9, color='green')
    plt.title(f'Bolide Altitude vs Time\n')
    plt.grid(True)
    # plt.legend(loc='upper right')
    plt.gcf().text(.10, .95, 'UNCLASSIFIED', fontsize=10, color='green')
    plt.gcf().text(.80, .05, 'UNCLASSIFIED', fontsize=10, color='green')
    pdf_pages.savefig(fig7)

    # velocity vs time
    fig8 = plt.figure(figsize=(9.0, 6.5))
    plt.plot(times, velocities/1000.)
    plt.xscale('linear')
    plt.yscale('linear')
    plt.xlabel('Time (seconds)')
    plt.ylabel('Velocity (km/s)')
    plt.gcf().text(.20, .89, header, fontsize=10, color='black')
    # plt.legend()
    plt.title(f'Bolide Velocity vs Time\n')
    if frag_times:
        plt.axvline(frag_times[0], linestyle='--', color='green', alpha=0.6)
        plt.gcf().text(.14, .86, f"Frag at {frag_alts[0]/1000.:.1f} km", fontsize=9, color='green')
    plt.grid(True)
    plt.gcf().text(.10, .95, 'UNCLASSIFIED', fontsize=10, color='green')
    plt.gcf().text(.80, .05, 'UNCLASSIFIED', fontsize=10, color='green')
    pdf_pages.savefig(fig8)

    # acceleration vs time
    fig9 = plt.figure(figsize=(9.0, 6.5))
    plt.plot(times, accels)
    plt.xscale('linear')
    plt.yscale('linear')
    plt.xlabel('Time (seconds)')
    plt.ylabel('Acceleration (m/s^2)')
    plt.gcf().text(.20, .89, header, fontsize=10, color='black')
    if frag_times:
        plt.axvline(frag_times[0], linestyle='--', color='green', alpha=0.6)
        plt.gcf().text(.14, .86, f"Frag at {frag_alts[0]/1000.:.1f} km", fontsize=9, color='green')
    plt.title(f'Bolide Acceleration vs Time\n')
    plt.grid(True)
    # plt.legend()
    plt.gcf().text(.10, .95, 'UNCLASSIFIED', fontsize=10, color='green')
    plt.gcf().text(.80, .05, 'UNCLASSIFIED', fontsize=10, color='green')
    pdf_pages.savefig(fig9)

    # dynamic pressure vs time
    fig10 = plt.figure(figsize=(9.0, 6.5))
    plt.plot(times, qs)
    plt.xscale('linear')
    plt.yscale('linear')
    plt.xlabel('Time (seconds)')
    plt.ylabel('Dynamic Pressure (Pa)')
    plt.gcf().text(.20, .89, header, fontsize=10, color='black')
    if frag_times:
        plt.axvline(frag_times[0], linestyle='--', color='green', alpha=0.6)
        plt.gcf().text(.14, .86, f"Frag at {frag_alts[0]/1000.:.1f} km", fontsize=9, color='green')
        plt.axhline(init_strength, linestyle='--', color='orange', alpha=0.6)
        plt.gcf().text(.14, .84, f"Init Strength =  {init_strength:.1e} Pa", fontsize=9, color='orange')
    plt.title(f'Bolide Dynamic Pressure vs Time\n')
    # plt.legend()
    plt.grid(True)
    plt.gcf().text(.10, .95, 'UNCLASSIFIED', fontsize=10, color='green')
    plt.gcf().text(.80, .05, 'UNCLASSIFIED', fontsize=10, color='green')
    pdf_pages.savefig(fig10)

    # dynamic pressure vs altitude
    fig11 = plt.figure(figsize=(9.0, 6.5))
    plt.plot(alts/1000., qs)
    plt.xscale('linear')
    plt.yscale('linear')
    plt.xlabel('Altitude (km)')
    plt.ylabel('Dynamic Pressure (Pa)')
    plt.gcf().text(.20, .89, header, fontsize=10, color='black')
    if frag_alts:
        plt.axhline(frag_alts[0], linestyle='--', color='green', alpha=0.6)
        plt.gcf().text(.14, .86, f"Frag at {frag_alts[0]/1000.:.1f} km", fontsize=9, color='green')
    plt.title(f'Bolide Dynamic Pressure vs Altitude\n')
    plt.grid(True)
    plt.gcf().text(.10, .95, 'UNCLASSIFIED', fontsize=10, color='green')
    plt.gcf().text(.80, .05, 'UNCLASSIFIED', fontsize=10, color='green')
    pdf_pages.savefig(fig11)

    # stagnation temperature vs time
    fig12 = plt.figure(figsize=(9.0, 6.5))
    plt.plot(times, T_stags)
    plt.xscale('linear')
    plt.yscale('linear')
    plt.xlabel('Time (seconds)')
    plt.ylabel('Stagnation Temperature (K)')
    plt.gcf().text(.20, .89, header, fontsize=10, color='black')
    if frag_times:
        plt.axvline(frag_times[0], linestyle='--', color='green', alpha=0.6)
        plt.gcf().text(.14, .86, f"Frag at {frag_alts[0]/1000.:.1f} km", fontsize=9, color='green')
    plt.title(f'Bolide Stagnation Temperature vs Time\n')
    plt.grid(True)
    plt.gcf().text(.10, .95, 'UNCLASSIFIED', fontsize=10, color='green')
    plt.gcf().text(.80, .05, 'UNCLASSIFIED', fontsize=10, color='green')
    pdf_pages.savefig(fig12)

    # stagnation pressure vs time
    fig13 = plt.figure(figsize=(9.0, 6.5))
    plt.plot(times, p_stags)
    plt.xscale('linear')
    plt.yscale('linear')
    plt.xlabel('Time (seconds)')
    plt.ylabel('Stagnation Pressure (Pa)')
    plt.gcf().text(.20, .89, header, fontsize=10, color='black')
    if frag_times:
        plt.axvline(frag_times[0], linestyle='--', color='green', alpha=0.6)
        plt.gcf().text(.14, .86, f"Frag at {frag_alts[0]/1000.:.1f} km", fontsize=9, color='green')
    plt.title(f'Bolide Stagnation Pressure vs Time\n')
    plt.grid(True)
    plt.gcf().text(.10, .95, 'UNCLASSIFIED', fontsize=10, color='green')
    plt.gcf().text(.80, .05, 'UNCLASSIFIED', fontsize=10, color='green')
    pdf_pages.savefig(fig13)

    # surface temperature vs time
    fig14 = plt.figure(figsize=(9.0, 6.5))
    plt.plot(times, T_surfs)
    plt.xscale('linear')
    plt.yscale('linear')
    plt.xlabel('Time (seconds)')
    plt.ylabel('Surface Temperature (K)')
    plt.gcf().text(.20, .89, header, fontsize=10, color='black')
    if frag_times:
        plt.axvline(frag_times[0], linestyle='--', color='green', alpha=0.6, label=f"Frag at {frag_alts[0]/1000.:.1f} km and {frag_times[0]:.2f} sec")
    plt.title(f'Bolide Surface Temperature vs Time\n')
    plt.grid(True)
    plt.gcf().text(.10, .95, 'UNCLASSIFIED', fontsize=10, color='green')
    plt.gcf().text(.80, .05, 'UNCLASSIFIED', fontsize=10, color='green')
    pdf_pages.savefig(fig14)

    # luminous efficiency (tau) vs time
    fig15 = plt.figure(figsize=(9.0, 6.5))
    plt.plot(times, taus, label=f'τ vs time')
    plt.xscale('linear')
    plt.yscale('linear')
    plt.xlabel('Time (seconds)')
    plt.ylabel('Luminous Efficiency (τ)')
    plt.gcf().text(.20, .89, header, fontsize=10, color='black')
    if frag_times:
        plt.axvline(frag_times[0], linestyle='--', color='green', alpha=0.6, label=f"Frag at {frag_alts[0]/1000.:.1f} km and {frag_times[0]:.2f} sec")
    plt.title(f'Bolide Luminous Efficiency vs Time\n')
    plt.grid(True)
    plt.legend(loc='upper left')
    plt.gcf().text(.10, .95, 'UNCLASSIFIED', fontsize=10, color='green')
    plt.gcf().text(.80, .05, 'UNCLASSIFIED', fontsize=10, color='green')
    pdf_pages.savefig(fig15)

    # luminous efficiency (tau) vs altitude
    fig16 = plt.figure(figsize=(9.0, 6.5))
    plt.plot(taus, alts / 1000., label=f'τ vs Altitude')
    plt.xscale('linear')
    plt.yscale('linear')
    plt.xlabel('Luminous Efficiency (τ)')
    plt.ylabel('Altitude (km)')
    plt.gcf().text(.20, .89, header, fontsize=10, color='black')
    if frag_alts:
        plt.axhline(frag_alts[0] / 1000., linestyle='--', color='green', alpha=0.6, label=f"Frag at {frag_alts[0]/1000.:.1f} km and {frag_times[0]:.2f} sec")
    plt.title(f'Bolide Luminous Efficiency vs Altitude\n')
    plt.grid(True)
    plt.legend(loc='upper right')
    plt.gcf().text(.10, .95, 'UNCLASSIFIED', fontsize=10, color='green')
    plt.gcf().text(.80, .05, 'UNCLASSIFIED', fontsize=10, color='green')
    pdf_pages.savefig(fig16)

    # total radiated energy vs time
    fig17 = plt.figure(figsize=(9.0, 6.5))
    plt.plot(times, E_rad_totals)
    plt.xscale('linear')
    plt.yscale('linear')
    plt.xlabel('Time (seconds)')
    plt.ylabel('Total Radiated Energy (J)')
    plt.gcf().text(.20, .89, header, fontsize=10, color='black')
    if frag_times:
        plt.axvline(frag_times[0], linestyle='--', color='green', alpha=0.6)
        plt.gcf().text(.14, .86, f"Frag at {frag_alts[0]/1000.:.1f} km", fontsize=9, color='green')
    plt.title(f'Bolide Total Radiated Energy vs Time\n')
    plt.grid(True)
    plt.gcf().text(.10, .95, 'UNCLASSIFIED', fontsize=10, color='green')
    plt.gcf().text(.80, .05, 'UNCLASSIFIED', fontsize=10, color='green')
    pdf_pages.savefig(fig17)

    pdf_pages.close()

    print("All post-processing complete")