#!/usr/bin/env python3
from pathlib import Path
import sys

project_root = Path(__file__).parent
src_path = project_root / "src"
sys.path.insert(0, str(src_path))

from bolum.core import bolide_luminosity_model
from bolum.io import rdinput, write_init_params
from bolum.plotting import plot_outputs

if __name__ == "__main__":

    #### File Path ####
    basename = f"/Users/julie/Desktop/Projects/bolum/Events/"

    # Main Program
    event_name = input('Enter event name: ')
    print(f'Event name input: {event_name}')
    inputpath  = f"{basename}{event_name}/bolum_in.dat"
    outputpath = f"{basename}{event_name}"

    # read inputs
    (diameter, velocity, theta_deg,
    init_strength, strength_secondary, strength_tertiary, strength_quaternary,
    density, porosity, Cd, L_ablation,
    n_fragments, n_frag_init, flare_duration, rho_tau_scale,
    zstart, xstart, tstart, dt, tstop) = rdinput(inputpath)

    # write out initial parameters
    write_init_params(outputpath, diameter, velocity, theta_deg,
    init_strength, strength_secondary, strength_tertiary, strength_quaternary,
    density, porosity, Cd, L_ablation,
    n_fragments, n_frag_init, flare_duration, rho_tau_scale,
    zstart, xstart, tstart, dt, tstop)

    # run simulation
    (header, strength, T_array, p_array, rho_array, a_array,
    tau_array, lum_array, alt_array, time_array, mass_array,
    dia_array, area_array, frag_times, frag_alts, KE_array,
    v_array, accel_array, heatflux_array, Fdrag_array, E_rad_ttl_array,
    q_array, Mach_array, T_stag_array, p_stag_array, T_surf_array,
    peak_power, peak_rad_i, E_rad_tot, E_event) = bolide_luminosity_model(outputpath, diameter, velocity, theta_deg,
                                                                        init_strength, strength_secondary, strength_tertiary, strength_quaternary,
                                                                        density, porosity, Cd, L_ablation,
                                                                        n_fragments, n_frag_init, flare_duration, rho_tau_scale,
                                                                        zstart, xstart, tstart, dt, tstop)
                                                                
    # end simulation and plot output
    # error checking
    print(f"total cycles:     {time_array.size}\n")

    print(f"STOP all done bolide simulation complete\n")
    plot_outputs(outputpath, header, time_array.flatten(), tau_array.flatten(), lum_array.flatten(),
                alt_array.flatten(),frag_times, frag_alts, v_array.flatten(), accel_array.flatten(), 
                q_array.flatten(), T_stag_array.flatten(), p_stag_array.flatten(), T_surf_array.flatten(),
                init_strength,
                E_rad_ttl_array.flatten(), peak_power, peak_rad_i, E_rad_tot, E_event)