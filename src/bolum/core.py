'''
bolum core file
'''

import numpy as np
from bolum.constants import *
from bolum.io import *
from bolum.physics import *
from bolum.trajectory import *
from bolum.fragmentation import *
from bolum.luminosity import *
from bolum.plotting import *

############ Main simulation function ############
def bolide_luminosity_model(outputpath, diameter, velocity, theta_deg,
                            init_strength, strength_secondary, strength_tertiary, strength_quaternary,
                            density, porosity, Cd, L_ablation,
                            n_fragments, n_frag_init, flare_duration, rho_tau_scale,
                            zstart, xstart, tstart, dt, tstop):

    # initialize timestep
    timestep = int(1)
    print(f"first timestep = {timestep}")

    # compute initial area and mass
    area = area_init_calc(diameter)
    mass = mass_init_calc(diameter, density)

    # compute initial material strength of bolide
    strength_initial = compute_strength(init_strength, porosity, density)
    current_strength = strength_initial

    # initialize fragmentation model
    fragmented     = False  # still useful as flag for first frag
    flare_end_time = None
    frag_stage     = 0
    frag_times     = []
    frag_alts      = []

    # define header
    header = (f"Dia: {diameter:.1f} m | Vel: {velocity/1000.:.1f} km/s | theta: {theta_deg:.1f} deg | "
    f"Strength: {strength_initial:.1e} Pa")

    # open output files
    # general output file
    bolide_outputs = open(f'{outputpath}/bolide_out.dat', 'w')
    bolide_outputs.write(f"Bolide Luminosity and Fragmentation Simulation — 8\n")
    bolide_outputs.write(f"{header}\n\n")
    bolide_outputs.write(f"Initial area     = {area:.4e} m²\n")
    bolide_outputs.write(f"Initial mass     = {mass:.4e} kg\n")
    bolide_outputs.write(f"Initial strength = {strength_initial:.4e} Pa\n\n")

    # power vs time output file
    bolide_pvt = open(f'{outputpath}/power_vs_time.dat', 'w')
    bolide_pvt.write(f'{header}\n')
    bolide_pvt.write(f'Power (W)    Time (sec)\n')

    # power vs altitude output file
    bolide_pva = open(f'{outputpath}/power_vs_alt.dat', 'w')
    bolide_pva.write(f'{header}\n')
    bolide_pva.write(f'Power (W)    Altitude (m)\n')

    # trajectory init (your original function — assuming it exists)
    (x, z, R, alt, v, t, theta, vx, vz, alpha, delta, gamma, phi, s,
    surfx, surfz, R_E) = trajectory_init(bolide_outputs, xstart, zstart, velocity, tstart, theta_deg)

    # initialize total radiated and deposited energy to zero
    E_rad_total = 0.0
    E_deposited = 0.0
    print(f"setting total radiated energy to {E_rad_total:.4e} J\n")
    bolide_outputs.write(f"Setting total radiated energy to {E_rad_total:.4e} J\n\n")

    # initialize arrays
    T_array           = np.zeros(1)
    p_array           = np.zeros(1)
    rho_array         = np.zeros(1)
    a_array           = np.zeros(1)
    tau_array         = np.zeros(1)
    lum_array         = np.zeros(1)
    alt_array         = np.zeros(1)
    time_array        = np.zeros(1)
    mass_array        = np.zeros(1)
    dia_array         = np.zeros(1)
    area_array        = np.zeros(1)
    KE_array          = np.zeros(1)
    v_array           = np.zeros(1)
    accel_array       = np.zeros(1)
    heatflux_array    = np.zeros(1)
    Fdrag_array       = np.zeros(1)
    E_rad_ttl_array   = np.zeros(1)
    E_deposited_array = np.zeros(1)
    q_array           = np.zeros(1)
    Mach_array        = np.zeros(1)
    T_stag_array      = np.zeros(1)
    p_stag_array      = np.zeros(1)
    T_surf_array      = np.zeros(1)

    # main loop
    while timestep < mxcycl:

        # calculate quantities
        T       = temp_altitude_calc(alt)
        p       = atmos_pressure_calc(alt, T)
        rho     = atmos_density_calc(alt, p, T)
        a       = sound_speed_calc(T, p, rho)
        M       = mach_num_calc(v, a)
        q       = dynamic_pressure_calc(rho, v)
        Fdrag   = Fdrag_calc(rho, v, area, Cd)
        KE      = KE_calc(mass, v)
        T_stag  = T_stag_calc(M, T)
        p_stag  = p_stag_calc(q, p)
        q_h     = heatflux_calc(Fdrag, v, area)
        T_surf  = T_surf_calc(epsilon, q_h)

        # Multi-stage fragmentation check
        (did_frag, area, mass, frag_stage, frag_times, frag_alts,
        n_fragments, flare_end_time) = frag_attempt(bolide_outputs, alt, t, frag_stage, frag_times, frag_alts,
                                    current_strength, n_fragments, density, mass, area, q, flare_end_time, flare_duration)

        if did_frag:
            fragmented = True
            # After first fragmentation → switch to weaker strength
            if frag_stage == 1:
                current_strength = strength_secondary
                print(f"Material strength updated to secondary value: {current_strength:.2e} Pa")
            if frag_stage == 2:
                current_strength = strength_tertiary
                print(f"Material strength updated to tertiary value: {current_strength:.2e} Pa")
            if frag_stage == 3:
                current_strength = strength_quaternary
                print(f"Material strength updated to quaternary value: {current_strength:.2e} Pa")
            if frag_stage > 3:
                current_strength = 3.e7
                print(f'Material strength updated to terminal value of {current_strength} Pa')

        # Luminosity (single term)
        tau        = tau_calc(v, mass, rho, rho_tau_scale, n_fragments if fragmented else 1)
        dmdt       = dmdt_calc(q_h, area, L_ablation)
        luminosity = lum_calc(tau, Fdrag, v, dmdt, n_frag_init, fragmented, t, flare_end_time)

        # Per-fragment reporting
        lum_per_fragment = luminosity / n_fragments if fragmented and n_fragments > 1 else None

        # trajectory update (your original)
        (x, z, R, alt, v, theta, vx, vz, acc, alpha, delta, gamma, phi, s,
        surfx, surfz, R_E) = trajectory_update(x, z, R, v, theta, vx, vz, alpha, delta, gamma, phi, s, Fdrag, dt, mass)

        # mass / size update
        mass     = mass_update_calc(mass, dmdt, dt)
        diameter = diameter_update_calc(mass, density, n_fragments if fragmented else 1, fragmented)
        area     = area_update_calc(diameter, n_fragments if fragmented else 1, fragmented)

        # energy accumulation
        E_rad_total = E_rad_total_calc(E_rad_total, luminosity, dt)
        E_deposited = E_deposited_calc(E_deposited, Fdrag, v, dmdt, dt)

        # store to arrays
        T_array[timestep-1]           = T
        p_array[timestep-1]           = p
        rho_array[timestep-1]         = rho
        a_array[timestep-1]           = a
        tau_array[timestep-1]         = tau
        lum_array[timestep-1]         = luminosity
        alt_array[timestep-1]         = alt
        time_array[timestep-1]        = t
        mass_array[timestep-1]        = mass
        dia_array[timestep-1]         = diameter
        area_array[timestep-1]        = area
        KE_array[timestep-1]          = KE
        v_array[timestep-1]           = v
        accel_array[timestep-1]       = acc
        heatflux_array[timestep-1]    = q_h
        Fdrag_array[timestep-1]       = Fdrag
        E_rad_ttl_array[timestep-1]   = E_rad_total
        E_deposited_array[timestep-1] = E_deposited
        q_array[timestep-1]           = q
        Mach_array[timestep-1]        = M
        T_stag_array[timestep-1]      = T_stag
        p_stag_array[timestep-1]      = p_stag
        T_surf_array[timestep-1]      = T_surf

        # advance time
        t = t + dt

        # write outputs
        wr_out(bolide_outputs, luminosity, alt, t, mass, diameter, area, KE, v, acc, q_h, Fdrag,
               E_rad_total, E_deposited, q, M, T_stag, p_stag, T_surf, T, p, rho, a, theta_deg, tau,
               frag_stage, lum_per_fragment)
        
        wr_power_vs_time(bolide_pvt, luminosity, t)

        wr_power_vs_alt(bolide_pva, luminosity, alt)

        # end conditions
        if t > tstop:
            print('STOP simulation: time exceeded end time\n\n')
            break
        if mass <= 0.:
            print('STOP simulation: mass <= 0 kg\n\n')
            break
        if alt <= 0.:
            print('STOP simulation: altitude <= 0 km\n\n')
            break
        if v <= 0.:
            print('STOP simulation: v <= 0 m/s\n\n')
            break

        # extend arrays
        T_array           = np.vstack([T_array, new_row])
        p_array           = np.vstack([p_array, new_row])
        rho_array         = np.vstack([rho_array, new_row])
        a_array           = np.vstack([a_array, new_row])
        tau_array         = np.vstack([tau_array, new_row])
        lum_array         = np.vstack([lum_array, new_row])
        alt_array         = np.vstack([alt_array, new_row])
        time_array        = np.vstack([time_array, new_row]) 
        mass_array        = np.vstack([mass_array, new_row])
        dia_array         = np.vstack([dia_array, new_row])
        area_array        = np.vstack([area_array, new_row])
        KE_array          = np.vstack([KE_array, new_row])
        v_array           = np.vstack([v_array, new_row])
        accel_array       = np.vstack([accel_array, new_row])
        heatflux_array    = np.vstack([heatflux_array, new_row])
        Fdrag_array       = np.vstack([Fdrag_array, new_row])
        E_rad_ttl_array   = np.vstack([E_rad_ttl_array, new_row])
        E_deposited_array = np.vstack([E_deposited_array, new_row])
        q_array           = np.vstack([q_array, new_row])
        Mach_array        = np.vstack([Mach_array, new_row])
        T_stag_array      = np.vstack([T_stag_array, new_row])
        p_stag_array      = np.vstack([p_stag_array, new_row])
        T_surf_array      = np.vstack([T_surf_array, new_row])

        # update timestep
        timestep = timestep + 1
        print(f'timestep updated: new time will be {t:.3e} sec\n\n')

    # final stats
    peak_power, peak_rad_i = peak_find(bolide_outputs, lum_array)
    E_event                = E_event_calc(bolide_outputs, E_rad_total, E_deposited) # your original

    # print out if bolide fragmented or not
    if fragmented == True:
        print(f'bolide fragmented at {frag_alts[0]/1000.:.1f} km and {frag_times[0]:.2f} sec')
    else:
        print(f'bolide did not fragment')

    bolide_outputs.close()
    bolide_pvt.close()
    bolide_pva.close()

    return (header, strength_initial, T_array, p_array, rho_array, a_array, tau_array,
            lum_array, alt_array, time_array, mass_array, dia_array, area_array,
            frag_times, frag_alts, KE_array, v_array, accel_array, heatflux_array,
            Fdrag_array, E_rad_ttl_array, q_array, Mach_array, T_stag_array,
            p_stag_array, T_surf_array, peak_power, peak_rad_i, E_rad_total, E_event)