# Simplified post-fragmentation evolution using arrays

# overall workflow and integration picture into bolum:
# Start simulation
#    ↓
# [Single bolide evolution]
#    ↓
# Check: q > strength?
#    ↓ No
# Keep evolving as one bolide
#    ↓
# Yes →
# Mark fragmentation
# Record state: t, h, v, m, area
#    ↓
# Break main loop
#    ↓
# Call evolve_fragments(...)
#    ↓
# [Simulate each fragment independently]
#    ↓
# Return full time histories of all fragments


def evolve_fragments(n_fragments, init_mass, init_v, init_h, area, density, angle,
                     t, dt, tstop, eta_0, v_ref, rho_air0, n, m, L_ablation):

    import numpy as np

    fragment_data = []

    # Initialize fragment arrays
    for i in range(n_fragments):
        mass = init_mass / n_fragments
        diameter = 2 * ((3 * mass) / (4 * np.pi * density))**(1/3)
        frag_area = np.pi * (diameter / 2)**2

        state = [mass, init_v, init_h, frag_area, diameter]
        arrays = {
            'time': [t], 'mass': [mass], 'v': [init_v], 'h': [init_h],
            'lum': [], 'area': [frag_area], 'KE': [], 'acc': [],
            'q': [], 'Mach': [], 'T_stag': [], 'p_stag': [], 'T_surf': []
        }
        fragment_data.append([state, arrays])

    # Main evolution loop
    while any(frag[0][1] > 0 and frag[0][2] > 0 and frag[0][0] > 0 for frag in fragment_data):
        for frag in fragment_data:
            state, arr = frag
            mass, v, h, area, diameter = state

            if v <= 0 or h <= 0 or mass <= 0:
                continue

            rho = atmos_density_calc(h)
            T = temp_altitude_calc(h)
            a = sound_speed_calc(T)
            M = mach_num_calc(v, a)
            q = dynamic_pressure(rho, v)
            drag = drag_calc(rho, v, area)
            acc = accel_calc(drag, mass)
            KE = KE_calc(mass, v)
            p = atmos_pressure_calc(h)
            T_stag = T_stag_calc(M, T)
            p_stag = p_stag_calc(q, p)
            q_h = heatflux_calc(drag, v, area)
            T_surf = T_surf_calc(q_h)
            dm_dt = dm_dt_calc(q_h, area, L_ablation)

            eta = luminous_efficiency_scaled(v, rho, eta_0, v_ref, rho_air0, n, m)
            luminosity = lum_calc(eta, drag, v, dt)

            # Update state
            v -= acc * dt
            h -= v * np.sin(angle) * dt
            mass -= dm_dt * dt
            diameter = 2 * ((3 * mass) / (4 * np.pi * density))**(1/3)
            area = np.pi * (diameter / 2)**2

            t += dt

            state[:] = [mass, v, h, area, diameter]
            arr['time'].append(t)
            arr['mass'].append(mass)
            arr['v'].append(v)
            arr['h'].append(h)
            arr['lum'].append(luminosity)
            arr['area'].append(area)
            arr['KE'].append(KE)
            arr['acc'].append(acc)
            arr['q'].append(q)
            arr['Mach'].append(M)
            arr['T_stag'].append(T_stag)
            arr['p_stag'].append(p_stag)
            arr['T_surf'].append(T_surf)

        if t > tstop:
            break

    return fragment_data
