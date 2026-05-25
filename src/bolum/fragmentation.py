'''
Bolum fragmentation model file
'''

from bolum.constants import *

# Main fragmentation model
def frag_attempt(bolide_outputs, alt, t, frag_stage, frag_times, frag_alts,
                 current_strength, n_fragments, density, mass, area, q,
                 flare_end_time, flare_duration):

    # Try to fragment if dynamic pressure > current strength
    if q > current_strength:
        frag_stage = frag_stage + 1
        frag_times.append(t)
        frag_alts.append(alt)

        n_fragments = n_fragments**frag_stage

        mass_fragment  = mass / n_fragments
        r_fragment     = ((3. * mass_fragment) / (4. * pi * density))**(1. / 3.)
        dia_fragment   = 2. * r_fragment
        area_fragment  = 4. * pi * (dia_fragment / 2.)**2.
        area_new       = n_fragments * area_fragment
        mass_new       = mass_fragment * n_fragments # conserved
        flare_end_time = t + flare_duration

        print(f"Fragmentation stage {frag_stage} at {alt/1000.:.1f} km, t={t:.3f} s (q={q:.2e} Pa)")
        print(f'flare end time = {flare_end_time:.4e} sec')
        # print(type(flare_end_time))
        bolide_outputs.write(f"Fragmentation stage {frag_stage} at {alt:.4e} m / {t:.4e} s\n\n")

        return True, area_new, mass_new, frag_stage, frag_times, frag_alts, n_fragments, flare_end_time
    else:   
        return False, area, mass, frag_stage, frag_times, frag_alts, n_fragments, flare_end_time