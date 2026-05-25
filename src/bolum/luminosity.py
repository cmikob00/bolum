'''
Bolum luminosity calculations file
'''

from bolum.constants import *

# function to calculate bolide lumonous efficiency
def tau_calc(v, mass, rho, rho_tau_scale, n_fragments):

    # This luminous efficiency function is designed to output
    # bolide luminous efficiency based on optical (broadband Si)
    # sensor response - not the entire spectrum

    v_kms = v / 1000.0
    if n_fragments > 1:
        m_eff = mass / n_fragments
    else:
        m_eff = mass

    ln_v = np.log(v_kms) if v_kms > 1e-5 else -10.0
    ln_m = np.log(m_eff) if m_eff > 1e-10 else -20.0

    if v_kms < 25.372:
        ln_tau = (0.567 - 10.307 * ln_v + 9.781 * ln_v**2. - 3.0414 * ln_v**3. + 0.3213 * ln_v**4. + 0.347 * np.tanh(0.38 * ln_m))
    else:
        ln_tau = -1.4286 + ln_v + 0.347 * np.tanh(0.38 * ln_m)

    # exponentiate
    tau = np.exp(ln_tau)

    # apply air density scaling for large bolides (small bolides don't need it)
    # this is not bulletproof
    # updated to 1.e-3 and 0.4 for v9 on 7 April 2026 cjm
    if rho_tau_scale == 1:
        tau = tau * (rho / 1.e-3)**0.45  # empirical term set to match Chelyabinsk luminosity

    # tau is now in percent, so divide by 100 to get true value
    tau = tau / 100.

    # Extra Stuff:
    # Bolometric τ (Borovička base) + airburst scaling
    # Small bolides: pure Borovička (realistic 2-7%)
    # Large (>~few m, airburst regime): empirical boost to match Chelyabinsk (12-17%)
    # if mass > 5.e5:  # Threshold: ~5m stony (~1e6 kg); below = fireballs
    #     tau_scale = 1.0 + 2.8 * (mass / 1.e7)**0.35  # ~1x small, 3.5x Chelyabinsk
    #     tau = tau * tau_scale
    # else: tau remains ~6% for v~19 km/s, m~1e3 kg (perfect for Benešov/Košice)

    return tau

def lum_calc(tau, Fdrag, v, dmdt, n_frag_init, fragmented, t, flare_end_time):

    # calculate luminosity (recall tau was set to match optical outputs)
    # based on ablation only - not modeling a shock front
    E_dot_total = (Fdrag * v) + (0.5 * v**2. * abs(dmdt))
    luminosity = tau * E_dot_total

    # apply "pancake" model to flare if fragmented
    if fragmented and t <= flare_end_time:
            luminosity = luminosity * n_frag_init**0.333

    return luminosity