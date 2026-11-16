# bolum v2.3 -- where each formulation comes from

Provenance: **lit** = taken from the literature, **adapt** = literature form adapted, **cal** = value calibrated here, **new** = construct introduced here (no direct published equivalent).
Citation details marked (verify) are from memory and were not confirmed against the paper.

| # | Formulation in the code | Provenance | Sources |
|---|---|---|---|
| 1 | Density vs altitude (US Std 1976, analytic to 86 km, table to 120 km) | lit | U.S. Standard Atmosphere 1976 (NOAA-S/T 76-1562; a standard, not a journal paper) |
| 2 | Drag force F = 1/2 Cd rho A v^2, Cd = 1 | lit | Ceplecha et al. 1998, Space Sci. Rev. 84, 327; Collins et al. 2005, MAPS 40, 817 |
| 3 | Ablation dm/dt = Lambda (1/2) rho v^3 A / Q, Lambda <= 1 on projected area | lit (form), cal (value) | Ceplecha et al. 1998; Collins et al. 2005 pancake equations; Johnston, Stern & Wheeler 2018, Icarus 309, 25 (canonical Lambda = 0.1; coupled radiation/ablation lowers heating). Default Lambda = 0.05 -> sigma = 0.005 s^2/km^2 is calibrated |
| 4 | Ablation onset when net flux exceeds radiative loss at T_ablate | adapt | Love & Brownlee 1991, Icarus 89, 26 |
| 5 | Breakup when stagnation pressure rho v^2 >= strength | lit | Chyba, Thomas & Zahnle 1993, Nature 361, 40; Hills & Goda 1993, AJ 105, 1114; Collins et al. 2005 (P_s = rho v^2 = Y); Wheeler, Register & Mathias 2017, Icarus 295, 149; Borovicka, Spurny & Shrbeny 2020, AJ 160, 42 |
| 6 | Child strength S = S_parent (m_parent/m_child)^alpha | lit (form), cal (alpha = 0.4) | Weibull 1951, J. Appl. Mech. 18, 293; Wheeler et al. 2017 (S2 = S1 (m1/m2)^alpha) |
| 7 | Pancake spreading dr/dt = c v sqrt(rho_air/rho_body); pancake limit f_max | adapt (form), cal (c = 0.5, f_max = 8) | Chyba et al. 1993; Hills & Goda 1993; Collins et al. 2005 (pancake factor ~7); Korycansky & Zahnle 2005, Planet. Space Sci. 53, 695 (verify pages; their dispersion coefficient is 0.41 under the square root, i.e. c ~ 0.64) |
| 8 | Progressive cascade: split into children (power-law mass shares), area kept continuous | new (concept: cloud/fragment models) | Concept from Ceplecha & ReVelle 2005, MAPS 40, 35; Wheeler et al. 2017; Register, Mathias & Wheeler 2017, Icarus 284, 157 (verify volume/pages); Borovicka et al. 2013 (Kosice), MAPS 48, 1757. The specific rules (4 children, k^-1 shares, f_split = 2) are not published values |
| 9 | Dust release at each breakup (20 %, 2 cm particles) | new, cal | Motivated by Borovicka et al. 2020 (first-phase fragmentation often loses >= 40 % of the mass), Borovicka et al. 2013 Nature 503, 235 (Chelyabinsk dust), Ceplecha & ReVelle 2005. Values calibrated to Chelyabinsk |
| 10 | Luminosity equation I = tau (v^2/2 dm/dt + m v dv/dt), summed over fragments | lit | Ceplecha et al. 1998; Ceplecha & ReVelle 2005; Borovicka et al. 2020 Eq. 2 |
| 11 | tau(v, m) polynomial/tanh fit | lit | Borovicka, Spurny & Shrbeny 2020, AJ 160, 42, Eqs. 3-4 (velocity form: Pecina & Ceplecha 1983; mass form adapted from ReVelle & Ceplecha 2001, as cited there). The published paper says tau may apply to the whole spectrum or a limited band - confirm which band applies to your use |
| 12 | tau clamp (v >= 2 km/s) and taper to zero by 4 km/s | new | -- |
| 13 | Thin/thick efficiency blend, eta = (1-w) tau + w eta_thick f_si, w = r/(r + r_half) | new, cal | Constrained by Brown et al. 2002, Nature 420, 294 (integral efficiency rises with energy) and the Chelyabinsk/TC3 radiated energies; physical motivation: radiation-dominated heating of large bodies (Johnston et al. 2018). eta_thick = 0.20, r_half = 6 m are calibrated |
| 14 | Si-band fraction f_si(T) = integral of Planck B over 0.4-1.1 um / (sigma T^4/pi), T = 6000 K | lit | Planck's law; 6000 K blackbody convention of Brown et al. 2002 (verify) |
| 15 | Brown population check E_total = 8.2508 E_rad^0.885 (kt) | lit (used as a check only) | Brown et al. 2002 |
| 16 | Energy ledger KE0 + W_grav = E_dep + KE_left | -- | Work-energy theorem; no citation needed |
| 17 | Validation targets | lit | Chelyabinsk: Brown et al. 2013 Nature 503, 238; Popova et al. 2013 Science 342, 1069; Borovicka et al. 2013 Nature 503, 235. TC3: Jenniskens et al. 2009 Nature 458, 485; Borovicka & Charvat 2009, A&A 507, 1015. Kosice: Borovicka et al. 2013, MAPS 48, 1757 |

## Strength convention
All strengths in bolum v2.3 are stagnation pressures rho v^2. Collins et al. 2005 state this directly (P_s = rho v^2 compared with the strength), the fragment-cloud model uses rho_air v^2 > strength, and Borovicka et al. 2020 use the same pressure. 1/2 rho v^2 is the incompressible dynamic pressure and understates the hypersonic load by about 2x (Rayleigh pitot with gamma = 1.4 gives 0.92 rho v^2). The drag force still uses 1/2 Cd rho A v^2.

## Added in v2.3: light-curve feature analysis
| # | Item | Provenance | Source |
|---|---|---|---|
| 18 | Feature labels (ablation onset, first yield, flares, peak, roll-off, tail) and their mapping to observed light-curve phases | new (labelling only; no physics) | Jenniskens, Pilorz, Robertson & Stern 2026, Meteorit. Planet. Sci. 61(9), 2351-2392, doi:10.1111/maps.70203. Only the abstract, section headings and press summaries were read; the full text was not accessible, so phase-by-phase details (heights, numbers) were not checked |
| 19 | Flare criterion: dlnL/dt > 3 x (v sin(theta) / density scale height), >= 2 per s, rise >= 25 % | new | Heuristic chosen here; the air-density tracking rate is the rate at which L would grow if it simply followed rho(h) |
