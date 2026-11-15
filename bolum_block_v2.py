#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
bolum_block_v2.py  (version 2.1)  --  reduced-order bolide entry / fragmentation / luminosity model

CHANGES IN 2.1 (vs 2.0)
  * Fragmentation loading is now ALWAYS the stagnation (ram) pressure  p = rho * v^2  [R5,R6,R7,R8,R2].
    q_factor is gone.  Strengths therefore follow the published convention directly.
    !!  If you carry over strengths from the ORIGINAL bolum_block.py (which used 0.5*rho*v^2),
    !!  MULTIPLY THEM BY 2 to get the same behaviour.
  * One standardized 6-panel summary figure per event (same layout for every event) in bolide_plots.pdf.
  * Log power-vs-time plots now extend ~10 decades below the peak so the early-time signal is visible.
  * Much richer command-line report (also written to run_report.txt).
  * Optional sensor-saturation overlay (--sat-W) for comparing against saturated data.
  * References / provenance for every formulation are listed below and tagged [R#] in the code.

WHAT THE MODEL DOES (details, equations and provenance in the REFERENCES block at the end of this docstring)
  Step 1  US Std 1976 atmosphere, Earth-centred trajectory, adaptive step, exact kinetic-energy ledger.
  Step 2  Ablation: dm/dt = Lambda * (1/2) rho v^3 * A_proj / Q  (Lambda <= 1, decoupled from Cd).
  Step 3  Fragment groups: yield at rho v^2 >= S, pancake spreading, Weibull-strengthened children, dust.
  Light   L = eta * (drag work + KE of ablated mass);  eta blends the empirical tau(v,m) (thin, meteor-like)
          with an optically-thick bolometric efficiency for large clouds.  Output: Si-band L_si and bolometric
          L_bol = L_si / f_si(T_rad),  f_si = in-band fraction of a blackbody at T_rad (Planck law).

INPUT DECK: unchanged format (first 2 lines skipped, 'key = value  # comment').  Legacy keys rho_tau_scale,
flare_duration, n_frag_init, q_factor are accepted and ignored.  strength_secondary/tertiary/quaternary are
optional overrides of the child strength at generations 1/2/3 (otherwise Weibull scaling).

USAGE
  python3 bolum_block_v2.py EVENT [--base DIR] [--no-plots] [--quiet] [--sat-W WATTS]

PROVENANCE KEY:  [lit] taken from the literature   [adapt] literature form, adapted   [cal] value calibrated here
                 [new] construct introduced here (no direct published equivalent)

REFERENCES
  [R1]  Ceplecha, Borovicka, Elford, ReVelle, Hawkes, Porubcan & Simek (1998), Meteor phenomena and bodies,
        Space Sci. Rev. 84, 327-471.   Drag, ablation, luminosity equation I = -tau (v^2/2 dm/dt + m v dv/dt).
  [R2]  Borovicka, Spurny & Shrbeny (2020), Two strengths of ordinary chondritic meteoroids as derived from their
        atmospheric fragmentation modeling, Astron. J. 160, 42 (doi:10.3847/1538-3881/ab9608).  tau(v,m) fit
        (their Eqs. 3-4; velocity dependence from Pecina & Ceplecha 1983, mass dependence adapted from ReVelle &
        Ceplecha 2001); fragmentation strengths expressed as p = rho v^2.
  [R4]  Ceplecha & ReVelle (2005), Fragmentation model of meteoroid motion, mass loss, and radiation in the
        atmosphere, Meteorit. Planet. Sci. 40, 35-54.   Sum of fragment luminosities; continuous + gross fragmentation.
  [R5]  Chyba, Thomas & Zahnle (1993), The 1908 Tunguska explosion: atmospheric disruption of a stony asteroid,
        Nature 361, 40-44.   Pancake model, stagnation-pressure breakup.
  [R6]  Hills & Goda (1993), The fragmentation of small asteroids in the atmosphere, Astron. J. 105, 1114-1144.
  [R7]  Collins, Melosh & Marcus (2005), Earth Impact Effects Program, Meteorit. Planet. Sci. 40, 817-840.
        Breakup when stagnation pressure rho v^2 equals strength; pancake factor ~7.
  [R8]  Wheeler, Register & Mathias (2017), A fragment-cloud model for asteroid breakup and atmospheric energy
        deposition, Icarus 295, 149-169.   Strength scaling S2 = S1 (m1/m2)^alpha, debris cloud + discrete fragments,
        lateral dispersion under a common bow shock.
  [R9]  Register, Mathias & Wheeler (2017), Asteroid fragmentation approaches for modeling atmospheric energy
        deposition, Icarus 284, 157-166 (volume/pages from memory - verify).
  [R10] Johnston, Stern & Wheeler (2018), Radiative heating of large meteoroids during atmospheric entry,
        Icarus 309, 25-44.   Canonical heat-transfer coefficient 0.1; radiative heating dominates for large bodies
        and is reduced by coupled ablation/radiation.
  [R11] Weibull (1951), A statistical distribution function of wide applicability, J. Appl. Mech. 18, 293-297.
  [R12] Brown, Spalding, ReVelle, Tagliaferri & Worden (2002), The flux of small near-Earth objects colliding with
        the Earth, Nature 420, 294-296.   E_total = 8.2508 E_rad^0.885 (kt); integral luminous efficiency rises with
        energy; 6000 K blackbody convention for the bolometric conversion.
  [R13] Brown et al. (2013), A 500-kiloton airburst over Chelyabinsk and an enhanced hazard from small impactors,
        Nature 503, 238-241.
  [R14] Popova et al. (2013), Chelyabinsk airburst, damage assessment, meteorite recovery, and characterization,
        Science 342, 1069-1073.
  [R15] Borovicka et al. (2013), The trajectory, structure and origin of the Chelyabinsk asteroidal impactor,
        Nature 503, 235-237.
  [R16] Borovicka et al. (2013), The Kosice meteorite fall: atmospheric trajectory, fragmentation, and orbit,
        Meteorit. Planet. Sci. 48, 1757-1779.
  [R17] Jenniskens et al. (2009), The impact and recovery of asteroid 2008 TC3, Nature 458, 485-488.
  [R18] Borovicka & Charvat (2009), Meteosat observation of the atmospheric entry of 2008 TC3 over Sudan and the
        associated dust cloud, Astron. Astrophys. 507, 1015-1022.
  [R19] U.S. Standard Atmosphere, 1976 (NOAA/NASA/USAF, NOAA-S/T 76-1562).  A government standard, not a journal paper.
  [R20] Love & Brownlee (1991), Heating and thermal transformation of micrometeoroids entering the Earth's
        atmosphere, Icarus 89, 26-43.   Radiative loss vs. heating at the surface (ablation-onset idea).
  [R21] Korycansky & Zahnle (2005), Modeling crater populations on Venus and Titan, Planet. Space Sci. 53, 695-710
        (pages from memory - verify).   Transverse dispersion velocity v_t = (0.41 rho/rho_body)^0.5 v.
  [T1]  Standard textbook physics: Planck's law (band fraction); Rayleigh pitot / Newtonian hypersonic stagnation
        pressure (p_stag = 0.92 rho v^2 for gamma = 1.4, -> rho v^2 for strong shocks).
"""
__version__ = "2.1"

import os
import sys
import time
import math
import argparse
from dataclasses import dataclass, fields, field
import numpy as np

# ----------------------------------------------------------------------------------------
# constants
# ----------------------------------------------------------------------------------------
PI = math.pi
RE = 6.371e6                 # mean Earth radius (m)
GM = 3.986004418e14          # Earth gravitational parameter (m^3/s^2)
G0 = 9.80665
JKT = 4.184e12               # J per kiloton TNT
SIGMA_SB = 5.670374419e-8
H_PL, C_LIGHT, K_B = 6.62607015e-34, 2.99792458e8, 1.380649e-23
DEFAULT_BASE = "/Users/julie/Desktop/Projects/bolum/Events/"

_trapz = getattr(np, "trapezoid", None) or np.trapz


# ----------------------------------------------------------------------------------------
# atmosphere  [lit R19]: US Standard 1976 (analytic layers to 86 km, tabulated density 86-120 km)
# ----------------------------------------------------------------------------------------
class Atmosphere:
    """Density vs geometric altitude.  rho(alt_m) is vectorised; extrapolates exponentially above 120 km."""

    _UPPER_KM = np.array([86., 90., 95., 100., 105., 110., 115., 120.])
    _UPPER_RHO = np.array([6.958e-6, 3.416e-6, 1.393e-6, 5.604e-7, 2.176e-7, 9.708e-8, 4.289e-8, 2.222e-8])

    def __init__(self):
        g0, Rg, r0 = 9.80665, 287.053, 6356766.0
        Hb = np.array([0, 11, 20, 32, 47, 51, 71, 84.852]) * 1e3          # geopotential base (m)
        Lb = np.array([-6.5, 0.0, 1.0, 2.8, 0.0, -2.8, -2.0]) * 1e-3     # K/m
        Tb, Pb = [288.15], [101325.0]
        for i in range(7):
            dH = Hb[i + 1] - Hb[i]
            if abs(Lb[i]) < 1e-12:
                Tb.append(Tb[-1])
                Pb.append(Pb[-1] * math.exp(-g0 * dH / (Rg * Tb[-1])))
            else:
                T1 = Tb[-1] + Lb[i] * dH
                Pb.append(Pb[-1] * (T1 / Tb[-1]) ** (-g0 / (Rg * Lb[i])))
                Tb.append(T1)
        h = np.arange(0.0, 86.0001e3, 100.0)
        Hg = r0 * h / (r0 + h)
        lay = np.clip(np.searchsorted(Hb, Hg, side="right") - 1, 0, 6)
        T = np.empty_like(h)
        P = np.empty_like(h)
        for i in range(7):
            s = lay == i
            dH = Hg[s] - Hb[i]
            if abs(Lb[i]) < 1e-12:
                T[s] = Tb[i]
                P[s] = Pb[i] * np.exp(-g0 * dH / (Rg * Tb[i]))
            else:
                T[s] = Tb[i] + Lb[i] * dH
                P[s] = Pb[i] * (T[s] / Tb[i]) ** (-g0 / (Rg * Lb[i]))
        rho_low = P / (Rg * T)
        hk = np.concatenate([h / 1e3, self._UPPER_KM[1:]])
        lr = np.concatenate([np.log(rho_low), np.log(self._UPPER_RHO[1:])])
        self._hk = np.arange(0.0, 120.0001, 0.1)
        self._lr = np.interp(self._hk, hk, lr)
        self._scale_h = 1e3 * (self._hk[-1] - self._hk[-51]) / (self._lr[-51] - self._lr[-1])

    def rho(self, alt_m):
        a = np.asarray(alt_m, dtype=float) / 1e3
        lr = np.interp(a, self._hk, self._lr)
        hi = a > 120.0
        if np.any(hi):
            lr = np.where(hi, self._lr[-1] - (a - 120.0) * 1e3 / self._scale_h, lr)
        return np.exp(lr)


ATM = Atmosphere()


# ----------------------------------------------------------------------------------------
# luminous efficiency, sensor band, bolometric <-> Si-band conversion
# ----------------------------------------------------------------------------------------
def tau_fit(v, m):
    """tau(v, m) as a fraction (the published fit is in percent).  [lit R2, Eqs. 3-4]

        ln tau = 0.567 - 10.307 ln v + 9.781 ln^2 v - 3.0414 ln^3 v + 0.3213 ln^4 v + 0.347 tanh(0.38 ln m)   v < 25.372 km/s
        ln tau = -1.4286 + ln v + 0.347 tanh(0.38 ln m)                                                         v > 25.372 km/s
        (v in km/s, m in kg).  Velocity form: Pecina & Ceplecha (1983); mass form adapted from ReVelle & Ceplecha (2001).

    [new] The fit is only meaningful for meteor-like speeds: v is clamped to >= 2 km/s inside the polynomial and
    tau is tapered smoothly to zero between 2 and 4 km/s (the raw polynomial turns upward below ~2 km/s)."""
    v = np.asarray(v, dtype=float)
    vk = np.maximum(v, 2000.0) / 1000.0
    lv = np.log(vk)
    th = 0.347 * np.tanh(0.38 * np.log(np.maximum(m, 1e-10)))
    lt_lo = 0.567 - 10.307 * lv + 9.781 * lv ** 2 - 3.0414 * lv ** 3 + 0.3213 * lv ** 4 + th
    lt_hi = -1.4286 + lv + th
    tau = np.exp(np.where(vk < 25.372, lt_lo, lt_hi)) / 100.0
    taper = np.clip((v - 2000.0) / 2000.0, 0.0, 1.0) ** 2
    return tau * taper


def band_fraction(T, lo_um=0.4, hi_um=1.1, lam_um=None, resp=None):
    """In-band fraction of blackbody exitance at temperature T [K].  [T1: Planck's law]

        f(T) = integral_band  R(lam) B_lam(T) dlam  /  (sigma T^4 / pi),     B_lam = 2hc^2 / (lam^5 (exp(hc/lam k T) - 1))

    Default response R = 1 inside [lo_um, hi_um] (broadband silicon, 0.4-1.1 um) and 0 outside.  Pass lam_um/resp
    (peak-normalised relative spectral response) to use a measured sensor curve instead."""
    T = np.atleast_1d(np.asarray(T, dtype=float))
    if resp is None:
        lam = np.linspace(lo_um, hi_um, int(np.clip((hi_um - lo_um) / 0.0015, 600, 200000))) * 1e-6
        w = np.ones_like(lam)
    else:
        lam = np.asarray(lam_um, dtype=float) * 1e-6
        w = np.asarray(resp, dtype=float)
    x = H_PL * C_LIGHT / (lam[None, :] * K_B * T[:, None])
    B = 2.0 * H_PL * C_LIGHT ** 2 / lam[None, :] ** 5 / np.expm1(x)
    return _trapz(B * w[None, :], lam, axis=1) / (SIGMA_SB * T ** 4 / PI)


def bolometric_from_si(L_si, T_rad=6000.0, lo_um=0.4, hi_um=1.1):
    """L_bol = L_si / f_si(T_rad).  Assumes the radiator is a blackbody at T_rad."""
    return np.asarray(L_si) / float(band_fraction(T_rad, lo_um, hi_um)[0])


def si_from_bolometric(L_bol, T_rad=6000.0, lo_um=0.4, hi_um=1.1):
    """L_si = L_bol * f_si(T_rad)."""
    return np.asarray(L_bol) * float(band_fraction(T_rad, lo_um, hi_um)[0])


# ----------------------------------------------------------------------------------------
# parameters
# ----------------------------------------------------------------------------------------
@dataclass
class BolumParams:
    # --- required ---
    diameter: float = 1.0            # m
    velocity: float = 20.0e3         # m/s (atmosphere-relative)
    theta_deg: float = 45.0          # entry angle to local horizontal at zstart
    init_strength: float = 2.0e6     # Pa: stagnation pressure rho*v^2 at which the intact body yields [lit R2,R5-R8]
    density: float = 3300.0          # bulk density, kg/m^3
    # --- material / aerodynamics ---
    porosity: float = 0.1            # only used if init_strength <= 0
    Cd: float = 1.0                  # drag coefficient [lit R1,R7]
    # --- ablation (step 2) ---
    L_ablation: float = 5.0e6        # J/kg, effective heat of ablation Q [lit R1]
    Lambda: float = 0.05             # heat-transfer coefficient (<=1) on projected area [cal; canonical 0.1, R10]
    T_ablate: float = 2300.0         # K, ablation onset: net flux must exceed radiative loss [adapt R20]
    emissivity: float = 0.9
    # --- fragmentation / pancake engine (step 3) ---
    n_fragments: int = 4             # children per breakup [cal]
    frag_index: float = 1.0          # child mass share ~ k^-frag_index [new]
    weibull_alpha: float = 0.40      # S_child = S_parent (m_parent/m_child)^alpha [lit form R8,R11; cal value]
    c_spread: float = 0.5            # lateral spreading coefficient (dr/dt = c v sqrt(rho_air/rho_body)) [adapt R5,R6,R21; cal]
    f_split: float = 2.0             # split after the group's linear size grew by this factor [new]
    f_max: float = 8.0               # pancake limit (linear); R7 uses ~7 [adapt]
    max_groups: int = 128
    m_min_frag: float = 1.0          # kg; groups below this do not split further
    dust_frac: float = 0.20          # mass fraction released as fine dust at each breakup [new; motivated by R2,R4,R15; cal]
    dust_radius: float = 2.0e-2      # m, characteristic dust particle radius [cal]
    strength_secondary: float = None
    strength_tertiary: float = None
    strength_quaternary: float = None
    # --- luminosity / sensor ---
    tau_basis: str = "si"            # 'si': tau fit is the Si-band efficiency; 'bolometric': it is bolometric
    T_rad: float = 6000.0            # K, blackbody temperature for band conversion [lit R12 convention]
    sensor_lo_um: float = 0.4
    sensor_hi_um: float = 1.1
    eta_thick_bol: float = 0.20      # bolometric radiated fraction of a large, optically thick cloud [new; cal to R13,R17; R12]
    lum_r_half: float = 6.0          # m, projected radius where thin/thick efficiencies blend 50/50 [new; cal]
    eta_max: float = 0.6             # hard cap on radiated fraction
    t_cool: float = 0.0              # s, optional first-order radiative relaxation of the light curve [new]
    # --- numerics ---
    zstart: float = 100.0e3
    xstart: float = 0.0              # legacy, unused
    tstart: float = 0.0
    dt: float = 0.003                # maximum time step (s)
    dt_min: float = 5.0e-5
    tstop: float = 120.0
    v_dark: float = 3.0e3            # groups slower than this leave luminous flight (m/s)

    def strength0(self):
        if self.init_strength > 0:
            return self.init_strength
        return 2.0e5 * (1.0 - self.porosity) * (self.density / 3000.0)


_LEGACY_IGNORED = {"rho_tau_scale", "flare_duration", "n_frag_init", "q_factor"}


def rdinput(inputpath):
    """Read a bolum_in.dat deck (same format as before)."""
    raw = {}
    with open(inputpath, "r") as fh:
        lines = fh.readlines()[2:]
    for line in lines:
        if "=" not in line:
            continue
        key, rest = line.split("=", 1)
        val = rest.split("#")[0].strip()
        try:
            raw[key.strip()] = float(val)
        except ValueError:
            if val.lower() in ("si", "bolometric"):
                raw[key.strip()] = val.lower()
    valid = {f.name for f in fields(BolumParams)}
    kw = {}
    for k, v in raw.items():
        if k in valid:
            kw[k] = v
        elif k == "q_factor":
            print("[bolum] note: 'q_factor' is no longer an input; the loading is always rho*v^2 (stagnation pressure).")
        elif k not in _LEGACY_IGNORED:
            print(f"[bolum] warning: unknown input key '{k}' ignored")
    for req in ("diameter", "velocity", "theta_deg", "init_strength", "density"):
        if req not in kw:
            raise ValueError(f"Missing required input parameter: {req}")
    for ik in ("n_fragments", "max_groups"):
        if ik in kw:
            kw[ik] = int(kw[ik])
    return BolumParams(**kw)


# ----------------------------------------------------------------------------------------
# result container
# ----------------------------------------------------------------------------------------
@dataclass
class BolumResult:
    params: BolumParams
    t: np.ndarray
    L_si: np.ndarray
    L_bol: np.ndarray
    P_dep: np.ndarray
    alt: np.ndarray            # energy-weighted mean altitude of the radiating groups (m)
    v: np.ndarray              # mass-weighted speed (m/s)
    mass: np.ndarray           # total mass still in flight (kg)
    n_groups: np.ndarray
    q_max: np.ndarray          # largest loading rho*v^2 among groups (Pa)
    eta_eff: np.ndarray        # L_si / P_dep
    events: list               # dicts: kind, t, alt, q, mass, strength, gen
    alt_bins: np.ndarray       # bin centres (m)
    dEdh_dep: np.ndarray       # J per m of altitude
    dEdh_rad_si: np.ndarray
    summary: dict = field(default_factory=dict)


# ----------------------------------------------------------------------------------------
# the engine
# ----------------------------------------------------------------------------------------
def run_bolide(p: BolumParams, verbose=False) -> BolumResult:
    wall0 = time.time()
    rho_b = p.density
    m0 = 4.0 / 3.0 * PI * (0.5 * p.diameter) ** 3 * rho_b
    NG = int(p.max_groups)
    n_child = int(p.n_fragments)
    w_share = np.arange(1, n_child + 1, dtype=float) ** (-p.frag_index)
    w_share /= w_share.sum()
    q_loss = p.emissivity * SIGMA_SB * p.T_ablate ** 4
    m_dead = 1.0e-9 * m0
    f_si = float(band_fraction(p.T_rad, p.sensor_lo_um, p.sensor_hi_um)[0])
    basis_si = p.tau_basis.lower() == "si"
    legacy_S = [p.strength_secondary, p.strength_tertiary, p.strength_quaternary]

    m = np.zeros(NG); x = np.zeros(NG); z = np.zeros(NG)
    vx = np.zeros(NG); vz = np.zeros(NG)
    f = np.ones(NG); fb = np.ones(NG); S = np.zeros(NG)
    gen = np.zeros(NG, dtype=int)
    alive = np.zeros(NG, dtype=bool)
    yielded = np.zeros(NG, dtype=bool)
    dust = np.zeros(NG, dtype=bool)

    th = math.radians(p.theta_deg)
    m[0], z[0] = m0, RE + p.zstart
    vx[0], vz[0] = p.velocity * math.cos(th), -p.velocity * math.sin(th)
    S[0] = p.strength0()
    alive[0] = True
    KE0 = 0.5 * m0 * p.velocity ** 2

    bin_w = 500.0
    nb = int(150e3 / bin_w)
    dep_bins = np.zeros(nb)
    rad_bins = np.zeros(nb)

    T_ = []; L_si = []; P_ = []; ALT = []; V_ = []; M_ = []; N_ = []; Q_ = []
    events = []
    E_dep_tot = E_rad_tot = W_g_tot = KE_retired = M_retired = 0.0
    retired = []
    t = p.tstart
    nsteps = 0
    split_refused = 0
    n_groups_max = 1
    a_max = 0.0
    q_peak, alt_q_peak = 0.0, float("nan")

    def spawn_children(i, t_now, alt_now, q_now):
        """Replace group i by n_child children plus a dust group.  [new; concept R4,R8,R16]
        Mass conserving; total projected area ~continuous; child strength by Weibull scaling [R8,R11]."""
        free = np.flatnonzero(~alive)
        need = n_child if p.dust_frac > 0 else n_child - 1
        mp, fp, Sp, gp = m[i], f[i], S[i], gen[i]
        m_solid = mp * (1.0 - p.dust_frac)
        if len(free) < need or m_solid * w_share[-1] < p.m_min_frag:
            return False
        slots = [i] + list(free[: n_child - 1])
        xi_, zi_, vxi_, vzi_ = x[i], z[i], vx[i], vz[i]
        for k, s_ in enumerate(slots):
            mk = m_solid * w_share[k]
            fk = max(1.0, fp * (mk / mp) ** (1.0 / 6.0))
            gnew = gp + 1
            if gnew <= 3 and legacy_S[gnew - 1] is not None:
                Sk = max(Sp, float(legacy_S[gnew - 1]))
            else:
                Sk = Sp * (mp / mk) ** p.weibull_alpha
            m[s_], x[s_], z[s_], vx[s_], vz[s_] = mk, xi_, zi_, vxi_, vzi_
            f[s_], fb[s_], S[s_], gen[s_] = fk, fk, Sk, gnew
            alive[s_] = True
            yielded[s_] = False
            dust[s_] = False
        if p.dust_frac > 0:
            d_ = free[n_child - 1]
            m[d_], x[d_], z[d_], vx[d_], vz[d_] = mp * p.dust_frac, xi_, zi_, vxi_, vzi_
            f[d_], fb[d_], S[d_], gen[d_] = 1.0, np.inf, np.inf, gp + 1
            alive[d_] = True
            yielded[d_] = False
            dust[d_] = True
        events.append(dict(kind="split", t=t_now, alt=alt_now, q=q_now, mass=mp, strength=Sp, gen=gp))
        return True

    while alive.any() and t < p.tstop:
        ia = np.flatnonzero(alive)
        mi, xi, zi = m[ia], x[ia], z[ia]
        du = dust[ia]
        vxi, vzi, fi, Si = vx[ia], vz[ia], f[ia], S[ia]
        ri = np.hypot(xi, zi)
        alt = ri - RE
        vi = np.hypot(vxi, vzi)
        rho = ATM.rho(alt)
        q = rho * vi * vi                                        # stagnation-pressure loading [lit R5-R8, R2]
        rs = (3.0 * mi / (4.0 * PI * rho_b)) ** (1.0 / 3.0)
        A = np.where(du, 3.0 * mi / (4.0 * rho_b * p.dust_radius), PI * rs * rs * fi * fi)   # projected area
        F = 0.5 * p.Cd * rho * vi * vi * A                       # drag [lit R1,R7]
        a_dr = F / mi
        # --- heating / ablation (step 2) [adapt R1, R10, R20] ---
        qh = p.Lambda * 0.5 * rho * vi ** 3
        mdot = np.maximum(qh - q_loss, 0.0) * A / p.L_ablation
        # --- yield + pancake spreading (step 3) [lit R5-R8] ---
        yl = q >= Si
        newly = yl & ~yielded[ia]
        if newly.any():
            for j in np.flatnonzero(newly):
                events.append(dict(kind="yield", t=t, alt=alt[j], q=q[j], mass=mi[j], strength=Si[j], gen=int(gen[ia[j]])))
            yielded[ia[newly]] = True
        dfdt = np.where(yl, p.c_spread * vi * np.sqrt(rho / rho_b) / rs, 0.0)
        # --- adaptive step ---
        with np.errstate(divide="ignore", invalid="ignore"):
            d1 = np.min(np.where(mdot > 0, np.where(du, 0.3, 0.03) * mi / mdot, np.inf))
            d2 = np.min(np.where(a_dr > 0, 0.02 * vi / a_dr, np.inf))
            d3 = np.min(np.where(dfdt > 0, 0.05 * fi / dfdt, np.inf))
        dt = max(min(p.dt, d1, d2, d3), p.dt_min)
        # --- update ---
        gr = GM / ri ** 2
        gxv, gzv = -gr * xi / ri, -gr * zi / ri
        dm = mi * (-np.expm1(-mdot * dt / mi))
        mnew = mi - dm
        dead = mnew < m_dead
        dm = np.where(dead, mi, dm)
        mnew = np.where(dead, 0.0, mnew)
        dvx_d, dvz_d = -a_dr * vxi / vi * dt, -a_dr * vzi / vi * dt
        vxn = vxi + dvx_d + gxv * dt
        vzn = vzi + dvz_d + gzv * dt
        vbx, vbz = 0.5 * (vxi + vxn), 0.5 * (vzi + vzn)
        E_drag = -mnew * (vbx * dvx_d + vbz * dvz_d)
        E_abl = 0.5 * dm * vi * vi
        W_g = mnew * (vbx * gxv + vbz * gzv) * dt
        E_dep = E_drag + E_abl
        xn, zn = xi + vbx * dt, zi + vbz * dt
        fn = np.minimum(fi + dfdt * dt, p.f_max)
        vn = np.hypot(vxn, vzn)
        altn = np.hypot(xn, zn) - RE
        # --- radiation: thin (tau fit) <-> thick (large cloud) blend [new; R2, R12, R13] ---
        tau = tau_fit(vi, mi)
        r_proj = np.sqrt(A / PI)
        w_thick = r_proj / (r_proj + p.lum_r_half)
        eta_thin = tau if basis_si else tau * f_si
        eta_thick = p.eta_thick_bol * f_si * np.clip((vi - 4.0e3) / 4.0e3, 0.0, 1.0)
        eta = np.minimum((1.0 - w_thick) * eta_thin + w_thick * eta_thick, p.eta_max * f_si)
        E_rad_si = eta * E_dep
        # --- ledgers / diagnostics ---
        E_dep_s, W_g_s = float(E_dep.sum()), float(W_g.sum())
        E_dep_tot += E_dep_s
        W_g_tot += W_g_s
        E_rad_tot += float(E_rad_si.sum())
        bi = np.clip((alt / bin_w).astype(int), 0, nb - 1)
        np.add.at(dep_bins, bi, E_dep)
        np.add.at(rad_bins, bi, E_rad_si)
        jq = int(np.argmax(q))
        if q[jq] > q_peak:
            q_peak, alt_q_peak = float(q[jq]), float(alt[jq])
        if (~du).any():
            a_max = max(a_max, float(a_dr[~du].max()))
        t += dt
        T_.append(t); L_si.append(E_rad_si.sum() / dt); P_.append(E_dep_s / dt)
        wts = E_dep if E_dep_s > 0 else mi
        ALT.append(float((wts * alt).sum() / wts.sum()))
        V_.append(float((mi * vi).sum() / mi.sum()) if mi.sum() > 0 else 0.0)
        M_.append(float(mnew.sum())); N_.append(int((mnew > 0).sum())); Q_.append(float(q.max()))
        n_groups_max = max(n_groups_max, int(alive.sum()))
        # --- commit state ---
        m[ia], x[ia], z[ia], vx[ia], vz[ia], f[ia] = mnew, xn, zn, vxn, vzn, fn
        gone = (mnew <= 0)
        alive[ia[gone]] = False
        out = (~gone) & ((altn <= 0.0) | (vn < p.v_dark))
        for j in np.flatnonzero(out):
            k = ia[j]
            ke = 0.5 * mnew[j] * vn[j] ** 2
            KE_retired += ke
            M_retired += mnew[j]
            retired.append(dict(mass=mnew[j], v=vn[j], alt=altn[j], t=t, ke=ke, ground=bool(altn[j] <= 0), dust=bool(du[j])))
            alive[k] = False
        # --- splitting ---
        cand = np.flatnonzero(alive & yielded & (f >= fb * p.f_split))
        for k in cand:
            altk = float(np.hypot(x[k], z[k]) - RE)
            vk = float(np.hypot(vx[k], vz[k]))
            qk = float(ATM.rho(altk)) * vk * vk
            if not spawn_children(k, t, altk, qk):
                fb[k] = np.inf                                   # cannot split: keep as a pancake
                split_refused += 1
        nsteps += 1

    # ------------------------------------------------------------------ assemble outputs
    T_ = np.array(T_); L_si = np.array(L_si); P_ = np.array(P_); ALT = np.array(ALT)
    if p.t_cool > 0.0 and len(T_) > 1:                          # optional energy-conserving relaxation
        dts = np.diff(np.concatenate([[p.tstart], T_]))
        Lr = np.empty_like(L_si); prev = 0.0
        for k in range(len(L_si)):
            prev += (L_si[k] - prev) * (1.0 - math.exp(-dts[k] / p.t_cool))
            Lr[k] = prev
        L_si = Lr
    L_bol = L_si / f_si                                         # bolometric = Si-band / in-band blackbody fraction
    dts = np.diff(np.concatenate([[p.tstart], T_])) if len(T_) else np.array([])
    E_rad_si_tot = float((L_si * dts).sum()) if len(T_) else 0.0
    ke_end_active = float(sum(0.5 * m[k] * (vx[k] ** 2 + vz[k] ** 2) for k in np.flatnonzero(alive)))
    closure = (KE0 + W_g_tot - E_dep_tot - KE_retired - ke_end_active) / KE0
    alt_c = (np.arange(nb) + 0.5) * bin_w
    dEdh = dep_bins / bin_w
    sm = np.convolve(dEdh, np.ones(3) / 3.0, mode="same")
    ipk = int(np.argmax(sm)) if sm.max() > 0 else 0
    half = sm >= 0.5 * sm[ipk]
    lo = hi = ipk
    while lo > 0 and half[lo - 1]:
        lo -= 1
    while hi < nb - 1 and half[hi + 1]:
        hi += 1
    cum_top = np.cumsum(dep_bins[::-1])[::-1]
    ih = int(np.flatnonzero(cum_top >= 0.5 * dep_bins.sum())[-1]) if dep_bins.sum() > 0 else 0
    ipl = int(np.argmax(L_si)) if len(L_si) else 0
    lit_end = np.flatnonzero(L_si > 0.01 * L_si.max()) if len(L_si) else np.array([], dtype=int)
    lit_start = np.flatnonzero(L_si > 1e-6 * L_si.max()) if len(L_si) else np.array([], dtype=int)
    first_yield = next((e for e in events if e["kind"] == "yield"), None)
    summary = dict(
        KE0_J=KE0, KE0_kt=KE0 / JKT, m0=m0,
        E_dep_J=E_dep_tot, E_dep_kt=E_dep_tot / JKT,
        E_rad_si_J=E_rad_si_tot, E_rad_bol_J=E_rad_si_tot / f_si,
        E_rad_bol_kt=E_rad_si_tot / f_si / JKT,
        band_fraction=f_si,
        eta_si_overall=E_rad_si_tot / E_dep_tot if E_dep_tot > 0 else 0.0,
        eta_bol_overall=E_rad_si_tot / f_si / E_dep_tot if E_dep_tot > 0 else 0.0,
        peak_L_si_W=float(L_si.max()) if len(L_si) else 0.0,
        peak_I_si_W_sr=float(L_si.max()) / (4 * PI) if len(L_si) else 0.0,
        peak_L_bol_W=float(L_bol.max()) if len(L_bol) else 0.0,
        alt_peak_L_km=float(ALT[ipl]) / 1e3 if len(ALT) else float("nan"),
        t_peak_L=float(T_[ipl]) if len(T_) else float("nan"),
        v_at_peak_L_kms=float(V_[ipl]) / 1e3 if len(V_) else float("nan"),
        alt_first_light_km=float(ALT[lit_start[0]]) / 1e3 if len(lit_start) else float("nan"),
        t_first_light=float(T_[lit_start[0]]) if len(lit_start) else float("nan"),
        alt_end_lum_km=float(ALT[lit_end[-1]]) / 1e3 if len(lit_end) else float("nan"),
        alt_peak_dEdh_km=alt_c[ipk] / 1e3,
        peak_dEdh_kt_per_km=sm[ipk] * 1e3 / JKT,
        half_range_km=(alt_c[lo] / 1e3, alt_c[hi] / 1e3),
        alt_half_energy_km=alt_c[ih] / 1e3,
        first_yield_alt_km=first_yield["alt"] / 1e3 if first_yield else float("nan"),
        first_yield_q_MPa=first_yield["q"] / 1e6 if first_yield else float("nan"),
        first_yield_t=first_yield["t"] if first_yield else float("nan"),
        n_split_events=sum(1 for e in events if e["kind"] == "split"),
        split_refused=split_refused, n_groups_max=n_groups_max,
        q_max_MPa=q_peak / 1e6, alt_q_max_km=alt_q_peak / 1e3, max_decel_g=a_max / G0,
        mass_survived_kg=M_retired, mass_survived_frac=M_retired / m0,
        mass_survived_solid_frac=sum(r["mass"] for r in retired if not r["dust"]) / m0,
        mass_survived_dust_frac=sum(r["mass"] for r in retired if r["dust"]) / m0,
        ground_ke_J=float(sum(r["ke"] for r in retired if r["ground"])),
        energy_closure=closure, reached_tstop=bool(alive.any()),
        steps=nsteps, wall_s=time.time() - wall0,
        duration_s=float(T_[-1] - p.tstart) if len(T_) else 0.0,
        retired=retired,
    )
    warns = []
    if abs(closure) > 1e-6:
        warns.append(f"energy ledger closure {closure:.1e} exceeds 1e-6")
    if split_refused:
        warns.append(f"{split_refused} breakup(s) refused (group cap max_groups={NG} or m_min_frag) - raise max_groups?")
    if summary["reached_tstop"]:
        warns.append(f"stopped at tstop={p.tstop:.0f} s with groups still in flight")
    if summary["E_dep_kt"] < 0.5 * summary["KE0_kt"]:
        warns.append("less than half of the initial kinetic energy was deposited (body survives/skips?)")
    summary["warnings"] = warns
    return BolumResult(p, T_, L_si, L_bol, P_, ALT, np.array(V_), np.array(M_), np.array(N_),
                       np.array(Q_), L_si / np.maximum(P_, 1e-30), events, alt_c, dEdh, rad_bins / bin_w, summary)


# ----------------------------------------------------------------------------------------
# command-line report
# ----------------------------------------------------------------------------------------
def format_report(res: BolumResult, max_events=12):
    p, s = res.params, res.summary
    L = []
    add = L.append
    add("=" * 100)
    add(f" bolum v{__version__} run report")
    add("=" * 100)
    add(" INPUT")
    add(f"   body            d = {p.diameter:.3f} m   m0 = {s['m0']:.4e} kg   rho = {p.density:.0f} kg/m3")
    add(f"   entry           v = {p.velocity/1e3:.3f} km/s   theta = {p.theta_deg:.2f} deg   start alt = {p.zstart/1e3:.0f} km")
    add(f"   KE0             {s['KE0_J']:.4e} J = {s['KE0_kt']:.4f} kt TNT")
    add(f"   strength        S0 = {p.strength0()/1e6:.4f} MPa  (loading = rho*v^2, stagnation pressure); Weibull alpha = {p.weibull_alpha}")
    add(f"   ablation        Lambda = {p.Lambda}  Q = {p.L_ablation:.2e} J/kg  -> sigma = Lambda/(2Q) = {p.Lambda/(2*p.L_ablation)*1e6:.4f} s^2/km^2")
    add(f"   fragmentation   {p.n_fragments} children/breakup, dust {p.dust_frac*100:.0f}% (r = {p.dust_radius*100:.1f} cm), spread c = {p.c_spread}, f_split = {p.f_split}, f_max = {p.f_max}")
    add(f"   sensor          band {p.sensor_lo_um}-{p.sensor_hi_um} um at T_rad = {p.T_rad:.0f} K -> f_si = {s['band_fraction']:.4f}   "
        f"(L_bol = L_si / f_si);  tau basis = {p.tau_basis}")
    add(f"   luminosity      eta = (1-w)*tau(v,m) + w*{p.eta_thick_bol:.2f}*f_si,  w = r/(r + {p.lum_r_half:g} m);  t_cool = {p.t_cool:g} s")
    add(" FRAGMENTATION")
    add(f"   first yield     alt {s['first_yield_alt_km']:.2f} km  at loading {s['first_yield_q_MPa']:.3f} MPa  (t = {s['first_yield_t']:.3f} s)")
    add(f"   breakups        {s['n_split_events']} (refused: {s['split_refused']});  max simultaneous groups {s['n_groups_max']}")
    add(f"   max loading     {s['q_max_MPa']:.2f} MPa at {s['alt_q_max_km']:.1f} km;   max deceleration of solid groups {s['max_decel_g']:.1f} g")
    if res.events:
        add("   first events:      kind    t(s)   alt(km)  loading(MPa)    mass(kg)  strength(MPa)  gen")
        for e in res.events[:max_events]:
            add(f"                    {e['kind']:>7s} {e['t']:7.3f} {e['alt']/1e3:9.2f} {e['q']/1e6:12.3f} {e['mass']:11.3e} {e['strength']/1e6:13.3f} {e['gen']:4d}")
        if len(res.events) > max_events:
            add(f"                    ... {len(res.events)-max_events} more in events.dat")
    add(" ENERGY")
    add(f"   deposited       {s['E_dep_J']:.4e} J = {s['E_dep_kt']:.3f} kt  ({100*s['E_dep_kt']/s['KE0_kt']:.2f}% of KE0; gravity work adds a little)")
    add(f"   dE/dh           peak {s['peak_dEdh_kt_per_km']:.1f} kt/km at {s['alt_peak_dEdh_km']:.1f} km;  half-peak range {s['half_range_km'][0]:.1f}-{s['half_range_km'][1]:.1f} km;  "
        f"half of the energy deposited above {s['alt_half_energy_km']:.1f} km")
    add(f"   ledger closure  (KE0 + W_grav - E_dep - KE_left)/KE0 = {s['energy_closure']:.2e}")
    add(" LIGHT")
    add(f"   first light     L_si > 1e-6 of peak at {s['alt_first_light_km']:.1f} km (t = {s['t_first_light']:.2f} s);  luminous end (L < 1% of peak) at {s['alt_end_lum_km']:.1f} km")
    add(f"   peak            L_si = {s['peak_L_si_W']:.3e} W  ({s['peak_I_si_W_sr']:.3e} W/sr)   L_bol = {s['peak_L_bol_W']:.3e} W")
    add(f"                   at t = {s['t_peak_L']:.2f} s, ~{s['alt_peak_L_km']:.1f} km, mass-weighted speed {s['v_at_peak_L_kms']:.1f} km/s")
    add(f"   radiated        Si band {s['E_rad_si_J']:.3e} J;  bolometric {s['E_rad_bol_J']:.3e} J = {s['E_rad_bol_kt']:.4f} kt")
    add(f"   efficiency      bolometric {100*s['eta_bol_overall']:.2f}% / Si-band {100*s['eta_si_overall']:.2f}% of deposited energy")
    kt_b = 8.2508 * max(s['E_rad_bol_kt'], 1e-12) ** 0.885
    add(f"   Brown+2002 check  E_total = 8.2508 E_rad^0.885 -> {kt_b:.3f} kt  (model KE0 = {s['KE0_kt']:.3f} kt, ratio {kt_b/s['KE0_kt']:.2f}; population scatter ~x2)")
    add(" SURVIVORS")
    add(f"   solid {100*s['mass_survived_solid_frac']:.2f}% and dust {100*s['mass_survived_dust_frac']:.2f}% of m0 left luminous flight (< {p.v_dark/1e3:.1f} km/s) or reached the ground;"
        f"  ground-impact KE {s['ground_ke_J']:.2e} J")
    add(" RUN")
    add(f"   {s['steps']} steps, simulated {s['duration_s']:.1f} s, wall time {s['wall_s']:.2f} s")
    if s["warnings"]:
        for w_ in s["warnings"]:
            add(f"   WARNING: {w_}")
    else:
        add("   no warnings")
    add("=" * 100)
    return "\n".join(L)


# ----------------------------------------------------------------------------------------
# output files
# ----------------------------------------------------------------------------------------
def make_header(p: BolumParams):
    return (f"v = {p.velocity/1e3:.2f} km/s   d = {p.diameter:.2f} m   theta = {p.theta_deg:.1f} deg   "
            f"S0 = {p.strength0()/1e6:.2f} MPa (rho v^2)   rho = {p.density:.0f} kg/m3")


def write_outputs(outdir, res: BolumResult):
    p, s = res.params, res.summary
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "init_params.dat"), "w") as fh:
        fh.write(f"Initialization Parameters for Bolide Luminosity and Fragmentation Simulation (bolum v{__version__})\n")
        for fl in fields(BolumParams):
            fh.write(f"{fl.name:22s} {getattr(p, fl.name)}\n")
    with open(os.path.join(outdir, "bolide_out.dat"), "w") as fh:
        fh.write("# time_s  alt_km  v_kms  mass_kg  n_groups  P_dep_W  L_si_W  L_bol_W  eta_si  loading_max_Pa\n")
        for k in range(len(res.t)):
            fh.write(f"{res.t[k]:.5f} {res.alt[k]/1e3:.4f} {res.v[k]/1e3:.4f} {res.mass[k]:.5e} {res.n_groups[k]:d} "
                     f"{res.P_dep[k]:.5e} {res.L_si[k]:.5e} {res.L_bol[k]:.5e} {res.eta_eff[k]:.5f} {res.q_max[k]:.5e}\n")
    with open(os.path.join(outdir, "power_vs_time.dat"), "w") as fh:
        fh.write("# L_si_W  time_s  L_bol_W\n")
        for k in range(len(res.t)):
            fh.write(f"{res.L_si[k]:.4e}   {res.t[k]:.4e}   {res.L_bol[k]:.4e}\n")
    with open(os.path.join(outdir, "power_vs_alt.dat"), "w") as fh:
        fh.write("# L_si_W  alt_m  L_bol_W\n")
        for k in range(len(res.t)):
            fh.write(f"{res.L_si[k]:.4e}   {res.alt[k]:.4e}   {res.L_bol[k]:.4e}\n")
    with open(os.path.join(outdir, "dEdh.dat"), "w") as fh:
        fh.write("# alt_km  dE_dep/dh_kt_per_km  dE_rad_si/dh_kt_per_km\n")
        for a, d, r in zip(res.alt_bins, res.dEdh_dep, res.dEdh_rad_si):
            fh.write(f"{a/1e3:.2f} {d*1e3/JKT:.5e} {r*1e3/JKT:.5e}\n")
    with open(os.path.join(outdir, "events.dat"), "w") as fh:
        fh.write("# kind  t_s  alt_km  loading_MPa  mass_kg  strength_MPa  generation\n")
        for e in res.events:
            fh.write(f"{e['kind']} {e['t']:.4f} {e['alt']/1e3:.3f} {e['q']/1e6:.4f} {e['mass']:.4e} {e['strength']/1e6:.4f} {e['gen']}\n")
    with open(os.path.join(outdir, "summary.dat"), "w") as fh:
        for k, v in s.items():
            if k not in ("retired", "warnings"):
                fh.write(f"{k:26s} {v}\n")
    with open(os.path.join(outdir, "run_report.txt"), "w") as fh:
        fh.write(format_report(res, max_events=40) + "\n")


# ----------------------------------------------------------------------------------------
# standardized plotting
# ----------------------------------------------------------------------------------------
C_SI, C_BOL, C_DEP, C_PUB, C_BRK, C_SAT = "#1f77b4", "#ff7f0e", "black", "crimson", "seagreen", "purple"


def _alt_window(res):
    """Altitude window common to all altitude panels: where L_si > 1e-6 of peak, padded."""
    lit = np.flatnonzero(res.L_si > 1e-6 * res.L_si.max()) if len(res.L_si) else []
    if len(lit) == 0:
        return 0.0, 100.0
    hi = res.alt[lit[0]] / 1e3 + 4.0
    lo = max(0.0, res.alt[lit[-1]] / 1e3 - 4.0)
    return lo, min(max(hi, lo + 10.0), 125.0)


def _time_log_window(res, decades=10.0):
    """(ymin, ymax, t_lo, t_hi) for log power-vs-time plots: ~10 decades below the peak so early-time signal shows."""
    pk = res.L_si.max()
    floor = max(pk * 10 ** (-decades), 1e-3)
    ok = np.flatnonzero(res.L_si > floor)
    t = res.t - res.params.tstart
    if len(ok) == 0:
        return floor, 3 * pk, t[0], t[-1]
    return floor, 3 * res.L_bol.max(), max(t[ok[0]] - 0.02 * (t[ok[-1]] - t[ok[0]]), 0.0), t[ok[-1]] + 0.02 * (t[ok[-1]] - t[ok[0]])


def model_metric_rows(res):
    s = res.summary
    return [
        ("KE0 (kt)", f"{s['KE0_kt']:.4g}", "", None),
        ("deposited energy (kt)", f"{s['E_dep_kt']:.4g}", "", None),
        ("first yield alt (km)", f"{s['first_yield_alt_km']:.1f}", "", None),
        ("breakups", f"{s['n_split_events']}", "", None),
        ("peak dE/dh (kt/km @ km)", f"{s['peak_dEdh_kt_per_km']:.1f} @ {s['alt_peak_dEdh_km']:.1f}", "", None),
        ("half-peak range (km)", f"{s['half_range_km'][0]:.1f}-{s['half_range_km'][1]:.1f}", "", None),
        ("peak L_Si (W)", f"{s['peak_L_si_W']:.2e}", "", None),
        ("peak I_Si (W/sr)", f"{s['peak_I_si_W_sr']:.2e}", "", None),
        ("alt of peak light (km)", f"{s['alt_peak_L_km']:.1f}", "", None),
        ("E_rad bolometric (kt)", f"{s['E_rad_bol_kt']:.4g}", "", None),
        ("bolometric efficiency", f"{100*s['eta_bol_overall']:.1f}%", "", None),
        ("luminous end alt (km)", f"{s['alt_end_lum_km']:.1f}", "", None),
        ("energy-ledger closure", f"{s['energy_closure']:.1e}", "", None),
    ]


def _banner(fig, header, title=None):
    fig.text(.01, .985, "UNCLASSIFIED", fontsize=9, color="green", va="top")
    fig.text(.99, .012, "UNCLASSIFIED", fontsize=9, color="green", ha="right")
    if title:
        fig.suptitle(title, fontsize=13, y=0.995)
    fig.text(.5, .955, header, fontsize=9, ha="center")


def plot_event_summary(res, ref=None, rows=None, sat_W=None, title=None):
    """Standardized 2x3 per-event figure (identical layout for every event).

       (a) power vs time, log y      (b) Si power vs altitude    (c) dE/dh vs altitude
       (d) speed & mass vs altitude  (e) loading vs altitude     (f) metrics table (model [vs published])
    ref (optional) keys: peak_L_W, peak_L_alt_km, peak_dEdh=(alt_km,kt/km), half_range_km=(lo,hi),
                         flares_km=[...], end_km, first_break_km=(lo,hi)
    """
    import matplotlib.pyplot as plt
    ref = ref or {}
    p, s = res.params, res.summary
    t = res.t - p.tstart
    alt_km = res.alt / 1e3
    lo, hi = _alt_window(res)
    ymin, ymax, tlo, thi = _time_log_window(res)
    fig, ax = plt.subplots(2, 3, figsize=(17, 9.6))
    fig.subplots_adjust(left=.055, right=.985, top=.90, bottom=.07, wspace=.28, hspace=.30)
    # (a) power vs time
    a = ax[0, 0]
    a.semilogy(t, np.maximum(res.L_si, 1e-30), color=C_SI, label="Si band")
    a.semilogy(t, np.maximum(res.L_bol, 1e-30), color=C_BOL, alpha=.75, label="bolometric")
    if sat_W:
        a.axhline(sat_W, color=C_SAT, ls=":", label=f"sensor saturation {sat_W:.1e} W")
        a.semilogy(t, np.maximum(np.minimum(res.L_si, sat_W), 1e-30), color=C_SAT, lw=.8, alpha=.8, label="Si band, clipped")
    if "peak_L_W" in ref:
        a.axhline(ref["peak_L_W"], color=C_PUB, ls="--", lw=1, label="published peak (Si band)")
    a.set_ylim(ymin, ymax); a.set_xlim(tlo, thi)
    a.set_xlabel("time since start (s)"); a.set_ylabel("power (W)"); a.set_title("(a) light curve"); a.legend(fontsize=8, loc="lower right"); a.grid(True, alpha=.3, which="both")
    # (b) L vs altitude
    a = ax[0, 1]
    a.semilogx(np.maximum(res.L_si, 1e-30), alt_km, color=C_SI, label="Si band")
    a.semilogx(np.maximum(res.L_bol, 1e-30), alt_km, color=C_BOL, alpha=.75, label="bolometric")
    if "peak_L_W" in ref and "peak_L_alt_km" in ref:
        a.plot([ref["peak_L_W"]], [ref["peak_L_alt_km"]], "*", color=C_PUB, ms=15, label="published peak")
    for h in ref.get("flares_km", []):
        a.axhline(h, color=C_PUB, ls="--", lw=.8, alpha=.7)
    if "end_km" in ref:
        a.axhline(ref["end_km"], color="gray", ls="--", lw=1, label="published luminous end")
    a.set_xlim(max(1e-8 * res.L_si.max(), 1e-3), 3 * res.L_bol.max()); a.set_ylim(lo, hi)
    a.set_xlabel("power (W)"); a.set_ylabel("altitude (km)"); a.set_title("(b) power vs altitude"); a.legend(fontsize=8); a.grid(True, alpha=.3, which="both")
    # (c) dE/dh
    a = ax[0, 2]
    a.plot(res.dEdh_dep * 1e3 / JKT, res.alt_bins / 1e3, color=C_DEP, label="deposited")
    a.plot(res.dEdh_rad_si * 1e3 / JKT, res.alt_bins / 1e3, color=C_SI, label="radiated (Si band)")
    if "half_range_km" in ref:
        a.axhspan(ref["half_range_km"][0], ref["half_range_km"][1], color="orange", alpha=.2, label="published half-energy range")
    if "peak_dEdh" in ref:
        a.plot([ref["peak_dEdh"][1]], [ref["peak_dEdh"][0]], "*", color=C_PUB, ms=15, label="published peak")
    for h in ref.get("flares_km", []):
        a.axhline(h, color=C_PUB, ls="--", lw=.8, alpha=.7)
    a.set_ylim(lo, hi); a.set_xlabel("dE/dh (kt TNT / km)"); a.set_ylabel("altitude (km)"); a.set_title("(c) energy deposition"); a.legend(fontsize=8); a.grid(True, alpha=.3)
    # (d) speed and mass
    a = ax[1, 0]
    a.plot(res.v / p.velocity, alt_km, color="firebrick", label="speed / v0 (mass-weighted)")
    a.plot(res.mass / s["m0"], alt_km, color="k", label="mass in flight / m0")
    a.set_xlim(0, 1.05); a.set_ylim(lo, hi); a.set_xlabel("fraction of initial value"); a.set_ylabel("altitude (km)")
    a.set_title("(d) speed and mass"); a.legend(fontsize=8); a.grid(True, alpha=.3)
    # (e) loading
    a = ax[1, 1]
    a.semilogx(res.q_max / 1e6, alt_km, color="slateblue", label="max loading, rho v^2")
    sp = [e for e in res.events if e["kind"] == "split"]
    if sp:
        a.scatter([e["q"] / 1e6 for e in sp], [e["alt"] / 1e3 for e in sp], s=[12 + 70 * (e["mass"] / s["m0"]) ** 0.4 for e in sp],
                  color=C_BRK, alpha=.6, label="breakups (size ~ parent mass)", zorder=3)
    a.axvline(p.strength0() / 1e6, color="gray", ls=":", label="initial strength S0")
    if "first_break_km" in ref:
        a.axhspan(ref["first_break_km"][0], ref["first_break_km"][1], color="orange", alpha=.2, label="published first breakup")
    a.set_ylim(lo, hi); a.set_xlabel("loading (MPa)"); a.set_ylabel("altitude (km)"); a.set_title("(e) loading and breakups"); a.legend(fontsize=8); a.grid(True, alpha=.3, which="both")
    # (f) table
    a = ax[1, 2]
    a.axis("off")
    rows = rows or model_metric_rows(res)
    has_pub = any(r[2] for r in rows)
    a.text(0.0, 1.0, "(f) metrics" + ("  (model | published | status)" if has_pub else "  (model)"), fontsize=11, va="top", weight="bold", transform=a.transAxes)
    y = 0.93
    for name, mod, pub, ok in rows:
        a.text(0.0, y, name, fontsize=8.8, va="top", family="monospace", transform=a.transAxes)
        a.text(0.46, y, mod, fontsize=8.8, va="top", family="monospace", transform=a.transAxes)
        if pub:
            a.text(0.70, y, pub, fontsize=8.8, va="top", family="monospace", transform=a.transAxes)
        if ok is not None:
            a.text(0.97, y, "ok" if ok else "MISS", fontsize=8.8, va="top", family="monospace", ha="right",
                   color=("green" if ok else "red"), weight="bold", transform=a.transAxes)
        y -= 0.062
    _banner(fig, make_header(p), title or "Bolide event summary")
    return fig


def plot_outputs(outdir, res: BolumResult, ref=None, rows=None, sat_W=None, title=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    p, s = res.params, res.summary
    header = make_header(p)
    t = res.t - p.tstart
    pdf = PdfPages(os.path.join(outdir, "bolide_plots.pdf"))
    fig = plot_event_summary(res, ref=ref, rows=rows, sat_W=sat_W, title=title)
    pdf.savefig(fig); plt.close(fig)

    def page(title_):
        fig_ = plt.figure(figsize=(9.5, 6.8))
        plt.title(title_ + "\n")
        plt.grid(True, alpha=.4, which="both")
        _banner(fig_, header)
        return fig_

    def save(fig_):
        pdf.savefig(fig_); plt.close(fig_)

    notes = (f"E_rad bolometric = {s['E_rad_bol_kt']:.3f} kt = {s['E_rad_bol_J']:.3e} J\n"
             f"E_rad Si band    = {s['E_rad_si_J']:.3e} J\n"
             f"E deposited      = {s['E_dep_kt']:.3f} kt  (KE0 = {s['KE0_kt']:.3f} kt)\n"
             f"peak L_Si        = {s['peak_L_si_W']:.3e} W  ({s['peak_I_si_W_sr']:.3e} W/sr)")
    for scale in ("linear", "log"):
        fig = page(f"Light curve ({scale})")
        plt.plot(t, res.L_si, color=C_SI, label="Si band")
        plt.plot(t, res.L_bol, color=C_BOL, alpha=.75, label="bolometric")
        if sat_W:
            plt.axhline(sat_W, color=C_SAT, ls=":", label=f"sensor saturation {sat_W:.1e} W")
        if scale == "log":
            ymin, ymax, tlo, thi = _time_log_window(res)
            plt.yscale("log"); plt.ylim(ymin, ymax); plt.xlim(tlo, thi)
        fy = next((e for e in res.events if e["kind"] == "yield"), None)
        if fy:
            plt.axvline(fy["t"] - p.tstart, ls="--", color=C_BRK, alpha=.6, label=f"first yield {fy['alt']/1e3:.1f} km")
        plt.xlabel("Time (s)"); plt.ylabel("Luminosity (W)"); plt.legend(loc="upper right")
        fig.text(.14, .66, notes, fontsize=8.5, va="top", family="monospace")
        save(fig)
    fig = page("Radiant intensity vs time (isotropic emitter)")
    plt.plot(t, res.L_si / (4 * PI), color=C_SI); plt.yscale("log")
    plt.ylim(max(1e-10 * res.L_si.max() / (4 * PI), 1e-3), 3 * res.L_si.max() / (4 * PI))
    plt.xlabel("Time (s)"); plt.ylabel("Si-band radiant intensity (W/sr)")
    save(fig)
    lo, hi = _alt_window(res)
    fig = page("Power vs altitude")
    plt.plot(res.L_si, res.alt / 1e3, color=C_SI, label="Si band")
    plt.plot(res.L_bol, res.alt / 1e3, color=C_BOL, alpha=.75, label="bolometric")
    plt.xscale("log"); plt.xlim(max(1e-8 * res.L_si.max(), 1e-3), 3 * res.L_bol.max()); plt.ylim(lo, hi)
    plt.xlabel("Luminosity (W)"); plt.ylabel("Altitude (km)"); plt.legend()
    save(fig)
    fig = page("Energy deposition vs altitude")
    plt.plot(res.dEdh_dep * 1e3 / JKT, res.alt_bins / 1e3, color=C_DEP, label="deposited")
    plt.plot(res.dEdh_rad_si * 1e3 / JKT, res.alt_bins / 1e3, color=C_SI, label="radiated (Si band)")
    plt.ylim(lo, hi); plt.xlabel("dE/dh (kt TNT per km)"); plt.ylabel("Altitude (km)"); plt.legend()
    save(fig)
    fig = page("Mass in flight and fragment-group count")
    ax = plt.gca(); ax.plot(t, res.mass, color="k"); ax.set_yscale("log")
    ax.set_xlabel("Time (s)"); ax.set_ylabel("Total mass in flight (kg)")
    ax2 = ax.twinx(); ax2.plot(t, res.n_groups, color="red", alpha=.7); ax2.set_ylabel("Fragment groups", color="red")
    save(fig)
    fig = page("Loading rho v^2 and speed vs altitude")
    plt.plot(res.alt / 1e3, res.q_max / 1e6, color="slateblue"); plt.yscale("log"); plt.xlim(hi, lo)
    plt.xlabel("Altitude (km)"); plt.ylabel("max loading (MPa)")
    ax2 = plt.gca().twinx(); ax2.plot(res.alt / 1e3, res.v / 1e3, color="firebrick", alpha=.8); ax2.set_ylabel("speed (km/s)", color="firebrick")
    save(fig)
    fig = page("Fragmentation history (marker size ~ parent mass)")
    sp = [e for e in res.events if e["kind"] == "split"]
    if sp:
        plt.scatter([e["t"] - p.tstart for e in sp], [e["alt"] / 1e3 for e in sp],
                    s=[10 + 60 * (e["mass"] / s["m0"]) ** 0.4 for e in sp], c=C_BRK, alpha=.6)
    plt.xlabel("Time (s)"); plt.ylabel("Breakup altitude (km)")
    save(fig)
    pdf.close()


# ----------------------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description=f"bolum v{__version__} bolide model")
    ap.add_argument("event", nargs="?", help="event folder name containing bolum_in.dat")
    ap.add_argument("--base", default=os.environ.get("BOLUM_BASE", DEFAULT_BASE), help="events base directory")
    ap.add_argument("--no-plots", action="store_true")
    ap.add_argument("--quiet", action="store_true", help="print only a one-line summary")
    ap.add_argument("--sat-W", type=float, default=None, help="sensor saturation level in W (draws it on the plots)")
    a = ap.parse_args(argv)
    event = a.event or input("Enter event name: ")
    outdir = os.path.join(a.base, event)
    p = rdinput(os.path.join(outdir, "bolum_in.dat"))
    res = run_bolide(p)
    write_outputs(outdir, res)
    if not a.no_plots:
        plot_outputs(outdir, res, sat_W=a.sat_W, title=f"Event: {event}")
    s = res.summary
    if a.quiet:
        print(f"bolum v{__version__}: {s['steps']} steps in {s['wall_s']:.2f} s | E_dep {s['E_dep_kt']:.3f} kt of KE0 {s['KE0_kt']:.3f} kt | "
              f"E_rad(bol) {s['E_rad_bol_kt']:.4f} kt | peak L_Si {s['peak_L_si_W']:.3e} W at {s['alt_peak_L_km']:.1f} km | closure {s['energy_closure']:.1e}")
    else:
        print(format_report(res))
        files = ["init_params.dat", "bolide_out.dat", "power_vs_time.dat", "power_vs_alt.dat", "dEdh.dat", "events.dat", "summary.dat", "run_report.txt"]
        if not a.no_plots:
            files.append("bolide_plots.pdf")
        print(f" outputs written to {outdir}/ : " + ", ".join(files))


if __name__ == "__main__":
    main()
