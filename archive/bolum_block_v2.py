#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
bolum_block_v2.py  (version 2.3)  --  reduced-order bolide entry / fragmentation / luminosity model

VERSION 2.3 ADDS LIGHT-CURVE FEATURE ANALYSIS (when and why the light curve changes character).  THE PHYSICS AND ALL
PREVIOUS NUMBERS ARE IDENTICAL TO 2.2 (verified bit-for-bit); new output only.
VERSION 2.2 WAS A STRUCTURAL REFACTOR OF 2.1: THE PHYSICS AND THE NUMBERS ARE IDENTICAL.
  The big time-stepping loop has been broken into small, named, commented functions so that each physical
  process can be read (and later modified) on its own.  Every function below says WHAT it computes, WHY, and
  which equation / reference it comes from.  Nothing about the model itself changed (verified bit-for-bit
  against 2.1 on 12 test cases; see the regression note at the bottom of the diagnostics script).

HOW THE MODEL WORKS  (read this first)
--------------------------------------
The bolide is represented as a set of "groups".  A group is a lump of material that moves as one body: it has a
mass, a position and velocity vector (Earth-centred x,z in metres), a strength, and an area-growth factor `f`
(f = 1 means a compact sphere; f > 1 means it has been flattened/spread into a "pancake").  At the start there is
ONE group: the intact body.  Breakups replace a group by several child groups plus a "dust" group.

Each time step does the following, in this order (this is exactly the body of the loop in run_bolide):

   1. LOCAL ENVIRONMENT   altitude, speed, air density, and the stagnation-pressure "loading" rho*v^2 of each group
   2. GEOMETRY            projected (cross-section) area of each group, including pancake growth / dust area
   3. FORCES & HEATING    drag deceleration, ablation (mass-loss) rate
   4. YIELD CHECK         a group whose loading >= its strength "yields" and starts to spread sideways (pancake)
   5. TIME STEP           choose dt small enough to resolve mass loss, deceleration and spreading
   6. ADVANCE             update mass, velocity, position, spread factor; book-keep the energy given to the air
   7. RADIATION           turn the deposited energy into light (Si-band power) with an efficiency eta
   8. RECORD              accumulate energy ledgers, dE/dh bins, light-curve samples
   9. RETIRE              remove groups that burned up, hit the ground, or slowed below luminous speeds
  10. SPLIT               groups that have spread far enough break into children + dust (progressive fragmentation)

ENERGY BOOK-KEEPING (why the ledger closes exactly)
  Each step the kinetic energy of a group changes by:  - E_drag (work done on the air by drag)
                                                       + W_grav (work done by gravity)
                                                       - E_abl  (kinetic energy carried away by ablated mass)
  E_dep = E_drag + E_abl is the energy deposited in the atmosphere; a fraction eta of it is radiated.

PROVENANCE KEY:  [lit] taken from the literature   [adapt] literature form, adapted   [cal] value calibrated here
                 [new] construct introduced here (no direct published equivalent)

USAGE
  python3 bolum_block_v2.py EVENT [--base DIR] [--no-plots] [--quiet] [--sat-W WATTS]

INPUT DECK: unchanged format (first 2 lines skipped, 'key = value  # comment').  Legacy keys rho_tau_scale,
flare_duration, n_frag_init, q_factor are accepted and ignored.  strength_secondary/tertiary/quaternary are optional
overrides of the child strength at generations 1/2/3 (otherwise Weibull scaling).  Strengths are stagnation
pressures rho*v^2 (the published convention); multiply strengths from the ORIGINAL bolum_block.py (0.5*rho*v^2) by 2.

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
        Icarus 309, 25-44.   Canonical heat-transfer coefficient 0.1; radiative heating dominates for large bodies.
  [R11] Weibull (1951), A statistical distribution function of wide applicability, J. Appl. Mech. 18, 293-297.
  [R12] Brown, Spalding, ReVelle, Tagliaferri & Worden (2002), The flux of small near-Earth objects colliding with
        the Earth, Nature 420, 294-296.   E_total = 8.2508 E_rad^0.885 (kt); integral luminous efficiency rises with
        energy; 6000 K blackbody convention for the bolometric conversion (verify).
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
  [R22] Jenniskens, Pilorz, Robertson & Stern (2026), Bolide light curve systematics from 75 recovered meteorites,
        Meteorit. Planet. Sci. 61(9), 2351-2392, doi:10.1111/maps.70203.  Seven-phase description of bolide light curves
        vs altitude (used only to LABEL features in the analysis; the model does not implement their melting, spin or
        plasma-instability physics).
  [T1]  Standard textbook physics: Planck's law (band fraction); Rayleigh pitot / Newtonian hypersonic stagnation
        pressure (p_stag = 0.92 rho v^2 for gamma = 1.4, -> rho v^2 for strong shocks).
"""
__version__ = "2.3"

import os
import sys
import time
import math
import argparse
import textwrap
from dataclasses import dataclass, fields, field
import numpy as np

# ========================================================================================
# 1. CONSTANTS
# ========================================================================================
PI = math.pi
RE = 6.371e6                 # mean Earth radius (m)
GM = 3.986004418e14          # Earth gravitational parameter (m^3/s^2)
G0 = 9.80665                 # standard gravity, used only to express decelerations in "g"
JKT = 4.184e12               # joules per kiloton of TNT
SIGMA_SB = 5.670374419e-8    # Stefan-Boltzmann constant
H_PL, C_LIGHT, K_B = 6.62607015e-34, 2.99792458e8, 1.380649e-23   # Planck, speed of light, Boltzmann
DEFAULT_BASE = "/Users/julie/Desktop/Projects/bolum/Events/"

_trapz = getattr(np, "trapezoid", None) or np.trapz   # numpy 1.x / 2.x compatibility


# ========================================================================================
# 2. ATMOSPHERE  [lit R19]
# ========================================================================================
class Atmosphere:
    """Air density as a function of geometric altitude: US Standard Atmosphere 1976.

    Below 86 km the standard's analytic layer equations (constant temperature or constant lapse rate, hydrostatic
    balance) are evaluated on a 100 m grid.  Between 86 and 120 km a tabulated density is used.  The result is
    stored as ln(rho) on a 100 m grid and linearly interpolated, which is fast and smooth enough for this model.
    Above 120 km the density is extrapolated exponentially using the local scale height."""

    _UPPER_KM = np.array([86., 90., 95., 100., 105., 110., 115., 120.])
    _UPPER_RHO = np.array([6.958e-6, 3.416e-6, 1.393e-6, 5.604e-7, 2.176e-7, 9.708e-8, 4.289e-8, 2.222e-8])

    def __init__(self):
        g0, Rg, r0 = 9.80665, 287.053, 6356766.0            # gravity, specific gas constant of air, standard Earth radius
        Hb = np.array([0, 11, 20, 32, 47, 51, 71, 84.852]) * 1e3          # geopotential altitudes of layer bases (m)
        Lb = np.array([-6.5, 0.0, 1.0, 2.8, 0.0, -2.8, -2.0]) * 1e-3     # temperature lapse rate in each layer (K/m)
        # --- temperature and pressure at the base of every layer (hydrostatic integration, layer by layer) ---
        Tb, Pb = [288.15], [101325.0]
        for i in range(7):
            dH = Hb[i + 1] - Hb[i]
            if abs(Lb[i]) < 1e-12:                           # isothermal layer
                Tb.append(Tb[-1])
                Pb.append(Pb[-1] * math.exp(-g0 * dH / (Rg * Tb[-1])))
            else:                                            # constant-lapse-rate layer
                T1 = Tb[-1] + Lb[i] * dH
                Pb.append(Pb[-1] * (T1 / Tb[-1]) ** (-g0 / (Rg * Lb[i])))
                Tb.append(T1)
        # --- temperature, pressure and density on a fine altitude grid up to 86 km ---
        h = np.arange(0.0, 86.0001e3, 100.0)
        Hg = r0 * h / (r0 + h)                               # geometric -> geopotential altitude
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
        rho_low = P / (Rg * T)                               # ideal gas law
        # --- splice the analytic part (<= 86 km) with the table (86-120 km), store ln(rho) on a 100 m grid ---
        hk = np.concatenate([h / 1e3, self._UPPER_KM[1:]])
        lr = np.concatenate([np.log(rho_low), np.log(self._UPPER_RHO[1:])])
        self._hk = np.arange(0.0, 120.0001, 0.1)
        self._lr = np.interp(self._hk, hk, lr)
        self._scale_h = 1e3 * (self._hk[-1] - self._hk[-51]) / (self._lr[-51] - self._lr[-1])   # scale height near 120 km (m)

    def rho(self, alt_m):
        """Air density (kg/m^3) at altitude alt_m (m); works on scalars or arrays."""
        a = np.asarray(alt_m, dtype=float) / 1e3
        lr = np.interp(a, self._hk, self._lr)
        hi = a > 120.0
        if np.any(hi):
            lr = np.where(hi, self._lr[-1] - (a - 120.0) * 1e3 / self._scale_h, lr)
        return np.exp(lr)


ATM = Atmosphere()


# ========================================================================================
# 3. LUMINOUS EFFICIENCY AND SENSOR-BAND CONVERSION
# ========================================================================================
def tau_fit(v, m):
    """Empirical luminous efficiency tau(v, m) as a FRACTION (the published fit is in percent).  [lit R2, Eqs. 3-4]

        ln tau = 0.567 - 10.307 ln v + 9.781 ln^2 v - 3.0414 ln^3 v + 0.3213 ln^4 v + 0.347 tanh(0.38 ln m)   v < 25.372 km/s
        ln tau = -1.4286 + ln v + 0.347 tanh(0.38 ln m)                                                         v > 25.372 km/s
        (v in km/s, m in kg).  Velocity form: Pecina & Ceplecha (1983); mass form adapted from ReVelle & Ceplecha (2001).

    The fit gives the fraction of the energy lost by a METEOR-SIZED fragment that comes out as light.
    [new] The fit is only meaningful for meteor-like speeds: v is clamped to >= 2 km/s inside the polynomial and
    tau is tapered smoothly to zero between 2 and 4 km/s (the raw polynomial turns upward below ~2 km/s)."""
    v = np.asarray(v, dtype=float)
    vk = np.maximum(v, 2000.0) / 1000.0                      # speed in km/s, clamped at the low end
    lv = np.log(vk)
    th = 0.347 * np.tanh(0.38 * np.log(np.maximum(m, 1e-10)))   # mass-dependent term
    lt_lo = 0.567 - 10.307 * lv + 9.781 * lv ** 2 - 3.0414 * lv ** 3 + 0.3213 * lv ** 4 + th   # v < 25.372 km/s branch
    lt_hi = -1.4286 + lv + th                                # v > 25.372 km/s branch
    tau = np.exp(np.where(vk < 25.372, lt_lo, lt_hi)) / 100.0
    taper = np.clip((v - 2000.0) / 2000.0, 0.0, 1.0) ** 2    # 0 at 2 km/s, 1 at 4 km/s  [new]
    return tau * taper


def band_fraction(T, lo_um=0.4, hi_um=1.1, lam_um=None, resp=None):
    """In-band fraction of blackbody exitance at temperature T [K].  [T1: Planck's law]

        f(T) = integral_band  R(lam) B_lam(T) dlam  /  (sigma T^4 / pi),     B_lam = 2hc^2 / (lam^5 (exp(hc/lam k T) - 1))

    This is the fraction of a blackbody's total radiated power that falls inside the sensor band, so
    L_si = f * L_bol  and  L_bol = L_si / f.
    Default response R = 1 inside [lo_um, hi_um] (broadband silicon, 0.4-1.1 um) and 0 outside.  Pass lam_um/resp
    (peak-normalised relative spectral response) to use a measured sensor curve instead."""
    T = np.atleast_1d(np.asarray(T, dtype=float))
    if resp is None:                                         # top-hat response: integrate Planck over [lo, hi]
        lam = np.linspace(lo_um, hi_um, int(np.clip((hi_um - lo_um) / 0.0015, 600, 200000))) * 1e-6
        w = np.ones_like(lam)
    else:                                                    # user-supplied sensor response curve
        lam = np.asarray(lam_um, dtype=float) * 1e-6
        w = np.asarray(resp, dtype=float)
    x = H_PL * C_LIGHT / (lam[None, :] * K_B * T[:, None])
    B = 2.0 * H_PL * C_LIGHT ** 2 / lam[None, :] ** 5 / np.expm1(x)      # Planck spectral radiance
    return _trapz(B * w[None, :], lam, axis=1) / (SIGMA_SB * T ** 4 / PI)


def bolometric_from_si(L_si, T_rad=6000.0, lo_um=0.4, hi_um=1.1):
    """L_bol = L_si / f_si(T_rad).  Assumes the radiator is a blackbody at T_rad."""
    return np.asarray(L_si) / float(band_fraction(T_rad, lo_um, hi_um)[0])


def si_from_bolometric(L_bol, T_rad=6000.0, lo_um=0.4, hi_um=1.1):
    """L_si = L_bol * f_si(T_rad)."""
    return np.asarray(L_bol) * float(band_fraction(T_rad, lo_um, hi_um)[0])


# ========================================================================================
# 4. PARAMETERS AND INPUT DECK
# ========================================================================================
@dataclass
class BolumParams:
    """Every tunable number in the model, with its default and provenance."""
    # --- required: the event ---
    diameter: float = 1.0            # m
    velocity: float = 20.0e3         # m/s (atmosphere-relative)
    theta_deg: float = 45.0          # entry angle to local horizontal at zstart
    init_strength: float = 2.0e6     # Pa: stagnation pressure rho*v^2 at which the intact body yields [lit R2,R5-R8]
    density: float = 3300.0          # bulk density, kg/m^3
    # --- material / aerodynamics ---
    porosity: float = 0.1            # only used if init_strength <= 0
    Cd: float = 1.0                  # drag coefficient [lit R1,R7]
    # --- ablation (step 2 of the original plan) ---
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
        """Initial strength of the intact body (Pa)."""
        if self.init_strength > 0:
            return self.init_strength
        return 2.0e5 * (1.0 - self.porosity) * (self.density / 3000.0)


_LEGACY_IGNORED = {"rho_tau_scale", "flare_duration", "n_frag_init", "q_factor"}


def rdinput(inputpath):
    """Read a bolum_in.dat deck (first two lines skipped; 'key = value  # comment') and return BolumParams."""
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


# ========================================================================================
# 5. RESULT CONTAINER
# ========================================================================================
@dataclass
class BolumResult:
    """Everything a run produces: time series, dE/dh profiles, breakup events, and a summary dictionary."""
    params: BolumParams
    t: np.ndarray
    L_si: np.ndarray           # Si-band luminosity (W)
    L_bol: np.ndarray          # bolometric luminosity (W) = L_si / f_si
    P_dep: np.ndarray          # power deposited in the atmosphere (W)
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
    P_drag: np.ndarray = None      # power given to the air by drag work (W), time series
    P_abl: np.ndarray = None       # power carried off as kinetic energy of ablated mass (W), time series
    L_si_dust: np.ndarray = None   # part of L_si emitted by dust groups (W)
    features: list = field(default_factory=list)   # light-curve features (see analyze_light_curve)


# ========================================================================================
# 6. FRAGMENT-GROUP STATE
# ========================================================================================
@dataclass
class GroupState:
    """All fragment groups, stored as arrays of length max_groups ("slots").

    A slot is in use when alive[i] is True.  Dead slots are recycled when a group splits.
        m        mass (kg)
        x, z     position in Earth-centred Cartesian coordinates (m); altitude = hypot(x, z) - RE
        vx, vz   velocity components (m/s)
        f        linear spread factor: projected radius = (solid-sphere radius) * f   (1 = compact, f_max = full pancake)
        fb       value of f at birth; the group splits when f >= fb * f_split
        S        strength: loading rho*v^2 (Pa) at which the group yields
        gen      generation (0 = original body, 1 = its children, ...)
        alive    slot in use
        yielded  group has already been recorded as yielding (used only to log the first yield once)
        dust     True for dust groups (fine particles with a large area-to-mass ratio; they never split)"""
    m: np.ndarray
    x: np.ndarray
    z: np.ndarray
    vx: np.ndarray
    vz: np.ndarray
    f: np.ndarray
    fb: np.ndarray
    S: np.ndarray
    gen: np.ndarray
    alive: np.ndarray
    yielded: np.ndarray
    dust: np.ndarray


def make_initial_state(p: BolumParams, m0: float) -> GroupState:
    """Create the empty slot arrays and put the intact body in slot 0.

    The body starts at altitude zstart above the surface at x = 0, moving at speed `velocity` in a straight line
    that makes angle theta with the local horizontal (down and to +x).  Earth-centred coordinates mean the local
    vertical rotates automatically as the body travels, so no flat-Earth approximation is needed."""
    NG = int(p.max_groups)
    st = GroupState(m=np.zeros(NG), x=np.zeros(NG), z=np.zeros(NG), vx=np.zeros(NG), vz=np.zeros(NG),
                    f=np.ones(NG), fb=np.ones(NG), S=np.zeros(NG), gen=np.zeros(NG, dtype=int),
                    alive=np.zeros(NG, dtype=bool), yielded=np.zeros(NG, dtype=bool), dust=np.zeros(NG, dtype=bool))
    th = math.radians(p.theta_deg)
    st.m[0], st.z[0] = m0, RE + p.zstart
    st.vx[0], st.vz[0] = p.velocity * math.cos(th), -p.velocity * math.sin(th)
    st.S[0] = p.strength0()
    st.alive[0] = True
    return st


# ========================================================================================
# 7. PHYSICS BUILDING BLOCKS  (each is a small pure function acting on arrays of groups)
# ========================================================================================
def local_environment(x, z, vx, vz):
    """Altitude, speed, air density and stagnation-pressure loading for each group.

    The loading is p = rho * v^2, the (hypersonic) stagnation pressure that squeezes the front of the body
    [lit R5-R8, R2].  (1/2 rho v^2 would be the incompressible dynamic pressure, which understates the load ~2x.)
    Returns r (distance from Earth's centre), alt, v, rho, q."""
    r = np.hypot(x, z)
    alt = r - RE
    v = np.hypot(vx, vz)
    rho = ATM.rho(alt)
    q = rho * v * v
    return r, alt, v, rho, q


def solid_radius(m, rho_body):
    """Radius of a compact sphere of mass m and bulk density rho_body (m)."""
    return (3.0 * m / (4.0 * PI * rho_body)) ** (1.0 / 3.0)


def projected_area(m, f, is_dust, rho_body, dust_radius):
    """Cross-sectional area each group presents to the airflow (m^2).  Also returns the solid radius.

    * Ordinary group:  A = pi (r_solid * f)^2.  f = 1 for a compact sphere; f grows as the group spreads into a
      pancake [lit R5, R6, R7].
    * Dust group:      many small spheres of radius `dust_radius` (default 2 cm) share the group's mass; their total
      cross-section is A = 3 m / (4 rho r_dust).  This large area-to-mass ratio makes dust decelerate and ablate
      quickly, which deposits energy and produces flares.  [new; cal]"""
    rs = solid_radius(m, rho_body)
    A = np.where(is_dust, 3.0 * m / (4.0 * rho_body * dust_radius), PI * rs * rs * f * f)
    return A, rs


def drag_deceleration(rho, v, A, m, Cd):
    """Drag force F = 1/2 Cd rho A v^2 and the resulting deceleration a = F/m.  [lit R1, R7]  Returns (F, a)."""
    F = 0.5 * Cd * rho * v * v * A
    return F, F / m


def ablation_mass_rate(rho, v, A, Lambda, q_loss, Q):
    """Mass-loss rate dm/dt (kg/s, positive) from aerodynamic heating.

    The free-stream kinetic-energy flux is (1/2) rho v^3.  A fraction Lambda (<= 1; canonical 0.1 [R10]) of it heats the
    surface:  q_heat = Lambda * (1/2) rho v^3.  Surface radiation at the ablation temperature carries away q_loss =
    emissivity*sigma*T_ablate^4, so only the NET flux beyond that ablates material (ablation onset idea [adapt R20]):
        dm/dt = max(q_heat - q_loss, 0) * A / Q,       Q = effective heat of ablation (J/kg)  [lit R1]
    The flux is applied to the PROJECTED area A (not the full sphere), so it can never exceed the kinetic-energy flux
    the airstream actually carries through that area."""
    qh = Lambda * 0.5 * rho * v ** 3
    return np.maximum(qh - q_loss, 0.0) * A / Q


def pancake_spreading_rate(yielding, v, rho, rho_body, rs, c_spread):
    """Rate of change of the spread factor, df/dt (1/s).  [adapt R5, R6, R21; cal c_spread]

    A group whose loading exceeds its strength ("yielding") is squeezed and spreads laterally, like a pancake,
    at a speed  dr/dt = c * v * sqrt(rho_air / rho_body).  Dividing by the solid radius turns this into growth of
    the linear spread factor f.  Groups that are not yielding do not spread."""
    return np.where(yielding, c_spread * v * np.sqrt(rho / rho_body) / rs, 0.0)


def gravity_acceleration(x, z, r):
    """Gravitational acceleration components (m/s^2) pointing toward Earth's centre."""
    gr = GM / r ** 2
    return -gr * x / r, -gr * z / r


def choose_time_step(p, mdot, m, is_dust, a_drag, v, dfdt, f):
    """Pick the time step (s) so that no group changes too much in one step.

    The step is the smallest of:  the user maximum p.dt;  3 % of the shortest mass-loss time (30 % for dust, which
    burns fast by design);  2 % of the shortest deceleration time v/a;  5 % of the shortest spreading time f/(df/dt).
    It is never smaller than p.dt_min."""
    with np.errstate(divide="ignore", invalid="ignore"):
        d1 = np.min(np.where(mdot > 0, np.where(is_dust, 0.3, 0.03) * m / mdot, np.inf))
        d2 = np.min(np.where(a_drag > 0, 0.02 * v / a_drag, np.inf))
        d3 = np.min(np.where(dfdt > 0, 0.05 * f / dfdt, np.inf))
    return max(min(p.dt, d1, d2, d3), p.dt_min)


@dataclass
class StepUpdate:
    """New state of the active groups after one step, plus the energy terms for the step."""
    mnew: np.ndarray
    xn: np.ndarray
    zn: np.ndarray
    vxn: np.ndarray
    vzn: np.ndarray
    fn: np.ndarray
    vn: np.ndarray
    altn: np.ndarray
    E_drag: np.ndarray   # work done on the air by drag during the step (J)
    E_abl: np.ndarray    # kinetic energy carried off by ablated mass (J)
    E_dep: np.ndarray    # = E_drag + E_abl: energy deposited in the atmosphere (J)
    W_g: np.ndarray      # work done on the groups by gravity (J)


def advance_groups(p, dt, m, x, z, vx, vz, f, r, v, a_drag, mdot, dfdt, m_dead):
    """Advance every active group by one step of length dt and book-keep the energy exchanged.  Returns StepUpdate.

    Mass:      exponential loss  m_new = m * exp(-mdot*dt/m)  (never negative, stable for fast ablation).  A group whose
               mass falls below m_dead is treated as completely ablated (all of its mass is removed).
    Velocity:  drag along the direction of motion plus gravity (explicit Euler).
    Position:  advanced with the MEAN of the old and new velocity (second-order accurate).
    Energy:    using the mean velocity makes the kinetic-energy bookkeeping exact:
                   KE_new - KE_old  =  -E_drag + W_grav - E_abl
               with E_drag = work against drag, E_abl = 1/2 dm v^2 carried away by ablated mass, W_grav = gravity work.
    Spread:    f grows by df/dt*dt, limited to the pancake maximum f_max."""
    gxv, gzv = gravity_acceleration(x, z, r)
    dm = m * (-np.expm1(-mdot * dt / m))                     # mass ablated this step (exponential law)
    mnew = m - dm
    dead = mnew < m_dead                                     # essentially nothing left: remove the whole group
    dm = np.where(dead, m, dm)
    mnew = np.where(dead, 0.0, mnew)
    dvx_d, dvz_d = -a_drag * vx / v * dt, -a_drag * vz / v * dt        # velocity change from drag (opposes motion)
    vxn = vx + dvx_d + gxv * dt
    vzn = vz + dvz_d + gzv * dt
    vbx, vbz = 0.5 * (vx + vxn), 0.5 * (vz + vzn)            # mean velocity over the step
    E_drag = -mnew * (vbx * dvx_d + vbz * dvz_d)             # work done against drag = energy given to the air
    E_abl = 0.5 * dm * v * v                                 # kinetic energy of the mass that was ablated
    W_g = mnew * (vbx * gxv + vbz * gzv) * dt                # work done by gravity
    E_dep = E_drag + E_abl
    xn, zn = x + vbx * dt, z + vbz * dt
    fn = np.minimum(f + dfdt * dt, p.f_max)
    vn = np.hypot(vxn, vzn)
    altn = np.hypot(xn, zn) - RE
    return StepUpdate(mnew, xn, zn, vxn, vzn, fn, vn, altn, E_drag, E_abl, E_dep, W_g)


def radiative_efficiency(p, f_si, basis_si, v, m, A):
    """Fraction eta of the deposited energy that is emitted as Si-band light, for each group.  [new blend; R2, R12, R13]

    Two limits are blended:
      * THIN  (small, meteor-like fragments):  eta_thin = tau(v, m), the published empirical fit [lit R2].  If the fit is
        interpreted as bolometric (tau_basis = 'bolometric') it is multiplied by f_si to give a Si-band value.
      * THICK (large, radiation-dominated clouds): eta_thick = eta_thick_bol * f_si, a bolometric efficiency that does
        not depend on speed (but is tapered to zero between 4 and 8 km/s so light ends when the body slows).
    The blend weight is  w = r/(r + r_half), r = projected radius of the group, so big clouds are 'thick' and small
    fragments are 'thin'.  [cal: eta_thick_bol = 0.20, r_half = 6 m, fitted to Chelyabinsk and TC3.]
    The result is capped at eta_max (bolometric) so the model can never radiate more than a sane fraction."""
    tau = tau_fit(v, m)
    r_proj = np.sqrt(A / PI)
    w_thick = r_proj / (r_proj + p.lum_r_half)
    eta_thin = tau if basis_si else tau * f_si
    eta_thick = p.eta_thick_bol * f_si * np.clip((v - 4.0e3) / 4.0e3, 0.0, 1.0)
    return np.minimum((1.0 - w_thick) * eta_thin + w_thick * eta_thick, p.eta_max * f_si)


# ========================================================================================
# 8. FRAGMENTATION BOOK-KEEPING
# ========================================================================================
def log_new_yields(st, ia, yl, t, alt, q, Si, events):
    """Record (once per group) the moment a group first yields, i.e. its loading first reaches its strength."""
    newly = yl & ~st.yielded[ia]
    if newly.any():
        for j in np.flatnonzero(newly):
            events.append(dict(kind="yield", t=t, alt=alt[j], q=q[j], mass=st.m[ia][j], strength=Si[j], gen=int(st.gen[ia[j]])))
        st.yielded[ia[newly]] = True


def child_mass_shares(p):
    """Mass fractions of the children of one breakup: share of child k (k = 1..n) is proportional to k^-frag_index,
    normalised to 1.  With the default index 1 and 4 children the shares are 0.48, 0.24, 0.16, 0.12.  [new]"""
    w = np.arange(1, int(p.n_fragments) + 1, dtype=float) ** (-p.frag_index)
    return w / w.sum()


def spawn_children(st, i, p, w_share, legacy_S, t_now, alt_now, q_now, events):
    """Break group i into n children plus a dust group (progressive fragmentation).  [new; concept R4, R8, R16]

    * Mass is conserved:  dust_frac of the parent goes to the dust group; the rest is shared among the children by
      `w_share`.
    * Children keep the parent's position and velocity (they separate later through differential drag).
    * STRENGTH: a child of mass m_k gets  S_k = S_parent * (m_parent/m_k)^alpha  (Weibull size scaling: smaller pieces
      are stronger because they contain fewer flaws) [lit form R8, R11].  If the deck gave strength_secondary/
      tertiary/quaternary these override generations 1, 2, 3 (never below the parent's strength).
    * AREA continuity: the spread factor of each child is chosen so that the children's total projected area is about
      equal to the parent's, so the drag does not jump at the instant of breakup.
    * Returns False (and does nothing) if there are not enough free slots or the smallest child would be lighter than
      p.m_min_frag; the caller then leaves the group as a pancake."""
    free = np.flatnonzero(~st.alive)
    n_child = int(p.n_fragments)
    need = n_child if p.dust_frac > 0 else n_child - 1
    mp, fp, Sp, gp = st.m[i], st.f[i], st.S[i], st.gen[i]     # parent mass, spread, strength, generation
    m_solid = mp * (1.0 - p.dust_frac)
    if len(free) < need or m_solid * w_share[-1] < p.m_min_frag:
        return False
    slots = [i] + list(free[: n_child - 1])                   # child 0 reuses the parent's slot
    xi_, zi_, vxi_, vzi_ = st.x[i], st.z[i], st.vx[i], st.vz[i]
    for k, s_ in enumerate(slots):
        mk = m_solid * w_share[k]
        fk = max(1.0, fp * (mk / mp) ** (1.0 / 6.0))          # area-continuity spread factor
        gnew = gp + 1
        if gnew <= 3 and legacy_S[gnew - 1] is not None:      # deck override of the stage strength
            Sk = max(Sp, float(legacy_S[gnew - 1]))
        else:                                                 # Weibull scaling
            Sk = Sp * (mp / mk) ** p.weibull_alpha
        st.m[s_], st.x[s_], st.z[s_], st.vx[s_], st.vz[s_] = mk, xi_, zi_, vxi_, vzi_
        st.f[s_], st.fb[s_], st.S[s_], st.gen[s_] = fk, fk, Sk, gnew
        st.alive[s_] = True
        st.yielded[s_] = False
        st.dust[s_] = False
    if p.dust_frac > 0:                                       # the dust group takes the next free slot
        d_ = free[n_child - 1]
        st.m[d_], st.x[d_], st.z[d_], st.vx[d_], st.vz[d_] = mp * p.dust_frac, xi_, zi_, vxi_, vzi_
        st.f[d_], st.fb[d_], st.S[d_], st.gen[d_] = 1.0, np.inf, np.inf, gp + 1   # infinite strength: never yields/splits
        st.alive[d_] = True
        st.yielded[d_] = False
        st.dust[d_] = True
    events.append(dict(kind="split", t=t_now, alt=alt_now, q=q_now, mass=mp, strength=Sp, gen=gp))
    return True


def split_spread_groups(st, p, w_share, legacy_S, t, events):
    """Split every group that has spread far enough (f >= fb * f_split) while yielding.

    Returns the number of splits that were REFUSED (no free slots / fragments too small).  A refused group is
    marked unsplittable (fb = inf): it simply stays a pancake limited by f_max."""
    refused = 0
    cand = np.flatnonzero(st.alive & st.yielded & (st.f >= st.fb * p.f_split))
    for k in cand:
        altk = float(np.hypot(st.x[k], st.z[k]) - RE)
        vk = float(np.hypot(st.vx[k], st.vz[k]))
        qk = float(ATM.rho(altk)) * vk * vk
        if not spawn_children(st, k, p, w_share, legacy_S, t, altk, qk, events):
            st.fb[k] = np.inf
            refused += 1
    return refused


def retire_groups(st, ia, up: StepUpdate, du, p, t, rec):
    """Remove groups that have left the simulation this step, updating the recorder's survivor ledgers.

    A group is retired if it has been completely ablated (mass 0 - nothing to keep), hit the ground (altitude <= 0),
    or slowed below p.v_dark (it is no longer luminous; it continues in 'dark flight').  Retired groups that still
    have mass are logged in rec.retired with their mass, speed, altitude and kinetic energy (these are the survivors),
    and their kinetic energy and mass are added to rec.KE_retired / rec.M_retired (needed for the energy closure)."""
    gone = (up.mnew <= 0)
    st.alive[ia[gone]] = False
    out = (~gone) & ((up.altn <= 0.0) | (up.vn < p.v_dark))
    for j in np.flatnonzero(out):
        ke = 0.5 * up.mnew[j] * up.vn[j] ** 2
        rec.KE_retired += ke
        rec.M_retired += up.mnew[j]
        rec.retired.append(dict(mass=up.mnew[j], v=up.vn[j], alt=up.altn[j], t=t, ke=ke, ground=bool(up.altn[j] <= 0), dust=bool(du[j])))
        st.alive[ia[j]] = False


# ========================================================================================
# 9. RUN RECORDER  (ledgers, light-curve samples, dE/dh bins)
# ========================================================================================
class RunRecorder:
    """Collects everything we want to report while the loop runs."""

    def __init__(self, p: BolumParams):
        self.bin_w = 500.0                                   # altitude-bin width for dE/dh (m)
        self.nb = int(150e3 / self.bin_w)
        self.dep_bins = np.zeros(self.nb)                    # energy deposited in each altitude bin (J)
        self.rad_bins = np.zeros(self.nb)                    # Si-band energy radiated in each altitude bin (J)
        self.T, self.L_si, self.P, self.ALT, self.V, self.M, self.N, self.Q = [], [], [], [], [], [], [], []
        self.P_drag, self.P_abl, self.L_dust = [], [], []   # which channel deposits the energy / how much light is from dust
        self.t_abl_onset, self.alt_abl_onset = None, float('nan')   # first time any group ablates
        self.events = []                                     # yield and split events
        self.retired = []                                    # groups that left the simulation with mass
        self.E_dep_tot = self.E_rad_tot = self.W_g_tot = 0.0  # running energy totals
        self.KE_retired = self.M_retired = 0.0
        self.split_refused = 0
        self.n_groups_max = 1
        self.a_max = 0.0                                     # largest deceleration of any non-dust group (m/s^2)
        self.q_peak, self.alt_q_peak = 0.0, float("nan")     # largest loading and where it occurred
        self.nsteps = 0

    def record_step(self, t, dt, alt, v, q, a_drag, du, m_before, up: StepUpdate, E_rad_si, n_alive_before, mdot):
        """Log one step.  `t` is the time at the END of the step."""
        E_dep_s, W_g_s = float(up.E_dep.sum()), float(up.W_g.sum())
        self.E_dep_tot += E_dep_s
        self.W_g_tot += W_g_s
        self.E_rad_tot += float(E_rad_si.sum())
        # dE/dh: add each group's deposited / radiated energy to the altitude bin it occupied at the start of the step
        bi = np.clip((alt / self.bin_w).astype(int), 0, self.nb - 1)
        np.add.at(self.dep_bins, bi, up.E_dep)
        np.add.at(self.rad_bins, bi, E_rad_si)
        jq = int(np.argmax(q))                               # track the most heavily loaded group
        if q[jq] > self.q_peak:
            self.q_peak, self.alt_q_peak = float(q[jq]), float(alt[jq])
        if (~du).any():
            self.a_max = max(self.a_max, float(a_drag[~du].max()))
        # time series (light curve etc.)
        self.T.append(t)
        self.L_si.append(E_rad_si.sum() / dt)
        self.P.append(E_dep_s / dt)
        wts = up.E_dep if E_dep_s > 0 else m_before          # altitude: weighted by deposited energy (or by mass)
        self.ALT.append(float((wts * alt).sum() / wts.sum()))
        self.V.append(float((m_before * v).sum() / m_before.sum()) if m_before.sum() > 0 else 0.0)
        self.M.append(float(up.mnew.sum()))
        self.N.append(int((up.mnew > 0).sum()))
        self.Q.append(float(q.max()))
        self.n_groups_max = max(self.n_groups_max, n_alive_before)
        # --- extra series for the light-curve analysis (do not affect any physics) ---
        self.P_drag.append(float(up.E_drag.sum()) / dt)
        self.P_abl.append(float(up.E_abl.sum()) / dt)
        self.L_dust.append(float(E_rad_si[du].sum()) / dt)
        if self.t_abl_onset is None and mdot.max() > 0:      # ablation begins: remember when and where
            self.t_abl_onset = t - dt
            self.alt_abl_onset = float(alt[int(np.argmax(mdot))])
        self.nsteps += 1


# ========================================================================================
# 10. POST-PROCESSING
# ========================================================================================
def relax_light_curve(L_si, T, tstart, t_cool):
    """Optional first-order radiative relaxation of the light curve (energy conserving).  [new]
    Each sample moves toward the instantaneous power with time constant t_cool, mimicking the cooling time of the
    hot gas.  The total radiated energy is unchanged; flares are smoothed and delayed.  t_cool = 0 disables it."""
    dts = np.diff(np.concatenate([[tstart], T]))
    Lr = np.empty_like(L_si)
    prev = 0.0
    for k in range(len(L_si)):
        prev += (L_si[k] - prev) * (1.0 - math.exp(-dts[k] / t_cool))
        Lr[k] = prev
    return Lr


def deposition_profile_metrics(dep_bins, bin_w):
    """Characterise the dE/dh profile.  Returns (altitudes, dE/dh smoothed, i_peak, i_lo, i_hi, i_half_energy).

    dE/dh is smoothed with a 3-bin boxcar.  i_lo/i_hi bound the contiguous region around the peak where dE/dh is at
    least half of the peak value (the 'half-energy range').  i_half_energy is the altitude bin above which half of the
    total deposited energy was released."""
    nb = len(dep_bins)
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
    cum_top = np.cumsum(dep_bins[::-1])[::-1]                # energy deposited above each bin
    ih = int(np.flatnonzero(cum_top >= 0.5 * dep_bins.sum())[-1]) if dep_bins.sum() > 0 else 0
    return alt_c, dEdh, sm, ipk, lo, hi, ih


def build_summary(p, st, rec: RunRecorder, T_, L_si, L_bol, ALT, V_, KE0, m0, f_si, wall0):
    """Collect headline numbers, energy-closure check and warnings into one dictionary."""
    dts = np.diff(np.concatenate([[p.tstart], T_])) if len(T_) else np.array([])
    E_rad_si_tot = float((L_si * dts).sum()) if len(T_) else 0.0
    # kinetic energy still carried by groups that are alive at the end (only non-zero if tstop was reached)
    ke_end_active = float(sum(0.5 * st.m[k] * (st.vx[k] ** 2 + st.vz[k] ** 2) for k in np.flatnonzero(st.alive)))
    # ENERGY CLOSURE: KE0 + gravity work = deposited + KE left behind (retired groups + still-flying groups)
    closure = (KE0 + rec.W_g_tot - rec.E_dep_tot - rec.KE_retired - ke_end_active) / KE0
    alt_c, dEdh, sm, ipk, lo, hi, ih = deposition_profile_metrics(rec.dep_bins, rec.bin_w)
    ipl = int(np.argmax(L_si)) if len(L_si) else 0
    lit_end = np.flatnonzero(L_si > 0.01 * L_si.max()) if len(L_si) else np.array([], dtype=int)     # last time above 1 % of peak
    lit_start = np.flatnonzero(L_si > 1e-6 * L_si.max()) if len(L_si) else np.array([], dtype=int)   # first time above 1e-6 of peak
    first_yield = next((e for e in rec.events if e["kind"] == "yield"), None)
    retired = rec.retired
    summary = dict(
        KE0_J=KE0, KE0_kt=KE0 / JKT, m0=m0,
        E_dep_J=rec.E_dep_tot, E_dep_kt=rec.E_dep_tot / JKT,
        E_rad_si_J=E_rad_si_tot, E_rad_bol_J=E_rad_si_tot / f_si,
        E_rad_bol_kt=E_rad_si_tot / f_si / JKT,
        band_fraction=f_si,
        eta_si_overall=E_rad_si_tot / rec.E_dep_tot if rec.E_dep_tot > 0 else 0.0,
        eta_bol_overall=E_rad_si_tot / f_si / rec.E_dep_tot if rec.E_dep_tot > 0 else 0.0,
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
        n_split_events=sum(1 for e in rec.events if e["kind"] == "split"),
        split_refused=rec.split_refused, n_groups_max=rec.n_groups_max,
        q_max_MPa=rec.q_peak / 1e6, alt_q_max_km=rec.alt_q_peak / 1e3, max_decel_g=rec.a_max / G0,
        mass_survived_kg=rec.M_retired, mass_survived_frac=rec.M_retired / m0,
        mass_survived_solid_frac=sum(r["mass"] for r in retired if not r["dust"]) / m0,
        mass_survived_dust_frac=sum(r["mass"] for r in retired if r["dust"]) / m0,
        ground_ke_J=float(sum(r["ke"] for r in retired if r["ground"])),
        energy_closure=closure, reached_tstop=bool(st.alive.any()),
        steps=rec.nsteps, wall_s=time.time() - wall0,
        duration_s=float(T_[-1] - p.tstart) if len(T_) else 0.0,
        retired=retired,
    )
    warns = []
    if abs(closure) > 1e-6:
        warns.append(f"energy ledger closure {closure:.1e} exceeds 1e-6")
    if rec.split_refused:
        warns.append(f"{rec.split_refused} breakup(s) refused (group cap max_groups={p.max_groups} or m_min_frag) - raise max_groups?")
    if summary["reached_tstop"]:
        warns.append(f"stopped at tstop={p.tstop:.0f} s with groups still in flight")
    if summary["E_dep_kt"] < 0.5 * summary["KE0_kt"]:
        warns.append("less than half of the initial kinetic energy was deposited (body survives/skips?)")
    summary["warnings"] = warns
    return summary, alt_c, dEdh


# ========================================================================================
# 11. THE DRIVER: run_bolide
# ========================================================================================
def run_bolide(p: BolumParams, verbose=False) -> BolumResult:
    """Simulate one entry event and return a BolumResult.  Reads as the 10-step recipe in the module docstring."""
    wall0 = time.time()
    # ---------- one-time set-up ----------
    rho_b = p.density
    m0 = 4.0 / 3.0 * PI * (0.5 * p.diameter) ** 3 * rho_b        # initial mass of the intact sphere
    w_share = child_mass_shares(p)                                # how a breakup divides the mass among children
    q_loss = p.emissivity * SIGMA_SB * p.T_ablate ** 4            # radiative loss at the ablation temperature (W/m^2)
    m_dead = 1.0e-9 * m0                                          # below this mass a group counts as fully ablated
    f_si = float(band_fraction(p.T_rad, p.sensor_lo_um, p.sensor_hi_um)[0])   # Si-band fraction of the blackbody
    basis_si = p.tau_basis.lower() == "si"
    legacy_S = [p.strength_secondary, p.strength_tertiary, p.strength_quaternary]
    st = make_initial_state(p, m0)
    rec = RunRecorder(p)
    KE0 = 0.5 * m0 * p.velocity ** 2                              # initial kinetic energy
    t = p.tstart

    # ---------- time loop: one pass = one time step for all groups that are alive ----------
    while st.alive.any() and t < p.tstop:
        ia = np.flatnonzero(st.alive)                             # indices of the active groups
        n_alive_before = int(st.alive.sum())
        mi, xi, zi = st.m[ia], st.x[ia], st.z[ia]
        vxi, vzi, fi, Si = st.vx[ia], st.vz[ia], st.f[ia], st.S[ia]
        du = st.dust[ia]

        # 1. local environment: altitude, speed, air density, stagnation-pressure loading
        ri, alt, vi, rho, q = local_environment(xi, zi, vxi, vzi)
        # 2. geometry: projected area (pancake growth / dust)
        A, rs = projected_area(mi, fi, du, rho_b, p.dust_radius)
        # 3. forces and heating: drag deceleration and ablation rate
        F, a_dr = drag_deceleration(rho, vi, A, mi, p.Cd)
        mdot = ablation_mass_rate(rho, vi, A, p.Lambda, q_loss, p.L_ablation)
        # 4. yield check: loading >= strength -> the group yields, logs the event, and spreads sideways
        yl = q >= Si
        log_new_yields(st, ia, yl, t, alt, q, Si, rec.events)
        dfdt = pancake_spreading_rate(yl, vi, rho, rho_b, rs, p.c_spread)
        # 5. time step
        dt = choose_time_step(p, mdot, mi, du, a_dr, vi, dfdt, fi)
        # 6. advance mass, velocity, position, spread; compute the energy given to the atmosphere
        up = advance_groups(p, dt, mi, xi, zi, vxi, vzi, fi, ri, vi, a_dr, mdot, dfdt, m_dead)
        # 7. radiation: Si-band energy emitted this step
        eta = radiative_efficiency(p, f_si, basis_si, vi, mi, A)
        E_rad_si = eta * up.E_dep
        # 8. record ledgers, dE/dh bins, light-curve samples
        t += dt
        rec.record_step(t, dt, alt, vi, q, a_dr, du, mi, up, E_rad_si, n_alive_before, mdot)
        # commit the new state of the active groups
        st.m[ia], st.x[ia], st.z[ia], st.vx[ia], st.vz[ia], st.f[ia] = up.mnew, up.xn, up.zn, up.vxn, up.vzn, up.fn
        # 9. retire groups that burned up, hit the ground, or dropped below luminous speed
        retire_groups(st, ia, up, du, p, t, rec)
        # 10. split groups that have spread enough (progressive fragmentation + dust release)
        rec.split_refused += split_spread_groups(st, p, w_share, legacy_S, t, rec.events)

    # ---------- assemble the outputs ----------
    T_ = np.array(rec.T); L_si = np.array(rec.L_si); P_ = np.array(rec.P); ALT = np.array(rec.ALT)
    if p.t_cool > 0.0 and len(T_) > 1:
        L_si = relax_light_curve(L_si, T_, p.tstart, p.t_cool)
    L_bol = L_si / f_si                                           # bolometric = Si-band / in-band blackbody fraction
    summary, alt_c, dEdh = build_summary(p, st, rec, T_, L_si, L_bol, ALT, np.array(rec.V), KE0, m0, f_si, wall0)
    summary["t_ablation_onset"] = rec.t_abl_onset if rec.t_abl_onset is not None else float("nan")
    summary["alt_ablation_onset_km"] = rec.alt_abl_onset / 1e3
    res = BolumResult(p, T_, L_si, L_bol, P_, ALT, np.array(rec.V), np.array(rec.M), np.array(rec.N),
                      np.array(rec.Q), L_si / np.maximum(P_, 1e-30), rec.events, alt_c, dEdh, rec.rad_bins / rec.bin_w, summary,
                      P_drag=np.array(rec.P_drag), P_abl=np.array(rec.P_abl), L_si_dust=np.array(rec.L_dust))
    res.features, lc_stats = analyze_light_curve(res)          # 15. light-curve feature analysis (new in 2.3)
    summary.update(lc_stats)
    return res



# ========================================================================================
# 11b. LIGHT-CURVE FEATURE ANALYSIS  (new in 2.3; reads the finished run, changes no physics)
# ========================================================================================
def _moving_average(y, n):
    """Boxcar smoothing with edge padding (n samples)."""
    k = np.ones(n) / n
    yp = np.concatenate([np.full(n, y[0]), y, np.full(n, y[-1])])
    return np.convolve(yp, k, mode="same")[n:-n]


def _moving_median(y, n):
    """Running median over n samples (edge padded); used as the 'local trend' when looking for small bumps."""
    n += (n + 1) % 2                                       # make the window odd
    h = n // 2
    yp = np.concatenate([np.full(h, y[0]), y, np.full(h, y[-1])])
    return np.median(np.lib.stride_tricks.sliding_window_view(yp, n), axis=1)


def analyze_light_curve(res, k_flare=3.0, min_rise=1.25, bump_min=0.04):
    """Find and classify the features of the light curve.  Returns (features, stats).

    WHY: a light curve changes character for several different reasons, and only some of them are real physics.
    Each feature gets a `flag`:
        physical        a real process that the model represents with reasonable fidelity
        stylized        a real process, but the SHAPE, SHARPNESS or NUMBER of the feature is set by the model's rules
                        (instant dust release, fixed 4-child cascade, sharp ablation threshold), so do not over-read it
        artifact        a numerical / bookkeeping effect (groups removed when they slow below v_dark)
        low-confidence  region controlled by artifacts or by assumptions outside the model's range
        info            bookkeeping marker (e.g. the arbitrary first-light threshold)

    THE REFERENCE RATE.  If the light simply tracked the air density along the path, ln L would grow at
        s_atm = (descent speed) / (density scale height) = v sin(theta) / H(alt)    (about 1-2 per second).
    A 'flare' is a rise in ln L steeper than k_flare * s_atm (and at least 2 per second) that lifts L by at least
    min_rise (default 25 %).  A 'bump' is a smaller local maximum (>= bump_min above the local trend) that is not part of
    a flare.  Both are matched to breakup events in time (flare: +-0.05 s, bump: +-0.1 s).

    LABELLING TO REAL LIGHT CURVES.  Jenniskens et al. (2026) [R22] describe seven phases in camera light curves of
    recovered meteorite falls.  The model has analogues of only some of them (see notes on each feature)."""
    p, s = res.params, res.summary
    t, L = res.t, res.L_si
    stats = {}
    if len(t) < 20 or not (L.max() > 0):
        return [], stats
    peak = float(L.max())
    ipk = int(np.argmax(L))
    sin0 = math.sin(math.radians(p.theta_deg))
    dts = np.diff(np.concatenate([[p.tstart], t]))
    feats = []

    def add(kind, tag, tt, flag, note, **extra):
        feats.append(dict(kind=kind, tag=tag, t=float(tt), alt_km=float(np.interp(tt, t, res.alt)) / 1e3,
                          L_si=float(np.interp(tt, t, L)), flag=flag, note=note, **extra))

    def ev_near(ta, tb, pad):
        return [e for e in res.events if ta - pad <= e["t"] <= tb + pad]

    # ---- 1. first light (arbitrary visibility threshold) ----
    t_first = s["t_first_light"] if np.isfinite(s["t_first_light"]) else t[0]
    add("first_light", "", t_first, "info", "L_si first exceeds 1e-6 of its peak. This is just a visibility threshold; the real value depends on "
        "your sensor floor and on where the simulation starts (zstart).")

    # ---- 2. ablation onset: where the model's light changes from drag-only to ablation-driven ----
    t_on = s.get("t_ablation_onset", float("nan"))
    if np.isfinite(t_on):
        add("ablation_onset", "A", t_on, "stylized",
            f"Net heating flux first exceeds radiative loss at T_ablate = {p.T_ablate:.0f} K, so ablation switches on (a kink in the slope). "
            f"Real onset is gradual over a few km. BELOW this point the model's light is drag work only (no ablated material): treat "
            f"light above ~{s['alt_ablation_onset_km']:.0f} km as low-confidence. Loosely like the end of Phase 1 / start of Phase 2 in [R22], "
            f"but the model has no spin, volatile (Na/K) emission or melt physics.")
    # ---- 3. first yield and first breakup (pancake phase) ----
    ys = next((e for e in res.events if e["kind"] == "yield"), None)
    sp1 = next((e for e in res.events if e["kind"] == "split"), None)
    if ys is not None:
        add("first_yield", "Y", ys["t"], "physical",
            f"Loading rho v^2 reaches S0 = {ys['strength']/1e6:.3g} MPa. The body starts to flatten (pancake): area, drag and heating all grow, so "
            f"the light rises faster. Analogue of the onset of fragmentation (Phase 5 / H_f in [R22]).")
    if sp1 is not None:
        extra = ""
        if ys is not None:
            ratio = float(np.interp(sp1["t"], t, L)) / max(float(np.interp(ys["t"], t, L)), 1e-30)
            extra = f" Light grew x{ratio:.1f} between first yield and this breakup (pancake spreading)."
        add("first_breakup", "B", sp1["t"], "stylized",
            f"Spread factor reached f_split = {p.f_split:g}: the body breaks into {p.n_fragments} children plus {p.dust_frac*100:.0f} % of its mass as dust. "
            f"The dust has a large area-to-mass ratio, decelerates and ablates within a fraction of a second, which makes a sharp flare. "
            f"The instantaneous release is a model simplification.{extra}")

    # ---- 4/5. flares and bumps from the log-slope of the light curve ----
    flares = []
    t0 = float(max(t_first, t[0])); t1 = float(t[-1])
    du = float(np.clip((t1 - t0) / 3000.0, 1e-3, 5e-3))
    grid = np.arange(t0, t1, du)
    if len(grid) > 50:
        floor = 1e-9 * peak
        yg = np.log(np.maximum(np.interp(grid, t, L), floor))
        nsm = max(3, int(round(0.03 / du)))
        ysm = _moving_average(yg, nsm)
        slope = np.gradient(ysm, du)
        alt_g, v_g = np.interp(grid, t, res.alt), np.interp(grid, t, res.v)
        invH = (np.log(ATM.rho(alt_g - 250.0)) - np.log(ATM.rho(alt_g + 250.0))) / 500.0   # 1 / density scale height
        s_atm = v_g * sin0 * invH
        hot = slope > np.maximum(k_flare * s_atm, 2.0)
        runs, i, n = [], 0, len(hot)
        while i < n:
            if hot[i]:
                j = i
                while j + 1 < n and hot[j + 1]:
                    j += 1
                if runs and (i - runs[-1][1]) * du < 0.03:
                    runs[-1][1] = j
                else:
                    runs.append([i, j])
                i = j + 1
            else:
                i += 1
        for i, j in runs:
            i0, j0 = max(i - nsm, 0), min(j + nsm, n - 1)
            if ysm[j0] - ysm[i0] < math.log(min_rise):
                continue
            ta, tb = float(grid[i0]), float(grid[j0])
            La, Lb = float(np.interp(ta, t, L)), float(np.interp(tb, t, L))
            near = ev_near(ta, tb, 0.05)
            dust_share = float("nan")
            if res.L_si_dust is not None and Lb > La:
                dust_share = (float(np.interp(tb, t, res.L_si_dust)) - float(np.interp(ta, t, res.L_si_dust))) / (Lb - La)
            nsplit = sum(1 for e in near if e["kind"] == "split")
            if near:
                what = ("; ".join(f"{e['kind']} gen {e['gen']} at {e['alt']/1e3:.1f} km ({e['mass']/s['m0']*100:.0f}% of m0)" for e in near[:3])
                        + (f" +{len(near)-3} more" if len(near) > 3 else ""))
            else:
                what = "no breakup within 0.05 s (growth of existing pancakes)"
            dtxt = f" Dust groups supply {100*dust_share:.0f}% of the rise." if np.isfinite(dust_share) else ""
            flares.append((ta, tb))
            add("flare", "F", 0.5 * (ta + tb), "stylized",
                f"Light rises x{Lb/max(La,1e-30):.1f} in {tb-ta:.2f} s (steeper than {k_flare:g}x the air-density tracking rate). Associated with: {what}.{dtxt} "
                f"A real fragmentation flare (Phase 5 in [R22]); its sharpness and timing here come from the model's deterministic cascade.",
                t_start=ta, t_end=tb, rise=Lb / max(La, 1e-30), n_splits=nsplit, dust_share=dust_share)
        # bumps: small local maxima not already part of a flare
        win = max(int(round(0.6 / du)), nsm * 3)
        detr = ysm - _moving_median(ysm, win)
        w = max(2, int(round(0.08 / du)))
        cand = []
        for i in range(w, n - w):
            if detr[i] > bump_min and detr[i] >= detr[i - w:i + w + 1].max():
                tt = float(grid[i])
                if any(ta - 0.1 <= tt <= tb + 0.1 for ta, tb in flares):
                    continue
                if cand and tt - cand[-1][0] < 0.1:
                    if detr[i] > cand[-1][1]:
                        cand[-1] = (tt, float(detr[i]))
                    continue
                cand.append((tt, float(detr[i])))
        for tt, h in cand:
            near = ev_near(tt, tt, 0.1)
            if near:
                e = min(near, key=lambda e: abs(e["t"] - tt))
                what = f"coincides with a {e['kind']} (gen {e['gen']}) at {e['alt']/1e3:.1f} km"
            else:
                what = "no breakup within 0.1 s (overlap of several groups' light / pancake growth)"
            add("bump", "b", tt, "stylized", f"Small fluctuation, {100*(math.exp(h)-1):.0f}% above the local trend; {what}. "
                f"The fine structure of this phase is set by the cascade rules, so its pattern is not predictive of a real event.", height=h)

    # ---- 6. peak, 7. roll-off ----
    add("peak", "P", t[ipk], "physical",
        f"Maximum light. It sits ~{s['alt_peak_L_km']-s['alt_peak_dEdh_km']:+.1f} km from the peak of energy deposition ({s['alt_peak_dEdh_km']:.1f} km) because the "
        f"luminous efficiency tau(v) falls as the fragments slow (mass-weighted speed here {s['v_at_peak_L_kms']:.1f} km/s).")
    after = np.flatnonzero(L[ipk:] < peak / math.e)
    if len(after):
        k = ipk + int(after[0])
        add("rolloff", "D", t[k], "physical",
            f"Light has fallen to 1/e of the peak {t[k]-t[ipk]:.2f} s after it. Cause: the kinetic energy is being used up (speed {res.v[k]/p.velocity*100:.0f}% of v0, "
            f"{res.mass[k]/s['m0']*100:.0f}% of the mass still in flight) and tau(v) drops with speed. Analogue of Phase 7 (decline until dark flight) in [R22].")
        stats["lc_decay_efold_s"] = float(t[k] - t[ipk])
    # drag overtakes ablation as the main energy channel
    if res.P_drag is not None:
        idx = np.flatnonzero((res.P_drag > res.P_abl) & (L > 1e-3 * peak) & (t > (ys["t"] if ys else t[0])))
        if len(idx):
            add("drag_overtakes_ablation", "", t[idx[0]], "stylized",
                "Drag work now exceeds the kinetic energy carried off by ablated mass. The ratio ablation/drag is about Lambda v^2/(2 Q Cd), so this happens "
                f"when the speed falls below ~{math.sqrt(2*p.L_ablation*p.Cd/p.Lambda)/1e3:.0f} km/s (from the model's heating law).")
    tail = np.flatnonzero(L[ipk:] < 1e-3 * peak)
    if len(tail):
        add("tail", "T", t[ipk + int(tail[0])], "low-confidence",
            f"Below 1e-3 of the peak the light is controlled by groups being retired when they slow below v_dark = {p.v_dark/1e3:.1f} km/s and by the "
            f"tau taper (2-4 km/s). Do not interpret structure below this level.")
    # ---- 10. retirement steps (artifact) ----
    last = -9.0
    nart = 0
    for r in s.get("retired", []):
        tr = r["t"]
        if tr - last < 0.02:
            continue
        Lb_, La_ = float(np.interp(tr - 0.015, t, L)), float(np.interp(tr + 0.015, t, L))
        if Lb_ > 1e-4 * peak and La_ < 0.85 * Lb_ and nart < 6:
            add("retirement_step", "R", tr, "artifact",
                f"A group is removed because it slowed below v_dark = {p.v_dark/1e3:.1f} km/s; the light drops {100*(1-La_/Lb_):.0f}% within 0.03 s. "
                f"A bookkeeping step, not a physical dimming.")
            nart += 1
            last = tr
    feats.sort(key=lambda f_: f_["t"])
    for k, f_ in enumerate(feats, 1):
        f_["id"] = k

    # ---- summary statistics ----
    E = np.cumsum(L * dts)
    Etot = float(E[-1]) if E[-1] > 0 else 1.0
    for frac in (0.1, 0.5, 0.9):
        tf = float(np.interp(frac * Etot, E, t))
        stats[f"lc_t{int(frac*100)}"] = tf
        stats[f"lc_alt{int(frac*100)}_km"] = float(np.interp(tf, t, res.alt)) / 1e3
    above = L >= 0.5 * peak
    stats["lc_fwhm_s"] = float((dts * above).sum())
    i10 = np.flatnonzero(L[:ipk + 1] >= 0.1 * peak)
    stats["lc_rise_time_s"] = float(t[ipk] - t[i10[0]]) if len(i10) else float("nan")
    stats["lc_E_prefrag_frac"] = float(np.interp(ys["t"], t, E)) / Etot if ys else float("nan")
    stats["lc_dust_E_frac"] = float((res.L_si_dust * dts).sum()) / Etot if res.L_si_dust is not None else float("nan")
    if res.P_drag is not None:
        Ed, Ea = float((res.P_drag * dts).sum()), float((res.P_abl * dts).sum())
        stats["lc_ablation_E_frac"] = Ea / (Ed + Ea) if (Ed + Ea) > 0 else float("nan")
    stats["lc_n_flares"] = sum(1 for f_ in feats if f_["kind"] == "flare")
    stats["lc_n_bumps"] = sum(1 for f_ in feats if f_["kind"] == "bump")
    return feats, stats


# ========================================================================================
# 12. COMMAND-LINE REPORT
# ========================================================================================
def format_report(res: BolumResult, max_events=12):
    """Human-readable run report (printed to the terminal and saved as run_report.txt)."""
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
    add(" LIGHT-CURVE FEATURES   (heuristic classifier; flags: physical = real process | stylized = real process but its shape/number is set by model rules |")
    add("                         artifact = bookkeeping effect | low-confidence = controlled by artifacts)")
    if res.features:
        add("   id tag  kind                       t(s)  alt(km)     L_si(W)  flag")
        shown = 0
        for f_ in res.features:
            if f_["kind"] == "bump" and shown >= 8:
                continue
            add(f"   {f_['id']:>2d}  {f_['tag'] or '-':<3s} {f_['kind']:<22s} {f_['t']-p.tstart:7.3f} {f_['alt_km']:8.2f} {f_['L_si']:11.3e}  {f_['flag']}")
            for ln in textwrap.wrap(f_["note"], width=92):
                add("          " + ln)
            shown += 1 if f_["kind"] == "bump" else 0
        nb_ = sum(1 for f_ in res.features if f_["kind"] == "bump")
        if nb_ > 8:
            add(f"   ... {nb_-8} more bumps in light_curve_features.dat")
        add(f"   timing          10/50/90 % of the radiated energy by t = {s['lc_t10']-p.tstart:.2f} / {s['lc_t50']-p.tstart:.2f} / {s['lc_t90']-p.tstart:.2f} s  "
            f"(~{s['lc_alt10_km']:.0f} / {s['lc_alt50_km']:.0f} / {s['lc_alt90_km']:.0f} km)")
        add(f"   shape           time above half-peak {s['lc_fwhm_s']:.2f} s;  rise (10% -> peak) {s['lc_rise_time_s']:.2f} s;  decay e-fold {s.get('lc_decay_efold_s', float('nan')):.2f} s")
        add(f"   budget          radiated before first yield {100*s['lc_E_prefrag_frac']:.1f}%;  from dust groups {100*s['lc_dust_E_frac']:.1f}%;  "
            f"ablation supplies {100*s.get('lc_ablation_E_frac', float('nan')):.0f}% of the deposited energy, drag the rest")
        add(f"   counts          {s['lc_n_flares']} flares, {s['lc_n_bumps']} small bumps")
    else:
        add("   (no features detected)")
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


# ========================================================================================
# 13. OUTPUT FILES
# ========================================================================================
def make_header(p: BolumParams):
    """One-line description of the event, printed at the top of every plot."""
    return (f"v = {p.velocity/1e3:.2f} km/s   d = {p.diameter:.2f} m   theta = {p.theta_deg:.1f} deg   "
            f"S0 = {p.strength0()/1e6:.2f} MPa (rho v^2)   rho = {p.density:.0f} kg/m3")


def write_outputs(outdir, res: BolumResult):
    """Write all data files into the event folder."""
    p, s = res.params, res.summary
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "init_params.dat"), "w") as fh:
        fh.write(f"Initialization Parameters for Bolide Luminosity and Fragmentation Simulation (bolum v{__version__})\n")
        for fl in fields(BolumParams):
            fh.write(f"{fl.name:22s} {getattr(p, fl.name)}\n")
    with open(os.path.join(outdir, "bolide_out.dat"), "w") as fh:
        fh.write("# time_s  alt_km  v_kms  mass_kg  n_groups  P_dep_W  L_si_W  L_bol_W  eta_si  loading_max_Pa  P_drag_W  P_abl_W  L_si_dust_W\n")
        for k in range(len(res.t)):
            fh.write(f"{res.t[k]:.5f} {res.alt[k]/1e3:.4f} {res.v[k]/1e3:.4f} {res.mass[k]:.5e} {res.n_groups[k]:d} "
                     f"{res.P_dep[k]:.5e} {res.L_si[k]:.5e} {res.L_bol[k]:.5e} {res.eta_eff[k]:.5f} {res.q_max[k]:.5e} "
                     f"{res.P_drag[k]:.5e} {res.P_abl[k]:.5e} {res.L_si_dust[k]:.5e}\n")
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
    with open(os.path.join(outdir, "light_curve_features.dat"), "w") as fh:
        fh.write("# id  tag  kind  t_s  alt_km  L_si_W  flag  note\n")
        for f_ in res.features:
            fh.write(f"{f_['id']} {f_['tag'] or '-'} {f_['kind']} {f_['t']:.4f} {f_['alt_km']:.3f} {f_['L_si']:.4e} {f_['flag']}  {f_['note']}\n")
    with open(os.path.join(outdir, "summary.dat"), "w") as fh:
        for k, v in s.items():
            if k not in ("retired", "warnings"):
                fh.write(f"{k:26s} {v}\n")
    with open(os.path.join(outdir, "run_report.txt"), "w") as fh:
        fh.write(format_report(res, max_events=40) + "\n")


# ========================================================================================
# 14. STANDARDIZED PLOTTING
# ========================================================================================
C_SI, C_BOL, C_DEP, C_PUB, C_BRK, C_SAT = "#1f77b4", "#ff7f0e", "black", "crimson", "seagreen", "purple"


def _alt_window(res):
    """Altitude window common to all altitude panels: where L_si > 1e-6 of peak, padded by 4 km."""
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
    """Rows (label, model value, published value, ok-flag) for the metrics table; model-only here."""
    s = res.summary
    return [
        ("KE0 (kt)", f"{s['KE0_kt']:.4g}", "", None),
        ("deposited energy (kt)", f"{s['E_dep_kt']:.4g}", "", None),
        ("first yield alt (km)", f"{s['first_yield_alt_km']:.1f}", "", None),
        ("breakups", f"{s['n_split_events']}", "", None),
        ("peak dE/dh (kt/km @ km)", f"{s['peak_dEdh_kt_per_km']:.3g} @ {s['alt_peak_dEdh_km']:.1f}", "", None),
        ("half-peak range (km)", f"{s['half_range_km'][0]:.1f}-{s['half_range_km'][1]:.1f}", "", None),
        ("peak L_Si (W)", f"{s['peak_L_si_W']:.2e}", "", None),
        ("peak I_Si (W/sr)", f"{s['peak_I_si_W_sr']:.2e}", "", None),
        ("alt of peak light (km)", f"{s['alt_peak_L_km']:.1f}", "", None),
        ("E_rad bolometric (kt)", f"{s['E_rad_bol_kt']:.4g}", "", None),
        ("bolometric efficiency", f"{100*s['eta_bol_overall']:.1f}%", "", None),
        ("luminous end alt (km)", f"{s['alt_end_lum_km']:.1f}", "", None),
        ("flares / bumps found", f"{s.get('lc_n_flares', 0)} / {s.get('lc_n_bumps', 0)}", "", None),
        ("time above half-peak (s)", f"{s.get('lc_fwhm_s', float('nan')):.2f}", "", None),
        ("radiated before 1st yield", f"{100*s.get('lc_E_prefrag_frac', float('nan')):.1f}%", "", None),
        ("energy-ledger closure", f"{s['energy_closure']:.1e}", "", None),
    ]


def _banner(fig, header, title=None):
    fig.text(.01, .985, "UNCLASSIFIED", fontsize=9, color="green", va="top")
    fig.text(.99, .012, "UNCLASSIFIED", fontsize=9, color="green", ha="right")
    if title:
        fig.suptitle(title, fontsize=13, y=0.995)
    fig.text(.5, .955, header, fontsize=9, ha="center")



_TAG_COLOR = {"A": "teal", "Y": "seagreen", "B": "seagreen", "F": "crimson", "b": "gray", "P": "black", "D": "black", "T": "gray", "R": "red"}


def _mark_features(ax, res, with_labels=True):
    """Thin vertical lines + a letter at the top of the axis for each detected light-curve feature.
       A ablation onset | Y first yield | B first breakup | F flare | b small bump | P peak | D 1/e decay | T low-confidence tail | R retirement artifact"""
    p = res.params
    for f_ in res.features:
        tag = f_["tag"]
        if not tag:
            continue
        x = f_["t"] - p.tstart
        ax.axvline(x, color=_TAG_COLOR.get(tag, "gray"), lw=0.6 if tag == "b" else 0.9, ls=":" if tag in ("b", "R", "T") else "--", alpha=.7)
        if with_labels and tag != "b":
            ax.text(x, 0.975, tag, transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=8.5,
                    color=_TAG_COLOR.get(tag, "gray"), weight="bold")


def plot_feature_guide(res, sat_W=None, title=None):
    """One page: the log light curve annotated with every detected feature, and a table saying what each one is."""
    import matplotlib.pyplot as plt
    p = res.params
    t = res.t - p.tstart
    ymin, ymax, tlo, thi = _time_log_window(res)
    fig = plt.figure(figsize=(14, 11))
    ax = fig.add_axes([0.07, 0.56, 0.90, 0.36])
    ax.semilogy(t, np.maximum(res.L_si, 1e-30), color=C_SI, label="Si band total")
    if res.L_si_dust is not None:
        ax.semilogy(t, np.maximum(res.L_si_dust, 1e-30), color=C_BRK, lw=0.9, alpha=.8, label="of which dust groups")
    if sat_W:
        ax.axhline(sat_W, color=C_SAT, ls=":", label=f"sensor saturation {sat_W:.1e} W")
    ax.set_ylim(ymin, ymax); ax.set_xlim(tlo, thi)
    _mark_features(ax, res)
    ax.set_xlabel("time since start (s)"); ax.set_ylabel("Si-band power (W)"); ax.grid(True, alpha=.3, which="both"); ax.legend(loc="lower right", fontsize=9)
    ax.set_title("Light-curve feature guide:  A ablation onset | Y first yield | B first breakup | F flare | b bump | P peak | D 1/e decay | T low-confidence tail | R retirement artifact", fontsize=9)
    lines = ["id tag kind                    t(s)  alt(km)   L_si(W)    flag"]
    for f_ in res.features:
        if f_["kind"] == "bump" and sum(1 for g in res.features if g["kind"] == "bump" and g["id"] <= f_["id"]) > 5:
            continue
        lines.append(f"{f_['id']:>2d}  {f_['tag'] or '-':<2s} {f_['kind']:<22s} {f_['t']-p.tstart:6.3f} {f_['alt_km']:7.2f} {f_['L_si']:10.2e}   {f_['flag']}")
        for ln in textwrap.wrap(f_["note"], width=135)[:3]:
            lines.append("      " + ln)
    nb_ = sum(1 for g in res.features if g["kind"] == "bump")
    if nb_ > 5:
        lines.append(f"   ... {nb_-5} more bumps listed in light_curve_features.dat")
    fig.text(0.05, 0.50, "\n".join(lines[:44]), fontsize=7.6, va="top", family="monospace")
    _banner(fig, make_header(p), title or "Light-curve feature guide")
    return fig


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
    _mark_features(a, res)
    a.set_xlabel("time since start (s)"); a.set_ylabel("power (W)"); a.set_title("(a) light curve  (letters: see feature guide)"); a.legend(fontsize=8, loc="lower right"); a.grid(True, alpha=.3, which="both")
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
        a.text(0.0, y, name, fontsize=8.5, va="top", family="monospace", transform=a.transAxes)
        a.text(0.46, y, mod, fontsize=8.5, va="top", family="monospace", transform=a.transAxes)
        if pub:
            a.text(0.70, y, pub, fontsize=8.5, va="top", family="monospace", transform=a.transAxes)
        if ok is not None:
            a.text(0.97, y, "ok" if ok else "MISS", fontsize=8.5, va="top", family="monospace", ha="right",
                   color=("green" if ok else "red"), weight="bold", transform=a.transAxes)
        y -= 0.055
    _banner(fig, make_header(p), title or "Bolide event summary")
    return fig


def plot_outputs(outdir, res: BolumResult, ref=None, rows=None, sat_W=None, title=None):
    """Write bolide_plots.pdf: the standardized summary page followed by single-quantity pages."""
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
    fig = plot_feature_guide(res, sat_W=sat_W, title=title)
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


# ========================================================================================
# 15. MAIN
# ========================================================================================
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
        files = ["init_params.dat", "bolide_out.dat", "power_vs_time.dat", "power_vs_alt.dat", "dEdh.dat", "events.dat", "light_curve_features.dat", "summary.dat", "run_report.txt"]
        if not a.no_plots:
            files.append("bolide_plots.pdf")
        print(f" outputs written to {outdir}/ : " + ", ".join(files))


if __name__ == "__main__":
    main()
