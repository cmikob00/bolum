#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
bolum_block_v2.py  --  reduced-order bolide entry / fragmentation / luminosity model

Rewrite of bolum_block.py (steps 1-3 of the improvement plan, plus Si-band output).

WHAT CHANGED vs. bolum_block.py
-------------------------------
Step 1 (bugs / hygiene)
  * n_fragments no longer compounds (3 -> 9 -> 729 ...).  Fragmentation is now a cascade of
    explicit fragment *groups*, each with its own mass, velocity, strength and area.
  * tau(v,m) fit is clamped to its calibration range and tapered to zero below ~4 km/s
    (the polynomial in ln v blows up below ~2 km/s).
  * Flare multiplier (n_frag * decay) removed.  It created energy and double counted area growth.
  * Atmosphere: US Standard 1976 below 86 km (analytic layers) + tabulated density 86-120 km.
    The NASA-GRC fit used before was 2x too dense at 80 km and ~10x too dense at 100 km.
  * Earth-centred Cartesian trajectory (no alpha/gamma/phi/law-of-cosines machinery, no
    theta = arctan(|vz|/vx) singularities), vectorised over groups, adaptive time step.
  * No per-step vstack / print / file I/O: preallocated lists, files written once at the end.
  * Exact kinetic-energy ledger: KE0 + gravity work = E_deposited + KE_left (to round-off).
Step 2 (ablation / heat transfer)
  * Heating is Lambda * (1/2) rho v^3 on the *projected* area (Lambda <= 1), decoupled from Cd.
    Previously the flux was applied to the full sphere area (4x the available kinetic-energy flux)
    with Lambda == Cd == 1, i.e. an effective sigma ~ 0.4 s^2/km^2 (typical stony: 0.01-0.05).
  * Ablation only above a net-flux threshold (radiative loss at T_ablate).
Step 3 (fragment engine)
  * Progressive fragmentation with Weibull strength scaling S ~ S_parent (m_parent/m)^alpha.
  * Pancake lateral spreading of any group loaded above its strength (Chyba/Hills-Goda/Collins
    style, d r/dt = c * v * sqrt(rho_air/rho_body)), splitting into child groups when the area has
    grown by f_split^2.  Area is continuous across a split.
  * Differential drag separates groups naturally -> multiple flares, no flare hack.
Luminosity
  * Per-group radiated power = eta * (drag work + kinetic energy of ablated mass), with
      eta = (1 - w) * tau(v, m) + w * eta_thick(v),   w = r_proj / (r_proj + lum_r_half)
    tau(v,m) is your original fit (taken as SILICON-BAND by default, tau_basis='si'), so small
    bolides reproduce the old light output.  eta_thick is a bounded, v-independent bolometric
    efficiency for large optically-thick clouds (replaces the (rho/1e-3)^0.45 tuning, which was
    compensating for tau(v) falling as the fragments decelerate).  PLACEHOLDER for step 4.
  * Output is Si-band power L_si AND bolometric L_bol = L_si / f_si(T_rad), with f_si the blackbody
    in-band fraction for the sensor band [sensor_lo_um, sensor_hi_um].

Input deck: same format as before (first 2 lines skipped, 'key = value  # comment').  Legacy keys
rho_tau_scale, flare_duration, n_frag_init are accepted and ignored.  strength_secondary/tertiary/
quaternary are optional overrides for the child strength at generations 1/2/3 (otherwise Weibull).
"""
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
JKT = 4.184e12               # J per kiloton TNT
SIGMA_SB = 5.670374419e-8
H_PL, C_LIGHT, K_B = 6.62607015e-34, 2.99792458e8, 1.380649e-23
DEFAULT_BASE = "/Users/julie/Desktop/Projects/bolum/Events/"

_trapz = getattr(np, "trapezoid", None) or np.trapz


# ----------------------------------------------------------------------------------------
# atmosphere: US Standard 1976 (analytic to 86 km, tabulated density 86-120 km)
# ----------------------------------------------------------------------------------------
class Atmosphere:
    """Density vs geometric altitude.  rho(alt_m) is vectorised; extrapolates exponentially >120 km."""

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
        # upper table (km); start just above 86 km so the analytic part owns <= 86 km
        hk = np.concatenate([h / 1e3, self._UPPER_KM[1:]])
        lr = np.concatenate([np.log(rho_low), np.log(self._UPPER_RHO[1:])])
        self._hk = np.arange(0.0, 120.0001, 0.1)
        self._lr = np.interp(self._hk, hk, lr)
        self._scale_h = 1e3 * (self._hk[-1] - self._hk[-51]) / (self._lr[-51] - self._lr[-1])  # m, near 120 km

    def rho(self, alt_m):
        a = np.asarray(alt_m, dtype=float) / 1e3
        lr = np.interp(a, self._hk, self._lr)
        hi = a > 120.0
        if np.any(hi):
            lr = np.where(hi, self._lr[-1] - (a - 120.0) * 1e3 / self._scale_h, lr)
        return np.exp(lr)


ATM = Atmosphere()


# ----------------------------------------------------------------------------------------
# luminous efficiency, sensor band, and radiation helpers
# ----------------------------------------------------------------------------------------
def tau_fit(v, m):
    """Original empirical tau(v, m) (fraction, not percent), vectorised, clamped and tapered.

    Valid only for meteor-like speeds.  Below ~4 km/s the luminous flight is ending, so the
    efficiency is tapered smoothly to zero (the raw polynomial in ln v rises again below ~2 km/s)."""
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
    """Fraction of blackbody exitance (temperature T [K]) falling in the sensor band.

    Default: top-hat [lo_um, hi_um] (broadband silicon, 0.4-1.1 um).  Pass lam_um/resp (peak-normalised
    relative spectral response) to use your own sensor curve."""
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


# ----------------------------------------------------------------------------------------
# parameters
# ----------------------------------------------------------------------------------------
@dataclass
class BolumParams:
    # --- required ---
    diameter: float = 1.0            # m
    velocity: float = 20.0e3         # m/s (atmosphere-relative)
    theta_deg: float = 45.0          # entry angle to local horizontal at zstart
    init_strength: float = 1.0e6     # Pa (loading at which the intact body yields; see q_factor)
    density: float = 3300.0          # bulk density, kg/m^3
    # --- material / aerodynamics ---
    porosity: float = 0.1            # only used if init_strength <= 0
    Cd: float = 1.0
    # --- ablation (step 2) ---
    L_ablation: float = 5.0e6        # J/kg, effective heat of ablation Q
    Lambda: float = 0.05             # heat-transfer coefficient (<=1), on projected area (sigma = Lambda/2Q ~ 0.005 s^2/km^2)
    T_ablate: float = 2300.0         # K, ablation onset (radiative-equilibrium threshold)
    emissivity: float = 0.9
    # --- fragmentation / pancake engine (step 3) ---
    n_fragments: int = 4             # children per breakup
    frag_index: float = 1.0          # child mass share ~ k^-frag_index
    weibull_alpha: float = 0.40      # S_child = S_parent (m_parent/m_child)^alpha
    q_factor: float = 0.5            # loading q = q_factor*rho*v^2.  0.5 = dynamic pressure (original code);
                                     # 1.0 = rho*v^2 convention used in Borovicka et al. strength tables
    c_spread: float = 0.5            # lateral spreading coefficient
    f_split: float = 2.0             # split after the group's linear size grew by this factor
    f_max: float = 8.0               # pancake limit (linear)
    max_groups: int = 128
    m_min_frag: float = 1.0          # kg; groups below this do not split further
    dust_frac: float = 0.20          # mass fraction released as fine dust at each breakup
    dust_radius: float = 2.0e-2      # m, characteristic dust particle radius (sets dust area/mass)
    strength_secondary: float = None
    strength_tertiary: float = None
    strength_quaternary: float = None
    # --- luminosity / sensor ---
    tau_basis: str = "si"            # 'si': tau fit is Si-band efficiency; 'bolometric'
    T_rad: float = 6000.0            # K, blackbody temperature for band conversion (Brown et al. convention)
    sensor_lo_um: float = 0.4
    sensor_hi_um: float = 1.1
    eta_thick_bol: float = 0.20      # bolometric radiated fraction of a large, optically thick cloud
    lum_r_half: float = 6.0          # m, projected radius at which thin/thick efficiencies are blended 50/50
    eta_max: float = 0.6             # hard cap on radiated fraction
    t_cool: float = 0.0              # s, optional first-order radiative relaxation of the light curve
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
        return 1.0e5 * (1.0 - self.porosity) * (self.density / 3000.0)


_LEGACY_IGNORED = {"rho_tau_scale", "flare_duration", "n_frag_init"}


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
    q_max: np.ndarray          # largest dynamic pressure among groups (Pa)
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

    # group state arrays
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

    # altitude bins for dE/dh
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
    dt_prev = p.dt

    def spawn_children(i, t_now, alt_now, q_now):
        """Replace group i by n_child children plus a dust group.  Mass conserving; area ~continuous."""
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
            fk = max(1.0, fp * (mk / mp) ** (1.0 / 6.0))      # keeps sum of areas ~ the parent's
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
        q = p.q_factor * rho * vi * vi
        rs = (3.0 * mi / (4.0 * PI * rho_b)) ** (1.0 / 3.0)
        A = np.where(du, 3.0 * mi / (4.0 * rho_b * p.dust_radius), PI * rs * rs * fi * fi)   # projected area
        F = 0.5 * p.Cd * rho * vi * vi * A
        a_dr = F / mi
        # --- heating / ablation (step 2): Lambda * (1/2) rho v^3 on the projected area ---
        qh = p.Lambda * 0.5 * rho * vi ** 3
        mdot = np.maximum(qh - q_loss, 0.0) * A / p.L_ablation
        # --- yield + pancake spreading (step 3) ---
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
        dm = mi * (-np.expm1(-mdot * dt / mi))                   # exponential (positive-definite) mass loss
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
        # --- radiation ---
        tau = tau_fit(vi, mi)
        r_proj = np.sqrt(A / PI)
        w_thick = r_proj / (r_proj + p.lum_r_half)               # 0: meteor-like (thin), 1: large cloud (thick)
        eta_thin = tau if basis_si else tau * f_si                # Si-band efficiency of the original fit
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
        t += dt
        T_.append(t); L_si.append(E_rad_si.sum() / dt); P_.append(E_dep_s / dt)
        wts = E_dep if E_dep_s > 0 else mi
        ALT.append(float((wts * alt).sum() / wts.sum()))
        V_.append(float((mi * vi).sum() / mi.sum()) if mi.sum() > 0 else 0.0)
        M_.append(float(mnew.sum())); N_.append(int((mnew > 0).sum())); Q_.append(float(q.max()))
        # --- commit state ---
        m[ia], x[ia], z[ia], vx[ia], vz[ia], f[ia] = mnew, xn, zn, vxn, vzn, fn
        gone = (mnew <= 0)
        alive[ia[gone]] = False
        # retire groups that hit the ground or leave luminous flight
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
            qk = p.q_factor * float(ATM.rho(altk)) * vk * vk
            if not spawn_children(k, t, altk, qk):
                fb[k] = np.inf                                   # cannot split: keep as a pancake
        nsteps += 1

    # ------------------------------------------------------------------ assemble outputs
    T_ = np.array(T_); L_si = np.array(L_si); P_ = np.array(P_)
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
    # smoothed dE/dh (3-bin boxcar) -> peak altitude and half-energy range
    dEdh = dep_bins / bin_w
    sm = np.convolve(dEdh, np.ones(3) / 3.0, mode="same")
    ipk = int(np.argmax(sm)) if sm.max() > 0 else 0
    half = sm >= 0.5 * sm[ipk]
    lo = hi = ipk
    while lo > 0 and half[lo - 1]:
        lo -= 1
    while hi < nb - 1 and half[hi + 1]:
        hi += 1
    ipl = int(np.argmax(L_si)) if len(L_si) else 0
    lit = np.flatnonzero(L_si > 0.01 * L_si.max()) if len(L_si) else np.array([], dtype=int)
    alt_end_lum = float(ALT[lit[-1]]) / 1e3 if len(lit) else float("nan")   # last altitude with L > 1% of peak
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
        alt_end_lum_km=alt_end_lum,
        alt_peak_dEdh_km=alt_c[ipk] / 1e3,
        peak_dEdh_kt_per_km=sm[ipk] * 1e3 / JKT,
        half_range_km=(alt_c[lo] / 1e3, alt_c[hi] / 1e3),
        first_yield_alt_km=first_yield["alt"] / 1e3 if first_yield else float("nan"),
        first_yield_q_MPa=first_yield["q"] / 1e6 if first_yield else float("nan"),
        n_split_events=sum(1 for e in events if e["kind"] == "split"),
        mass_survived_kg=M_retired, mass_survived_frac=M_retired / m0,
        mass_survived_solid_frac=sum(r['mass'] for r in retired if not r['dust']) / m0,
        mass_survived_dust_frac=sum(r['mass'] for r in retired if r['dust']) / m0,
        ground_ke_J=float(sum(r["ke"] for r in retired if r["ground"])),
        energy_closure=closure,
        steps=nsteps, wall_s=time.time() - wall0,
        duration_s=float(T_[-1] - p.tstart) if len(T_) else 0.0,
        retired=retired,
    )
    return BolumResult(p, T_, L_si, L_bol, P_, np.array(ALT), np.array(V_), np.array(M_), np.array(N_),
                       np.array(Q_), L_si / np.maximum(P_, 1e-30), events, alt_c, dEdh, rad_bins / bin_w, summary)


# ----------------------------------------------------------------------------------------
# output files and plots
# ----------------------------------------------------------------------------------------
def make_header(p: BolumParams):
    return (f"v = {p.velocity/1e3:.2f} km/s   d = {p.diameter:.2f} m   theta = {p.theta_deg:.1f} deg   "
            f"S0 = {p.init_strength:.2e} Pa   rho = {p.density:.0f} kg/m3")


def write_outputs(outdir, res: BolumResult):
    p, s = res.params, res.summary
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "init_params.dat"), "w") as fh:
        fh.write("Initialization Parameters for Bolide Luminosity and Fragmentation Simulation (bolum v2)\n")
        for fl in fields(BolumParams):
            fh.write(f"{fl.name:22s} {getattr(p, fl.name)}\n")
    with open(os.path.join(outdir, "bolide_out.dat"), "w") as fh:
        fh.write("# time_s  alt_km  v_kms  mass_kg  n_groups  P_dep_W  L_si_W  L_bol_W  eta_si  q_max_Pa\n")
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
        fh.write("# kind  t_s  alt_km  q_MPa  mass_kg  strength_MPa  generation\n")
        for e in res.events:
            fh.write(f"{e['kind']} {e['t']:.4f} {e['alt']/1e3:.3f} {e['q']/1e6:.4f} {e['mass']:.4e} {e['strength']/1e6:.4f} {e['gen']}\n")
    with open(os.path.join(outdir, "summary.dat"), "w") as fh:
        for k, v in s.items():
            if k != "retired":
                fh.write(f"{k:26s} {v}\n")


def plot_outputs(outdir, res: BolumResult):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    p, s = res.params, res.summary
    header = make_header(p)
    t = res.t - p.tstart
    pdf = PdfPages(os.path.join(outdir, "bolide_plots.pdf"))

    def finish(fig, title):
        plt.title(title + "\n")
        plt.grid(True)
        fig.text(.10, .95, "UNCLASSIFIED", fontsize=10, color="green")
        fig.text(.80, .05, "UNCLASSIFIED", fontsize=10, color="green")
        fig.text(.12, .91, header, fontsize=9)
        pdf.savefig(fig)
        plt.close(fig)

    notes = (f"Rad E (bol) = {s['E_rad_bol_kt']:.3f} kt = {s['E_rad_bol_J']:.3e} J\n"
             f"Rad E (Si band) = {s['E_rad_si_J']:.3e} J\n"
             f"Deposited E = {s['E_dep_kt']:.3f} kt (KE0 = {s['KE0_kt']:.3f} kt)\n"
             f"Peak L_Si = {s['peak_L_si_W']:.3e} W  ({s['peak_I_si_W_sr']:.3e} W/sr)")
    for scale in ("linear", "log"):
        fig = plt.figure(figsize=(9, 6.5))
        plt.plot(t, res.L_si, color="blue", label="Si-band L")
        plt.plot(t, res.L_bol, color="darkorange", alpha=0.7, label="Bolometric L")
        if scale == "log":
            plt.xscale("log"); plt.yscale("log")
            plt.xlim(max(t[0], 1e-3), 1.05 * t[-1]); plt.ylim(max(1e-3 * res.L_si.max(), 1.0), 3 * res.L_bol.max())
        for e in res.events[:1]:
            plt.axvline(e["t"] - p.tstart, ls="--", color="green", alpha=.6,
                        label=f"first yield {e['alt']/1e3:.1f} km")
        plt.xlabel("Time (s)"); plt.ylabel("Luminosity (W)")
        fig.text(.14, .66, notes, fontsize=9, va="top")
        plt.legend(loc="upper right")
        finish(fig, f"Bolide Light Curve ({scale})")
    fig = plt.figure(figsize=(9, 6.5))
    plt.plot(t, res.L_si / (4 * PI), color="blue")
    plt.xlabel("Time (s)"); plt.ylabel("Si-band radiant intensity (W/sr)")
    finish(fig, "Radiant Intensity vs Time (isotropic)")
    fig = plt.figure(figsize=(9, 6.5))
    plt.plot(res.L_si, res.alt / 1e3, color="blue", label="Si band")
    plt.plot(res.L_bol, res.alt / 1e3, color="darkorange", alpha=.7, label="Bolometric")
    plt.xscale("log"); plt.xlim(max(1e-3 * res.L_si.max(), 1.0), 3 * res.L_bol.max())
    plt.xlabel("Luminosity (W)"); plt.ylabel("Altitude (km)"); plt.legend()
    finish(fig, "Luminosity vs Altitude")
    fig = plt.figure(figsize=(9, 6.5))
    plt.plot(res.dEdh_dep * 1e3 / JKT, res.alt_bins / 1e3, label="Deposited")
    plt.plot(res.dEdh_rad_si * 1e3 / JKT, res.alt_bins / 1e3, label="Radiated (Si band)")
    plt.ylim(0, min(120, res.alt_bins[-1] / 1e3)); plt.xlabel("dE/dh (kt TNT per km)"); plt.ylabel("Altitude (km)")
    plt.legend()
    finish(fig, "Energy Deposition vs Altitude")
    fig = plt.figure(figsize=(9, 6.5))
    ax = plt.gca(); ax.plot(t, res.mass, color="k"); ax.set_yscale("log")
    ax.set_xlabel("Time (s)"); ax.set_ylabel("Total mass in flight (kg)")
    ax2 = ax.twinx(); ax2.plot(t, res.n_groups, color="red", alpha=.7); ax2.set_ylabel("Fragment groups", color="red")
    finish(fig, "Mass and Fragment-Group Count")
    fig = plt.figure(figsize=(9, 6.5))
    plt.plot(res.alt / 1e3, res.q_max / 1e6, label="max dynamic pressure")
    plt.yscale("log"); plt.xlabel("Altitude (km)"); plt.ylabel("q (MPa)"); plt.gca().invert_xaxis()
    ax2 = plt.gca().twinx(); ax2.plot(res.alt / 1e3, res.v / 1e3, color="red", alpha=.7); ax2.set_ylabel("v (km/s)", color="red")
    finish(fig, "Dynamic Pressure and Speed vs Altitude")
    fig = plt.figure(figsize=(9, 6.5))
    sp = [e for e in res.events if e["kind"] == "split"]
    if sp:
        plt.scatter([e["t"] - p.tstart for e in sp], [e["alt"] / 1e3 for e in sp],
                    s=[10 + 60 * (e["mass"] / s["m0"]) ** 0.4 for e in sp], c="green", alpha=.6)
    plt.xlabel("Time (s)"); plt.ylabel("Breakup altitude (km)")
    finish(fig, "Fragmentation History (marker size ~ parent mass)")
    pdf.close()


# ----------------------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description="bolum v2 bolide model")
    ap.add_argument("event", nargs="?", help="event folder name containing bolum_in.dat")
    ap.add_argument("--base", default=os.environ.get("BOLUM_BASE", DEFAULT_BASE), help="events base directory")
    ap.add_argument("--no-plots", action="store_true")
    a = ap.parse_args(argv)
    event = a.event or input("Enter event name: ")
    outdir = os.path.join(a.base, event)
    p = rdinput(os.path.join(outdir, "bolum_in.dat"))
    res = run_bolide(p)
    write_outputs(outdir, res)
    if not a.no_plots:
        plot_outputs(outdir, res)
    s = res.summary
    print(f"bolum v2: {s['steps']} steps in {s['wall_s']:.2f} s | E_dep {s['E_dep_kt']:.3f} kt of KE0 {s['KE0_kt']:.3f} kt | "
          f"E_rad(bol) {s['E_rad_bol_kt']:.4f} kt | peak L_Si {s['peak_L_si_W']:.3e} W at {s['alt_peak_L_km']:.1f} km | "
          f"energy closure {s['energy_closure']:.2e}")


if __name__ == "__main__":
    main()
