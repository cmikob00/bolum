#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
bolum_diagnostics.py -- verification, convergence and validation suite for bolum_block_v2.py

  A. component checks      atmosphere vs US Std 1976, tau(v,m) behaviour, Si band fraction
  B. verification          engine vs an independent RK45 integration of the same single-body ODEs
  C. conservation          kinetic-energy ledger closure on every run in this suite
  D. convergence           time step, group cap, children per breakup
  E. validation            Chelyabinsk 2013, 2008 TC3, Kosice 2010 vs published values
  F. size / geometry sweep 0.5-25 m, 3 strengths; trend checks; radiated fraction vs Brown et al. (2002)
  G. run time

NOTE: Chelyabinsk and TC3 were used to calibrate the *global* parameters (c_spread, f_split, alpha,
dust_frac/radius, eta_thick/lum_r_half).  They are calibration checks, not independent validation.
Kosice and the size/geometry sweeps are not used for calibration.
Usage:  python bolum_diagnostics.py [outdir]
"""
import sys
import os
import math
import time
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from bolum_block_v2 import (BolumParams, run_bolide, ATM, tau_fit, band_fraction, RE, GM, JKT, PI, SIGMA_SB)

OUT = sys.argv[1] if len(sys.argv) > 1 else "."
os.makedirs(OUT, exist_ok=True)
REPORT, RUNS = [], []


def check(name, ok, detail):
    REPORT.append((name, bool(ok), detail))
    print(("PASS  " if ok else "FAIL  ") + f"{name}: {detail}")


def run(tag, **kw):
    r = run_bolide(BolumParams(**kw))
    RUNS.append((tag, r))
    return r


# =====================================================================  A. components
print("\n=== A. component checks ===")
std = {0: 1.225, 10: 0.4135, 20: 0.08891, 30: 0.01841, 40: 3.996e-3, 50: 1.027e-3,
       60: 3.097e-4, 70: 8.283e-5, 80: 1.846e-5}
worst = max(abs(float(ATM.rho(h * 1e3)) / r - 1) for h, r in std.items())
check("atmosphere vs US Std 1976 (0-80 km)", worst < 0.02, f"worst error {worst*100:.2f}%")
r86 = float(ATM.rho(86e3))
check("atmosphere continuity at 86 km", abs(r86 / 6.958e-6 - 1) < 0.02, f"rho(86 km) = {r86:.3e} vs 6.958e-6")
vv = np.linspace(500, 40e3, 400)
t_hi = tau_fit(vv, 1e3)
check("tau(v) finite and non-negative 0.5-40 km/s", np.all(np.isfinite(t_hi)) and np.all(t_hi >= 0),
      f"max {t_hi.max()*100:.2f}%, tau(1 km/s) = {float(tau_fit(1e3, 1e3)):.1e} (raw polynomial would give ~1.8%)")
t20 = float(tau_fit(20e3, 1e3))
check("tau(20 km/s, 1 t) matches original fit", abs(t20 - 0.0682) < 0.004, f"tau = {t20*100:.2f}% (original code: ~6.8%)")
f6 = float(band_fraction(6000.)[0])
ftot = float(band_fraction(6000., 0.05, 100.)[0])
check("Si-band blackbody fraction (0.4-1.1 um, 6000 K)", abs(f6 - 0.643) < 0.01,
      f"f_si(6000 K) = {f6:.3f}; f_si(4500 K) = {float(band_fraction(4500.)[0]):.3f}; f_si(3500 K) = {float(band_fraction(3500.)[0]):.3f}")
check("band integral sanity (0.05-100 um -> 1)", abs(ftot - 1) < 0.02, f"{ftot:.3f}")

# =====================================================================  B. verification vs independent ODE
print("\n=== B. verification against independent RK45 integration (single intact body) ===")


def reference(p, alt_stop_m):
    rho_b = p.density
    m0 = 4 / 3 * PI * (p.diameter / 2) ** 3 * rho_b
    q_loss = p.emissivity * SIGMA_SB * p.T_ablate ** 4
    th = math.radians(p.theta_deg)

    def rhs(t, y):
        x, z, vx, vz, m = y
        r = math.hypot(x, z)
        v = math.hypot(vx, vz)
        rho = float(ATM.rho(r - RE))
        rs = (3 * m / (4 * PI * rho_b)) ** (1 / 3)
        A = PI * rs * rs
        a = 0.5 * p.Cd * rho * v * A / m
        md = max(p.Lambda * 0.5 * rho * v ** 3 - q_loss, 0.0) * A / p.L_ablation
        g = GM / r ** 2
        return [vx, vz, -a * vx - g * x / r, -a * vz - g * z / r, -md]

    def ev(t, y):
        return math.hypot(y[0], y[1]) - RE - alt_stop_m
    ev.terminal = True
    y0 = [0.0, RE + p.zstart, p.velocity * math.cos(th), -p.velocity * math.sin(th), m0]
    s = solve_ivp(rhs, [0, 300], y0, events=ev, rtol=1e-10, atol=[1e-3, 1e-3, 1e-6, 1e-6, 1e-9], method="DOP853")
    y = s.y[:, -1]
    return y[4], math.hypot(y[2], y[3])


for tag, kw in {"2 m, 15 km/s, 30 deg": dict(diameter=2.0, velocity=15e3, theta_deg=30., density=3300.),
                "0.5 m, 20 km/s, 60 deg": dict(diameter=0.5, velocity=20e3, theta_deg=60., density=3300.),
                "5 m, 25 km/s, 45 deg": dict(diameter=5.0, velocity=25e3, theta_deg=45., density=3300.)}.items():
    p = BolumParams(init_strength=1e12, **kw)
    res = run("verify " + tag, init_strength=1e12, **kw)
    alt_stop = 35e3 if kw["diameter"] > 1 else 45e3
    m_ref, v_ref = reference(p, alt_stop)
    a_ = res.alt[::-1]
    m_eng = float(np.interp(alt_stop, a_, res.mass[::-1]))
    v_eng = float(np.interp(alt_stop, a_, res.v[::-1]))
    em, ev_ = abs(m_eng / m_ref - 1), abs(v_eng / v_ref - 1)
    check(f"engine vs RK45 [{tag}] at {alt_stop/1e3:.0f} km", em < 0.01 and ev_ < 0.01,
          f"mass err {em*100:.3f}% (m_ref/m0 = {m_ref/res.summary['m0']:.4f}), speed err {ev_*100:.3f}%")

# =====================================================================  small-bolide recovery of original Si tau
print("\n=== B2. small bolides recover the original Si-band efficiency ===")
res = run("0.3 m intact", diameter=0.3, velocity=20e3, theta_deg=45., init_strength=1e9, density=3300.)
dts = np.diff(np.concatenate([[0.0], res.t]))
E_tau = float(np.sum(tau_fit(res.v, np.maximum(res.mass, 1e-3)) * res.P_dep * dts))
ratio = res.summary["E_rad_si_J"] / E_tau
check("0.3 m: integrated L_Si / (tau_fit x deposited energy)", 0.93 < ratio < 1.12, f"ratio {ratio:.3f}")

# =====================================================================  D. convergence
print("\n=== D. convergence (Chelyabinsk-like) ===")
ch = dict(diameter=19.6, velocity=19.16e3, theta_deg=18.3, init_strength=0.7e6, density=3300.)
base = run("chely", **ch)
rows = {}
for dt in (0.012, 0.006, 0.0015):
    rows[dt] = run(f"chely dt={dt}", dt=dt, **ch).summary
for dt, s in rows.items():
    d1 = abs(s["E_rad_bol_kt"] / base.summary["E_rad_bol_kt"] - 1)
    d2 = abs(s["peak_L_si_W"] / base.summary["peak_L_si_W"] - 1)
    print(f"     dt_max={dt:<7} steps={s['steps']:6d}  dE_rad={d1*100:.2f}%  dL_peak={d2*100:.2f}%")
s = rows[0.0015]
check("time-step convergence (3 ms vs 1.5 ms)",
      abs(s["E_rad_bol_kt"] / base.summary["E_rad_bol_kt"] - 1) < 0.01 and abs(s["peak_L_si_W"] / base.summary["peak_L_si_W"] - 1) < 0.01,
      "E_rad and peak L change < 1%")
s2 = run("chely NG=512", max_groups=512, **ch).summary
check("group cap not binding (128 vs 512)", abs(s2["E_rad_bol_kt"] / base.summary["E_rad_bol_kt"] - 1) < 0.005,
      f"E_rad change {abs(s2['E_rad_bol_kt']/base.summary['E_rad_bol_kt']-1)*100:.2f}%")
for n in (3, 6):
    sn = run(f"chely n={n}", n_fragments=n, **ch).summary
    print(f"     n_fragments={n}: E_rad(bol) {sn['E_rad_bol_kt']:.1f} kt, peak dE/dh {sn['alt_peak_dEdh_km']:.1f} km")

# =====================================================================  E. validation
print("\n=== E. event comparison ===")
cb = base.summary
lo, hi = cb["half_range_km"]
KE_pub = 590.0
print(f"     Chelyabinsk: KE0 (model, 1.3e7 kg) = {cb['KE0_kt']:.0f} kt; published 430-640 kt (Brown+2013), 590+-50 (Popova+2013)")
check("Chelyabinsk: deposited energy ~ all of KE0", cb["E_dep_kt"] / cb["KE0_kt"] > 0.98, f"E_dep/KE0 = {cb['E_dep_kt']/cb['KE0_kt']:.4f}")
check("Chelyabinsk: peak dE/dh altitude (pub. ~29.5 km)", abs(cb["alt_peak_dEdh_km"] - 29.5) < 1.5, f"{cb['alt_peak_dEdh_km']:.1f} km")
check("Chelyabinsk: half-energy range (pub. 27-33 km)", abs(lo - 27) < 1.5 and abs(hi - 33) < 1.5, f"{lo:.1f}-{hi:.1f} km")
check("Chelyabinsk: peak dE/dh (pub. ~80 kt/km)", 0.6 < cb["peak_dEdh_kt_per_km"] / 80 < 1.4, f"{cb['peak_dEdh_kt_per_km']:.0f} kt/km")
check("Chelyabinsk: radiated energy (pub. 3.75e14 J bolometric)", abs(math.log(cb["E_rad_bol_J"] / 3.75e14)) < math.log(1.25),
      f"{cb['E_rad_bol_J']:.2e} J ({cb['eta_bol_overall']*100:.1f}% of deposited)")
check("Chelyabinsk: peak Si-band intensity (pub. 2.7e13 W/sr)", abs(math.log(cb["peak_I_si_W_sr"] / 2.7e13)) < math.log(1.4),
      f"{cb['peak_I_si_W_sr']:.2e} W/sr")
check("Chelyabinsk: altitude of peak light (pub. 29.7 km)", abs(cb["alt_peak_L_km"] - 29.7) < 2.5, f"{cb['alt_peak_L_km']:.1f} km")

tc = run("TC3", diameter=4.1, velocity=12.8e3, theta_deg=19.5, init_strength=0.15e6, density=2300.)
ct = tc.summary
y1 = next(e for e in tc.events if e["kind"] == "yield")
print(f"     TC3: KE0 = {ct['KE0_kt']:.2f} kt (pub. 1.1-2.1 kt)")
check("TC3: first fragmentation flare (pub. 45-53 km)", 43 < y1["alt"] / 1e3 < 54, f"first breakup {y1['alt']/1e3:.1f} km at q = {y1['q']/1e6:.2f} MPa")
check("TC3: peak light altitude (pub. 37.5 km)", abs(ct["alt_peak_L_km"] - 37.5) < 2.5, f"{ct['alt_peak_L_km']:.1f} km")
check("TC3: peak energy deposition (pub. burst ~37 km)", abs(ct["alt_peak_dEdh_km"] - 37.0) < 3.0, f"{ct['alt_peak_dEdh_km']:.1f} km")
check("TC3: radiated energy (pub. ~0.1 kt, factor-1.6)", 1 / 1.6 < ct["E_rad_bol_kt"] / 0.1 < 1.6, f"{ct['E_rad_bol_kt']:.3f} kt")

ko = run("Kosice", diameter=1.25, velocity=15e3, theta_deg=60., init_strength=0.1e6, strength_secondary=1.0e6, density=3400., q_factor=1.0)
ck = ko.summary
print(f"     Kosice (NOT used in calibration; published strengths 0.1 then ~1 MPa, rho*v^2 convention -> q_factor=1): "
      f"peak L {ck['alt_peak_L_km']:.1f} km (pub. 36), end of luminosity {ck['alt_end_lum_km']:.1f} km (pub. 17), "
      f"first breakup {ck['first_yield_alt_km']:.1f} km (pub. 55-57)")
check("Kosice: first breakup altitude (pub. 55-57 km)", abs(ck["first_yield_alt_km"] - 56) < 6, f"{ck['first_yield_alt_km']:.1f} km")
check("Kosice: altitude of max brightness (pub. 36 km)", abs(ck["alt_peak_L_km"] - 36) < 5, f"{ck['alt_peak_L_km']:.1f} km")
check("Kosice: luminous end height (pub. 17 km)", abs(ck["alt_end_lum_km"] - 17) < 6, f"{ck['alt_end_lum_km']:.1f} km")

# =====================================================================  F. sweeps
print("\n=== F. size sweep (v = 17 km/s, 45 deg, 3300 kg/m^3) ===")
diams = [0.5, 1, 2, 4, 8, 12, 20, 25]
strengths = [0.15e6, 0.7e6, 2.0e6]
sweep = {}
print(f"     {'d(m)':>5} {'S0(MPa)':>8} {'KE0(kt)':>9} {'Edep/KE0':>9} {'1st brk km':>10} {'dEdh pk km':>10} {'Erad_bol/Edep':>13} {'Brown-corridor':>14} {'peak L_Si (W)':>13} {'wall s':>7}")
for S0 in strengths:
    for d in diams:
        r = run(f"sweep d={d} S0={S0:.0e}", diameter=d, velocity=17e3, theta_deg=45., init_strength=S0, density=3300.)
        s = r.summary
        sweep[(d, S0)] = s
        e_brown = 8.2508 * max(s["E_rad_bol_kt"], 1e-12) ** 0.885          # Brown+2002 total energy from radiated energy
        corr = e_brown / s["KE0_kt"]
        sweep[(d, S0)]["brown_ratio"] = corr
        print(f"     {d:5.1f} {S0/1e6:8.2f} {s['KE0_kt']:9.3f} {s['E_dep_kt']/s['KE0_kt']:9.3f} {s['first_yield_alt_km']:10.1f} {s['alt_peak_dEdh_km']:10.1f} "
              f"{s['eta_bol_overall']:13.3f} {corr:14.2f} {s['peak_L_si_W']:13.2e} {s['wall_s']:7.2f}")
corr_all = np.array([v["brown_ratio"] for v in sweep.values()])
check("radiated fraction vs Brown+2002 relation (E_total = 8.25 E_rad^0.885)",
      np.all((corr_all > 0.5) & (corr_all < 2.0)), f"Brown-implied / model KE0 in [{corr_all.min():.2f}, {corr_all.max():.2f}] (population scatter ~x2)")
eta = np.array([v["eta_bol_overall"] for v in sweep.values()])
check("bolometric radiated fraction stays physical (<30%)", eta.max() < 0.30, f"range {eta.min()*100:.1f}-{eta.max()*100:.1f}%")
for S0 in strengths:   # larger body of the same strength: more energy deposited, never above KE0
    ok = all(sweep[(d, S0)]["E_dep_kt"] <= 1.01 * sweep[(d, S0)]["KE0_kt"] for d in diams)
    check(f"E_dep <= KE0 (+gravity work) for all sizes, S0={S0/1e6:.2f} MPa", ok, "")
# trends
alt_by_theta = {}
for th in (15, 30, 45, 75):
    alt_by_theta[th] = run(f"5 m, theta={th}", diameter=5, velocity=17e3, theta_deg=float(th), init_strength=0.7e6, density=3300.).summary["alt_peak_dEdh_km"]
check("steeper entry -> deeper peak deposition (5 m, 17 km/s)", all(alt_by_theta[a] >= alt_by_theta[b] for a, b in ((15, 30), (30, 45), (45, 75))),
      ", ".join(f"{k} deg: {v:.1f} km" for k, v in alt_by_theta.items()))
alt_by_v = {}
for vk in (12, 17, 25, 35):
    alt_by_v[vk] = run(f"5 m, v={vk}", diameter=5, velocity=vk * 1e3, theta_deg=45., init_strength=0.7e6, density=3300.).summary["first_yield_alt_km"]
check("faster entry -> higher first breakup (q = rho v^2/2)", all(alt_by_v[a] < alt_by_v[b] for a, b in ((12, 17), (17, 25), (25, 35))),
      ", ".join(f"{k} km/s: {v:.1f} km" for k, v in alt_by_v.items()))
for S0 in strengths[:2]:
    ys = [sweep[(d, S0)]["first_yield_alt_km"] for d in diams]
    check(f"first-breakup altitude set by q = S0 (size independent), S0={S0/1e6:.2f} MPa", max(ys) - min(ys) < 4.0,
          f"spread {max(ys)-min(ys):.1f} km")
for d in (2, 8, 20):
    pk = [sweep[(d, S0)]["alt_peak_dEdh_km"] for S0 in strengths]
    check(f"stronger body -> deeper peak deposition (d={d} m)", pk[0] > pk[1] > pk[2] - 0.6, ", ".join(f"{a:.1f}" for a in pk) + " km")

# =====================================================================  C/G. ledger and run time
print("\n=== C/G. conservation and run time (all %d runs) ===" % len(RUNS))
clo = np.array([abs(r.summary["energy_closure"]) for _, r in RUNS])
wall = np.array([r.summary["wall_s"] for _, r in RUNS])
check("kinetic-energy ledger closure on every run", clo.max() < 1e-9, f"max |closure| = {clo.max():.1e} of KE0")
check("run time well under 1 minute", wall.max() < 20, f"mean {wall.mean():.2f} s, max {wall.max():.2f} s over {len(RUNS)} runs")

# =====================================================================  figures
print("\n=== writing diagnostics figures ===")
pdf = PdfPages(os.path.join(OUT, "bolum_v2_diagnostics.pdf"))
fig, ax = plt.subplots(2, 3, figsize=(16, 9))
a = ax[0, 0]
a.plot(base.dEdh_dep * 1e3 / JKT, base.alt_bins / 1e3, label="model E deposited")
a.axhspan(27, 33, color="orange", alpha=.2, label="pub. half-energy range 27-33 km")
a.plot([80], [29.5], "r*", ms=14, label="pub. peak ~80 kt/km @ 29.5 km")
a.set_ylim(10, 70); a.set_xlabel("dE/dh (kt/km)"); a.set_ylabel("altitude (km)"); a.set_title("Chelyabinsk: energy deposition"); a.legend(fontsize=8)
a = ax[0, 1]
a.plot(base.L_si, base.alt / 1e3, label="model L_Si")
a.plot([3.4e14], [29.7], "r*", ms=14, label="pub. peak 2.7e13 W/sr x 4pi @ 29.7 km")
a.set_xscale("log"); a.set_xlim(1e11, 1e15); a.set_ylim(10, 70); a.set_xlabel("Si-band luminosity (W)"); a.set_ylabel("altitude (km)")
a.set_title("Chelyabinsk: Si-band light"); a.legend(fontsize=8)
a = ax[0, 2]
a.plot(tc.L_si, tc.alt / 1e3, label="model L_Si (TC3)")
for h, c in ((53, "gray"), (45, "gray"), (37.5, "r")):
    a.axhline(h, color=c, ls="--", lw=1)
a.text(a.get_xlim()[0] if a.get_xlim()[0] > 0 else 1e8, 53.5, "pub. flares 53 / 45 / 37.5 km", fontsize=8)
a.set_xscale("log"); a.set_xlim(1e8, 5e11); a.set_ylim(20, 70); a.set_xlabel("Si-band luminosity (W)"); a.set_ylabel("altitude (km)"); a.set_title("2008 TC3"); a.legend(fontsize=8)
a = ax[1, 0]
a.plot(ko.L_si, ko.alt / 1e3, label="model L_Si (Kosice)")
a.axhline(36, color="r", ls="--", lw=1, label="pub. max brightness 36 km"); a.axhline(17, color="gray", ls="--", lw=1, label="pub. end 17 km")
a.set_xscale("log"); a.set_xlim(1e7, 1e11); a.set_ylim(10, 70); a.set_xlabel("Si-band luminosity (W)"); a.set_ylabel("altitude (km)"); a.set_title("Kosice (not calibrated)"); a.legend(fontsize=8)
a = ax[1, 1]
ek = np.logspace(-3.5, 3.2, 100)
a.plot(ek, (8.2508 * ek ** 0.885)[:] * 0 + (ek / (8.2508 * ek ** 0.885)), "k-", label="Brown+2002 relation")
a.fill_between(ek, 0.5 * ek / (8.2508 * ek ** 0.885), 2 * ek / (8.2508 * ek ** 0.885), color="gray", alpha=.2, label="x2 band")
for S0, mk in zip(strengths, "osd"):
    xs = [sweep[(d, S0)]["E_dep_kt"] for d in diams]; ys = [sweep[(d, S0)]["eta_bol_overall"] for d in diams]
    a.plot(xs, ys, mk + "-", label=f"model S0={S0/1e6:.2f} MPa")
a.plot([cb["E_dep_kt"], ct["E_dep_kt"]], [cb["eta_bol_overall"], ct["eta_bol_overall"]], "r*", ms=12, label="Chelyabinsk, TC3 runs")
a.set_xscale("log"); a.set_xlabel("deposited energy (kt)"); a.set_ylabel("bolometric radiated fraction"); a.set_title("Radiated fraction vs energy"); a.legend(fontsize=8)
a = ax[1, 2]
for S0 in strengths:
    a.plot(diams, [sweep[(d, S0)]["alt_peak_dEdh_km"] for d in diams], "o-", label=f"peak dE/dh, S0={S0/1e6:.2f} MPa")
a.set_xscale("log"); a.set_xlabel("diameter (m)"); a.set_ylabel("altitude of peak deposition (km)"); a.set_title("Sweep: 17 km/s, 45 deg"); a.legend(fontsize=8)
plt.tight_layout(); pdf.savefig(fig); fig.savefig(os.path.join(OUT, "diag_events.png"), dpi=90); plt.close(fig)

fig, ax = plt.subplots(1, 3, figsize=(16, 4.5))
h = np.linspace(0, 110, 300)
ax[0].semilogy(h, ATM.rho(h * 1e3), label="bolum v2")
ax[0].semilogy(list(std) + [90, 100, 110], list(std.values()) + [3.416e-6, 5.604e-7, 9.708e-8], "ro", label="US Std 1976")
ax[0].set_xlabel("altitude (km)"); ax[0].set_ylabel("density (kg/m^3)"); ax[0].set_title("Atmosphere"); ax[0].legend()
v = np.linspace(0.5e3, 40e3, 300)
raw = lambda vk: np.exp(0.567 - 10.307 * np.log(vk) + 9.781 * np.log(vk) ** 2 - 3.0414 * np.log(vk) ** 3 + 0.3213 * np.log(vk) ** 4 + 0.347 * np.tanh(0.38 * np.log(1e3))) / 100
ax[1].semilogy(v / 1e3, raw(v / 1e3), "r--", label="original (raw polynomial)")
ax[1].semilogy(v / 1e3, np.maximum(tau_fit(v, 1e3), 1e-6), label="v2 (tapered)")
ax[1].set_xlabel("speed (km/s)"); ax[1].set_ylabel("tau (1 t body)"); ax[1].set_title("Luminous efficiency fit"); ax[1].legend()
for dt, s in rows.items():
    pass
dts_ = [0.012, 0.006, 0.003, 0.0015]
ER = [rows[0.012]["E_rad_bol_kt"], rows[0.006]["E_rad_bol_kt"], cb["E_rad_bol_kt"], rows[0.0015]["E_rad_bol_kt"]]
ax[2].semilogx(dts_, ER, "o-"); ax[2].set_xlabel("max time step (s)"); ax[2].set_ylabel("E_rad bolometric (kt)"); ax[2].set_title("Time-step convergence (Chelyabinsk)")
plt.tight_layout(); pdf.savefig(fig); fig.savefig(os.path.join(OUT, "diag_components.png"), dpi=90); plt.close(fig)
pdf.close()

npass = sum(1 for _, ok, _ in REPORT if ok)
print(f"\nSUMMARY: {npass}/{len(REPORT)} checks passed")
for n, ok, d in REPORT:
    if not ok:
        print(f"  FAILED: {n}: {d}")
with open(os.path.join(OUT, "bolum_v2_diagnostics_report.txt"), "w") as fh:
    for n, ok, d in REPORT:
        fh.write(f"{'PASS' if ok else 'FAIL'}  {n}: {d}\n")
    fh.write(f"\n{npass}/{len(REPORT)} checks passed\n")
