#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
bolum_diagnostics.py (for bolum_block_v2.py v2.3) -- verification, convergence and validation suite

  A. components        atmosphere vs US Std 1976, tau(v,m) behaviour, Si-band fraction, loading convention
  B. verification      engine vs an independent RK45 integration of the same single-body ODEs
  C. conservation      kinetic-energy ledger closure on every run
  D. convergence       time step, group cap, children per breakup
  E. events            Chelyabinsk 2013, 2008 TC3, Kosice 2010 vs published values  (one STANDARDIZED figure per event)
  F. sweeps            0.5-25 m, three strengths, entry-angle and speed trends, radiated fraction vs Brown et al. (2002)
  G. regression/CLI    v2.3 reproduces v2.0 numbers (after the strength-convention change); command-line smoke test
  H. run time
  I. light-curve feature analysis (new in 2.3): ablation-onset altitude vs analytic, attribution, robustness

Chelyabinsk and TC3 were used to calibrate the GLOBAL parameters (c_spread, f_split, alpha, dust, eta_thick, lum_r_half),
so they are calibration checks.  Kosice and the sweeps were not used for calibration.
All strengths are stagnation pressures rho*v^2 (the published convention); the strengths calibrated in v2.0 under
0.5*rho*v^2 appear here doubled.

Usage:  python bolum_diagnostics.py [outdir]
"""
import sys
import os
import math
import tempfile
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import bolum_block_v2 as B
from bolum_block_v2 import (BolumParams, run_bolide, ATM, tau_fit, band_fraction, RE, GM, JKT, PI, SIGMA_SB,
                            plot_event_summary, plot_feature_guide, model_metric_rows)

OUT = sys.argv[1] if len(sys.argv) > 1 else "."
os.makedirs(OUT, exist_ok=True)
REPORT, RUNS = [], []


def check(name, ok, detail):
    REPORT.append((name, bool(ok), detail))
    print(("PASS  " if ok else "FAIL  ") + f"{name}: {detail}")
    return bool(ok)


def run(tag, **kw):
    r = run_bolide(BolumParams(**kw))
    RUNS.append((tag, r))
    return r


def rowcheck(rows, label, check_name, model_str, pub_str, ok, detail=None):
    """Record a check AND add the same line (model | published | status) to the event's metrics table."""
    check(check_name, ok, detail or f"{model_str} vs {pub_str}")
    rows.append((label, model_str, pub_str, bool(ok)))


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
      f"max {t_hi.max()*100:.2f}%, tau(1 km/s) = {float(tau_fit(1e3, 1e3)):.1e} (raw polynomial: ~1.8%)")
t20 = float(tau_fit(20e3, 1e3))
check("tau(20 km/s, 1 t) matches the published fit (Borovicka+2020 Eq. 3)", abs(t20 - 0.0682) < 0.004, f"tau = {t20*100:.2f}%")
f6 = float(band_fraction(6000.)[0])
ftot = float(band_fraction(6000., 0.05, 100.)[0])
check("Si-band blackbody fraction (0.4-1.1 um, 6000 K)", abs(f6 - 0.643) < 0.01,
      f"f_si(6000 K) = {f6:.3f}; f_si(4500 K) = {float(band_fraction(4500.)[0]):.3f}; f_si(3500 K) = {float(band_fraction(3500.)[0]):.3f}")
check("band integral sanity (0.05-100 um -> 1)", abs(ftot - 1) < 0.02, f"{ftot:.3f}")
x_ = np.array([1.0, 7.3e5, 3.3e4])
check("Si <-> bolometric conversion round trip", np.allclose(B.si_from_bolometric(B.bolometric_from_si(x_)), x_), "L_si -> L_bol -> L_si")

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

print("\n=== B2. small bolides recover the original Si-band efficiency ===")
res = run("0.3 m intact", diameter=0.3, velocity=20e3, theta_deg=45., init_strength=1e9, density=3300.)
dts = np.diff(np.concatenate([[0.0], res.t]))
E_tau = float(np.sum(tau_fit(res.v, np.maximum(res.mass, 1e-3)) * res.P_dep * dts))
ratio = res.summary["E_rad_si_J"] / E_tau
check("0.3 m: integrated L_Si / (tau_fit x deposited energy)", 0.93 < ratio < 1.12, f"ratio {ratio:.3f}")

# =====================================================================  D. convergence
print("\n=== D. convergence (Chelyabinsk-like) ===")
CH = dict(diameter=19.6, velocity=19.16e3, theta_deg=18.3, init_strength=1.4e6, density=3300.)
base = run("chely", **CH)
rows_dt = {}
for dt in (0.012, 0.006, 0.0015):
    rows_dt[dt] = run(f"chely dt={dt}", dt=dt, **CH).summary
for dt, s in rows_dt.items():
    d1 = abs(s["E_rad_bol_kt"] / base.summary["E_rad_bol_kt"] - 1)
    d2 = abs(s["peak_L_si_W"] / base.summary["peak_L_si_W"] - 1)
    print(f"     dt_max={dt:<7} steps={s['steps']:6d}  dE_rad={d1*100:.2f}%  dL_peak={d2*100:.2f}%")
s = rows_dt[0.0015]
check("time-step convergence (3 ms vs 1.5 ms)",
      abs(s["E_rad_bol_kt"] / base.summary["E_rad_bol_kt"] - 1) < 0.01 and abs(s["peak_L_si_W"] / base.summary["peak_L_si_W"] - 1) < 0.01,
      "E_rad and peak L change < 1%")
s2 = run("chely NG=512", max_groups=512, **CH).summary
check("group cap not binding (128 vs 512)", abs(s2["E_rad_bol_kt"] / base.summary["E_rad_bol_kt"] - 1) < 0.005,
      f"E_rad change {abs(s2['E_rad_bol_kt']/base.summary['E_rad_bol_kt']-1)*100:.2f}%")
for n in (3, 6):
    sn = run(f"chely n={n}", n_fragments=n, **CH).summary
    print(f"     n_fragments={n}: E_rad(bol) {sn['E_rad_bol_kt']:.1f} kt, peak dE/dh {sn['alt_peak_dEdh_km']:.1f} km")

# =====================================================================  E. events (standardized figures)
print("\n=== E. events ===")
EVENT_FIGS = []

# ---- Chelyabinsk
cb = base.summary
rows_ch = []
lo, hi = cb["half_range_km"]
rowcheck(rows_ch, "KE0 (kt)", "Chelyabinsk: KE0 within published range", f"{cb['KE0_kt']:.0f}", "430-640", 430 <= cb["KE0_kt"] <= 640,
         f"model KE0 = {cb['KE0_kt']:.0f} kt for 1.3e7 kg (published 430-640 kt, Brown+2013; Popova+2013)")
rowcheck(rows_ch, "E_dep / KE0", "Chelyabinsk: deposited energy ~ all of KE0", f"{cb['E_dep_kt']/cb['KE0_kt']:.3f}", ">0.98", cb["E_dep_kt"] / cb["KE0_kt"] > 0.98)
rowcheck(rows_ch, "peak dE/dh alt (km)", "Chelyabinsk: peak dE/dh altitude", f"{cb['alt_peak_dEdh_km']:.1f}", "29.5", abs(cb["alt_peak_dEdh_km"] - 29.5) < 1.5)
rowcheck(rows_ch, "half-peak range (km)", "Chelyabinsk: half-energy range", f"{lo:.1f}-{hi:.1f}", "27-33", abs(lo - 27) < 1.5 and abs(hi - 33) < 1.5)
rowcheck(rows_ch, "peak dE/dh (kt/km)", "Chelyabinsk: peak dE/dh", f"{cb['peak_dEdh_kt_per_km']:.0f}", "~80", 0.6 < cb["peak_dEdh_kt_per_km"] / 80 < 1.4)
rowcheck(rows_ch, "E_rad bol (J)", "Chelyabinsk: radiated energy", f"{cb['E_rad_bol_J']:.2e}", "3.75e14", abs(math.log(cb["E_rad_bol_J"] / 3.75e14)) < math.log(1.25))
rowcheck(rows_ch, "peak I_Si (W/sr)", "Chelyabinsk: peak Si-band intensity", f"{cb['peak_I_si_W_sr']:.2e}", "2.7e13", abs(math.log(cb["peak_I_si_W_sr"] / 2.7e13)) < math.log(1.4))
rowcheck(rows_ch, "alt of peak light (km)", "Chelyabinsk: altitude of peak light", f"{cb['alt_peak_L_km']:.1f}", "29.7", abs(cb["alt_peak_L_km"] - 29.7) < 2.5)
rows_ch += [r for r in model_metric_rows(base) if r[0] in ("first yield alt (km)", "breakups", "luminous end alt (km)", "energy-ledger closure")]
REF_CH = dict(peak_L_W=2.7e13 * 4 * PI, peak_L_alt_km=29.7, peak_dEdh=(29.5, 80.0), half_range_km=(27., 33.))
EVENT_FIGS.append((base, REF_CH, rows_ch, "Chelyabinsk 2013  (calibration event)"))
check("first-yield loading equals S0 (rho v^2 convention)", abs(cb["first_yield_q_MPa"] / 1.4 - 1) < 0.05, f"first yield at {cb['first_yield_q_MPa']:.3f} MPa for S0 = 1.4 MPa")

# ---- TC3
tc = run("TC3", diameter=4.1, velocity=12.8e3, theta_deg=19.5, init_strength=0.3e6, density=2300.)
ct = tc.summary
y1 = next(e for e in tc.events if e["kind"] == "yield")
rows_tc = []
rowcheck(rows_tc, "KE0 (kt)", "TC3: KE0 within published range", f"{ct['KE0_kt']:.2f}", "1.1-2.1", 1.0 <= ct["KE0_kt"] <= 2.2)
rowcheck(rows_tc, "first breakup alt (km)", "TC3: first fragmentation flare", f"{y1['alt']/1e3:.1f}", "45-53", 43 < y1["alt"] / 1e3 < 54,
         f"first breakup {y1['alt']/1e3:.1f} km at loading {y1['q']/1e6:.2f} MPa (published flares 53/45 km)")
rowcheck(rows_tc, "alt of peak light (km)", "TC3: peak light altitude", f"{ct['alt_peak_L_km']:.1f}", "37.5", abs(ct["alt_peak_L_km"] - 37.5) < 2.5)
rowcheck(rows_tc, "peak dE/dh alt (km)", "TC3: peak energy deposition", f"{ct['alt_peak_dEdh_km']:.1f}", "~37", abs(ct["alt_peak_dEdh_km"] - 37.0) < 3.0)
rowcheck(rows_tc, "E_rad bol (kt)", "TC3: radiated energy (factor 1.6)", f"{ct['E_rad_bol_kt']:.3f}", "~0.1", 1 / 1.6 < ct["E_rad_bol_kt"] / 0.1 < 1.6)
rows_tc += [r for r in model_metric_rows(tc) if r[0] in ("peak dE/dh (kt/km @ km)", "peak L_Si (W)", "luminous end alt (km)", "energy-ledger closure")]
REF_TC = dict(flares_km=[53, 45, 37.5], peak_L_alt_km=37.5)
EVENT_FIGS.append((tc, REF_TC, rows_tc, "2008 TC3  (calibration event)"))

# ---- Kosice
ko = run("Kosice", diameter=1.25, velocity=15e3, theta_deg=60., init_strength=0.1e6, strength_secondary=1.0e6, density=3400.)
ck = ko.summary
rows_ko = []
rowcheck(rows_ko, "first breakup alt (km)", "Kosice: first breakup altitude", f"{ck['first_yield_alt_km']:.1f}", "55-57", abs(ck["first_yield_alt_km"] - 56) < 6)
rowcheck(rows_ko, "alt of peak light (km)", "Kosice: altitude of max brightness", f"{ck['alt_peak_L_km']:.1f}", "36", abs(ck["alt_peak_L_km"] - 36) < 5)
rowcheck(rows_ko, "luminous end alt (km)", "Kosice: luminous end height", f"{ck['alt_end_lum_km']:.1f}", "17", abs(ck["alt_end_lum_km"] - 17) < 6,
         f"{ck['alt_end_lum_km']:.1f} km (published 17 km) - KNOWN WEAKNESS: small strong fragments are over-decelerated by the pancake")
rows_ko += [r for r in model_metric_rows(ko) if r[0] in ("KE0 (kt)", "peak L_Si (W)", "E_rad bolometric (kt)", "energy-ledger closure")]
REF_KO = dict(peak_L_alt_km=36.0, end_km=17.0, first_break_km=(55., 57.))
EVENT_FIGS.append((ko, REF_KO, rows_ko, "Kosice 2010  (NOT used in calibration)"))

# =====================================================================  F. sweeps
print("\n=== F. size sweep (v = 17 km/s, 45 deg, 3300 kg/m^3; strengths are rho*v^2) ===")
diams = [0.5, 1, 2, 4, 8, 12, 20, 25]
strengths = [0.3e6, 1.4e6, 4.0e6]
sweep = {}
print(f"     {'d(m)':>5} {'S0(MPa)':>8} {'KE0(kt)':>9} {'Edep/KE0':>9} {'1st brk km':>10} {'dEdh pk km':>10} {'Erad_bol/Edep':>13} {'Brown-ratio':>11} {'peak L_Si (W)':>13} {'wall s':>7}")
for S0 in strengths:
    for d in diams:
        r = run(f"sweep d={d} S0={S0:.0e}", diameter=d, velocity=17e3, theta_deg=45., init_strength=S0, density=3300.)
        s = r.summary
        sweep[(d, S0)] = s
        corr = 8.2508 * max(s["E_rad_bol_kt"], 1e-12) ** 0.885 / s["KE0_kt"]
        s["brown_ratio"] = corr
        print(f"     {d:5.1f} {S0/1e6:8.2f} {s['KE0_kt']:9.3f} {s['E_dep_kt']/s['KE0_kt']:9.3f} {s['first_yield_alt_km']:10.1f} {s['alt_peak_dEdh_km']:10.1f} "
              f"{s['eta_bol_overall']:13.3f} {corr:11.2f} {s['peak_L_si_W']:13.2e} {s['wall_s']:7.2f}")
corr_all = np.array([v["brown_ratio"] for v in sweep.values()])
check("radiated fraction vs Brown+2002 relation (E_total = 8.25 E_rad^0.885)",
      np.all((corr_all > 0.5) & (corr_all < 2.0)), f"Brown-implied / model KE0 in [{corr_all.min():.2f}, {corr_all.max():.2f}] (population scatter ~x2)")
eta = np.array([v["eta_bol_overall"] for v in sweep.values()])
check("bolometric radiated fraction stays physical (<30%)", eta.max() < 0.30, f"range {eta.min()*100:.1f}-{eta.max()*100:.1f}%")
for S0 in strengths:
    ok = all(sweep[(d, S0)]["E_dep_kt"] <= 1.01 * sweep[(d, S0)]["KE0_kt"] for d in diams)
    check(f"E_dep <= KE0 (+gravity work) for all sizes, S0={S0/1e6:.2f} MPa", ok, "")
alt_by_theta = {}
for th in (15, 30, 45, 75):
    alt_by_theta[th] = run(f"5 m, theta={th}", diameter=5, velocity=17e3, theta_deg=float(th), init_strength=1.4e6, density=3300.).summary["alt_peak_dEdh_km"]
check("steeper entry -> deeper peak deposition (5 m, 17 km/s)", all(alt_by_theta[a] >= alt_by_theta[b] for a, b in ((15, 30), (30, 45), (45, 75))),
      ", ".join(f"{k} deg: {v:.1f} km" for k, v in alt_by_theta.items()))
alt_by_v = {}
for vk in (12, 17, 25, 35):
    alt_by_v[vk] = run(f"5 m, v={vk}", diameter=5, velocity=vk * 1e3, theta_deg=45., init_strength=1.4e6, density=3300.).summary["first_yield_alt_km"]
check("faster entry -> higher first breakup (loading = rho v^2)", all(alt_by_v[a] < alt_by_v[b] for a, b in ((12, 17), (17, 25), (25, 35))),
      ", ".join(f"{k} km/s: {v:.1f} km" for k, v in alt_by_v.items()))
for S0 in strengths[:2]:
    ys = [sweep[(d, S0)]["first_yield_alt_km"] for d in diams]
    check(f"first-breakup altitude set by rho v^2 = S0 (size independent), S0={S0/1e6:.2f} MPa", max(ys) - min(ys) < 4.0, f"spread {max(ys)-min(ys):.1f} km")
for d in (2, 8, 20):
    pk = [sweep[(d, S0)]["alt_peak_dEdh_km"] for S0 in strengths]
    check(f"stronger body -> deeper peak deposition (d={d} m)", pk[0] > pk[1] > pk[2] - 0.6, ", ".join(f"{a:.1f}" for a in pk) + " km")

# =====================================================================  G. regression and CLI smoke test
print("\n=== G. regression vs v2.0 (strengths doubled for the convention change) and CLI smoke test ===")
V20 = {"Chelyabinsk": (93.04, 2.749e14, cb), "TC3": (0.1393, 1.844e11, ct), "Kosice": (None, 2.831e10, ck)}
for name, (er, pk, sm) in V20.items():
    ok = abs(sm["peak_L_si_W"] / pk - 1) < 0.005 and (er is None or abs(sm["E_rad_bol_kt"] / er - 1) < 0.005)
    check(f"v2.3 reproduces v2.0: {name}", ok, f"peak L {sm['peak_L_si_W']:.3e} (v2.0 {pk:.3e})" + (f", E_rad {sm['E_rad_bol_kt']:.3f} kt (v2.0 {er})" if er else ""))
with tempfile.TemporaryDirectory() as td:
    d_ = os.path.join(td, "smoke")
    os.makedirs(d_)
    with open(os.path.join(d_, "bolum_in.dat"), "w") as fh:
        fh.write("smoke test\nformat\ndiameter = 3.0\nvelocity = 18000\ntheta_deg = 40\ninit_strength = 1.0e6\ndensity = 3000\nq_factor = 0.5\n")
    import io, contextlib
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        B.main(["smoke", "--base", td])
    files = set(os.listdir(d_))
    need = {"run_report.txt", "bolide_plots.pdf", "summary.dat", "dEdh.dat", "events.dat", "power_vs_time.dat", "power_vs_alt.dat", "bolide_out.dat", "init_params.dat"}
    check("CLI run writes all outputs and a full report", need <= files and "ledger closure" in buf.getvalue() and "q_factor" in buf.getvalue(),
          f"{len(files)} files; legacy q_factor key handled")


# =====================================================================  I. light-curve feature analysis
print("\n=== I. light-curve feature analysis ===")
CL = run("Cleveland-like", diameter=1.8, velocity=17.5e3, theta_deg=31.9, init_strength=0.33e6, density=3300.)
kinds = lambda r_: [f_["kind"] for f_ in r_.features]
for nm, r_ in (("Cleveland-like", CL), ("Chelyabinsk", base), ("TC3", tc), ("Kosice", ko)):
    p_ = r_.params
    q_loss = p_.emissivity * SIGMA_SB * p_.T_ablate ** 4
    rho_on = q_loss / (p_.Lambda * 0.5 * p_.velocity ** 3)
    hs = np.linspace(40e3, 115e3, 7501)
    alt_on = hs[np.argmin(abs(ATM.rho(hs) - rho_on))] / 1e3
    got = r_.summary["alt_ablation_onset_km"]
    check(f"ablation-onset altitude equals the analytic threshold [{nm}]", abs(got - alt_on) < 1.0, f"found {got:.1f} km, analytic {alt_on:.1f} km")
    ok = all(k_ in kinds(r_) for k_ in ("ablation_onset", "first_yield", "first_breakup", "peak", "rolloff"))
    check(f"light-curve features found [{nm}]", ok, ", ".join(sorted(set(kinds(r_)))))
    fl = [f_ for f_ in r_.features if f_["kind"] == "flare"]
    check(f"first flare is attributed to a breakup/yield [{nm}]", len(fl) > 0 and ("split" in fl[0]["note"] or "yield" in fl[0]["note"]),
          f"{len(fl)} flares; first: t = {fl[0]['t']-r_.params.tstart:.2f} s" if fl else "no flare")
    st_ = r_.summary
    check(f"radiated-energy timing is ordered [{nm}]", st_["lc_t10"] < st_["lc_t50"] < st_["lc_t90"] and 0 <= st_["lc_E_prefrag_frac"] <= 1,
          f"t10/50/90 = {st_['lc_t10']-r_.params.tstart:.2f}/{st_['lc_t50']-r_.params.tstart:.2f}/{st_['lc_t90']-r_.params.tstart:.2f} s, pre-yield fraction {st_['lc_E_prefrag_frac']*100:.1f}%")
    pk_f = next(f_ for f_ in r_.features if f_["kind"] == "peak")
    check(f"'peak' feature is the light-curve maximum [{nm}]", abs(pk_f["L_si"] / r_.L_si.max() - 1) < 1e-6, "")
tail_ok = all(f_["t"] > next(g["t"] for g in CL.features if g["kind"] == "peak") for f_ in CL.features if f_["kind"] in ("rolloff", "tail"))
check("roll-off and tail features come after the peak", tail_ok, "")
intact = run("intact 0.5 m (no breakup)", diameter=0.5, velocity=20e3, theta_deg=45., init_strength=1e9, density=3300.)
ki = kinds(intact)
check("intact body: no yield/flare features, still gets onset/peak/roll-off", ("first_yield" not in ki and "flare" not in ki and "ablation_onset" in ki and "peak" in ki),
      ", ".join(sorted(set(ki))))
check("drag+ablation channels sum to the deposited power", np.allclose(CL.P_drag + CL.P_abl, CL.P_dep, rtol=1e-9), "P_drag + P_abl = P_dep")

# =====================================================================  C/H. ledger and run time
print("\n=== C/H. conservation and run time (all %d runs) ===" % len(RUNS))
clo = np.array([abs(r.summary["energy_closure"]) for _, r in RUNS])
wall = np.array([r.summary["wall_s"] for _, r in RUNS])
check("kinetic-energy ledger closure on every run", clo.max() < 1e-9, f"max |closure| = {clo.max():.1e} of KE0")
check("run time well under 1 minute", wall.max() < 20, f"mean {wall.mean():.2f} s, max {wall.max():.2f} s over {len(RUNS)} runs")

# =====================================================================  figures
print("\n=== writing diagnostics figures ===")
pdf = PdfPages(os.path.join(OUT, "bolum_v2_diagnostics.pdf"))
for k, (res_, ref_, rows_, title_) in enumerate(EVENT_FIGS):
    fig = plot_event_summary(res_, ref=ref_, rows=rows_, title=title_)
    pdf.savefig(fig)
    fig.savefig(os.path.join(OUT, f"diag_event_{k+1}_{title_.split()[0].lower()}.png"), dpi=80)
    plt.close(fig)
    fig = plot_feature_guide(res_, title=title_ + " - light-curve feature guide")
    pdf.savefig(fig)
    fig.savefig(os.path.join(OUT, f"diag_guide_{k+1}_{title_.split()[0].lower()}.png"), dpi=80)
    plt.close(fig)
fig = plot_event_summary(CL, title="Cleveland-like (illustrative inputs, no published comparison)")
pdf.savefig(fig); fig.savefig(os.path.join(OUT, "diag_event_cleveland.png"), dpi=80); plt.close(fig)
fig = plot_feature_guide(CL, title="Cleveland-like - light-curve feature guide")
pdf.savefig(fig); fig.savefig(os.path.join(OUT, "diag_guide_cleveland.png"), dpi=80); plt.close(fig)

fig, ax = plt.subplots(2, 2, figsize=(15, 10))
fig.subplots_adjust(left=.07, right=.98, top=.91, bottom=.07, wspace=.25, hspace=.30)
a = ax[0, 0]
ek = np.logspace(-3.5, 3.2, 100)
brown = ek / (8.2508 * ek ** 0.885)
a.plot(ek, brown, "k-", label="Brown+2002 relation"); a.fill_between(ek, 0.5 * brown, 2 * brown, color="gray", alpha=.2, label="x2 band")
for S0, mk in zip(strengths, "osd"):
    a.plot([sweep[(d, S0)]["E_dep_kt"] for d in diams], [sweep[(d, S0)]["eta_bol_overall"] for d in diams], mk + "-", label=f"model, S0 = {S0/1e6:.1f} MPa")
a.plot([cb["E_dep_kt"], ct["E_dep_kt"]], [cb["eta_bol_overall"], ct["eta_bol_overall"]], "*", color="crimson", ms=13, label="Chelyabinsk, TC3 runs")
a.set_xscale("log"); a.set_xlabel("deposited energy (kt)"); a.set_ylabel("bolometric radiated fraction"); a.set_title("(a) radiated fraction vs energy"); a.legend(fontsize=8); a.grid(True, alpha=.3)
a = ax[0, 1]
for S0 in strengths:
    a.plot(diams, [sweep[(d, S0)]["alt_peak_dEdh_km"] for d in diams], "o-", label=f"S0 = {S0/1e6:.1f} MPa")
a.set_xscale("log"); a.set_xlabel("diameter (m)"); a.set_ylabel("altitude of peak deposition (km)"); a.set_title("(b) sweep: 17 km/s, 45 deg"); a.legend(fontsize=8); a.grid(True, alpha=.3)
a = ax[1, 0]
Tg = np.linspace(2500, 10000, 100)
a.plot(Tg, band_fraction(Tg), color=B.C_SI)
a.axvline(6000, color="gray", ls=":"); a.set_xlabel("blackbody temperature (K)"); a.set_ylabel("f_si: fraction of power in 0.4-1.1 um")
a.set_title("(c) Si-band fraction (L_bol = L_si / f_si)"); a.grid(True, alpha=.3)
a = ax[1, 1]
dts_ = [0.012, 0.006, 0.003, 0.0015]
ER = [rows_dt[0.012]["E_rad_bol_kt"], rows_dt[0.006]["E_rad_bol_kt"], cb["E_rad_bol_kt"], rows_dt[0.0015]["E_rad_bol_kt"]]
a.semilogx(dts_, ER, "o-"); a.set_xlabel("max time step (s)"); a.set_ylabel("E_rad bolometric (kt)"); a.set_title("(d) time-step convergence (Chelyabinsk)"); a.grid(True, alpha=.3)
fig.text(.01, .985, "UNCLASSIFIED", fontsize=9, color="green", va="top"); fig.text(.99, .012, "UNCLASSIFIED", fontsize=9, color="green", ha="right")
fig.suptitle("Sweeps and components", fontsize=13)
pdf.savefig(fig); fig.savefig(os.path.join(OUT, "diag_sweeps.png"), dpi=80); plt.close(fig)

fig, ax = plt.subplots(1, 3, figsize=(17, 4.8))
h = np.linspace(0, 110, 300)
ax[0].semilogy(h, ATM.rho(h * 1e3), label="bolum v2.3.1")
ax[0].semilogy(list(std) + [90, 100, 110], list(std.values()) + [3.416e-6, 5.604e-7, 9.708e-8], "ro", label="US Std 1976")
ax[0].set_xlabel("altitude (km)"); ax[0].set_ylabel("density (kg/m^3)"); ax[0].set_title("Atmosphere"); ax[0].legend(); ax[0].grid(True, alpha=.3)
v = np.linspace(0.5e3, 40e3, 300)
raw = lambda vk: np.exp(0.567 - 10.307 * np.log(vk) + 9.781 * np.log(vk) ** 2 - 3.0414 * np.log(vk) ** 3 + 0.3213 * np.log(vk) ** 4 + 0.347 * np.tanh(0.38 * np.log(1e3))) / 100
ax[1].semilogy(v / 1e3, raw(v / 1e3), "r--", label="raw published polynomial")
ax[1].semilogy(v / 1e3, np.maximum(tau_fit(v, 1e3), 1e-6), label="v2.3 (clamped + tapered)")
ax[1].set_xlabel("speed (km/s)"); ax[1].set_ylabel("tau (1 t body)"); ax[1].set_title("Luminous efficiency fit"); ax[1].legend(); ax[1].grid(True, alpha=.3)
ax[2].plot(alt_by_theta.keys(), alt_by_theta.values(), "o-"); ax[2].set_xlabel("entry angle (deg)"); ax[2].set_ylabel("altitude of peak deposition (km)")
ax[2].set_title("Entry-angle trend (5 m, 17 km/s)"); ax[2].grid(True, alpha=.3)
fig.text(.01, .985, "UNCLASSIFIED", fontsize=9, color="green", va="top"); fig.text(.99, .012, "UNCLASSIFIED", fontsize=9, color="green", ha="right")
plt.tight_layout(); pdf.savefig(fig); fig.savefig(os.path.join(OUT, "diag_components.png"), dpi=80); plt.close(fig)
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
