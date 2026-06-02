import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.integrate import trapz

# ========================== PARAMETERS ==========================
# Bolide / geometry assumptions (CUSTOMIZE THESE)
R_GLM = 3.6e7          # meters, approximate GEO distance (35,000–42,000 km)
APERTURE_AREA = 0.0097 # m², GLM effective collecting area (approx)
GLM_BAND_CENTER = 777.4e-9  # meters
GLM_BAND_WIDTH = 1.1e-9     # FWHM approx
T_BB = 6000.0               # K, blackbody temperature for continuum (common assumption)
NARROW_FRACTION = None      # Set to None to compute; or hardcode e.g. 0.001–0.003 (~1/1018 is the lit value)

# Atmospheric model (very simple)
# need to improve this with a continuous model of atmospheric attenuation
def atmospheric_transmittance(h_bolide_km):
    """Rough vertical transmittance from h to space at ~777 nm"""
    if h_bolide_km >= 100:
        return 0.999
    elif h_bolide_km >= 60:
        tau = 0.002
    elif h_bolide_km >= 40:
        tau = 0.008
    elif h_bolide_km >= 30:
        tau = 0.018
    else:
        tau = 0.025
    # Simple slant (assume μ ≈ 0.8 average)
    return np.exp(-tau / 0.8)

# Blackbody narrowband fraction (numerical)
def bb_narrow_fraction(T=6000.0, center_nm=777.4, width_nm=1.1):
    from scipy.integrate import quad
    lam_c = center_nm * 1e-9
    dlam = width_nm * 1e-9 / 2
    h, c, k = 6.62607015e-34, 2.99792458e8, 1.380649e-23

    def B(lam, T):
        return (2 * h * c**2 / lam**5) / (np.exp(h*c/(lam*k*T)) - 1)

    integral, _ = quad(lambda lam: B(lam, T), lam_c - dlam, lam_c + dlam)
    total = (5.670374419e-8 * T**4) / np.pi   # W/m²/sr
    return integral / total

if NARROW_FRACTION is None:
    NARROW_FRACTION = bb_narrow_fraction(T_BB)

print(f"Using narrowband fraction: {NARROW_FRACTION:.5f}")

# ========================== LOAD DATA ==========================
# Your bolide output: CSV with columns ['time_s', 'broadband_power_W_sr']
bolide_df = pd.read_csv('your_bolide_output.csv')   # ← CHANGE FILENAME
t_bol = bolide_df['time_s'].values
P_broadband_sr = bolide_df['broadband_power_W_sr'].values   # W/sr (isotropic assumed)

# GLM data: assume CSV with ['time_s', 'joules_per_sample']
glm_df = pd.read_csv('glm_lightcurve.csv')   # ← CHANGE
t_glm = glm_df['time_s'].values
E_glm = glm_df['joules_per_sample'].values   # Joules received at aperture per ~2-3 ms sample

# ========================== CONVERSIONS ==========================
dt_glm = np.median(np.diff(t_glm))          # should be ~0.002–0.003 s
P_glm_received = E_glm / dt_glm             # Watts received at aperture

# Irradiance at GLM (W/m²)
irradiance_glm = P_glm_received / APERTURE_AREA

# ========================== SIMULATE GLM FROM BOLIDE ==========================
# Assume same time grid or interpolate
# For simplicity, use bolide time grid
h_bolide = 60.0  # km — change or make array if you have altitude vs time
T_atm = atmospheric_transmittance(h_bolide)

P_narrow_sr = P_broadband_sr * NARROW_FRACTION          # W/sr narrowband
P_received = (P_narrow_sr * T_atm) / (R_GLM ** 2)       # W at sensor (isotropic)
E_simulated = P_received * dt_glm                       # Joules per sample (to match GLM)

# Convert GLM to irradiance for overlay on your broadband plot if desired
# (your code outputs power at source; here we show received irradiance)

# ========================== PLOTTING ==========================
fig, ax1 = plt.subplots(figsize=(10, 6))

# Left: GLM-style (received Joules or irradiance)
ax1.plot(t_glm, E_glm, 'b-', label='GLM Observed (Joules/sample)', linewidth=2)
ax1.plot(t_bol, E_simulated, 'r--', label='Simulated from your model', linewidth=2)
ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Energy per sample (Joules) at GLM aperture', color='b')
ax1.tick_params(axis='y', labelcolor='b')

# Right: Irradiance at sensor (W/m²) to compare with your power curve style
ax2 = ax1.twinx()
ax2.plot(t_glm, irradiance_glm, 'b:', label='GLM Irradiance (W/m²)', alpha=0.7)
ax2.set_ylabel('Irradiance at GLM (W/m²)', color='blue')

ax1.legend(loc='upper left')
ax2.legend(loc='upper right')
plt.title('GLM Bolide Light Curve vs Your Model (Narrowband + Geometry)')
plt.grid(True)
plt.tight_layout()
plt.show()

# Optional: save simulated curve
sim_df = pd.DataFrame({'time_s': t_bol, 'E_sim_J': E_simulated, 'irradiance_W_m2': P_received / APERTURE_AREA})
sim_df.to_csv('simulated_glm.csv', index=False)