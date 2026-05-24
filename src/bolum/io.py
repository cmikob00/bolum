'''
bolum I/O file
'''

# reading in input deck
def rdinput(inputpath):

    params = {}

    with open(inputpath, "r") as f:
        lines = f.readlines()[2:]

    for line in lines:

        # skip blank/comment lines
        if "=" not in line:
            continue

        key, rest = line.split("=", 1)

        key = key.strip()
        value_str = rest.split("#")[0].strip()

        try:
            value = float(value_str)
        except ValueError:
            continue

        params[key] = value

    #
    # REQUIRED PARAMETERS
    #
    try:
        diameter      = params["diameter"]
        velocity      = params["velocity"]
        theta_deg     = params["theta_deg"]
        init_strength = params["init_strength"]
        density       = params["density"]

    except KeyError as e:
        raise ValueError(f"Missing required input parameter: {e}")

    #
    # OPTIONAL PARAMETERS WITH DEFAULTS
    #
    strength_secondary  = params.get("strength_secondary", 1.e7)
    strength_tertiary   = params.get("strength_tertiary", 2.e7)
    strength_quaternary = params.get("strength_quaternary", 3.e7)
    porosity            = params.get("porosity", 0.1)
    Cd                  = params.get("Cd", 1.0)
    L_ablation          = params.get("L_ablation", 5.e6)
    n_fragments         = params.get("n_fragments", 3)
    flare_duration      = params.get("flare_duration", 0.1)
    rho_tau_scale       = params.get("rho_tau_scale", 0.0)
    zstart              = params.get("zstart", 100000.0)
    xstart              = params.get("xstart", 0.0)
    tstart              = params.get("tstart", 0.0)
    dt                  = params.get("dt", 0.003)
    tstop               = params.get("tstop", 60.0)

    n_frag_init = n_fragments

    return (
        diameter,
        velocity,
        theta_deg,
        init_strength,
        strength_secondary,
        strength_tertiary,
        strength_quaternary,
        density,
        porosity,
        Cd,
        L_ablation,
        n_fragments,
        n_frag_init,
        flare_duration,
        rho_tau_scale,
        zstart,
        xstart,
        tstart,
        dt,
        tstop
    )

# writing out initial parameters to output file
def write_init_params(outputpath, diameter, velocity, theta_deg,
                      init_strength, strength_secondary, strength_tertiary, strength_quaternary,
                      density, porosity, Cd, L_ablation,
                      n_fragments, n_frag_init, flare_duration, rho_tau_scale,
                      zstart, xstart, tstart, dt, tstop):

    init_params = open(f'{outputpath}/init_params.dat', 'w')

    init_params.write(f"Initialization Parameters for Bolide Luminosity and Fragmentation Simulation\n")
    init_params.write(f"Diameter (m):                      {diameter:.3f}\n")
    init_params.write(f"velocity (m/s):                    {velocity:.3f}\n")
    init_params.write(f"Entry Angle to Horizon (deg):      {theta_deg:.3f}\n")
    init_params.write(f"Initial Strength (Pa):             {init_strength:.3e}\n")
    init_params.write(f"Secondary Strength (Pa):           {strength_secondary:.3e}\n")
    init_params.write(f"Tertiary Strength (Pa):            {strength_tertiary:.3e}\n")
    init_params.write(f"Quaternary Strength (Pa):          {strength_quaternary:.3e}\n")
    init_params.write(f"Density (kg/m^3):                  {density:.3f}\n")
    init_params.write(f"Porosity:                          {porosity:.3f}\n")
    init_params.write(f"Drag Coefficient:                  {Cd:.3f}\n")
    init_params.write(f"Heat of Ablation (J/kg):           {L_ablation:.3e}\n")
    init_params.write(f"Number of Fragments per stage:     {n_fragments}\n")
    init_params.write(f"Initially set number of fragments: {n_frag_init}\n")
    init_params.write(f"Flare Duration (sec):              {flare_duration:.3f}\n")
    init_params.write(f"Luminous Efficiency scaling term:  {rho_tau_scale}\n")
    init_params.write(f"Initial Altitude (m):              {zstart:.3f}\n")
    init_params.write(f"Initial Position (m):              {xstart:.3f}\n")
    init_params.write(f"Starting Time (s):                 {tstart:.3f}\n")
    init_params.write(f"Timestep (s):                      {dt:.3f}\n")
    init_params.write(f"Maximum Time (s):                  {tstop:.3f}\n")
    
    init_params.close()

# define general output file
def wr_out(bolide_outputs, luminosity, alt, t, mass, diameter, area, KE, v, acc, q_h, Fdrag,
           E_rad_total, E_deposited, q, M, T_stag, p_stag, T_surf, T, p, rho, a, theta_deg, tau,
           frag_stage, lum_per_fragment):

    bolide_outputs.write(f"time = {t:.4e}       lum = {luminosity:.4e}     tau = {tau:.4e}   dia = {diameter:.4e}    alt = {alt:.4e}\n")
    bolide_outputs.write(f"E_rad_ttl = {E_rad_total:.4e}  E_dep = {E_deposited:.4e}   KE = {KE:.4e}    Fdrag = {Fdrag:.4e}  frag_stage = {frag_stage}\n")
    bolide_outputs.write(f"v = {v:.4e}          acc = {acc:.4e}     theta_deg = {theta_deg:.3f} mass = {mass:.4e}   area = {area:.4e}\n")
    bolide_outputs.write(f"T_stag = {T_stag:.4e}     T_surf = {T_surf:.4e}  q_h = {q_h:.4e}   q = {q:.4e}      p_stag = {p_stag:.4e}\n")
    bolide_outputs.write(f"T = {T:.4e}          p = {p:.4e}       rho = {rho:.4e}   a = {a:.4e}      M = {M:.4e}\n\n")
    if lum_per_fragment is not None:
        bolide_outputs.write(f"lum_per_fragment = {lum_per_fragment:.4e}\n")
        bolide_outputs.write(f"\n")

# define power vs time output file
def wr_power_vs_time(bolide_pvt, luminosity, t):

    bolide_pvt.write(f'{luminosity:.4e}   {t:.4e}\n')

def wr_power_vs_alt(bolide_pva, luminosity, alt):

    bolide_pva.write(f'{luminosity:.4e}   {alt:.4e}\n')