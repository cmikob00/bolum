Bolum Readme

   March 2026 - v8:
   Moved drag coefficient (Cd) and heat of ablation (L_ablation) to input deck
   Removed broadband Si flag since luminous efficiency τ term assumes a broadband optical sensor
   Renamed bolum_in.txt to bolum_in.dat
   Set some default input parameters so the user need not define everything in the input deck (unless desired)
   Added inout deck and event paths for user ease.
   Added tertiary and quaternary material strength input parameters for even more fragmentation options.

   February 2026 - v7:
   Single luminous efficiency τ from Borovička et al. (2020)
   One total luminosity term = τ × total KE dissipation rate
   Per-fragment luminosity reporting after fragmentation
   Multi-stage fragmentation: initial strength → secondary (weaker) strength after first break
   Introduced "pancake" flare feature to model dust cloud radiating post-fragmentation
   Removed separate ablation/shock terms and related arrays
   Kept original structure, output format, and plotting as much as possible

   February 2026 - v6:
   Changed function lum_drag_calc to lum_ablation_calc for clarity
   Eliminated SISF (shock intensity scaling factor) - now using more physically
   realistic ablation and shock luminosity terms.
   Introduced new ablation and shock luminosity terms, including scaled luminous efficiency
   term from Drolshagen et al 2021 for ablation and dynamic emissivity term from Revelle
   and Ceplecha 2005, 2001-2002 for shock.
   Outputs are applicable for optical (~360 - 700 nm) emission.  Can scale up by factor of 1.15
   to 1.5 for broadband silicon - need to integrate broadband si response into code
   Need to improve fragment tracking to account for fragments of different masses
   
   October 2025 - v5:
   improved plotting and visualization
   added altitude-dependent atmospheric temperature, pressure, and density model from NASA GRC
   added writing out atmospheric temperature, pressure, and density to output file
   fixed Mach scaling term in shock luminosity calculation (lum_shock_calc)
   combined Mach scaling and density scaling term into one term called Shock Intensity Scaling Factor (SISF)
   added total radiated energy in kilotons to output file
   output quantities now written directly to numpy arrays instead of lists: started with np.zeros arrays,
   then append with a new row after each timestep
   added timestep counter (helps with array indexing)
   added total event energy calculation from Brown et al, 2002
   added max number of cycles parameter (100,000)
   improved physics and logic for fragmentation: now accounting for total surface area of n_fragments
   and no artificially-modeled "flare"
   updates to diameter and area update functions to account for tracking multiple fragments
   updates to drag force calculations and heat flux to account for cross-sectional area
   consolidated fragmentation logic and math into one self-contained function
   added ballistic trajectory functionality - seems to work fine

   July 2025 - v4:
   added scaled luminous efficiency function based on velocity and density (Ceplecha et al 1998)
   added print statement to tell when the bolide fragments
   added plotting output vs time in W/sr and output vs altitude in W/sr
   changed angle variable in input deck to be angle_deg for consistency
   added in x and y components for position, velocity, and acceleration due to gravity
   added ballistic coefficient (default 1.5) for drag calculcation
   added shock front luminosity calculations - includes scaling for density and stagnation pressure
   made reading in input deck function
   added peak power, radiant intensity, and total energy radiated finding functions
