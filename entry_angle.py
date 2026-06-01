import numpy as np

def meteor_entry_angle(lat_deg: float, lon_deg: float, vx: float, vy: float, vz: float):
    """
    Calculate meteor entry angle from horizontal (flight path angle).

    Parameters:
        lat_deg : float - Latitude in degrees (e.g. 41.2)
        lon_deg : float - Longitude in degrees (e.g. -82.0)
        vx, vy, vz : float - ECEF velocity components in km/s

    Returns:
        total_speed : float - Speed magnitude (km/s)
        angle_from_horizontal : float - Degrees from horizontal (0° = grazing, 90° = straight down)
        v_radial : float - Radial velocity component (negative = inbound)
    """
    # Convert to radians
    phi = np.deg2rad(lat_deg)
    lon = np.deg2rad(lon_deg)

    # Local radial (up/outward) unit vector in ECEF frame
    radial = np.array([
        np.cos(phi) * np.cos(lon),   # x
        np.cos(phi) * np.sin(lon),   # y
        np.sin(phi)                 # z
    ])

    # Velocity vector and magnitude
    v = np.array([vx, vy, vz])
    v_mag = np.linalg.norm(v)

    # Project onto local vertical
    v_radial = np.dot(v, radial)

    # Angle from horizontal
    if v_mag > 0:
        angle_deg = np.rad2deg(np.arcsin(abs(v_radial) / v_mag))
    else:
        angle_deg = 0.0

    return v_mag, angle_deg, v_radial


# ========================
# Main
# ========================
if __name__ == "__main__":
    lat = 41.2
    lon = -82.0
    vx = 2.3
    vy = -2.0
    vz = -14.6

    speed, angle, vrad = meteor_entry_angle(lat, lon, vx, vy, vz)

    print(f"Total speed:      {speed:.2f} km/s")
    print(f"Angle from horizontal: {angle:.1f}°")
    print(f"v_radial (inbound):   {vrad:.2f} km/s")