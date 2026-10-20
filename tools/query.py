import requests
import json
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
import os
from datetime import datetime, timezone
import numpy as np

# API endpoints
EVENT_LIST_URL = "https://neo-bolide.ndc.nasa.gov/service/event/public"
EVENT_URL_BASE = "https://neo-bolide.ndc.nasa.gov/service/event/"

def fetch_event_list(limit=500):  # Increase limit as needed (API usually allows hundreds)
    """Fetch bolide events (larger batch for filtering)."""
    params = {"limit": limit}
    response = requests.get(EVENT_LIST_URL, params=params)
    response.raise_for_status()
    data = response.json()

    df = pd.DataFrame(data.get('data', []))

    if not df.empty and 'datetime' in df.columns:
        # Handle both epoch-ms and ISO string formats
        if pd.api.types.is_numeric_dtype(df['datetime']):
            df['datetime'] = pd.to_datetime(df['datetime'], unit='ms', utc=True)
        else:
            df['datetime'] = pd.to_datetime(df['datetime'], utc=True)

    print(f"Fetched {len(df)} total events.")
    return df

def filter_by_date(df, start_date=None, end_date=None):
    """Filter DataFrame by datetime range (ISO format or datetime objects)."""
    if start_date:
        if isinstance(start_date, str):
            start_date = pd.to_datetime(start_date, utc=True)
        df = df[df['datetime'] >= start_date]

    if end_date:
        if isinstance(end_date, str):
            end_date = pd.to_datetime(end_date, utc=True)
        df = df[df['datetime'] <= end_date]

    print(f"After date filtering: {len(df)} events.")
    return df

def fetch_light_curve(event_id):
    """Fetch light curve + metadata for one event."""
    url = f"{EVENT_URL_BASE}{event_id}"
    response = requests.get(url)
    response.raise_for_status()
    data = response.json()

    attachments = data.get('data', [{}])[0].get('attachments', [])
    light_curves = {}

    for att in attachments:
        platform = att.get('platformId')
        geodata = att.get('geoData', [])
        if geodata:
            times = [point['time'] / 1000 for point in geodata]   # seconds since epoch
            energies = [point['energy'] for point in geodata]
            light_curves[platform] = {'time': times, 'energy': energies}

    return light_curves, data.get('data', [{}])[0]

def download_light_curves(event_ids, output_dir="bolide_light_curves"):
    """Download and save data."""
    os.makedirs(output_dir, exist_ok=True)

    for eid in tqdm(event_ids, desc="Downloading"):
        try:
            lcs, metadata = fetch_light_curve(eid)

            # Save metadata
            with open(os.path.join(output_dir, f"{eid}_metadata.json"), "w") as f:
                json.dump(metadata, f, indent=2, default=str)

            # Save light curves
            for sat, lc in lcs.items():
                df_lc = pd.DataFrame({
                    'time_seconds': lc['time'],
                    'energy': lc['energy']
                })
                df_lc.to_csv(os.path.join(output_dir, f"{eid}_{sat}.csv"), index=False)

                # Plot energy vs time
                plt.figure(figsize=(9, 5))
                plt.plot(lc['time'], lc['energy'], 'b.-', markersize=4)
                plt.xlabel('Time (seconds)')
                plt.ylabel('Energy (J)')
                plt.title(f'{eid} — {sat}')
                plt.grid(True, alpha=0.3)
                plt.savefig(os.path.join(output_dir, f"E_vs_t_{eid}_{sat}.png"), dpi=150)

                # Plot power vs time
                energy = lc['energy']
                power = np.array(energy) / 0.002
                plt.figure(figsize=(9, 5))
                plt.plot(lc['time'], power, 'b.-', markersize=4)
                plt.xlabel('Time (seconds)')
                plt.ylabel('Power (W)')
                plt.title(f'{eid} — {sat}')
                plt.grid(True, alpha=0.3)
                plt.savefig(os.path.join(output_dir, f"P_vs_t_{eid}_{sat}.png"), dpi=150)
                plt.close()

        except Exception as e:
            print(f"⚠️ Error on {eid}: {e}")

# ====================== EXAMPLE USAGE ======================
if __name__ == "__main__":
    # 1. Fetch a large batch
    events = fetch_event_list(limit=1000)

    # 2. Filter by date/time (adjust as needed)
    events = filter_by_date(
        events,
        start_date="2026-03-16",      # YYYY-MM-DD or full ISO
        end_date="2026-03-18"         # or datetime(2025, 7, 1, tzinfo=timezone.utc)
    )

    # Optional: further filters
    # events = events[events['confidenceRating'] == 'high']

    # 3. Download light curves
    sample_ids = events['_id'].tolist()   # or .head(20).tolist() for testing
    download_light_curves(sample_ids, output_dir="my_bolide_data")

    print("✅ Finished! Check the output folder.")