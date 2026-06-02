'''
Bolide Query and Plotting Tool
'''

import matplotlib.pyplot as plt
import pkg_resources
from bolides import BolideDataFrame
import pandas as pd

# Load the latest GLM bolides from NASA
bdf = BolideDataFrame(source='glm')  # This fetches/parses from [neo-bolide.ndc.nasa.gov](http://neo-bolide.ndc.nasa.gov)

# Filter around the Cleveland event (use datetime strings or pandas Timestamp)
event_date = '2026-03-17'
subset = bdf[bdf['time'].dt.strftime('%Y-%m-%d') == event_date]  # or broader range

# Or more flexibly:
subset = bdf[(bdf['time'] >= '2026-03-17 12:00') & (bdf['time'] < '2026-03-18')]

# Search by location (Cleveland approx 41.5°N, 81.7°W)
cleveland_events = bdf[bdf['lat'].between(40, 43) & bdf['lon'].between(-83, -80)]

print(subset)  # or subset.head(), etc.

# For a specific event (once you identify the row/index or ID)
event = subset.iloc[0]  # or whichever matches
event.plot_lightcurve()  # pulls and plots from the NASA site

glm_lc = event.get_lightcurve()  # or however the method works; check docs

# Plot
# plt.figure(figsize=(10, 6))
# plt.plot(bolum['time'], bolum['power'], label='Your Model', linewidth=2)
# plt.plot(glm_lc['time'], glm_lc['radiance'], label='GLM Observed', linestyle='--')
# plt.xlabel('Time (s relative to peak or entry)')
# plt.ylabel('Power / Radiance')
# plt.title('Boston Bolide: GLM vs Your Model')
# plt.legend()
# plt.grid(True)
# plt.show()