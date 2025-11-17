# Run plots for all Heatwave files

# Write a function to plot a heatwave plot for a given file. I would like to loop through a number of filenames and and plottitles, and call the function that we just have designed.

from typing import Any


import numpy as np
import rioxarray
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm, Normalize
from matplotlib.colors import ListedColormap
from matplotlib import rc_context
from matplotlib.ticker import FixedLocator, MultipleLocator
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import os
import geopandas as gpd

# Get all files in folder "Data/bgd"
folder = "data/bgd"
filenames = os.listdir(folder)
# ['3deg_Tw_mean_.tiff', '3deg_Tw_3hrmax_.tiff', '2deg_Tw_mean_.tiff', '4deg_Tw_3hrmax_.tiff', '3deg_tas_mean_.tiff', '2deg_tas_mean_.tiff', '4deg_tas_mean_.tiff', '2deg_tas_3hrmax_.tiff', '3deg_tas_3hrmax_.tiff', '4deg_tas_3hrmax_.tiff', '4deg_Tw_mean_.tiff', '2deg_Tw_3hrmax_.tiff']
# Also include the fiules ending with tif and tiff
filenames = [f for f in filenames if f.endswith(".tif") or f.endswith(".tiff")]

# Remove file extensions from filenames
filenames = [os.path.splitext(f)[0] for f in filenames]
# ['3deg_Tw_mean_', '3deg_Tw_3hrmax_', '2deg_Tw_mean_', '4deg_Tw_3hrmax_', '3deg_tas_mean_', '2deg_tas_mean_', '4deg_tas_mean_', '2deg_tas_3hrmax_', '3deg_tas_3hrmax_', '4deg_tas_3hrmax_', '4deg_Tw_mean_', '2deg_Tw_3hrmax_']

degrees = ['2', '3', '4']
vars = ['Tw', 'tas']
freq = ['mean', '3hrmax']

# Option 1: Same  color bounds for both mean / 3hrmax:
# Color boundaries used for mean(Tw):
color_boundaries_mean_Tw = np.arange(14, 34.5, .5)
cbar_ticks_mean_Tw = np.arange(14, 35, 2)
# Color boundaries used for max(Tw):
color_boundaries_max_Tw=color_boundaries_mean_Tw
cbar_ticks_max_Tw=cbar_ticks_mean_Tw
# Color boundaries used for mean(T):
color_boundaries_mean_T = np.arange(16, 46.5, .5)
cbar_ticks_mean_T = np.arange(16, 47, 2)
# Color boundaries used for max(T):
color_boundaries_max_T=color_boundaries_mean_T
cbar_ticks_max_T=cbar_ticks_mean_T

# # Option 2: Different color bounds for mean / 3hrmax:
# # Color boundaries used for mean(Tw):
# color_boundaries_mean_Tw = np.arange(14, 28.5, .5)
# cbar_ticks_mean_Tw = np.arange(14, 29, 2)
# # Color boundaries used for max(Tw):
# color_boundaries_max_Tw = np.arange(23, 34.5, .5)
# cbar_ticks_max_Tw = np.arange(20, 35, 2)

# # Color boundaries used for mean(T):
# color_boundaries_mean_T = np.arange(14, 30.5, .5)
# cbar_ticks_mean_T = np.arange(14, 31, 2)
# # Color boundaries used for max(T):
# color_boundaries_max_T = np.arange(26, 46.5, .5)
# cbar_ticks_max_T = np.arange(26, 47, 2)


filename_list = []
plottitles = []
color_boundaries_list = []
cbar_ticks_list = []
for freq in freq:
    for degree in degrees:
        for var in vars:
            filename = f"{degree}deg_{var}_{freq}_"
            if freq == 'mean':
                if var == 'Tw':
                    title = f"Mean wet-bulb temp. @ {degree}°C global warming."
                    color_boundaries = color_boundaries_mean_Tw
                    cbar_ticks = cbar_ticks_mean_Tw
                elif var == 'tas':
                    title = f"Mean air temperature @ {degree}°C global warming."
                    color_boundaries = color_boundaries_mean_T
                    cbar_ticks = cbar_ticks_mean_T
            elif freq == '3hrmax':
                if var == 'Tw':
                    title = f"Max wet-bulb temp. @ {degree}°C global warming."
                    color_boundaries = color_boundaries_max_Tw
                    cbar_ticks = cbar_ticks_max_Tw
                elif var == 'tas':
                    title = f"Max air temperature @ {degree}°C global warming."
                    color_boundaries = color_boundaries_max_T
                    cbar_ticks = cbar_ticks_max_T
            filename_list.append(filename)
            plottitles.append(title)
            color_boundaries_list.append(color_boundaries)
            cbar_ticks_list.append(cbar_ticks)

# ['2deg_Tw_mean_', '2deg_Tw_3hrmax_', '3deg_Tw_mean_', '3deg_Tw_3hrmax_', '4deg_Tw_mean_', '4deg_Tw_3hrmax_', '2deg_tas_mean_', '2deg_tas_3hrmax_', '3deg_tas_mean_', '3deg_tas_3hrmax_', '4deg_tas_mean_', '4deg_tas_3hrmax_']
print(plottitles)
print(color_boundaries_list)
print(cbar_ticks_list)

# Drop into interactive Python session to explore variables
# All variables (folder, filenames, etc.) will be available
# import code
# code.interact(local=locals())

# Now call the file "plotting_BGD.py" and pass the filenames and plottitles to it.
from plotting_bgd import plot_heatwave


# for filename, plottitle in zip(filenames, plottitles):
for filename, plottitle, color_boundaries, cbar_ticks in zip(filename_list, plottitles, color_boundaries_list, cbar_ticks_list):
    print(filename)
    print(plottitle)
    plot_heatwave(filename, plottitle, color_boundaries, cbar_ticks)
