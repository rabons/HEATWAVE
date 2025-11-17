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


def plot_heatwave(filename, plottitle, color_boundaries, cbar_ticks):
    """
    Plot a heatwave map for Bangladesh.
    
    Parameters:
    -----------
    filename : str
        Path to the raster file (without extension), e.g., "Data/BGD/3deg_Tw_mean_"
    plottitle : str
        Title for the plot
    """
    # Construct file paths
    # filename comes as "Data/BGD/3deg_Tw_mean_" (without extension)
    # Try both .tiff and .tif extensions
    rasterfile = f"data/bgd/{filename}.tiff"
    
    # Output file: use the base filename without folder path
    base_filename = os.path.basename(filename)
    outputfile = f"figures/bgd/{base_filename}_cartopy.png"
    
    # Open the GeoTIFF file using rioxarray
    data_array = rioxarray.open_rasterio(rasterfile)
    print(f"Processing: {base_filename}")
    print("Data shape:", data_array.shape)
    print("Data dims:", data_array.dims)

    # Read boundary file
    boundaryfile_a0 = "data/boundaries/WB_gpkg/World Bank Official Boundaries - Admin 0 BGD.gpkg"
    boundaryfile_a1 = "data/boundaries/WB_gpkg/World Bank Official Boundaries - Admin 1 BGD.gpkg"
    boundaryfile_a2 = "data/boundaries/WB_gpkg/World Bank Official Boundaries - Admin 2 BGD.gpkg"
    # Load both boundary files
    boundaries_a0 = gpd.read_file(boundaryfile_a0)
    boundaries_a1 = gpd.read_file(boundaryfile_a1)
    boundaries_a2 = gpd.read_file(boundaryfile_a2)
    # Color management
    #  colors = ["#FEFECB", "#FBEB99", "#F4CC68", "#EBA754", "#E48650", "#D1624C", "#A44642", "#72372E", "#422818", "#191900"] IPCC sequential
    # colors = ["#93C0DB", "#D7E6EE", "#FEFECB","#FBEB99","#F4CC68","#EBA754","#E48650","#D1624C","#B8524A","#A44642","#8A3D3D","#72372E","#5A2E38","#40204D","#352045","#2B1B3D"] #,"#21102A"] #,"#1A0D20","#191900"]
    # colors = ["#93C0DB","#A8CEE0","#BDD9E6","#D7E6EE","#E5F0F5","#F2F5E8","#FEFEE0","#FEFECB","#FBEB99","#F4CC68","#EBA754","#E48650","#D1624C","#B8524A","#A44642","#8A3D3D","#72372E","#5A2E38","#40204D","#352045","#2B1B3D","#21102A"]
    # colors = ["#93C0DB","#A8CEE0","#BDD9E6","#D7E6EE","#E5F0F5","#F2F5E8","#FEFEE0","#FEFECB","#FBEB99","#F4CC68","#EBA754","#E48650","#D1624C","#B8524A","#A44642","#8A3D3D","#72372E","#5A2E38","#40204D","#524261","#6B6370","#8C8A8E","#B5B5B8","#E8E8EA"]
    # colors = ["#FEFECB","#FBEB99","#F4CC68","#EBA754","#E48650","#D1624C","#B8524A","#A44642","#8A3D3D","#72372E","#5A2E38","#40204D","#4A2B55","#4E3A5C","#524261","#6B6370","#8C8A8E","#B5B5B8"] # Fewer purple steps
    # colors = ["#FEFECB","#FBEB99","#F4CC68","#EBA754","#E48650","#D1624C","#B8524A","#A44642","#8A3D3D","#72372E","#5A2E38","#40204D","#452854","#4D365B","#524261","#6B6370","#8C8A8E","#B5B5B8"] # More purple steps. Claudes' proposal
    # colors = ["#FEFECB","#FBEB99","#F4CC68","#EBA754","#E48650","#D1624C","#B8524A","#A44642","#8A3D3D","#72372E","#5A2E38","#40204D","#4D365B","#524261","#6B6370","#8C8A8E","#B5B5B8"] # More purple steps. Claudes' proposal with one less purple step; #452854 removed
     # colors = ["#FEFECB","#FBEB99","#F4CC68","#EBA754","#E48650","#DB7360","#C77062","#B55D58","#A45A56","#8B5548","#6F4852","#573A63","#654A73","#6B5672","#7D7380","#9A989C","#B5B5B8"] # Same than above, but lighter
    # colors = ["#FEFECB","#FBEB99","#F8DE85","#F4CC68","#F0BB5E","#EBA754","#E48650","#D1624C","#B8524A","#A44642","#8A3D3D","#72372E","#5A2E38","#40204D","#4D365B","#524261","#6B6370","#8C8A8E","#B5B5B8"] # Same than above, but more yellow steps
    colors = ["#FEFECB","#FBEB99","#F8DE85","#F4CC68","#F0BB5E","#EBA754","#E48650","#DB7360","#C77062","#B55D58","#A45A56","#8B5548","#6F4852","#573A63","#654A73","#6B5672","#7D7380","#9A989C","#B5B5B8"] # Same than above, but somewhat lighter
   
    ncolors = 256
    cmap = LinearSegmentedColormap.from_list("Tw", colors, N=ncolors)
    cmap.set_bad(color="white")
    # Color boundaries will be defined in the main loop
    # color_boundaries = np.arange(16, 34.5, .5)  # Include 30 by going slightly beyond
    # cbar_ticks = np.arange(16, 35, 2)  # Includes 30
    # norm = BoundaryNorm(color_boundaries, ncolors=cmap.N, clip=True)
    norm = BoundaryNorm(color_boundaries, ncolors=ncolors, clip=True)

    ################################################################################
    # Squeeze the band dimension for single-band rasters
    img = data_array.squeeze()

    # Set figure width in cm and convert to inches
    ny, nx = img.shape[-2], img.shape[-1]
    aspect = ny / nx
    width_cm = 8 #
    width_in = width_cm / 2.54
    height_in = width_in * aspect

    # Start context manager for fonts
    with rc_context({'font.family': 'sans-serif', 'font.size': 8}): # For ppt use Tahoma
        ################################################################################
        # Plot with +proj=tmerc +lon_0=90.4 +lat_0=23.8
        fig = plt.figure(figsize=(width_in, height_in))
        ax = fig.add_subplot(
            1, 1, 1,
            # Projection for Dhaka
            # projection=ccrs.TransverseMercator(
            #     central_longitude=90.4,
            #     central_latitude=23.8)
            projection=ccrs.PlateCarree()
        )
        # Gridlines
        gl = ax.gridlines(
            draw_labels=True,  # Enable labels
            linewidth=0.0,     # Very thin line (needed for ticks to appear)
            color='#333333',
            linestyle='--',
            alpha=0.3
        )
        
        # Style the labels
        gl.xlabel_style = {'size': 4,'color': '#333333'}
        gl.ylabel_style = {'size': 4,'color': '#333333'}
        
        # Rotate labels for better readability
        gl.rotate_labels = False
        
        # Control label format (degrees with formatters)
        gl.xformatter = LongitudeFormatter()
        gl.yformatter = LatitudeFormatter()

        # Control which sides have labels (left, top, bottom)
        gl.left_labels = True
        gl.right_labels = False
        gl.top_labels = False
        gl.bottom_labels = True
        

        # Reproject and plot the raster data
        # img=img.rio.reproject(ccrs.TransverseMercator(
        #         central_longitude=90.4,
        #         central_latitude=23.8
        #     ))
        im = img.plot.imshow(
            ax=ax,
            interpolation = "gaussian",
            #alpha=0.5,
            alpha=0.95,
            cmap=cmap,
            norm=norm,
            add_colorbar=True,
            cbar_kwargs={
                'boundaries': color_boundaries,
                'ticks': cbar_ticks,
                'label': f'T_{{w}} [°C]',
                'shrink': 0.65,  # Reduce height (0.8 = 80% of default height)
                'aspect': 30,    # Make it thinner (higher aspect ratio = thinner)
                'location': 'right'
            }
        )
        # A light gray, nearly white:
        '#F0F0F0'
        # Plot Admin0 boundaries (country level) - thicker lines
        for geometry in boundaries_a0.geometry:
            ax.add_geometries(
                [geometry],
                crs=ccrs.PlateCarree(),  # EPSG:4326
                facecolor='none',
                edgecolor='#333333',
                linewidth=0.2,  # Thicker line for country boundaries
                linestyle='-',
                zorder=10  # Draw on top of raster
            )
        
        # # Plot Admin2 boundaries (subnational level) - thinner lines
        for geometry in boundaries_a1.geometry:
            ax.add_geometries(
                [geometry],
                crs=ccrs.PlateCarree(),
                facecolor='none',
                edgecolor='#333333',
                linewidth=0.1,
                linestyle='-',
                zorder=11
            )

        for geometry in boundaries_a2.geometry:
            ax.add_geometries(
                [geometry],
                crs=ccrs.PlateCarree(),
                facecolor='none',
                edgecolor='#333333',
                linewidth=0.05,
                linestyle='--',
                zorder=12
            )    

        ax.set_title(plottitle, fontsize=8)

        # Style the plot frame (spines) with gray50 and thin linewidth
        for spine in ax.spines.values():
            spine.set_color('#808080')
            spine.set_linewidth(0.1)
        
        # Style the axis labels and tick labels to gray20
        # ax.set_xlabel('Longitude', color='#333333')
        # ax.set_ylabel('Latitude', color='#333333')
        # Style tick marks on the axis spines
        # Note: In cartopy, tick marks are controlled by gridlines, but we can also style them here
            # Ensure ticks are visible by accessing the underlying axes
        # Get the geographic extent to set ticks properly
        # ax.set_xticks(np.arange(87, 95, 1), crs=ccrs.PlateCarree())
        # ax.set_yticks(np.arange(19, 27, 1), crs=ccrs.PlateCarree())

        # # Set the tick formatters
        # ax.xaxis.set_major_formatter(LongitudeFormatter())
        # ax.yaxis.set_major_formatter(LatitudeFormatter())

        # ax.tick_params(
        #     axis='both',
        #     # which='both',
        #     which='major',
        #     direction='out',  # Outward-facing ticks
        #     length=4,         # Length of tick marks
        #     width=0.5,        # Width of tick marks
        #     color='#333333',
        #     labelsize=6,
        #     bottom=True,      # Show ticks on bottom
        #     top=False,        # Hide ticks on top (labels are shown via gridlines)
        #     left=True,        # Show ticks on left
        #     right=False       # Hide ticks on right
        # )
        
        
        # Style the colorbar frame and move label to top
        cbar = im.colorbar  # Get the colorbar object directly
        cbar_ax = cbar.ax  # Get the colorbar axes
        
        # Style the colorbar frame
        if hasattr(cbar_ax, 'spines'):
            for spine in cbar_ax.spines.values():
                spine.set_color('#333333')
                spine.set_linewidth(0.1)
        
        # Style the colorbar labels and position label at top
        cbar.set_label(fr'[°C]', color='#333333', labelpad=10, rotation=0)  # LaTeX with subscript
        # Position label at top by setting y-coordinate to 1.02 (above the colorbar)
        cbar_ax.yaxis.set_label_coords(1.05, 1.1)
        
        # Configure major and minor ticks
        # Major ticks: every 2 degrees (already set via cbar_ticks)
        # Minor ticks: every 1 degree
        cbar_ax.yaxis.set_minor_locator(MultipleLocator(1))  # Minor ticks every 1 degree

        
        # Style the ticks
        cbar_ax.tick_params(colors='#333333', which='major', labelsize=8)
        cbar_ax.tick_params(colors='#333333', which='minor', length=3, width=0.5)  # Smaller minor ticks

    
        plt.savefig(outputfile, dpi=300, bbox_inches="tight", pad_inches=0.01)
        plt.close()
    
    print(f"Saved successfully: {outputfile}")