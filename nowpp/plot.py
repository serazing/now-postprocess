import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER


def format_lon_lat(ax, proj, lon_min, lon_max, lat_min, lat_max, title=''):
    from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
    #from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
    
    ax.set_extent([lon_min, lon_max, 
                   lat_min, lat_max], proj)
    # Add coastline
    ax.coastlines('50m')
    
    # Modify the title
    ax.set_title(title)
            
    # Set lon labels
    lon_labels = np.arange(lon_min, lon_max + 1, 10)
    lon_labels[lon_labels > 180] -= 360
    ax.set_xticks(lon_labels, crs=proj)
    ax.set_xticklabels(lon_labels, rotation=45)
    ax.xaxis.set_major_formatter(LongitudeFormatter())
    ax.set_xlabel('')
            
    # Set lat labels
    lat_labels = np.arange(lat_min, lat_max + 1, 10)
    ax.set_yticks(lat_labels, crs=proj)
    ax.set_yticklabels(lat_labels)
    ax.yaxis.tick_left()
    ax.yaxis.set_major_formatter(LatitudeFormatter())
    ax.set_ylabel('')
                        
    # Plot the grid
    ax.grid(linewidth=2, color='black', alpha=0.5, linestyle='--')

            
def add_map(lon_min=-180, lon_max=180, lat_min=-90, lat_max=90,
            central_longitude=0., scale='auto', ax=None):
    """
    Add the map to the existing plot using cartopy

    Parameters
    ----------
    lon_min : float, optional
        Western boundary, default is -180
    lon_max : float, optional
        Eastern boundary, default is 180
    lat_min : float, optional
        Southern boundary, default is -90
    lat_max : float, optional
        Northern boundary, default is 90
    central_longitude : float, optional
        Central longitude, default is 180
    scale : {?auto?, ?coarse?, ?low?, ?intermediate?, ?high, ?full?}, optional
        The map scale, default is 'auto'
    ax : GeoAxes, optional
        A new GeoAxes will be created if None

    Returns
    -------
    ax : GeoAxes
    Return the current GeoAxes instance
    """

    extent = (lon_min, lon_max, lat_min, lat_max)
    if ax is None:
        ax = plt.subplot(1, 1, 1,
                         projection=ccrs.PlateCarree(                                  central_longitude=central_longitude))
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    land = cfeature.GSHHSFeature(scale=scale,
                                 levels=[1],
                                 facecolor=cfeature.COLORS['land'])
    ax.add_feature(land)
    gl = ax.gridlines(draw_labels=True, linestyle=':', color='black',
                      alpha=0.5)
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    return ax


def plot_one_map(data, ax=None, lon_min=110, lon_max=180, 
                 lat_min=-45, lat_max=0, **kwargs):
    data['nav_lon'] = data['nav_lon'] % 360
    add_map(ax=ax, lon_min=lon_min, lon_max=lon_max, 
            lat_min=lat_min, lat_max=lat_max)
    mapped_data = data.plot(x='nav_lon', y='nav_lat', **kwargs)
    return mapped_data
    
