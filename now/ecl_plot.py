import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

def add_map(lon_min=-180, lon_max=180, lat_min=-90, lat_max=90,
            central_longitude=0., scale='auto', ax=None,
			facecolor=cfeature.COLORS['land']):
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
						 projection=ccrs.PlateCarree(central_longitude=central_longitude))
	ax.set_extent(extent)
    #land = cfeature.GSHHSFeature(scale=scale,
                                 #levels=[1],
                                 #facecolor=facecolor)
    #ax.add_feature(land)
	ax.coastlines()
	gl = ax.gridlines(draw_labels=True, linestyle=':', color='black', alpha=0.5)
	gl.xlabels_top = False
	gl.ylabels_right = False
	gl.xformatter = LONGITUDE_FORMATTER
	gl.yformatter = LATITUDE_FORMATTER
	return ax


def plot_ecl_trajectories(ecl_data, season=None, variable=None, lw=1, alpha=0.5, **kwargs):
    for n in ecl_data['event']:
        ecl_event = ecl_data.sel(event=n)
        if season is 'DJF':
            ecl_event = ecl_event.where((ecl_event['month'] % 12) <= 2)
        elif season is 'JJA':
            ecl_event = ecl_event.where((ecl_event['month'] >= 6) & (ecl_event['month'] <= 8))
        plt.plot(ecl_event['lon'], ecl_event['lat'], color='black', lw=lw, alpha=alpha)
        if variable is not None:
            plt.scatter(ecl_event['lon'], ecl_event['lat'], c=ecl_event[variable], **kwargs)
            
            
def bin_data(data, lon_res=1., lat_res=1., 
                   lon_min=0., lon_max=360., 
                   lat_min=-80., lat_max=80):
    # Define the latitudinal and longitudinal binning
    #lon_bins = np.arange(data['lon'].min().data
    #                     data['lon'].max().data, lon_res)
    lon_bins = np.arange(lon_min, lon_max, lon_res)
    lon_labels = lon_bins[:-1] - np.diff(lon_bins) / 2
    #lat_bins = np.arange(data['lat'].min().data,
    #                     data['lat'].max().data, lat_res)
    lat_bins = np.arange(lat_min, lat_max, lat_res)
    lat_labels = lat_bins[:-1] - np.diff(lat_bins) / 2
    mean_data = []
    total_nobs = []
    lat_values = []
    for i, ds in list(data.groupby_bins('lat', lat_bins, 
                                        labels=lat_labels, 
                                        include_lowest=True)):
        try:
            group = ds.groupby_bins('lon', lon_bins, 
                                    labels=lon_labels, 
                                    include_lowest=True)
            bins = group.median().sortby('lon_bins')
            nobs = group.count().sortby('lon_bins')
            mean_data.append(bins)
            total_nobs.append(nobs)
            lat_values.append(i)
        except (ValueError, StopIteration):
            dummy_array = xr.DataArray(np.full(len(lon_labels), 
                                               np.nan), 
                                       dims='lon_bins', 
                                       coords={'lon_bins': 
                                               ('lon_bins', lon_labels)
                                              }
                                      )
            mean_data.append(dummy_array)
            lat_values.append(i)
    res_bins = (xr.concat(mean_data, dim='lat')
                  .assign_coords(lat=lat_values)
                  .rename({'lon_bins': 'lon'})
                  .sortby('lat')
               )
    res_obs = (xr.concat(total_nobs, dim='lat')
                 .assign_coords(lat=lat_values)
                 .rename({'lon_bins': 'lon'})
                 .sortby('lat')
              )
    return xr.Dataset({data.name:res_bins, 'nobs':res_obs})