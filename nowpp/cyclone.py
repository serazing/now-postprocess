import numpy as np
import xarray as xr
import pandas as pd
import wrf
from . import geometry
import matplotlib.pyplot as plt

EARTH_RADIUS = 6371 * 1e3


def get_cyclone_properties(cyclone_tracking):
    """
    Parameters
    ----------
    cyclone_tracking : xarray.Dataset
        The raw dataset containing the tracking of cyclone from Di Luca's
        code

    Returns
    -------
    cyclone_data : xarray.Dataset
        The tracking dataset reorganized and ready to process
    """
    def get_attribute(i):
        return cyclone_tracking.sel(Attributes=i)['ECL']
    lat = get_attribute(6)  # LAT_CEN
    lon = get_attribute(7)  # LON_CEN
    # size = get_attribute(32)   # SIZE_MEAN
    pressure = get_attribute(8)  # P_CEN
    laplacian = get_attribute(40)  # LAP_MEAN
    size = get_attribute(32)  # SIZE_MEAN
    year = get_attribute(2)  # YEAR
    month = get_attribute(3)  # MONTH
    day = get_attribute(4)  # DAY
    hour = get_attribute(5)  # HOUR
    dt = float(cyclone_tracking.temporal_scale_hours)
    duration = get_attribute(25).isel(TimeSteps=0) * dt
    
    cyclone_km = xr.IndexVariable('distance', [50, 100, 150, 200, 300,
                                               400, 500, 600])
    index_attributes = range(13, 21)
    pressure_gradient = xr.concat([get_attribute(ind)
                                   for ind in index_attributes],
                                  dim=cyclone_km)  #PG_300
    # Create a Dataset gathering the different variable
    cyclone_data = xr.Dataset({'lat': lat, 'lon': lon,
                               'pressure_gradient': pressure_gradient,
                               'pressure': pressure,
                               'laplacian': laplacian,
                               'size': size,
                               'duration': duration,
                               'year': year,
                               'month': month,
                               'day': day,
                               'hour': hour})
    cyclone_data.attrs['dt'] = dt
    # Cleaning data and renaming dimensions and
    cyclone_data = cyclone_data.where((np.abs(lat) < 90) &
                                      (laplacian < 1e10) &
                                      (size < 1e10)
                                      )
    cyclone_data = cyclone_data.rename({'ECLevent': 'event',
                                        'TimeSteps': 'time'})
    cyclone_data = cyclone_data.set_coords(('lat', 'lon'))
    cyclone_data = cyclone_data.assign(event=range(cyclone_data.sizes['event']),
                                       time=range(cyclone_data.sizes['time']))
    return cyclone_data


def assign_cyclone_trajectory(cyclone_data, dim='time'):
    """
    Compute and assign the trajectory of the cyclone in a new dataset. The
    trajectory of the cyclone consists in a speed translation and an angle.

    Paraneters
    ----------
    cyclone_data: xr.Dataset
        The cyclone dataset
    dim : str, optional
        The name of the time dimension

    Returns
    -------
    res : xr.Dataset
        The cyclone dataset with information about the trajectories added
    """
    angle = geometry.latlon2heading(cyclone_data['lat'], cyclone_data['lon'],
                                    dim)
    dy, dx = geometry.latlon2dydx(cyclone_data['lat'], cyclone_data['lon'],
                                  dim)
    dt = cyclone_data.attrs['dt'] * 3600
    u, v = dx / dt, dy / dt
    speed = np.sqrt(u ** 2 + v ** 2)
    return cyclone_data.assign(speed=speed, angle=angle)


def assign_time_and_location(cyclone_data, wrf_grid):
    ecl_location = (cyclone_data.drop_dims('distance')
                                .stack(disk=('event', 'time'))
                                .dropna('disk'))
    i_centre, j_centre = wrf.ll_to_xy(wrf_grid,
                                      ecl_location['lat'],
                                      ecl_location['lon'], meta=False)
    df = pd.DataFrame({'year': ecl_location['year'],
                       'month': ecl_location['month'],
                       'day': ecl_location['day'],
                       'hour': ecl_location['hour']})
    datetime = pd.to_datetime(df)
    coords = ecl_location.coords
    ecl_location = ecl_location.assign_coords(x=xr.DataArray(i_centre,
                                                             coords=coords),
                                              y=xr.DataArray(j_centre,
                                                             coords=coords),
                                              date=xr.DataArray(datetime,
                                                                coords=coords))
    res = (ecl_location.unstack('disk')
                       .assign({'pressure_gradient':
                                cyclone_data['pressure_gradient']})
           )
    return res


def match_cyclone_events(cyclone_test, cyclone_ref,
                         delta_hours=24., delta_x=600):
    """
    Match similar cyclone events between two simulations
    
    Parameters
    ----------
    cyclone_test : xr.Dataset
        The cyclone tested dataset to be tested against the reference dataset
    cyclone_ref : 
        The cyclone reference dataset
    delta_hours : float, optional
        Default is 24 hr.
    delta_x : float, optional
        Default is 600 km.
        
    Returns
    -------
    event_test : 
        The indices of events in the tested dataset
    event_ref : 
        The indices of events in the reference dataset    
    """
    ds_ref = cyclone_ref.rename({'event': 'event_ref', 'time': 'time_ref'})
    event_ref = xr.IndexVariable('event_ref',
                                 range(ds_ref.sizes['event_ref']))
    ds_ref = ds_ref.assign_coords(event_ref=event_ref)
    ds_ref_stacked = (ds_ref.stack(pass_ref=('event_ref', 'time_ref'))
                            .dropna('pass_ref').chunk({'pass_ref': 5e3}))

    ds_test = cyclone_test.rename({'event': 'event_test', 'time': 'time_test'})
    event_test = xr.IndexVariable('event_test',
                                  range(ds_test.sizes['event_test']))
    ds_test = ds_test.assign_coords(event_test=event_test)
    ds_test_stacked = (ds_test.stack(pass_test=('event_test', 'time_test'))
                              .dropna('pass_test').chunk({'pass_test': 5e3}))
    # Compute time lag between cyclones
    dt = np.abs(ds_test_stacked['date'] - ds_ref_stacked['date'])
    dt = dt.astype('f4') * 1e-9 / 3600.
    # Compute distance between cyclone centers
    lon_ref, lat_ref = ds_ref_stacked['lon'], ds_ref_stacked['lat']
    dlon = ds_test_stacked['lon'] - lon_ref
    dlat = ds_test_stacked['lat'] - lat_ref
    dx = np.cos(np.pi / 180. * lat_ref) * np.pi / 180. * EARTH_RADIUS * dlon    
    dy = np.pi / 180. * EARTH_RADIUS * dlat
    distance = 1e-3 * np.sqrt(dx ** 2 + dy ** 2)
    # Elaborate a score based on space and time separation
    score = np.sqrt((distance / delta_x) ** 2 + (dt / delta_hours) ** 2)
    score_valid = score.where((dt <= delta_hours) & (distance <= delta_x),
                              drop=True)
    # Get the corresponding events based on the best score
    best_score = (score_valid.unstack('pass_test').min('time_test')
                             .unstack('pass_ref').min('time_ref')
                             .compute())
    event_index = best_score.groupby('event_ref').argmin('event_test')
    event_test = best_score.event_test[event_index].data
    event_ref = best_score.event_ref.data
    return event_test, event_ref


def get_cyclone_imprint(cyclone_location, wrfout, n=0):
    """
    Parameters
    ----------
    cyclone_location : xr.Dataset
        The dataset containing the time and location of the cyclones
    wrfout : xr.Dataset
        The dataset containing the WRF outputs
    n : int
        Number of the cyclone to extract
    Returns
    -------

    """
    cyclone_event = cyclone_location.isel(event=n).dropna('time')
    gdata_event = wrfout.sel(time=cyclone_event['date'].data)
    list_of_pass = []
    # Loop over ECL timesteps
    for t in range(cyclone_event.sizes['time']):
        ecl_pass = cyclone_event.isel(time=t)
        i_centre, j_centre = int(ecl_pass['x'].data), int(ecl_pass['y'].data)
        gdata_pass = gdata_event.isel(time=t)\
                                .isel(west_east=slice(i_centre - 50,
                                                      i_centre + 50),
                                      south_north=slice(j_centre - 50,
                                                        j_centre + 50))
        list_of_pass.append(gdata_pass)
    gdata_cyclone = xr.concat(list_of_pass, dim='time')
    return gdata_cyclone


def get_data_by_season(ecl_data, season='Warm'):
    """
    Get the data corresponding to a particular season.

    Parameters
    ----------
    ecl_data : xarray.Dataset
        The ECL data corresponding to one simulation
    season : {'Warm', 'Cool'}, optional
        The seasonal months to select a particular season. Default is 'DJF'

    Returns
    -------
    res : xarray.Dataset
        A Dataset including only the season precised
    """
    if season == 'Warm':
        months = [11, 12, 1, 2, 3, 4]
    elif season == 'Cool':
        months = [5, 6, 7, 8, 9, 10]
    else:
        raise ValueError
    cond = ecl_data.month.isin(months)
    ecl_data_by_season = ecl_data.where(cond, drop=True)
    return ecl_data_by_season


def get_number_of_events(ecl_data, season=None):
    """
    Get the number of events by returning the size of the event dimension.
    If the season is precised, return only the number of event for this season.

    Parameters
    ----------
    ecl_data : xarray.DataSet
        The ECL data corresponding to one simulation
    season : {'Warm', 'Cool'}, optional
        The seasonal months to select a particular season

    Returns
    -------
    n_events : int
        The number of ECL events
    """
    if season is None:
        n_events = ecl_data.sizes['event']
    else:
        n_events = get_data_by_season(ecl_data, season=season).sizes['event']
    return n_events


def filter_events(cyclone_dataset, season='Cool',
                  lat_min=-50, lat_max=0, min_duration=18):
    filtered_dataset = {}
    for sim in cyclone_dataset:
        cyclone_data = get_data_by_season(cyclone_dataset[sim], season=season)
        lat = cyclone_data['lat']
        cyclone_data = cyclone_data.where((lat >= lat_min) & (lat <= lat_max))
        cyclone_data = cyclone_data.where(cyclone_data.duration
                                          >= min_duration, drop=True)
        filtered_dataset[sim] = cyclone_data
    return filtered_dataset


def align_cyclone_with_max_pressure_gradient(cyclone_dataset, lag=6, dt=6):
    aligned_dataset = {}
    for sim in cyclone_dataset:
        cyclone_data = cyclone_dataset[sim]
        list_of_cyclone = []
        lag_coord = np.arange(-lag * dt, (lag + 1) * dt, dt)
        mpi = cyclone_data.pressure_gradient.sel(distance=300).argmax('time')
        event_index = 0
        for ind in mpi.data:
            timesteps = np.arange(-lag + ind, lag + ind + 1)
            cond = np.where(timesteps > 0)
            try:
                cyclone_mpi = (cyclone_data.isel(time=timesteps[cond])
                                           .isel(event=event_index))
                cyclone_mpi = cyclone_mpi.rename({'time': 'lag'})
                cyclone_mpi['lag'] = lag_coord[cond]
                list_of_cyclone.append(cyclone_mpi)
            except IndexError:
                pass
            event_index += 1
        cyclone_mpi = xr.concat(list_of_cyclone, dim='event')
        aligned_dataset[sim] = cyclone_mpi
    return aligned_dataset


def bin_data(data,
             lon_min=0., lon_max=360., lon_res=1.,
             lat_min=-80., lat_max=80, lat_res=1.):
    # Define the latitudinal and longitudinal binning
    lon_bins = np.arange(lon_min, lon_max, lon_res)
    lon_labels = lon_bins[:-1] - np.diff(lon_bins) / 2
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
    res_bins = xr.concat(mean_data, dim='lat')
    res_bins = (res_bins.assign_coords(lat=lat_values)
                        .rename({'lon_bins': 'lon'})
                        .sortby('lat'))
    res_obs = xr.concat(total_nobs, dim='lat')
    res_obs = (res_obs.assign_coords(lat=lat_values)
                      .rename({'lon_bins': 'lon'})
                      .sortby('lat'))
    return xr.Dataset({data.name: res_bins, 'nobs': res_obs})

