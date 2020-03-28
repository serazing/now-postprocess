import xarray as xr
import os
from oocgcm.oceanmodels.nemo import grids
import pandas as pd
#------------------------------------------------------------
# Config file sectopm
#------------------------------------------------------------
def read_config_file(config_file):
    from configparser import ConfigParser
    cfg = ConfigParser()
    cfg.read(config_file)
    return cfg


def _check_simulations(config_file, simulations):
    cfg = read_config_file(config_file)
    if simulations is None:
        simulations = [key for key in cfg 
                       if key not in ['DEFAULT', 'GENERAL', 'WRF', 'NEMO']]
    elif not isinstance(simulations, (list, tuple)):
        simulations = [simulations,]
    return simulations


#------------------------------------------------------------
# NEMO section
#------------------------------------------------------------
def get_nemo_grid_from_config(config_file):
    cfg = read_config_file(config_file)
    nemo_grid = grids.nemo_2d_grid(nemo_coordinate_file=cfg['NEMO']['MeshGridFile'],
                                   nemo_byte_mask_file=cfg['NEMO']['MaskFile'])
    return nemo_grid

    
def get_nemo_meshgrid(config_file):
    cfg = read_config_file(config_file)
    meshgrid_file = cfg['NEMO']['MeshGridFile']
    meshgrid = xr.open_dataset(meshgrid_file).squeeze()
    return meshgrid


def get_nemo_filenames_from_config(config_file, simulation, 
                                   grid='T'):
    """
    Get the full path of NEMO filenames from the configuration file
    for a particular simulation and a particular grid
    
    Parameters
    ----------
    config_file : str
        The path of the configuration file
    simlation : str
        The name 
    grid : {'T', 'U', 'V', 'W'}
        Name of the grid
    
    Returns
    -------
    full_filenames : str
        The full filenames using wildcards    
    """
    cfg = read_config_file(config_file)
    try:
        nemo_path = cfg['NEMO']['NetCDFPath']
    except KeyError:
        nemo_path = './'
    try:
        nemo_dim = cfg['NEMO']['Dimensions']
    except KeyError:
        nemo_dim = '2D'
    try:
        nemo_freq = cfg['NEMO']['Frequency']
    except KeyError:
        nemo_freq = '*'        
    for section in cfg:
        if section == simulation: 
            main_path = cfg[simulation]['Path']
            try:
                prefix = cfg[simulation]['Prefix']  
            except KeyError:
                prefix = '*' 
            filenames = '%s_%s_*_*_grid_%s_%s.nc' %(prefix,
                                                    nemo_freq,
                                                    grid,
                                                    nemo_dim)
            full_filenames = main_path + nemo_path + filenames
            return full_filenames    

        
#------------------------------------------------------------
# WRF section
#------------------------------------------------------------ 
def get_wrf_filenames_from_config(config_file, simulation, prefix='wrfout'):
 
    """
    Get the full path of WRF filenames from the configuration file
    for a particular simulation and a particular grid
    
    Parameters
    ----------
    config_file : str
        The path of the configuration file
    simlation : str
        The name of the simulation as writtent in the configuration file
    
    Returns
    -------
    full_filenames : str
        The full filenames using wildcards    
    """
    cfg = read_config_file(config_file)
    try:
        wrf_path = cfg['WRF']['NetCDFPath']
    except KeyError:
        wrf_path = './'
    try:
        wrf_freq = cfg['WRF']['Frequency']
    except:
        wrf_freq = '*'
    try:
        wrf_prefix = cfg['WRF']['Prefix']
    except:
        prefix = '*'
    for section in cfg:
        if section == simulation: 
            main_path = cfg[simulation]['Path']                     
            filenames = '%s_*_%s_*.nc' %(prefix,
                                         wrf_freq)
            full_filenames = main_path + wrf_path + filenames
            return full_filenames
        
        
def get_nemo_zarr_folder_from_config(config_file, simulation, grid='U'):
    cfg = read_config_file(config_file)
    year_start = cfg['GENERAL']['YearStart']
    year_stop = cfg['GENERAL']['YearStop']
    try:
        nemo_path = cfg['NEMO']['ZarrPath']
    except KeyError:
        nemo_path = ''
    try:
        nemo_dim = cfg['NEMO']['Dimensions']
    except KeyError:
        nemo_dim = '2D'
    try:
        nemo_freq = cfg['NEMO']['Frequency']
    except KeyError:
        nemo_freq = '*'
    main_path = cfg[simulation]['Path']
    try:
        prefix = cfg[simulation]['Prefix']  
    except KeyError:
            prefix = '*'
    name = '%s_%s_%s_%s_grid_%s_%s.zarr'%(prefix, nemo_freq,
                                          year_start, year_stop,
                                          grid, nemo_dim)
    zarr_path = main_path + nemo_path + name
    return zarr_path        
        
        
        
# NEMO Section
# ------------------------------------------------------------
def get_nemo_mask(config_file, grid='U'):
    cfg = read_config_file(config_file)
    mask_file = cfg['NEMO']['MaskFile']
    dims = cfg['NEMO']['Dimensions']
    mask_array = xr.open_dataset(mask_file, chunks={}).squeeze()
    if grid is 'U':
        mask = mask_array['umask']
    elif grid is 'V':
        mask = mask_array['vmask']
    elif grid is 'T':
        mask = mask_array['tmask']
    if dims == '2D':
        mask = mask.isel(z=0)
    return mask


def open_nemo_griddata_from_netcdf(config_file, simulations=None, 
                                   grid='T', **kwargs):
    """
    Get the full path of NEMO filenames from the configuration file
    for a particular simulation and a particular grid
    
    Parameters
    ----------
    config_file : str
        The path of the configuration file
    simlation : str or list of str
        Names of simulations to open as written in the configuration file 
    grid : {'T', 'U', 'V', 'W'}, optional
        Type of 'C' grid to look for
    **kwargs : optional
        Additional arguments passed on to `xarray.open_mfdataset`
    
    Returns
    -------
    griddata : xarray.Dataset
        The NEMO outputs represented by a Dataset
    """
    depth_dims = {'U': 'depthu', 'V': 'depthv', 'T': 'deptht'} 
    cfg = read_config_file(config_file)
    year_start = cfg['GENERAL']['YearStart']
    year_stop = cfg['GENERAL']['YearStop']
    nemo_dim = cfg['NEMO']['Dimensions']
    simulations = _check_simulations(config_file, simulations)
    list_of_gdata = []
    # Get the corresponding mask
    mask = get_nemo_mask(config_file, grid=grid)
    for sim in simulations:
        filenames = get_nemo_filenames_from_config(config_file, sim, grid=grid)
        gdata = xr.open_mfdataset(filenames, **kwargs)
        # Rename time_counter
        gdata = gdata.set_index(time_counter='time_average_1d')
        # Rename depth coordinates if 3D fields
        if nemo_dim == '3D':
            gdata = gdata.rename({depth_dims[grid]: 'z'})
        # Remove duplicated values
        #gdata = gdata.sel(time_counter=\
        #                        ~gdata.indexes['time_counter'].duplicated())
        # Select a particular time range
        gdata = gdata.sel(time_counter=slice(year_start, year_stop))
        # Insure the data is masked
        gdata = gdata.where(mask == 1)
        list_of_gdata.append(gdata)
    griddata = xr.concat(list_of_gdata, dim='simulation')
    griddata['simulation'] = simulations
    return griddata


def open_nemo_griddata_from_zarr(config_file, simulations=None,
                                 grid='U', chunks=None):
    simulations = _check_simulations(config_file, simulations)
    list_of_griddata = []
    for sim in simulations:
        zarr_path = get_nemo_zarr_folder_from_config(config_file, sim, grid=grid)
        list_of_griddata.append(xr.open_zarr(zarr_path))
    griddata = xr.concat(list_of_griddata, dim='simulation')
    griddata['simulation'] = simulations
    return griddata.chunk(chunks=chunks)


# WRF Section
# ------------------------------------------------------------

def get_wrf_meshgrid(config_file):
    cfg = read_config_file(config_file)
    meshgrid_file = cfg['WRF']['MeshGridFile']
    meshgrid = xr.open_dataset(meshgrid_file).squeeze()
    return meshgrid


def open_wrf_griddata_from_netcdf(config_file, simulations=None,
                                  prefix='wrfout', domain='d01', 
                                  variables=None, year=None,
                                  grid=True, **kwargs):
    """
    Open WRF netCDF from the configuration file
    
    Parameters
    ----------
    config_file : str
        The path of the configuration file
    simulation : str or list of str, optional
        Names of simulations to open as written in the configuration file
    prefix : str, optional
        The prefix of WRF output files (usually 'wrfout' or 'wrfhrly'). 
    domain : str, optional
        The name of the WRF domain ('d01' is the parent domain by default).
    variables : list of str, optional
        A list of variables to read
    year : int, optional
        A specific year to read
    grid : bool, optional
        If True, metric from the grid file is added to the Dataset
    **kwargs : optional
        Additional arguments passed on to `xarray.open_mfdataset`
    
    Returns
    -------
    griddata : xarray.Dataset
        The WRF outputs represented by a Dataset
    """
    cfg = read_config_file(config_file)
    if simulations is None:
        simulations = [key for key in cfg 
                       if key not in ['DEFAULT', 'GENERAL', 'WRF', 'NEMO'] ]
    elif not isinstance(simulations, (list, tuple)):
        simulations = [simulations,]
    list_of_gdata = []
    for sim in simulations:        
        gdata = _open_wrf_griddata_from_netcdf(config_file, 
                                               sim, 
                                               prefix=prefix, 
                                               domain=domain, 
                                               variables=variables, 
                                               year=year,
                                               grid=grid, 
                                               **kwargs)
        list_of_gdata.append(gdata)
    griddata = xr.concat(list_of_gdata, dim='simulation')
    griddata['simulation'] = simulations
    return griddata


def _open_wrf_griddata_from_netcdf(config_file, simulation,
                                   prefix='wrfout', domain='d01', 
                                   variables=None, year=None,
                                   grid=True, **kwargs):
    cfg = read_config_file(config_file)
    main_path = cfg[simulation]['Path']    
    year_start = cfg['GENERAL']['YearStart']
    year_stop = cfg['GENERAL']['YearStop']    
    try:
        wrf_path = cfg['WRF']['NetCDFPath']
    except KeyError:
        wrf_path = './'
    try:
        wrf_freq = cfg['WRF']['Frequency']
    except:
        wrf_freq = '*'
    try:
        wrf_prefix = cfg['WRF']['Prefix']
    except:
        wrf_prefix = prefix
    if year is None:
        year = '*'
    filenames = '%s_%s_%s-*' %(wrf_prefix, domain, year)
    full_filenames = main_path + wrf_path + filenames
    ds = xr.open_mfdataset(full_filenames, concat_dim='Time',  
                           decode_coords=False, decode_times=False,
                           decode_cf=False, 
                           drop_variables=['XLON', 'XLAT'], 
                           **kwargs)
    ds = xr.decode_cf(ds)
    time = pd.to_datetime(ds.Times.load().astype('str'), format='%Y-%m-%d_%H:%M:%S')
    ds = ds.assign_coords(Time=time).rename({'Time': 'time'})
    if variables is not None:
        ds = ds[variables]
    if grid:
        meshgrid = get_wrf_meshgrid(config_file)
        ds = ds.assign_coords(XLAT=meshgrid.XLAT_M, 
                              XLONG=meshgrid.XLONG_M, 
                              LANDMASK=meshgrid.LANDMASK)
    return ds


def _open_wrfout_plev_from_netcdf(config_file, simulation, plev, **kwargs):
    cfg = read_config_file(config_file)
    main_path = cfg[simulation]['Path']    
    year_start = cfg['GENERAL']['YearStart']
    year_stop = cfg['GENERAL']['YearStop']    
    try:
        wrf_path = cfg['WRF']['NetCDFPath']
    except KeyError:
        wrf_path = './'
    try:
        wrf_freq = cfg['WRF']['Frequency']
    except:
        wrf_freq = '*'
    try:
        wrf_prefix = cfg['WRF']['Prefix']
    except:
        prefix = '*'     
    list_of_griddata = []  
    # First level is at the 10 meter high
    filenames = '%s_%s_%s.nc' %(wrf_prefix, wrf_freq, plev)
    full_filenames = main_path + wrf_path + filenames
    griddata_plev = xr.open_mfdataset(full_filenames, 
                                      concat_dim='Time', 
                                      **kwargs)
    griddata_plev = griddata_plev.sel(Time=slice(year_start, year_stop))
    return griddata_plev


# NETCDF to ZARR section
#------------------------
def netcdf_to_zarr(config_file, simulations=None,
                   variables=None,
                   nemo=True, wrf=True, 
                   nemo_grids=['U', 'V', 'T'],
                   nemo_chunks={'time_counter': 500},
                   wrf_chunks={'time': 500}, 
                   overwrite=False, compute=True):
    from zarr import Blosc
    compressor = Blosc(cname='zstd', clevel=3, shuffle=0)
    cfg = read_config_file(config_file)
    simulations = _check_simulations(config_file, simulations)
    list_of_res = []
    for sim in simulations:
        if nemo:
            for grid in nemo_grids:
                griddata = open_nemo_griddata_from_netcdf(config_file, simulations=sim,
                                                          grid=grid, parallel=True)
                                                          #chunks={'time_counter': 1})
                if variables is not None:
                    griddata = griddata[variables]
                zarr_path = get_nemo_zarr_folder_from_config(config_file, sim, grid=grid)
                #encoding = {var: {'compressor': compressor} for var in griddata.variables}
                encoding = {var: {'compressor': None} for var in griddata.variables}
                print(zarr_path)
                # Rechunk and save to the zarr format
                if not os.path.isdir(zarr_path):
                    res = griddata.chunk(nemo_chunks).to_zarr(zarr_path, mode='w-', 
                                                              encoding=encoding, 
                                                              compute=compute)
                elif os.path.isdir(zarr_path) and overwrite:
                    res = griddata.chunk(nemo_chunks).to_zarr(zarr_path, mode='w',
                                                              encoding=encoding,
                                                              compute=compute)
                else:
                    print("Skipping the conversion to zarr of grid %s of" 
                          " simulation %s because the folder aready "
                          "exits." %(grid, sim))
                list_of_res.append(res)
        if wrf:
            raise NotImplementedError
    return list_of_res


def _rename_wrfout_plev_variables(griddata, 
                                  list_of_variables=['U', 'V', 'W', 
                                                     'T', 'Vort', 'MoCo', 
                                                     'The', 'RH']):
    for var_origin in griddata.variables:
        for var in list_of_variables:
            if '%s_'%var in var_origin:
                griddata.rename({var_origin: var}, inplace=True)

                
def to_postprocess(griddata, config_file, output_name, type='netcdf'):
    cfg = read_config_file(config_file)
    postprocess_path = cfg['GENERAL']['PostProcessPath']
    if type == 'netcdf':
        griddata.to_netcdf(postprocess_path + output_name)
    elif type == 'zarr':
        griddata.to_zarr(postprocess_path + output_name)

		
def open_postprocess(config_file, name, chunks=None):
	cfg = read_config_file(config_file)
	postprocess_path = cfg['GENERAL']['PostProcessPath']
	try:
		griddata = xr.open_dataset(postprocess_path + name, chunks=chunks)
	except IOError:
		griddata = xr.open_zarr(postprocess_path + name).chunk(chunks)
	return griddata
