import xarray as xr
import os
from oocgcm.oceanmodels.nemo import grids
import pandas as pd
#------------------------------------------------------------
# Config file sectopm
#------------------------------------------------------------
def read_config_file(config_file):
    import configparser
    from configparser import ConfigParser
    cfg = ConfigParser(interpolation=configparser.ExtendedInterpolation())
    cfg.read(config_file)
    return cfg


def get_simulations(config_file):
    return _check_simulations(config_file, None)


def get_nemo_grid_from_config(config_file):
    cfg = read_config_file(config_file)
    nemo_grid = grids.nemo_2d_grid(
        nemo_coordinate_file=cfg['NEMO']['MeshGridFile'],
        nemo_byte_mask_file=cfg['NEMO']['MaskFile'])
    return nemo_grid


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


def get_nemo_meshgrid(config_file, grid='U'):
    cfg = read_config_file(config_file)
    meshgrid_file = cfg['NEMO']['MeshGridFile']
    dims = cfg['NEMO']['Dimensions']
    meshgrid = xr.open_dataset(meshgrid_file).squeeze()
    if grid is 'U':
        meshgrid = meshgrid[['e1u', 'e2u', 'e3u']]
    elif grid is 'V':
        meshgrid = meshgrid[['e1v', 'e2v', 'e3v']]
    elif grid is 'T':
        meshgrid = meshgrid[['e1t', 'e2t', 'e3t']]
    if dims == '2D':
        meshgrid = meshgrid.isel(z=0)
    return meshgrid


def _check_simulations(config_file, simulations):
    cfg = read_config_file(config_file)
    if simulations is None:
        simulations = [key for key in cfg 
                       if key not in ['DEFAULT', 'GENERAL', 'WRF', 'NEMO']]
    elif not isinstance(simulations, (list, tuple)):
        simulations = [simulations,]
    return simulations


PATH_DICT = {'raw': 'RawPath',
             'tmp': 'TmpPath',
             'climatology': 'ClimPath'}
OA_DICT = {'nemo': 'NEMO',
           'wrf': 'WRF'}


class Cursor:

    def __init__(self, config_file, which='nemo',
                 simulation=None, grid='T',
                 ystart=None, ystop=None, where='raw'):
        self.cfg = read_config_file(config_file)
        self.config_file = config_file
        if simulation is None:
            self.simulation = get_simulations(config_file)[0]
        else:
            self.simulation = simulation
        if ystart is None:
            self.ystart = self.cfg['GENERAL']['YearStart']
        else:
            self.ystart = ystart
        if ystop is None:
            self.ystop = self.cfg['GENERAL']['YearStop']
        else:
            self.ystop = ystop
        self.which = which
        self.grid = grid
        self.where = where
        self.nemo_mask = get_nemo_mask(config_file, grid=self.grid)
        self.nemo_grid = get_nemo_meshgrid(config_file, grid=self.grid)
        self.wrf_grid = get_wrf_meshgrid(config_file)
        self.set_path()
        self.set_basename()

    def set_nemo_grid(self):
        self.nemo_mask = get_nemo_mask(config_file, grid=grid)
        self.nemo_grid =

    def set_path(self):
        self.path = '%s/%s/%s/' % (self.cfg['GENERAL']['WorkDir'],
                                   self.cfg[self.simulation]['SimName'],
                                   self.cfg[OA_DICT[self.which]][PATH_DICT[self.where]])

    def set_basename(self):
        if self.which == 'nemo':
            prefix = self.cfg[self.simulation]['NemoPrefix']
            freq = self.cfg['NEMO']['Frequency']
            dim = self.cfg['NEMO']['Dimensions']
            if self.where == 'raw':
                self.basename =  '%s_%s_*_*_grid_%s_%s.nc' %(prefix, freq,
                                                          self.grid, dim)
            else:
                self.basename =  '%s_%s_%s_%s_grid_%s_%s' %(prefix, freq,
                                                            self.ystart,
                                                            self.ystop,
                                                            self.grid,
                                                            dim)
        elif self.which == 'wrf':
            prefix = self.cfg['WRF']['Prefix']
            dim = self.cfg['NEMO']['Dimensions']
            domain = 1
            if self.where == 'raw':
                self.basename = '%s_d%02i_????-??-??_00:00:00' % (prefix,
                                                                  domain)
            else:
                self.basename = '%s_d%02i_%s-%s' % (prefix, domain,
                                                    self.ystart, self.ystop,)

    def set(self, **kwargs):
        if 'which' in kwargs:
            self.which = kwargs['which']
        if 'simulation' in kwargs['which']:
            self.simulation = kwargs['simulation']
        if 'grid' in kwargs['grid']:
            self.grid = kwargs['grid']
        if 'ystart' in kwargs['ystart']:
            self.ystart = kwargs['ystart']
        if 'ystop' in kwargs['ystop']:
            self.ystop = kwargs['ystop']
        if 'where' in kwargs['where']:
            self.where = kwargs['where']

    def read(self, **kwargs):
        if self.where == 'raw':
            if self.which == 'nemo':
                depth_dims = {'U': 'depthu', 'V': 'depthv', 'T': 'deptht'}
                gdata = xr.open_mfdataset(self.path + self.basename, **kwargs)
                gdata = gdata.set_index(time_counter='time_average_1d')
                if self.cfg['NEMO']['Dimensions'] == '3D':
                    gdata = gdata.rename({depth_dims[self.grid]: 'z'})
                gdata = gdata.sel(time_counter=slice(self.ystart, self.ystop))
                gdata = gdata.where(self.nemo_mask == 1)
            elif self.which == 'wrf':
                pass
        elif self.where == 'tmp':
            gdata = xr.open_zarr(self.path + self.basename + '.zarr', **kwargs)
        elif self.where == 'climatology':
            pass
        return gdata

    def write(self, gdata, analysis='analysis', where=None, overwrite=False):
        if where is None:
            where = self.where
        if where == 'raw':
            raise ValueError('Cannot write in the raw directory')
        elif where == 'tmp':
            gdata.to_zarr(self.path + self.basename + '_%s.zarr' % analysis)
        elif where == 'climatology':
            if self.which == 'nemo':
                filename = self.cfg['NemoPrefix']
            elif self.which == 'wrf':
                pass


class DataBase:

    def __init__(self, config_file):
        self.cfg = read_config_file(config_file)
        self.config_file = config_file
        self.cursor = Cursor(config_file)

    def open(self, simulations=None, grids=None, which='nemo'):
        list_of_gdata = []
        for sim in simulations:
            self.cursor.set(simulation=sim)
            ds = self.cursor.read()

    def netcdf_to_zarr(self, simulations=None, grid=None, which='nemo',
                       chunks={}):
        from zarr import Blosc
        compressor = Blosc(cname='zstd', clevel=3, shuffle=0)
        simulations = _check_simulations(self.config_file, simulations)
        self.cursor.set(which=which)
        for sim in simulations:
            for grid in grids:
                self.cursor.set(simulation=sim, grid=grid)
                gdata = self.cursor.read(parallel=True)
                encoding = {var: {'compressor': compressor}
                            for var in gdata.variables}
                self.cursor.set(where='tmp')
                self.cursor.write(gdata, encoding=encoding)

def open_nemo_griddata_from_netcdf(simulations=None,
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
    list_of_gdata = []
    # Get the corresponding mask
    mask = get_nemo_mask(config_file, grid=grid)
    for sim in simulations:
        filenames = get_nemo_filenames_from_config(config_file, sim, grid=grid)
        print(filenames)
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


def netcdf_to_zarr(config_file, simulations=None,
                   variables=None,
                   nemo=True, wrf=True, 
                   nemo_grids=['U', 'V', 'T'],
                   nemo_chunks={'time_counter': 300},
                   wrf_chunks={'time': 300}, 
                   overwrite=False, compute=True):
    from zarr import Blosc
    compressor = Blosc(cname='zstd', clevel=3, shuffle=0)
    cfg = read_config_file(config_file)
    simulations = _check_simulations(config_file, simulations)
    list_of_res = []
    for sim in simulations:
        if nemo:
            for grid in nemo_grids:
                griddata = open_nemo_griddata_from_netcdf(config_file,
                                                          simulations=sim,
                                                          grid=grid,
                                                          parallel=True)
                                                          #chunks={'time_counter': 1})
                if variables is not None:
                    griddata = griddata[variables]
                zarr_path = get_nemo_zarr_folder_from_config(config_file,
                                                             sim, grid=grid)
                encoding = {var: {'compressor': compressor}
                            for var in griddata.variables}
                # Rechunk and save to the zarr format
                if not os.path.isdir(zarr_path):
                    res = griddata.chunk(nemo_chunks).to_zarr(zarr_path,
                                                              mode='w-',
                                                              encoding=encoding, 
                                                              compute=compute)
                elif os.path.isdir(zarr_path) and overwrite:
                    res = griddata.chunk(nemo_chunks).to_zarr(zarr_path,
                                                              mode='w',
                                                              encoding=encoding,
                                                              compute=compute)
                else:
                    print("Skipping the conversion to zarr of grid %s of" 
                          " simulation %s because the folder aready "
                          "exits." % (grid, sim))
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


def save_to(config_file, griddata, analysis,
            where='climatology', which='nemo', engine='netcdf'):
    sim = griddata.simulation.data
    cfg = read_config_file(config_file)
    path = '%s/%s/' % (cfg['GENERAL']['WorkDir'], cfg[sim]['SimName'])

    if which == 'nemo':
        path += cfg['NEMO']['ClimatoPath']
    elif which == 'wrf':
    #+ cfg['NEMO']['ClimatoPath'] + output_name + '.nc'
    if engine == 'netcdf':

    elif engine == 'zarr':
        griddata.to_zarr()






#def use_config(config_file):
#    def decorator(func):
#        def call (**args, **kwargs):
