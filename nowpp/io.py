import xarray as xr
import os
#from oocgcm.oceanmodels.nemo import grids
import pandas as pd


NEMO_NEW_DIMS = {'T': {'x': 'x_c', 'y': 'y_c', 'z': 'z_c'},
                 'U': {'x': 'x_r', 'y': 'y_c', 'z': 'z_c'},
                 'V': {'x': 'x_c', 'y': 'y_r', 'z': 'z_c'},
                 'F': {'x': 'x_r', 'y': 'y_r', 'z': 'z_c'},
                 'W': {'x': 'x_c', 'y': 'y_c', 'z': 'z_r'}
                }

NEMO_NEW_COORDS = {'T': {'nav_lon': 'lon_T', 'nav_lat': 'lat_T', 'depth': 'depth_T'},
                   'U': {'nav_lon': 'lon_U', 'nav_lat': 'lat_U', 'depth': 'depth_U'},
                   'V': {'nav_lon': 'lon_V', 'nav_lat': 'lat_V', 'depth': 'depth_V'},
                   'W': {'nav_lon': 'lon_T', 'nav_lat': 'lat_T', 'depth': 'depth_W'}
                  }

PATH_DICT = {'raw': 'RawPath',
             'tmp': 'TmpPath',
             'climatology': 'ClimPath'}

OA_DICT = {'nemo': 'NEMO',
           'wrf': 'WRF'}


WRF_NEW_DIMENSIONS = {}

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


def get_nemo_meshgrid(config_file):
    cfg = read_config_file(config_file)
    meshgrid_file = cfg['NEMO']['MeshGridFile']
    dims = cfg['NEMO']['Dimensions']
    original_meshgrid = xr.open_dataset(meshgrid_file).squeeze()
    list_of_meshgrid = []
    for grid in ['U', 'V', 'T', 'F', 'W']:
        meshgrid = original_meshgrid[list(NEMO_NEW_GRID[grid])].rename(NEMO_NEW_GRID[grid])
        coords = original_meshgrid[list(NEMO_NEW_COORDS[grid])].rename(NEMO_NEW_COORDS[grid])
        meshgrid = meshgrid.assign_coords(coords).rename_dims(NEMO_NEW_DIMS[grid])
        list_of_meshgrid.append(meshgrid)
    new_meshgrid = xr.merge(meshgrid)
    if dims == '2D':
        new_meshgrid = new_meshgrid.isel(z_c=0, z_r=0)
    return new_meshgrid

def _check_simulations(config_file, simulations):
    cfg = read_config_file(config_file)
    if simulations is None:
        simulations = [key for key in cfg 
                       if key not in ['DEFAULT', 'GENERAL', 'WRF', 'NEMO']]
    elif not isinstance(simulations, (list, tuple)):
        simulations = [simulations,]
    return simulations


class Cursor:

    def __init__(self, config_file, model='nemo',
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
        self.model = model
        self.grid = grid
        self.where = where
        self.nemo_mask = get_nemo_mask(config_file, grid=self.grid)
        #self.nemo_grid = get_nemo_meshgrid(config_file, grid=self.grid)
        #self.wrf_grid = get_wrf_meshgrid(config_file)
        self._update_path()
        self._update_basename()

    def _update_path(self):
        self.path = '%s/%s/%s/' % (self.cfg['GENERAL']['WorkDir'],
                                   self.cfg[self.simulation]['SimName'],
                                   self.cfg[OA_DICT[self.model]][PATH_DICT[self.where]])

    def _update_basename(self):
        if self.model == 'nemo':
            prefix = self.cfg[self.simulation]['NemoPrefix']
            freq = self.cfg['NEMO']['Frequency']
            dim = self.cfg['NEMO']['Dimensions']
            if self.where == 'raw':
                self.basename =  '%s_%s_*_*_grid_%s_%s.nc' % (prefix, freq,
                                                          self.grid, dim)
            else:
                self.basename =  '%s_%s_%s_%s_%s' % (prefix, freq,
                                                     self.ystart, self.ystop,
                                                     dim)
        elif self.model == 'wrf':
            prefix = self.cfg['WRF']['Prefix']
            dim = self.cfg['NEMO']['Dimensions']
            domain = 1
            if self.where == 'raw':
                self.basename = '%s_d%02i_????-??-??_00:00:00' % (prefix,
                                                                  domain)
            else:
                self.basename = '%s_d%02i_%s-%s' % (prefix, domain,
                                                    self.ystart, self.ystop,)

    def sel(self, **kwargs):
        if 'model' in kwargs:
            self.model = kwargs['model']
        if 'simulation' in kwargs:
            self.simulation = kwargs['simulation']
        if 'grid' in kwargs:
            self.grid = kwargs['grid']
        if 'ystart' in kwargs:
            self.ystart = kwargs['ystart']
        if 'ystop' in kwargs:
            self.ystop = kwargs['ystop']
        if 'where' in kwargs:
            self.where = kwargs['where']
        self._update_path()
        self._update_basename()
        self.nemo_mask = get_nemo_mask(self.config_file, grid=self.grid)

    def read(self, **kwargs):
        if self.where == 'raw':

            # Case for opening raw NEMO outputs
            if self.model == 'nemo':
                new_dims = NEMO_NEW_DIMS[self.grid]
                new_coords = NEMO_NEW_COORDS[self.grid]
                sim_coord = xr.DataArray([self.simulation, ], dims='simulation')
                if self.cfg['NEMO']['Dimensions'] == '2D':
                    new_dims = {dim: new_dims[dim] for dim in new_dims if dim != 'z'}
                    new_coords = {coord: new_coords[coord] for coord in new_coords if coord != 'depth'}
                gdata = xr.open_mfdataset(self.path + self.basename, **kwargs)
                gdata = gdata.where(self.nemo_mask == 1)
                gdata = gdata.set_index(time_counter='time_average_1d')
                gdata = gdata.rename({'time_counter': 'time'})
                gdata = gdata.rename_dims(new_dims)
                gdata = gdata.rename(new_coords)
                gdata = gdata.sel(time=slice(self.ystart, self.ystop))
                gdata = gdata.expand_dims('simulation')
                gdata = gdata.assign(simulation=sim_coord)
                
            # Case for opening raw WRF outputs
            elif self.model == 'wrf':
                #raise NotImplementedError
                gdata = xr.open_mfdataset(self.path + self.basename,
                                          concat_dim='Time',
                                          decode_coords=False,
                                          decode_times=False,
                                          decode_cf=False,
                                          drop_variables=['XLON', 'XLAT'],
                                          **kwargs)
                gdata = xr.decode_cf(gdata)
                time = pd.to_datetime(gdata.Times.load().astype('str'),
                                      format='%Y-%m-%d_%H:%M:%S')
                gdata = gdata.assign_coords(Time=time).rename({'Time': 'time'})
                gdata = gdata.assign_coords(XLAT=self.wrf_grid.XLAT_M,
                                            XLONG=self.wrf_grid.XLONG_M,
                                            LANDMASK=self.wrf_grid.LANDMASK)

        elif self.where == 'tmp':
            # Case for opening raw NEMO outputs
            gdata = xr.open_zarr(self.path + self.basename + '.zarr', **kwargs)
        elif self.where == 'climatology':
            raise NotImplementedError

        return gdata

    def write(self, gdata, extension='', where=None, chunks={'time': 300},
              engine='zarr', **kwargs):
        gdata = gdata.chunk(chunks)
        if where is None:
            where = self.where
        if where == 'raw':
            raise ValueError('Cannot write in the raw directory')
        else:
            if engine == 'zarr':
                gdata.to_zarr(self.path + self.basename + '%s.zarr' %
                              extension, **kwargs)
            elif engine == 'netcdf':
                gdata.to_netcdf(self.path + self.basename + '%s.nc' %
                                extension, **kwargs)


class DataBase:

    def __init__(self, config_file):
        self.config_file = config_file
        self.cfg = read_config_file(config_file)
        self.cs = Cursor(config_file)
        self.simulations = get_simulations(config_file)
        self.grids = ['U', 'V', 'T', 'W']
        self.models = ['NEMO', 'WRF']
        self.xgrids = None


    def _combine_grids(self, grids,  **kwargs):
        list_of_gdata = []
        for grid in grids:
            self.cs.sel(grid=grid)
            gdata_grid = self.cs.read(**kwargs)
            list_of_gdata.append(gdata_grid)
        gdata = xr.merge(list_of_gdata)
        return gdata

    def open(self, simulations=None, grids=None, model='nemo', where='raw',  **kwargs):
        list_of_gdata = []
        for sim in simulations:
            self.cs.sel(model=model, simulation=sim, where=where)
            if model == 'nemo' and where == 'raw':
                gdata = self._combine_grids(grids, **kwargs)
            else:
                gdata = self.cs.read(**kwargs)
            list_of_gdata.append(gdata)
        return xr.concat(list_of_gdata, dim='simulation')

    def netcdf_to_zarr(self, simulations=None, grids=['T'], model='nemo',
                       chunks={'time': 300},
                       read_kwargs={}, write_kwargs={}):
        if simulations is None:
            simulations = self.simulations
        else:
            simulations = _check_simulations(self.config_file, simulations)
        if grids is None:
            grids = self.grids
        if 'compressor' not in write_kwargs:
            from zarr import Blosc                                            
            compressor = Blosc(cname='zstd', clevel=3, shuffle=0)
        self.cs.sel(model=model)
        for sim in simulations:
            self.cs.sel(model=model, simulation=sim, where='raw')
            if model == 'nemo':
                gdata = self._combine_grids(grids, **read_kwargs)
            else:
                gdata = self.cs.read(**read_kwargs)
            encoding = {var: {'compressor': compressor}
                             for var in gdata.variables
                        }
            self.cs.sel(where='tmp')
            self.cs.write(gdata, chunks=chunks, encoding=encoding, **write_kwargs)

    def apply(self, func, simulations=None, model=None):
        apply_to_database(self, simulations=simulations, model=model)(func)

    def get_model_xgrid(self, model=None):
        """
        Return the XGCM grid operator for computing on the gridded dataset

        Parameters
        ----------
        model : basestring
            The model for which the grid is needed

        Returns
        -------
        xgrid  : xgcm.Grid
            The XGCM grid object for performing computation on the grid
        """
        if model is None:
            model = self.model[0]
        xgrid = self.xgrids[model]
        return xgrid


def apply_to_database(db, **options):
    def decorator(func):
        def call(*args, **kwargs):
            output_name = 'analysis'
            # Checking kwargs
            if 'simulations' in options:
                simulations = options['simulations']
            else:
                simulations = db.simulations
            if 'model' in options:
                model = options['model']
            else:
                model = db.model
            if 'save_engine' in options:
                save_engine = options['save_engine']
            else:
                save_engine = 'netcdf'
            # Loop over simulations
            for sim in simulations:
                db.cursor.sel(simulation=sim, model=model,
                              where='tmp', grid=grid)
                print("Applying %s on %s outputs (%s)" % (func.__name__, model, sim))
                gdata = db.cursor.read()
                res = func(gdata, **kwargs)
                db.cursor.sel(where='climatology')
                db.cursor.write(res, engine=save_engine,
                                extension=func.__name__)
        return call
    return decorator