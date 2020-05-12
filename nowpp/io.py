import xarray as xr
import os
import glob
import sh
import pandas as pd
from .models import nemo, wrf


def read_config_file(config_file):
    from yaml import load, Loader
    cfg = load(open(config_file, 'r'), Loader)
    return cfg


def read_config_file_old(config_file):
    import configparser
    from configparser import ConfigParser
    cfg = ConfigParser(interpolation=configparser.ExtendedInterpolation())
    cfg.read(config_file)
    return cfg


class Config:

    def __init__(self, file):
        self.conf = read_config_file(file)

    def get_start_date(self):
        return pd.to_datetime(str(self.conf['date']['start']))

    def get_end_date(self):
        return pd.to_datetime(str(self.conf['date']['end']),
                              format="%Y-%m-%d")

    def get_models(self):
        return list(self.conf['models'])

    def get_simulations(self):
        simulations = {sim['name']: sim for sim in self.conf['simulations']}
        return simulations

    def get_model_simulations(self, model):
        simulations = {sim['name']: sim for sim in self.conf['simulations']
                       if model in sim['models']}
        return simulations

    def get_simulation_models(self, simulation):
        simulations = self.get_simulations()
        return simulations[simulation]['models']

    def check_simulations(self, simulations, model):
        config_simulations = self.get_model_simulations(model)
        if simulations is None:
            simulations = config_simulations
        else:
            simulations = [sim for sim in config_simulations
                           if sim in simulations]
        return simulations

    def get_model_grids(self, model):
        grids = self.conf['models'][model]['grids']
        return grids

    def get_model_directories(self, model):
        directories = self.conf['models'][model]['directories']
        return directories

    def get_model_mask(self, model, grid):
        file = self.conf['models'][model]['files']['mask']
        if model == 'nemo':
            mask = nemo.read_mask(file, grids=[grid, ])
        elif model == 'wrf':
            mask = wrf.read_mask(file)
        else:
            raise ValueError
        return mask

    def get_model_mesh(self, model, grid):
        file = self.conf['models'][model]['files']['mesh']
        if model == 'nemo':
            mesh = nemo.read_mesh(file, grids=[grid, ])
        elif model == 'wrf':
            mesh = wrf.read_mesh(file)
        else:
            raise ValueError
        return mesh

    def get_model_mesh_file(self, model):
        file = self.conf['models'][model]['files']['mesh']
        return file

    def get_workdir(self):
        return self.conf['general']['directories']['work']

    def get_simulation_folder(self, simulation):
        simulations = self.get_simulations()
        return simulations[simulation]['folder']

    def get_model_directory(self, model, where):
        return self.conf['models'][model]['directories'][where]

    def get_path(self, model, simulation, where):
        path = os.path.join(self.get_workdir(),
                            self.get_simulation_folder(simulation),
                            model,
                            self.get_model_directory(model, where))
        return path

    def get_mdss_path(self, simulation):
        simulations = self.get_simulations()
        return simulations[simulation]['mdss']['directory']

    def get_mdss_project(self, simulation):
        simulations = self.get_simulations()
        return simulations[simulation]['mdss']['project']

    def get_basename(self, model, simulation, where, grid):
        simulations = self.get_model_simulations(model)
        ystart = self.get_start_date().year
        ystop = self.get_end_date().year
        if model == 'nemo':
            prefix = simulations[simulation]['prefix']
            freq = self.conf['models'][model]['frequency']
            dim = self.conf['models'][model]['dimensions']
            if where == 'raw':
                basename = '%s_%s_*_*_grid_%s_%s.nc' % (prefix, freq, grid,
                                                        dim)
            else:
                basename = '%s_%s_%s_%s_%s' % (prefix, freq, ystart,
                                               ystop, dim)
        elif model == 'wrf':
            domain = 1
            prefix = self.conf['models'][model]['prefix']
            if where == 'raw':
                basename = '%s_d%02i_????-??-??_00:00:00' % (prefix, domain)
            else:
                basename = '%s_d%02i_%s-%s' % (prefix, domain, ystart,
                                               ystop,)
        else:
            raise ValueError
        return basename


class Cursor:

    def __init__(self, config_file, model=None, simulation=None, grid=None,
                 where='raw', time_slice=None):
        self.config_file = config_file
        self.cfg = Config(self.config_file)
        self.model = None
        self.simulation = None
        self.grid = None
        self.where = where
        self.time_slice = time_slice
        self.path = None
        self.basename = None
        self.mesh_file = None
        # Apply methods to update previous values
        self._update_model(model)
        self._update_simulation(simulation)
        self._update_grid(grid)
        self._full_update()

    def _update_model(self, model):
        models = self.cfg.get_models()
        if model is None:
            self.model = models[0]
        elif model in models:
            self.model = model
        else:
            raise ValueError("Cannot recognise this type of model")

    def _update_simulation(self, simulation):
        simulations = list(self.cfg.get_model_simulations(self.model))
        if simulation is None:
            self.simulation = simulations[0]
        elif simulation in simulations:
            self.simulation = simulation
        else:
            raise ValueError("Cannot use this simulation name")

    def _update_grid(self, grid):
        grids = self.cfg.get_model_grids(self.model)
        if grid is None:
            self.grid = grids[0]
        elif grid in grids:
            self.grid = grid
        else:
            raise ValueError("Cannot use this grid type. Available values "
                             "are: %s " % grids)

    def _update_where(self, where):
        directories = self.cfg.get_model_directories(self.model)
        if where is None:
            self.directories = directories[0]
        elif where in directories:
            self.where = where
        else:
            raise ValueError("Cannot use this directory")

    def _update_path(self):
        self.path = self.cfg.get_path(self.model, self.simulation, self.where)

    def _update_basename(self):
        self.basename = self.cfg.get_basename(self.model, self.simulation,
                                              self.where, self.grid)

    def _update_mask(self):
        try:
            self.mask = self.cfg.get_model_mask(self.model, self.grid)
        except FileNotFoundError:
            self.mask = xr.Dataset()

    def _update_mesh(self):
        self.mesh_file = self.cfg.get_model_mesh_file(self.model)
        # try:
        #    self.mesh = self.cfg.get_model_mesh(self.model, self.grid)
        # except FileNotFoundError:
        #    self.mesh = xr.Dataset()

    def _full_update(self):
        self._update_path()
        self._update_basename()
        self._update_mesh()
        # self._update_mask()

    def sel(self, **kwargs):
        if 'model' in kwargs:
            self._update_model(kwargs['model'])
        if 'simulation' in kwargs:
            self._update_simulation(kwargs['simulation'])
        if 'grid' in kwargs:
            self._update_grid(kwargs['grid'])
        if 'where' in kwargs:
            self._update_where(kwargs['where'])
        self._full_update()

    def read(self, extension='', engine='zarr', variables=None, **kwargs):
        if self.where == 'raw':
            filenames = sorted(glob.glob(os.path.join(self.path,
                                                      self.basename)))
            # Case for opening raw NEMO outputs
            if self.model == 'nemo':
                gdata = nemo.open_netcdf_dataset(filenames, 
                                                 mesh_file=self.mesh_file, 
                                                 grid=self.grid, **kwargs)
            # Case for opening raw WRF outputs
            elif self.model == 'wrf':
                gdata = wrf.open_netcdf_dataset(filenames,
                                                mesh_file=self.mesh_file,
                                                variables=variables,
                                                **kwargs)
            else:
                raise ValueError("Cannot recognise this type of model")
            # Assign a new dimension corresponding to the simulations
            sim_coord = xr.DataArray([self.simulation, ], dims='simulation')
            gdata = gdata.assign(simulation=sim_coord)
        else:
            if engine == 'zarr':
                zarr_folder = self.basename + '%s.zarr' % extension
                zarr_path = os.path.join(self.path, zarr_folder)
                gdata = xr.open_zarr(zarr_path, **kwargs)
            elif engine == 'netcdf':
                filename = self.basename + '%s.nc' % extension
                netcdf_path = os.path.join(self.path, filename)
                gdata = xr.open_dataset(netcdf_path, **kwargs)
            else:
                raise ValueError("This engine is not supported")
        return gdata

    def write(self, gdata, extension='', where=None, chunks=None,
              engine='zarr', **kwargs):
        gdata = gdata.chunk(chunks)
        if where is None:
            where = self.where
        elif where == 'raw':
            raise ValueError('Cannot write in the raw directory')
        self.sel(where=where)
        sh.mkdir('-p', self.path)
        if engine == 'zarr':
            zarr_folder = self.basename + '%s.zarr' % extension
            zarr_path = os.path.join(self.path, zarr_folder)
            gdata.to_zarr(zarr_path, **kwargs)
        elif engine == 'netcdf':
            filename = self.basename + '%s.nc' % extension
            netcdf_path = os.path.join(self.path, filename)
            gdata.to_netcdf(netcdf_path, **kwargs)
        else:
            raise ValueError("This engine is not supported")

    def __repr__(self):
        message = ("Cursor \n" 
                   "====== \n"
                   "|-->SIMULATION: %s \n"
                   "|------->MODEL: %s \n" 
                   "|------->WHERE: %s \n"
                   "|-------->GRID: %s \n" 
                   "\n" 
                   "Filenames to read/write:\n" 
                   "%s" % (self.model, self.simulation, self.where,
                           self.grid, os.path.join(self.path, self.basename))
                   )
        return message


class DataBase:

    def __init__(self, config_file, model=None, simulations=None,
                 grids=None, decode_xgrids=False):
        self.cfg = Config(config_file)
        self.cs = Cursor(config_file)
        if model is None:
            self.model = self.cfg.get_models()[0]
        else:
            self.model = model
        if simulations is None:
            self.simulations = self.cfg.get_model_simulations(self.model)
        else:
            self.simulations = simulations
        if grids is None:
            self.grids = self.cfg.get_model_grids(self.model)
        else:
            self.grids = grids
        if decode_xgrids:
            raise NotImplementedError

    def __repr__(self):
        simulations = self.cfg.get_simulations()
        models = self.cfg.get_models()
        nb_simulations = len(simulations)
        nb_models = len(models)
        message = ("Database (simulations: %i, models : %i) \n \n" 
                   "Simulations: %s \n"
                   "Models: %s \n \n"
                   % (nb_simulations, nb_models, list(simulations), models)
                   )
        wkdir = self.cfg.get_workdir()
        wkdir_message = "%s / \n" % wkdir
        message += wkdir_message
        for sim in self.cfg.get_simulations():
            sim_folder = self.cfg.get_simulation_folder(sim)
            sim_message = "    | -- %s (%s)/ \n" % (sim_folder, sim)
            message += sim_message
            for model in self.cfg.get_simulation_models(sim):
                model_message = "    | | -- %s / \n" % model
                message += model_message
                for where in self.cfg.get_model_directories(model):
                    if (self.cs.simulation == sim and self.cs.model == model
                            and self.cs.where == where):
                        where_message = "--->| | | -- %s / \n" % where
                    else:
                        where_message = "    | | | -- %s / \n" % where
                    message += where_message

                    path = self.cfg.get_path(model, sim, where)
                    try:
                        filenames = os.listdir(path)
                    except FileNotFoundError:
                        filenames = []
                    if 0 < len(filenames) <= 10:
                        for fname in filenames:
                            file_message = "    | | | | -- %s  \n" % fname
                            message += file_message
                    elif len(filenames) > 10:
                        file_message = ("    | | | | -- %s  \n"
                                        "    | | | | -- %s  \n"
                                        "    | | | | -- ...  \n"
                                        "    | | | | -- %s  \n"
                                        "    | | | | -- %s  \n"
                                        % (filenames[0],  filenames[1],
                                           filenames[-2],  filenames[-1])
                                        )
                        message += file_message
                    else:
                        pass
        return message

    def _combine_grids(self, grids,  **kwargs):
        list_of_gdata = []
        for grid in grids:
            self.cs.sel(grid=grid)
            gdata_grid = self.cs.read(**kwargs)
            list_of_gdata.append(gdata_grid)
        gdata = xr.merge(list_of_gdata)
        return gdata

    def open(self, simulations=None, grids=None, model='nemo', where='raw',
             extension='', variables=None, **kwargs):
        if model is None:
            model = self.model
        if simulations is None:
            simulations = self.cfg.get_model_simulations(model)
        list_of_gdata = []
        for sim in simulations:
            self.cs.sel(model=model, simulation=sim, where=where)
            if model == 'nemo' and where == 'raw':
                gdata = self._combine_grids(grids, **kwargs)
            else:
                gdata = self.cs.read(extension=extension,
                                     variables=variables, **kwargs)
            list_of_gdata.append(gdata)
        return xr.concat(list_of_gdata, dim='simulation')

    def netcdf_to_zarr(self, simulations=None, grids=None, model='nemo',
                       chunks=None, read_kwargs=None, write_kwargs=None):
        if write_kwargs is None:
            write_kwargs = {}
        if read_kwargs is None:
            read_kwargs = {}
        if grids is None:
            grids = self.grids
        if simulations is None:
            simulations = self.simulations
        if grids is None:
            grids = self.grids
        if 'compressor' not in write_kwargs:
            from zarr import Blosc
            compressor = Blosc(cname='zstd', clevel=3, shuffle=0)
        else:
            compressor = None
        self.cs.sel(model=model)
        for sim in simulations:
            self.cs.sel(model=model, simulation=sim, where='raw')
            if model == 'nemo':
                gdata = self._combine_grids(grids, **read_kwargs)
            else:
                gdata = self.cs.read(**read_kwargs)
            encoding = {var: {'compressor': compressor}
                        for var in gdata.variables}
            print(gdata)
            self.cs.sel(where='tmp')
            self.cs.write(gdata, chunks=chunks, encoding=encoding,
                          **write_kwargs)

    def mean(self, dim=None, write=False, engine='netcdf'):
        @apply_to_database(self, engine=engine,
                           extension='_mean', write=write)
        def func(ds, **kwargs):
            return ds.mean(**kwargs)
        return func(dim=dim)

    def std(self, dim=None, write=False, engine='netcdf'):
        @apply_to_database(self, engine=engine,
                           extension='_std', write=write)
        def func(ds, **kwargs):
            return ds.std(**kwargs)
        return func(dim=dim)

    def seasonal_cycle(self, freq='dayofyear', write=False):
        @apply_to_database(self, engine='zarr',
                           extension='_seasonal_cycle', write=write)
        def func(ds):
            return ds.groupby('time.%s' % freq).mean('time')
        return func()

    def apply(self, func, simulations=None, model=None, **kwargs):
        if model is None:
            self.model = self.cfg.get_models()[0]
        else:
            self.model = model
        if simulations is None:
            self.simulations = self.cfg.get_model_simulations(self.model)
        else:
            self.simulations = simulations
        db_func = apply_to_database(self,
                                    simulations=simulations,
                                    model=model)(func)
        return db_func(**kwargs)

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
        raise NotImplementedError
        # if model is None:
        #    model = self.model[0]
        # xgrid = self.xgrids[model]
        # return xgrid

    def write(self, ds, extension='', engine='netcdf', **kwargs):
        ds_dict = dict(ds.groupby('simulation'))
        for sim in ds_dict:
            self.cs.sel(simulation=sim, where='tmp')
            self.cs.write(ds_dict[sim], extension=extension, engine=engine,
                          **kwargs)


def apply_to_database(db, **options):
    def decorator(func):
        def call(*args, **kwargs):
            # Checking kwargs
            if 'model' in options:
                model = options['model']
            else:
                model = db.model
            if 'simulations' in options:
                simulations = options['simulations']
            else:
                simulations = db.simulations
            if 'engine' in options:
                engine = options['engine']
            else:
                save_engine = 'netcdf'
            if 'extension' in options:
                extension = options['extension']
            else:
                extension = '_%s' % func.__name__
            if 'write' in options:
                write = options['write']
            else:
                write = False
            # Loop over simulations
            gdata = db.open(model=model, simulations=simulations, where='tmp')
            print("Applying %s on %s outputs" % (extension, model))
            res = func(gdata, **kwargs)
            if write:
                db.write(res, engine=engine, extension=extension)
            return res
        return call
    return decorator
