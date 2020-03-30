from . import io
import xarray as xr

def loop_over_dataset(config_file, **options):
    def decorator(func):
        def call(*args, **kwargs):
            output_name = 'analysis'
            cfg = io.read_config_file(config_file)
            # Checking kwargs
            if 'simulations' in options:
                simulations = io._check_simulations(config_file, options['simulations'])
            else:
                simulations = io.get_simulations(config_file)
            if 'model' in options:
                model = options['model']
            else:
                model = 'nemo'
            if 'grids' in options:
                grids = options['grids']
            else:
                grids = ['T']
            if 'save_engine' in options:
                save_engine = options['save_engine']
            else:
                save_engine='netcdf'
            # Loop over simulations
            for sim in simulations:
                for grid in grids:
                    db.cursor.set(simulation=sim, model=model, grid=grid)
                    print("Applying %s on %s outputs (%s, grid %s)" % (func.__name__, model, sim, grid))
                    gdata = db.cursor.read()
                    res = func(griddata, **kwargs)
                    db.cursor.set(where='climatology')
        return call
    return decorator




def open(config_file):
    griddata = io.open_nemo_griddata_from_zarr(config_file,
                                               simulations=[sim],
                                               grid=grid)
    return griddata

def save(config_file):
    pass

def decorate(func):
    griddata = open()
    def call(**kwargs):
        res = func(griddata, **kwargs)
        save(res, func)
    return call


mean_method = getattr(xr.Dataset, 'mean')
decorate(mean_method, dim='time')

def xarray_analysis(ds):
    ds = xr.Dataset.isel(ds, x=215)
    ds = xr.Dataset.mean(ds, dim='time')
    return ds

def get_griddata(config_file, )

    return griddata

def wrap(func):
    simulations = io.get_simulations(config_file)
    cfg = io.read_config_file(config_file)
    for sim in simulations:
        if nemo:
            for grid in nemo_grids:
                griddata = io.open_nemo_griddata_from_zarr(config_file,
                                                           simulations=[sim],
                                                           grid=grid)
                griddata.func(**kwargs)
                griddata_mean = griddata.func('time_counter')
                output_name = 'nemo_1990-2009_grid_%s_2D_mean.nc' % grid
                dest = cfg['NEMO']['ClimPath']
        if wrf:
            raise NotImplementedError
     pass


def apply_func_to_config(func_name):



def apply_to_config(config_file, which='nemo', nemo_grids=['U', 'V', 'T']):
    def decorate(func):
        def call(**kwargs):
            simulations = io.get_simulations(config_file)
            cfg = io.read_config_file(config_file)
            for sim in simulations:
                if which == 'nemo':
                    for grid in nemo_grids:
                        griddata = io.open_nemo_griddata_from_zarr(config_file,
                                                                   simulations=[sim],
                                                                   grid=grid)
                        res = func(griddata, **kwargs)
                elif which == 'wrf':
                    raise NotImplementedError
                elif which == 'both':
                    raise NotImplementedError
            return res
        return call
    return decorate


def nemo_mean(config_file, nemo_grids=['U', 'V', 'T'], **kwargs):
    @apply_to_config(config_file, which='nemo', nemo_grids=nemo_grids,
                     output_type='netcdf')
    def dataset_mean(**kwargs):
        mean_method = getattr(xr.Dataset, 'mean')
        return mean_method(**kwargs)
    return dataset_mean(**kwargs)

    def dataset_mean(**kwargs):
        mean_method = getattr(xr.Dataset, 'mean')
        return mean_method(**kwargs)
    return dataset_mean(**kwargs)

output_name = 'nemo_1990-2009_grid_%s_2D_%s.nc' % grid
path = cfg['NEMO']['ClimPath']
griddata.to_netcdf(path + output_name)

def trend(config_file, ):
    for sim in sim_dict:
        for grid in ['T', 'U', 'V']:
            griddata = nowpp.io.open_nemo_griddata_from_zarr(config_file, simulations=[sim], grid=grid)
            if grid is 'T':
                griddata = griddata.drop(('time_maximum_1d', 'time_minimum_1d'))
            griddata_trend, _ = xfit.linreg(griddata, dim='time_counter')
            output_name = 'nemo_1990-2008_grid_%s_2D_trend_slope.nc' % grid
            griddata_trend.to_netcdf('%s/%s/NEMO/SURFACE_CLIMATOLOGY/%s' % (output_path, sim_dict[sim] , output_name))


def seasonal_cycle(config_file, ):
    pass


def smooth_seasonal_cycle(config_file,):
    pass


def variability(config_file, ):
    pass
