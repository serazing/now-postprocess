import now.config as cfg
import now.io as io
import os

if not os.path.isdir(cfg.POSTPROCESS_PATH + 'MEAN'):
    os.mkdir(cfg.POSTPROCESS_PATH + 'MEAN')
    
# Open grids
for grid in cfg.GRIDS:
    griddata = io.open_nemo_griddata_from_zarr(grid=grid)
    griddata_mean = griddata.mean('time_counter')
    griddata_mean.to_netcdf('%sMEAN/%s_MEAN_%s-%s_grid_%s_%s.nc' \
                            %(cfg.POSTPROCESS_PATH,
                              cfg.POSTPROCESS_PREFIX,
                              cfg.YEAR_START, cfg.YEAR_STOP,
                              grid, cfg.DIMS)
                            )