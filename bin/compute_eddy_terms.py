import xarray as xr
import os
import now.config as cfg
import now.io as io
import now.physics as phy


if not os.path.isdir(cfg.POSTPROCESS_PATH + 'EDDY'):
    os.mkdir(cfg.POSTPROCESS_PATH + 'EDDY')

# Open datasets from zarr folders
gridU = io.open_nemo_griddata_from_zarr(grid='U')
gridV = io.open_nemo_griddata_from_zarr(grid='V')
gridU_mean = xr.open_dataset('%sMEAN/%s_MEAN_%s-%s_grid_%s_%s.nc' \
                            %(cfg.POSTPROCESS_PATH,
                              cfg.POSTPROCESS_PREFIX,
                              cfg.YEAR_START, cfg.YEAR_STOP,
                              'U', cfg.DIMS),
                             chunks={'simulation':1}
                            )
gridV_mean = xr.open_dataset('%sMEAN/%s_MEAN_%s-%s_grid_%s_%s.nc' \
                            %(cfg.POSTPROCESS_PATH,
                              cfg.POSTPROCESS_PREFIX,
                              cfg.YEAR_START, cfg.YEAR_STOP,
                              'V', cfg.DIMS),
                             chunks={'simulation':1})

ux_eddy = gridU['uos'] - gridU_mean['uos']
vy_eddy = gridV['vos'] - gridV_mean['vos']
taux_eddy = gridU['tauuo'] - gridU_mean['tauuo']
tauy_eddy = gridV['tauvo'] - gridV_mean['tauvo']

# Computation of Eddy Kinetic Energy
eke = phy.compute_ke(ux_eddy, vy_eddy).mean(dim='time_counter')
eke.attrs['long_name'] = "Eddy Kinetic Energy"
# Computation of Eddy Wind Work
eww = phy.compute_wind_work(ux_eddy, vy_eddy, 
                            taux_eddy, tauy_eddy).mean(dim='time_counter')
eww.attrs['long_name'] = "Eddy Wind Work"
#Create a DataSet and save the resuls to a netCDF
eddy_terms = xr.Dataset({'eke': eke, 'eww': eww})
eddy_terms.to_netcdf('%sEDDY/%s_EDDY_%s_%s-%s_grid_%s_%s.nc' \
                      %(cfg.POSTPROCESS_PATH,
                        cfg.POSTPROCESS_PREFIX,
                        cfg.FREQ, cfg.YEAR_START, 
                        cfg.YEAR_STOP, 'T', cfg.DIMS)
                    )