import xarray as xr
from now import io 

# NAMELIST
NEMO_PATH = '/g/data1/e14/gs9353/NOW_OUTPUTS/'
SIMULATIONS = ['Reference', 'Future', 'No current feedback',
               'Present high-resolution', 'Future high-resolution']
GRIDS = ['T']


gridT = now.io.open_nemo_griddata_from_zarr(cfg, grid='T')

# Loop on simulations and grids
for sim in SIMULATIONS:
    for grid in GRIDS:
        print("Saving grid %s of the simulation %s to a zarr folder..." %(grid, sim))    
        # Open grids from netCDF files
        griddata = io.open_nemo_griddata_from_zarr(cfg, grid='T')
        
        io.open_nemo_griddata_from_netcdf(NEMO_PATH, simulation=sim, 
                                                     grid=grid, parallel=True,
                                                     autoclose=True)
        # Save the xarray.Dataset to a zarr folder
        io.save_nemo_griddata_to_zarr(griddata, NEMO_PATH, simulation=sim, 
                                      grid=grid)