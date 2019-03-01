import xarray as xr
from now import io


# Define the parser for input parameters
import argparse
description = ("Convert NEMO outputs and post-processed WRF outptus"
               "to zarr format"
              ) 
parser = argparse.ArgumentParser(description=description)
parser.add_argument('path', nargs='+', help='Input files')
parser.add_argument('--wrf', action='store_true', default=False,
                    help='Convert or not WRF files')
parser.add_argument('--nemo', action='store_true', default=False,
                    help='Convert or not NEMO files')
parser.add_argument('--grids', nargs='+', default=['U', 'V', 'T'],
                    help='Convert or not NEMO files')
args = parser.parse_args()

cfg = io.read_config_file(args.cfg)
if simulations is None:
        simulations = [key for key in cfg 
                       if key not in ['DEFAULT', 'GENERAL', 'WRF', 'NEMO']]
    elif not isinstance(simulations, (list, tuple)):
        simulations = [simulations,]
# Loop on simulations and grids
for sim in SIMULATIONS:
    for grid in args.grid:
        print("Saving grid %s of the simulation %s to a zarr folder..." %(grid, sim))    
        # Open grids from netCDF files
        griddata = io.open_nemo_griddata_from_netcdf(NEMO_PATH, simulation=sim, 
                                                     grid=grid, parallel=True,
                                                     autoclose=True)
        # Save the xarray.Dataset to a zarr folder
        io.save_nemo_griddata_to_zarr(griddata, NEMO_PATH, simulation=sim, 
                                      grid=grid)