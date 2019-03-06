from .. import io
import glob
import pkg_resources
import pytest

CONFIG_FILE = pkg_resources.resource_filename('now',
                                              'cfg/config_test.ini')
 
def test_read_config_file():
    cfg = io.read_config_file(CONFIG_FILE)
    assert('GENERAL' in cfg)
 

def test_get_nemo_filenames_from_config():
    nemo_filenames = io.get_nemo_filenames_from_config(CONFIG_FILE, "Present")
    assert(glob.glob(nemo_filenames))
    assert(io.get_nemo_filenames_from_config(CONFIG_FILE, "Wrong") is None)
   

def test_get_wrf_filenames_from_config():
    wrf_filenames = io.get_wrf_filenames_from_config(CONFIG_FILE, "Present")
    assert(glob.glob(wrf_filenames))
    assert(io.get_wrf_filenames_from_config(CONFIG_FILE, "Wrong") is None)

    
def test_open_nemo_griddata_from_netcdf():
    io.open_nemo_griddata_from_netcdf(CONFIG_FILE, simulations="Present", 
                                      parallel=True, lock=True)
    io.open_nemo_griddata_from_netcdf(CONFIG_FILE, parallel=True, lock=True)

    
def test_open_nemo_griddata_from_zarr(): 
    io.open_nemo_griddata_from_zarr(CONFIG_FILE, grid='U')

    
def test_open_wrfout_from_netcdf():
    pressure_levels = [700, 750, 800, 850]
    io._open_wrfout_from_netcdf(CONFIG_FILE, pressure_levels, "Present",
                                parallel=True, lock=True)
    

def test_open_wrfout_griddata_from_netcdf(): 
    pressure_levels = [700, 750, 800, 850]
    io.open_wrfout_griddata_from_netcdf(CONFIG_FILE, pressure_levels,
                                        simulations=["Present", "Future"],
                                        parallel=True, lock=True)

    