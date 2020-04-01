from .. import io
import glob
import pkg_resources
import pytest

CONFIG_TEST = pkg_resources.resource_filename('nowpp',
                                              'cfg/config_example.yml')


def test_read_config_file():
    cfg = io.read_config_file(CONFIG_TEST)
    assert('general' in cfg)


def test_class_config():
    cfg = io.Config(CONFIG_TEST)
    assert(1990 == cfg.get_start_date().year)
    assert(2009 == cfg.get_end_date().year)
    assert(['nemo', 'wrf'] == cfg.get_models())
    assert(list(cfg.get_model_simulations('nemo')) ==
           ['NOW-CTRL', 'NOW-NoCFB', 'NOW-NoTFB'])
    assert(list(cfg.check_simulations(None, 'wrf')) ==
           ['NOW-CTRL', 'NOW-NoCFB', 'NOW-NoTFB', 'WRF-ONLY'])
    assert(list(cfg.check_simulations(['NOW-CTRL', 'NOW-NoCFB',
                                       'NOW-NoTFB', 'WRF-ONLY'], 'nemo')) ==
           ['NOW-CTRL', 'NOW-NoCFB', 'NOW-NoTFB'])
    assert(cfg.get_model_grids('nemo') == ['U', 'V', 'T'])
    assert(list(cfg.get_model_directories('nemo')) == ['raw', 'tmp', 'res'])


def test_class_cursor():
    cs = io.Cursor(CONFIG_TEST)
    #assert(csfg.check_simulations)
    #print(cfg)


def test_class_database():
    cs = io.Database(CONFIG_TEST)

#def test_class_cursor():
#    cfg = io.Config(CONFIG_TEST)
#    print(cfg)

#def test_get_nemo_filenames_from_config():
#    nemo_filenames = io.get_nemo_filenames_from_config(CONFIG_FILE, "Present")
#    assert(glob.glob(nemo_filenames))
#    assert(io.get_nemo_filenames_from_config(CONFIG_FILE, "Wrong") is None)
   

#def test_get_wrf_filenames_from_config():
#    wrf_filenames = io.get_wrf_filenames_from_config(CONFIG_FILE, "Present")
#    assert(glob.glob(wrf_filenames))
#    assert(io.get_wrf_filenames_from_config(CONFIG_FILE, "Wrong") is None)

    
#def test_open_nemo_griddata_from_netcdf():
#    io.open_nemo_griddata_from_netcdf(CONFIG_FILE, simulations="Present",
#                                      parallel=True, lock=True)
#    io.open_nemo_griddata_from_netcdf(CONFIG_FILE, parallel=True, lock=True)

    
#def test_open_nemo_griddata_from_zarr():
#    io.open_nemo_griddata_from_zarr(CONFIG_FILE, grid='U')

    
#def test_open_wrfout_from_netcdf():
#    pressure_levels = [700, 750, 800, 850]
#    io._open_wrfout_from_netcdf(CONFIG_FILE, pressure_levels, "Present",
#                                parallel=True, lock=True)
    

#def test_open_wrfout_griddata_from_netcdf():
#    pressure_levels = [700, 750, 800, 850]
#    io.open_wrfout_griddata_from_netcdf(CONFIG_FILE, pressure_levels,
#                                        simulations=["Present", "Future"],
#                                        parallel=True, lock=True)

    