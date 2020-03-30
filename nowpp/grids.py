from . import io
from xgcm import Grid
import xarray as xr


# NEMO section
###############################################################################
NEMO_NEW_DIMS = {'T': {'x': 'x_c', 'y': 'y_c', 'z': 'z_c'},
                 'U': {'x': 'x_r', 'y': 'y_c', 'z': 'z_c'},
                 'V': {'x': 'x_c', 'y': 'y_r', 'z': 'z_c'},
                 'F': {'x': 'x_r', 'y': 'y_r', 'z': 'z_c'},
                 'W': {'x': 'x_c', 'y': 'y_c', 'z': 'z_r'}
                 }

NEMO_NEW_COORDS = {'T': {'glamt': 'lon_T',
                         'gphit': 'lat_T',
                         'gdept': 'depth_T'},
                   'U': {'glamu': 'lon_U',
                         'gphiu': 'lat_U',
                         'gdepu': 'depth_U'},
                   'V': {'glamv': 'lon_V',
                         'gphiv': 'lat_V',
                         'gdepv': 'depth_V'},
                   'F': {'glamf': 'lon_F',
                         'gphif': 'lat_F'},
                   'W': {'gdept': 'depth_W'}
                   }

NEMO_NEW_GRID = {'T': {'tmask': 'mask_T',
                       'e1t': 'dx_T',
                       'e2t': 'dy_T',
                       'e3t': 'dz_T'},
                 'U': {'umask': 'mask_U',
                       'e1u': 'dx_U',
                       'e2u': 'dy_U',
                       'e3u': 'dz_U'},
                 'V': {'vmask': 'mask_V',
                       'e1v': 'dx_V',
                       'e2v': 'dy_V',
                       'e3v': 'dz_V'},
                 'F': {'fmask': 'mask_F',
                       'e1f': 'dx_F',
                       'e2f': 'dy_F', 
                       'ff': 'f_F'},
                 'W': {'e3w': 'dz_W'}
                 }


def read_nemo_grid(config_file):
    cfg = io.read_config_file(config_file)
    meshgrid_file = cfg['NEMO']['MeshGridFile']
    list_of_meshgrid = []
    original_meshgrid = xr.open_dataset(meshgrid_file).squeeze()
    for grid in ['U', 'V', 'T', 'F', 'W']:
        meshgrid = original_meshgrid[list(NEMO_NEW_GRID[grid])].rename(NEMO_NEW_GRID[grid])
        coords = original_meshgrid[list(NEMO_NEW_COORDS[grid])].rename(NEMO_NEW_COORDS[grid])
        meshgrid = meshgrid.assign_coords(coords).rename_dims(NEMO_NEW_DIMS[grid])
        list_of_meshgrid.append(meshgrid)
    ds_meshgrid = xr.merge(list_of_meshgrid)
    ds_meshgrid = update_nemo_grid(ds_meshgrid)
    return ds_meshgrid


def build_nemo_xgrid(ds, periodic=False):
    xgrid = Grid(ds, periodic=periodic,
                 coords={'X': {'center': 'x_c', 'right': 'x_r'},
                         'Y': {'center': 'y_c', 'right': 'y_r'},
                         'Z': {'center': 'z_c', 'right': 'z_r'},
                         },
                 metrics={('X',): ['dx_T', 'dx_U', 'dx_V', 'dx_F'],
                          ('Y',): ['dy_T', 'dy_U', 'dy_V', 'dy_F'],
                          ('Z',): ['dz_T', 'dz_U', 'dz_V'],
                          }
                 )
    return xgrid


def update_nemo_grid(ds):
    xgrid = build_nemo_xgrid(ds)
    metrics = {'dxdy_T': xgrid.get_metric(ds.lon_T, ('X', 'Y')),
               'dxdy_U': xgrid.get_metric(ds.lon_U, ('X', 'Y')),
               'dxdy_V': xgrid.get_metric(ds.lon_V, ('X', 'Y')),
               'dxdydz_T': xgrid.get_metric(ds.depth_T, ('X', 'Y', 'Z')),
               'dxdydz_U': xgrid.get_metric(ds.depth_U, ('X', 'Y', 'Z')),
               'dxdydz_V': xgrid.get_metric(ds.depth_V, ('X', 'Y',  'Z'))
               }
    ds = ds.assign(metrics)
    return ds


# WRF section
###############################################################################

WRF_NEW_DIMS = {'south_north': 'y_c', 'south_north_stag': 'y_g',
                'west_east': 'x_c', 'west_east_stag': 'x_g'}

WRF_NEW_COORDS = {'XLONG_M': 'lon_T', 'XLAT_M': 'lat_T',
                  'XLONG_U': 'lon_U', 'XLAT_U': 'lat_U',
                  'XLONG_V': 'lon_V', 'XLAT_V': 'lat_V'}

WRF_NEW_GRID = {'LANDMASK': 'mask_T', 'F': 'f_T'}


def read_wrf_grid(config_file):
    cfg = io.read_config_file(config_file)
    meshgrid_file = cfg['WRF']['MeshGridFile']
    meshgrid = xr.open_dataset(meshgrid_file).squeeze()
    return meshgrid


def build_wrf_xgrid(ds, periodic=False):
    xgrid = Grid(ds, periodic=periodic,
                 coords={'X': {'center': 'x_c', 'left': 'x_g'},
                         'Y': {'center': 'y_c', 'left': 'y_g'}})
    return xgrid

