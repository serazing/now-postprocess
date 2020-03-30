from . import io
from xgcm import Grid
import xarray as xr

NEMO_NEW_DIMS = {'T': {'x': 'x_c', 'y': 'y_c', 'z': 'z_c'},
                 'U': {'x': 'x_r', 'y': 'y_c', 'z': 'z_c'},
                 'V': {'x': 'x_c', 'y': 'y_r', 'z': 'z_c'},
                 'F': {'x': 'x_r', 'y': 'y_r', 'z': 'z_c'},
                 'W': {'x': 'x_c', 'y': 'y_c', 'z': 'z_r'}
                }

NEMO_NEW_COORDS = {'T': {'glamt': 'lon_T', 'gphit': 'lat_T', 'gdept': 'depth_T'},
                   'U': {'glamu': 'lon_U', 'gphiu': 'lat_U', 'gdepu': 'depth_U'},
                   'V': {'glamv': 'lon_V', 'gphiv': 'lat_V', 'gdepv': 'depth_V'},
                   'F': {'glamf': 'lon_F', 'gphif': 'lat_F'},
                   'W': {'gdept': 'depth_W'}
                  }

NEMO_NEW_GRID = {'T': {'tmask': 'mask_T', 'e1t': 'dx_T', 'e2t': 'dy_T', 'e3t': 'dz_T'},
                 'U': {'umask': 'mask_U', 'e1u': 'dx_U', 'e2u': 'dy_U', 'e3u': 'dz_U'},
                 'V': {'vmask': 'mask_V', 'e1v': 'dx_V', 'e2v': 'dy_V', 'e3v': 'dz_V'},
                 'F': {'fmask': 'mask_F', 'e1f': 'dx_F', 'e2f': 'dy_F', 'ff': 'f_F'},
                 'W': {'e3w': 'dz_W'}
                 }


def read_nemo_grid(config_file):
    cfg = io.read_config_file(config_file)
    meshgrid_file = cfg['NEMO']['MeshGridFile']
    dims = cfg['NEMO']['Dimensions']
    list_of_meshgrid = []
    original_meshgrid = xr.open_dataset(meshgrid_file).squeeze()
    for grid in ['U', 'V', 'T', 'F', 'W']:
        meshgrid = original_meshgrid[list(NEMO_NEW_GRID[grid])].rename(NEMO_NEW_GRID[grid])
        coords = original_meshgrid[list(NEMO_NEW_COORDS[grid])].rename(NEMO_NEW_COORDS[grid])
        meshgrid = meshgrid.assign_coords(coords).rename_dims(NEMO_NEW_DIMS[grid])
        list_of_meshgrid.append(meshgrid)
    ds_meshgrid = xr.merge(list_of_meshgrid)
    return ds_meshgrid


def get_nemo_xgrid(ds, periodic=False):
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
    xgrid = get_nemo_xgrid(ds)
    metrics = {'dxdy_T': xgrid.get_metric(ds.lon_T, ('X', 'Y')),
               'dxdy_U': xgrid.get_metric(ds.lon_U, ('X', 'Y')),
               'dxdy_V': xgrid.get_metric(ds.lon_V, ('X', 'Y')),
               'dxdydz_T': xgrid.get_metric(ds.depth_T, ('X', 'Y', 'Z')),
               'dxdydz_U': xgrid.get_metric(ds.depth_U, ('X', 'Y', 'Z')),
               'dxdydz_V': xgrid.get_metric(ds.depth_V, ('X', 'Y',  'Z'))
               }
    ds = ds.assign(metrics)
    return ds
