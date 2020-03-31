from . import io
from xgcm import Grid
import xarray as xr

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


def build_wrf_xgrid(ds, periodic=False):
    xgrid = Grid(ds, periodic=periodic,
                 coords={'X': {'center': 'x_c', 'left': 'x_g'},
                         'Y': {'center': 'y_c', 'left': 'y_g'}})
    return xgrid

