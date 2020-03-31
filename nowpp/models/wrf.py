import xarray as xr
import pandas as pd
from xgcm import Grid

WRF_NEW_DIMS = {'south_north': 'y_c', 'south_north_stag': 'y_g',
                'west_east': 'x_c', 'west_east_stag': 'x_g'}

WRF_NEW_COORDS = {'XLONG_M': 'lon_T', 'XLAT_M': 'lat_T',
                  'XLONG_U': 'lon_U', 'XLAT_U': 'lat_U',
                  'XLONG_V': 'lon_V', 'XLAT_V': 'lat_V'}

WRF_NEW_MESH = {'F': 'f_T'}

WRF_NEW_MASK = {'LANDMASK': 'mask_T'}


def read_mask(file):
    original_mesh = xr.open_dataset(file).squeeze()
    mask = original_mesh[list(WRF_NEW_MASK)].rename(WRF_NEW_MASK)
    coords = original_mesh[list(WRF_NEW_COORDS)].rename(WRF_NEW_COORDS)
    mask = mask.assign_coords(coords).rename_dims(WRF_NEW_DIMS)
    return mask


def read_mesh(file):
    original_mesh = xr.open_dataset(file).squeeze()
    mesh = original_mesh[list(WRF_NEW_MESH)].rename(WRF_NEW_MESH)
    coords = original_mesh[list(WRF_NEW_COORDS)].rename(WRF_NEW_COORDS)
    mesh = mesh.assign_coords(coords).rename_dims(WRF_NEW_DIMS)
    return mesh


def rename_coords_and_dims(ds, grid='T'):
    ds = ds.rename(WRF_NEW_COORDS).rename_dims(WRF_NEW_DIMS)
    return ds


def open_netcdf_dataset(files, grid='T',  **kwargs):
    ds = xr.open_mfdataset(files, concat_dim='Time',
                           decode_coords=False, decode_times=False,
                           decode_cf=False, drop_variables=['XLON', 'XLAT'],
                           **kwargs)
    ds = xr.decode_cf(ds)
    time = pd.to_datetime(ds.Times.load().astype('str'),
                          format='%Y-%m-%d_%H:%M:%S')
    ds = ds.assign_coords(Time=time).rename({'Time': 'time'})
    ds = rename_coords_and_dims(ds)
    return ds


def build_xgrid(ds, periodic=False):
    xgrid = Grid(ds, periodic=periodic,
                 coords={'X': {'center': 'x_c', 'left': 'x_g'},
                         'Y': {'center': 'y_c', 'left': 'y_g'}})
    return xgrid

