import numpy as np
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

def dlondlat_dxdy(dlon, dlat, lon, lat):
    EARTH_RADIUS = 6371 * 1000
    dx = np.cos(np.pi / 180. * lat) * np.pi / 180. * EARTH_RADIUS * dlon
    dy = (lon * 0 + 1) * np.pi / 180. * EARTH_RADIUS * dlat
    return dx, dy

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
    mesh = mesh.isel(x_g=slice(None, -1), y_g=slice(None, -1))
    xgrid = build_xgrid(mesh)
    lon_F = xgrid.interp(mesh.lon_U, 'Y', boundary='fill', fill_value=np.nan)
    lat_F = xgrid.interp(mesh.lat_V, 'X', boundary='fill', fill_value=np.nan)
    mesh = mesh.assign_coords(lon_F=lon_F, lat_F=lat_F)
    dlon_T = xgrid.diff(mesh.lon_U, 'X', boundary='fill', fill_value=np.nan)
    dlon_U = xgrid.diff(mesh.lon_T, 'X', boundary='fill', fill_value=np.nan)
    dlon_V = xgrid.diff(mesh.lon_F, 'X', boundary='fill', fill_value=np.nan)
    dlon_F = xgrid.diff(mesh.lon_V, 'X', boundary='fill', fill_value=np.nan)
    dlat_T = xgrid.diff(mesh.lat_V, 'Y', boundary='fill', fill_value=np.nan)
    dlat_U = xgrid.diff(mesh.lat_F, 'Y', boundary='fill', fill_value=np.nan)
    dlat_V = xgrid.diff(mesh.lat_T, 'Y', boundary='fill', fill_value=np.nan)
    dlat_F = xgrid.diff(mesh.lat_U, 'Y', boundary='fill', fill_value=np.nan)
    dx_T, dy_T = dlondlat_dxdy(dlon_T, dlat_T, mesh.lon_T, mesh.lat_T)
    dx_U, dy_U = dlondlat_dxdy(dlon_U, dlat_U, mesh.lon_U, mesh.lat_U)
    dx_V, dy_V = dlondlat_dxdy(dlon_V, dlat_V, mesh.lon_V, mesh.lat_V)
    dx_F, dy_F = dlondlat_dxdy(dlon_F, dlat_F, mesh.lon_F, mesh.lat_F)
    mesh = mesh.assign_coords(dx_T=dx_T, dy_T=dy_T, dx_U=dx_U, dy_U=dy_U,
                              dx_V=dx_V, dy_V=dy_V, dx_F=dx_F, dy_F=dy_F)
    return mesh


def rename_dims(ds):
    ds = ds.rename_dims(WRF_NEW_DIMS)
    return ds


def rename_coords_and_dims(ds, grid='T'):
    ds = ds.rename(WRF_NEW_COORDS).rename(WRF_NEW_DIMS)
    return ds


def open_netcdf_dataset(files, mesh_file=None, grid='T', variables=None,
                        **kwargs):
    ds = xr.open_mfdataset(files, concat_dim='Time', combine='nested',
                           drop_variables=['XLON', 'XLAT'], **kwargs)
    ds = ds.rename_dims(WRF_NEW_DIMS)
    ds = ds.isel(x_g=slice(None, -1), y_g=slice(None, -1))
    time = pd.to_datetime(ds.Times.load().astype('str').data,
                          format='%Y-%m-%d_%H:%M:%S')
    ds = ds.assign_coords(Time=time).rename({'Time': 'time'})
    mesh = read_mesh(mesh_file)
    ds = ds.assign_coords(mesh.variables)
    if variables is not None:
        ds = ds[variables]
    return ds


def build_xgrid(ds, periodic=False, metric=None):
    xgrid = Grid(ds, periodic=periodic,
                 coords={'X': {'center': 'x_c', 'left': 'x_g'},
                         'Y': {'center': 'y_c', 'left': 'y_g'}})
    return xgrid