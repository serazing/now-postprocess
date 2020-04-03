import xarray as xr
from xgcm import Grid

NEMO_NEW_DIMS = {'T': {'x': 'x_c', 'y': 'y_c', 'z': 'z_c'},
                 'U': {'x': 'x_r', 'y': 'y_c', 'z': 'z_c'},
                 'V': {'x': 'x_c', 'y': 'y_r', 'z': 'z_c'},
                 'F': {'x': 'x_r', 'y': 'y_r', 'z': 'z_c'},
                 'W': {'x': 'x_c', 'y': 'y_c', 'z': 'z_r'}
                 }

NEMO_NEW_DIMS_2D = {'T': {'x': 'x_c', 'y': 'y_c'},
                    'U': {'x': 'x_r', 'y': 'y_c'},
                    'V': {'x': 'x_c', 'y': 'y_r'},
                    'F': {'x': 'x_r', 'y': 'y_r'},
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

NEMO_NEW_COORDS_2D = {'T': {'glamt': 'lon_T',
                            'gphit': 'lat_T'},
                      'U': {'glamu': 'lon_U',
                            'gphiu': 'lat_U'},
                      'V': {'glamv': 'lon_V',
                            'gphiv': 'lat_V'},
                      'F': {'glamf': 'lon_F',
                            'gphif': 'lat_F'},
                      }

NEMO_NEW_MASK = {'T': {'tmask': 'mask_T'},
                 'U': {'umask': 'mask_U'},
                 'V': {'vmask': 'mask_V'},
                 'F': {'fmask': 'mask_F'}}


NEMO_NEW_MESH = {'T': {'e1t': 'dx_T',
                       'e2t': 'dy_T',
                       'e3t': 'dz_T'},
                 'U': {'e1u': 'dx_U',
                       'e2u': 'dy_U',
                       'e3u': 'dz_U'},
                 'V': {'e1v': 'dx_V',
                       'e2v': 'dy_V',
                       'e3v': 'dz_V'},
                 'F': {'e1f': 'dx_F',
                       'e2f': 'dy_F',
                       'ff': 'f_F'},
                 'W': {'e3w': 'dz_W'}
                 }


def read_mask(file, grids=None):
    if grids is None:
        grids = ['U', 'V', 'T', 'F']
    original_mask = xr.open_dataset(file, chunks={}).squeeze()
    list_of_mask = []
    for grid in grids:
        var_dict = NEMO_NEW_MASK[grid]
        coord_dict = NEMO_NEW_COORDS[grid]
        sim_dict = NEMO_NEW_DIMS[grid]
        mask = original_mask[list(var_dict)].rename(var_dict)
        coords = original_mask[list(coord_dict)].rename(coord_dict)
        mask = mask.assign_coords(coords).rename_dims(sim_dict)
        list_of_mask.append(mask)
    return xr.merge(list_of_mask)


def read_mesh(file, grids=None):
    if grids is None:
        grids = ['U', 'V', 'T', 'F', 'W']
    original_mesh = xr.open_dataset(file, chunks={}).squeeze()
    list_of_mesh = []
    for grid in grids:
        mesh_dict = NEMO_NEW_MESH[grid]
        coord_dict = NEMO_NEW_COORDS[grid]
        sim_dict = NEMO_NEW_DIMS[grid]
        mesh = original_mesh[list(mesh_dict)].rename(mesh_dict)
        coords = original_mesh[list(coord_dict)].rename(coord_dict)
        mesh = mesh.assign_coords(coords).rename_dims(sim_dict)
        list_of_mesh.append(mesh)
    return xr.merge(list_of_mesh)


def rename_dims(ds, grid='T'):
    if 'z' in ds.dims:
        dims_dict = NEMO_NEW_DIMS[grid]
    else:
        dims_dict = NEMO_NEW_DIMS_2D[grid]
    ds = ds.rename_dims(dims_dict)
    try:
        ds = ds.set_index(time_counter='time_average_1d')
    except ValueError:
        pass
    try:
        ds = ds.rename({'time_counter': 'time'})
    except ValueError:
        pass
    return ds


def assign_coords(ds, mesh, grid='T'):
    if grid == 'W':
        ds = ds.drop(('nav_lon', 'nav_lat', 'nav_lev'))
        ds = ds.assign_coords(mesh.variables)
    elif grid == ['U', 'V', 'T', 'F']:
        if 'z_c' in ds.dims:
            ds = ds.drop(('nav_lon', 'nav_lat', 'nav_lev'))
            ds = ds.assign_coords(mesh.variables)
        else:
            ds = ds.drop(('nav_lon', 'nav_lat'))
            ds = ds.assign_coords(mesh.drop('z_c').variables)
    else:
        raise ValueError
    return ds


def open_netcdf_dataset(files, mesh_file=None, grid='T',  **kwargs):
    ds = xr.open_mfdataset(files, **kwargs)
    ds = rename_dims(ds, grid=grid)
    if mesh_file is not None:
        mesh = read_mesh(mesh_file, grids=[grid])
        ds = assign_coords(ds, mesh, grid=grid)
        mask = ds['mask_%s'] % grid
        ds = ds.where(mask == 1)
    return ds


def build_xgrid(ds, periodic=False):
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


def update_grid(ds):
    xgrid = build_xgrid(ds)
    metrics = {'dxdy_T': xgrid.get_metric(ds.lon_T, ('X', 'Y')),
               'dxdy_U': xgrid.get_metric(ds.lon_U, ('X', 'Y')),
               'dxdy_V': xgrid.get_metric(ds.lon_V, ('X', 'Y')),
               'dxdydz_T': xgrid.get_metric(ds.depth_T, ('X', 'Y', 'Z')),
               'dxdydz_U': xgrid.get_metric(ds.depth_U, ('X', 'Y', 'Z')),
               'dxdydz_V': xgrid.get_metric(ds.depth_V, ('X', 'Y',  'Z'))
               }
    ds = ds.assign(metrics)
    return ds
