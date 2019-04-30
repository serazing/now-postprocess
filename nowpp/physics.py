from oocgcm.core.grids import VectorField2d
from . import io


def compute_wind_curl(config_file, taux, tauy):
	meshgrid = io.get_nemo_grid_from_config(config_file)
	taux.attrs['grid_location'] = 'u'
	tauy.attrs['grid_location'] = 'v'
	tau = VectorField2d(taux, tauy,
						x_component_grid_location='u', 
						y_component_grid_location='v')
	wind_curl = meshgrid.vertical_component_of_curl(tau)
	return wind_curl.squeeze()


def compute_geostrophic_velocities(config_file, sla):
	meshgrid = io.get_nemo_grid_from_config(config_file)
	sla.attrs['grid_location'] = 't'
	ug = meshgrid.geostrophic_current_from_sea_surface_height(sla)
	ux_geo, vy_geo = ug.x_component.squeeze(), ug.y_component.squeeze() 
	return ux_geo, vy_geo


def compute_ke(config_file, ux, vy):
	meshgrid = io.get_nemo_grid_from_config(config_file)
	ux.attrs['grid_location'] = 'u'
	vy.attrs['grid_location'] = 'v'
	u = VectorField2d(ux, vy, 
					  x_component_grid_location='u', 
					  y_component_grid_location='v')
	ke = 0.5 * meshgrid.scalar_product(u, u)
	return ke.squeeze()


def compute_wind_work(config_file, ux, vy, taux, tauy):
	meshgrid = io.get_nemo_grid_from_config(config_file)
	ux.attrs['grid_location'] = 'u'
	vy.attrs['grid_location'] = 'v'
	taux.attrs['grid_location'] = 'u'
	tauy.attrs['grid_location'] = 'v'    
	u = VectorField2d(ux, vy, 
					  x_component_grid_location='u', 
					  y_component_grid_location='v')
	tau = VectorField2d(taux, tauy, 
						x_component_grid_location='u', 
						y_component_grid_location='v')
	wind_work = meshgrid.scalar_product(u, tau)
	wind_work.attrs['long_name'] = 'Wind work'
	wind_work.attrs['units'] = 'm3/s3'
	return wind_work.squeeze()