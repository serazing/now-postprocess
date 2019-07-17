from oocgcm.core.grids import VectorField2d
from . import io


def wind_curl(config_file, taux, tauy):
	"""
	Compute the wind curl on NEMO C-grid evaluated at the 'F' point.

	Parameters
	----------
	config_file : str
		The configuration file including the path of the grid file
	taux : xr.DataArray
		The zonal wind stress evaluated at the 'U' point
	tauy : xr.DataArray
		The meridional wind stress evaluated at the 'V' point

	Returns
	-------
	wind_curl : xr.DataArray
		The wind curl evaluated at the 'F' point
	"""
	meshgrid = io.get_nemo_grid_from_config(config_file)
	taux.attrs['grid_location'] = 'u'
	tauy.attrs['grid_location'] = 'v'
	tau = VectorField2d(taux, tauy,
						x_component_grid_location='u', 
						y_component_grid_location='v')
	wind_curl = meshgrid.vertical_component_of_curl(tau)
	return wind_curl.squeeze()


def geostrophic_velocities(config_file, ssh):
	"""
	Compute the geostrophic velocities from the sea level field on NEMO C-grid.

	Parameters
	----------
	config_file : str
		The configuration file including the path of the grid file
	ssh : xr.DataArray
		The sea surface height evaluated at the 'T' grid

	Returns
	-------
	ux_geo : xr.DataArray
		The zonal geostrophic velocity evaluated at the 'U' point
	vy_geo : xr.DataArray
		The meridional geostrophic velocity evaluated at the 'V' point
	"""
	meshgrid = io.get_nemo_grid_from_config(config_file)
	ssh.attrs['grid_location'] = 't'
	ug = meshgrid.geostrophic_current_from_sea_surface_height(ssh)
	ux_geo, vy_geo = ug.x_component.squeeze(), ug.y_component.squeeze() 
	return ux_geo, vy_geo


def kinetic_energy(config_file, ux, vy):
	"""
	Compute the kinetic energy from the zonal and meridional velocities
	on NEMO C-grid.

	Parameters
	----------
	config_file : str
		The configuration file including the path of the grid file
	ux : xr.DataArray
		The zonal velocity evaluated at the 'U' point
	vy : xr.DataArray
		The meridional velocity evaluated at the 'V' point

	Returns
	-------
	ke : xr.DataArray
		The kinetic energy evaluated at the 'T' point
	"""
	meshgrid = io.get_nemo_grid_from_config(config_file)
	ux.attrs['grid_location'] = 'u'
	vy.attrs['grid_location'] = 'v'
	u = VectorField2d(ux, vy, 
					  x_component_grid_location='u', 
					  y_component_grid_location='v')
	ke = 0.5 * meshgrid.scalar_product(u, u)
	return ke.squeeze()


def wind_energy_transfer(config_file, ux, vy, taux, tauy):
	"""
	Compute the wind work energy from the zonal and meridional velocities
	on NEMO C-grid.

	Parameters
	----------
	config_file : str
		The configuration file including the path of the grid file
	ux : xr.DataArray
		The zonal velocity evaluated at the 'U' point
	vy : xr.DataArray
		The meridional velocity evaluated at the 'V' point
	taux : xr.DataArray
		The zonal wind stress evaluated at the 'U' point
	tauy : xr.DataArray
		The meridional wind stress evaluated at the 'V' point

	Returns
	-------
	wind_energy : xr.DataArray
		The transfer of energy from the wind evaluated at the 'T' point
	"""
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
	wind_energy = meshgrid.scalar_product(u, tau)
	wind_energy.attrs['long_name'] = 'Wind work'
	wind_energy.attrs['units'] = 'm3/s3'
	return wind_energy.squeeze()


def tracer_flux(config_file, ux, vy, tracer):
	"""
	Compute a flux of tracer using the zonal and meridional velocities
	on NEMO C-grid.

	Parameters
	----------
	config_file : str
		The configuration file including the path of the grid file
	ux : xr.DataArray
		The zonal velocity evaluated at the 'U' point
	vy : xr.DataArray
		The meridional velocity evaluated at the 'V' point
	tracer : xr.DataArray
		The tracer evaluated at the 'T' point
	
	Returns
	-------		
	ugradT : xr.DataArray
		The tracer flux evaluated at the 'T' point
	"""
	meshgrid = io.get_nemo_grid_from_config(config_file)
	ux.attrs['grid_location'] = 'u'
	vy.attrs['grid_location'] = 'v'
	tracer.attrs['grid_location'] = 't'
	u = VectorField2d(ux, vy,
					  x_component_grid_location='u',
					  y_component_grid_location='v')
	gradT = meshgrid.horizontal_gradient(tracer)
	ugradT = meshgrid.scalar_product(u, gradT)
	return ugradT.squeeze()
