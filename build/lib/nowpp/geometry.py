import numpy as np
import pandas as pd
import xarray as xr


EARTH_RADIUS = 6371 * 1000


def latlon2yx(lat, lon):
    """
    Convert latitude and longitude arrays to y and x arrays in m
    Parameters
    ----------
    lat : array-like
        Latitudinal spherical coordinates
    lon : array-like
        Longitudinal spherical coordinates
    Returns
    -------
    y : array-like
        Zonal cartesian coordinates
    x : array-like
        Meridional cartesian coordinates
    """
    y = np.pi / 180. * EARTH_RADIUS * lat
    x = np.cos(np.pi / 180. * lat) * np.pi / 180. * EARTH_RADIUS * lon
    return y, x


def latlon2dydx(lat, lon, dim, label='upper'):
    """
    Convert latitude and longitude arrays to elementary displacements in dy
    and dx
    Parameters
    ----------
    lat : array-like
        Latitudinal spherical coordinates
    lon : array-like
        Longitudinal spherical coordinates
    dim : str
        Dimension along which the differentiation is performed, generally
        associated with the time dimension.
    label : {'upper', 'lower'}, optional
        The new coordinate in dimension dim will have the values of
        either the minuend’s or subtrahend’s coordinate for values ‘upper’
        and ‘lower’, respectively.
    Returns
    -------
    dy : array-like
        Zonal elementary displacement in cartesian coordinates
    dx : array-like
        Meridional elementary displacement in cartesian coordinates
    """
    dlat = lat.diff(dim, label=label)
    dlon = lon.diff(dim, label=label)
    dy = np.pi / 180. * EARTH_RADIUS * dlat
    dx = np.cos(np.pi / 180. * lat) * np.pi / 180. * EARTH_RADIUS * dlon
    return dy, dx


def latlon2heading(lat, lon, dim, label='upper'):
    """
    Calculates the bearing between two points.
    The formulae used is the following:
        θ = atan2(sin(Δlong).cos(lat2),
                  cos(lat1).sin(lat2) − sin(lat1).cos(lat2).cos(Δlong))
    Parameters
    ----------
    Returns
    -------
      The bearing in degrees
    """
    dy, dx = latlon2dydx(lat, lon, dim, label=label)
    initial_heading = np.arctan2(dx, dy) * 180. / np.pi
    # Normalize the initial heading
    compass_heading = (initial_heading + 360) % 360
    return compass_heading
