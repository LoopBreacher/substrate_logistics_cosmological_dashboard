"""
HEALPix compatibility layer for Windows / astropy-healpix.
Provides drop-in replacements for standard healpy functions used in cosmological mapping:
- nside2npix
- npix2nside
- pix2vec
- pix2ang
- ang2pix
- vec2pix
- mollview
- graticule
- projscatter
"""

import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy_healpix import HEALPix

def nside2npix(nside):
    return 12 * nside * nside

def npix2nside(npix):
    return int(np.sqrt(npix // 12))

def pix2ang(nside, ipix, nest=False):
    hp = HEALPix(nside=nside, order='nested' if nest else 'ring')
    lon, lat = hp.healpix_to_lonlat(ipix)
    phi = lon.to_value(u.rad)
    theta = (90.0 * u.deg - lat).to_value(u.rad)
    return theta, phi

def pix2vec(nside, ipix, nest=False):
    theta, phi = pix2ang(nside, ipix, nest=nest)
    x = np.sin(theta) * np.cos(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(theta)
    return np.array([x, y, z])

def ang2pix(nside, theta, phi, nest=False):
    hp = HEALPix(nside=nside, order='nested' if nest else 'ring')
    lat = 90.0 * u.deg - (theta * u.rad)
    lon = phi * u.rad
    return hp.lonlat_to_healpix(lon, lat)

def vec2pix(nside, x, y, z, nest=False):
    r = np.sqrt(x**2 + y**2 + z**2)
    theta = np.arccos(np.clip(z / r, -1.0, 1.0))
    phi = np.arctan2(y, x) % (2.0 * np.pi)
    return ang2pix(nside, theta, phi, nest=nest)

def mollview(map_data, fig=None, sub=None, title="", unit="", cmap="plasma", 
             min=None, max=None, coord=None, cbar=True, flip="astro", **kwargs):
    """
    Renders a HEALPix map in Mollweide projection with default flip='astro' orientation.
    """
    if sub is not None:
        if isinstance(sub, tuple):
            nrows, ncols, index = sub
        elif isinstance(sub, int):
            nrows = sub // 100
            ncols = (sub % 100) // 10
            index = sub % 10
        else:
            nrows, ncols, index = 1, 1, 1
    else:
        nrows, ncols, index = 1, 1, 1

    if fig is not None:
        fig_obj = plt.figure(fig) if isinstance(fig, int) else fig
        ax = fig_obj.add_subplot(nrows, ncols, index, projection="mollweide")
    else:
        ax = plt.subplot(nrows, ncols, index, projection="mollweide")

    npix = len(map_data)
    nside = npix2nside(npix)
    
    # Create lon/lat grid for Mollweide interpolation
    n_lon = 240
    n_lat = 120
    lon_grid = np.linspace(-np.pi, np.pi, n_lon)
    lat_grid = np.linspace(-np.pi / 2.0, np.pi / 2.0, n_lat)
    Lon, Lat = np.meshgrid(lon_grid, lat_grid)

    theta_grid = np.pi / 2.0 - Lat
    
    # Astronomical convention: East to left (l increases leftward)
    if flip == "geo":
        phi_grid = Lon % (2.0 * np.pi)
    else:
        phi_grid = (-Lon) % (2.0 * np.pi)

    pix_indices = ang2pix(nside, theta_grid, phi_grid)
    grid_vals = map_data[pix_indices]

    vmin = min if min is not None else np.nanmin(map_data)
    vmax = max if max is not None else np.nanmax(map_data)

    mesh = ax.pcolormesh(Lon, Lat, grid_vals, cmap=cmap, vmin=vmin, vmax=vmax, shading='auto')
    if title:
        ax.set_title(title, pad=12, fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3, linestyle=':')
    
    if cbar:
        cb = plt.colorbar(mesh, ax=ax, orientation='horizontal', pad=0.08, shrink=0.7)
        if unit:
            cb.set_label(unit, fontsize=10)
    return ax

def graticule():
    ax = plt.gca()
    ax.grid(True, alpha=0.4, linestyle='--')

def projscatter(lon, lat, lonlat=True, coord=None, **kwargs):
    ax = plt.gca()
    if lonlat:
        lon_rad = np.radians(lon)
        lon_mapped = (lon_rad + np.pi) % (2.0 * np.pi) - np.pi
        x_plot = -lon_mapped
        lat_rad = np.radians(lat)
    else:
        theta, phi = lon, lat
        lat_rad = np.pi / 2.0 - theta
        phi_mapped = (phi + np.pi) % (2.0 * np.pi) - np.pi
        x_plot = -phi_mapped
    return ax.scatter(x_plot, lat_rad, **kwargs)