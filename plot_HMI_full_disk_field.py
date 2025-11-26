import numpy as np
import os
import re
import astropy.units as u
from astropy.coordinates import SkyCoord
import sunpy.map
from sunpy.coordinates import frames
from sunpy.coordinates.ephemeris import get_earth
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy.ndimage import median_filter

NOAA_num = 2796
save_or_not = 1

lon_ft = 101.4
lat_ft = -27.2
carr_lon = lon_ft * u.deg
carr_lat = lat_ft * u.deg 

Br_dir = 'E:/Research/Data/SDO/HMI/FullDisk/Magnetogram/4096/20210116-20210117/'
save_dir = 'E:/Research/Work/HCS_multifold_flapping_by_PSP_SDO/20210117/magnetogram/'

# (hmi.m_720s.yyyymmdd_HHMMSS_TAI.3.magnetogram.fits)
pattern = re.compile(
    rf'hmi\.m_720s\.(\d{{8}})_(\d{{6}})_TAI\.3\.magnetogram+\.fits'
)

Br_files = [f for f in os.listdir(Br_dir) if f.endswith('.fits')]

for Br_fn in Br_files:
    match = pattern.match(Br_fn)
    if not match:
        continue
    
    REC_date = match.group(1)
    REC_time = match.group(2)
    
    if REC_time[2:] != '0000':
        continue
    
    Br_map = sunpy.map.Map(os.path.join(Br_dir, Br_fn))
    
    Br_data = Br_map.data
    # Br_data_smoothed = median_filter(Br_data, size=17)
    # Br_data_smoothed = median_filter(Br_data, size=3)
    
    obs_time = Br_map.date
    earth_coord = get_earth(obs_time)
    
    carr_coord = SkyCoord(
        lon=carr_lon,
        lat=carr_lat,
        frame=frames.HeliographicCarrington,
        obstime=obs_time,
        observer=earth_coord
    )

    hpc_coord = carr_coord.transform_to(Br_map.coordinate_frame)
    x_pix, y_pix = Br_map.world_to_pixel(hpc_coord)
    
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection=Br_map)
    
    norm = colors.SymLogNorm(linthresh=5, vmin=-500, vmax=500, base=10)
    # im = plt.imshow(Br_data_smoothed, cmap='bwr', norm=norm)
    im = plt.imshow(Br_data, cmap='bwr', norm=norm)
    cbar = plt.colorbar(im, ax=ax, shrink=0.9, pad=0.05, label='Br (G)', extend='both')
    
    # contour = ax.contour(Br_data_smoothed, levels=[0], colors='k', linewidths=1, alpha=0.8)
    
    ax.scatter(x_pix, y_pix, color='lime', edgecolors='k', s=100, linewidths=2, 
               label=f'Footpoint in Carr.Coord. \n ({carr_lon.value} deg., {carr_lat.value} deg.)')
    ax.legend(loc='upper right', fontsize=10)
    
    arcsec_per_pixel = Br_map.scale.axis1 
    fov_arcsec = 200 * u.arcsec
    fov_pixels = (fov_arcsec / arcsec_per_pixel).value
    half_fov = fov_pixels / 2 
    
    ax.set_xlim(x_pix.value - half_fov, x_pix.value + half_fov)
    ax.set_ylim(y_pix.value - half_fov, y_pix.value + half_fov)
    
    ax.set_aspect('equal')
    ax.set_title(REC_date + '_' + REC_time + '_magnetogram')
    ax.invert_xaxis()
    ax.invert_yaxis()
    plt.tight_layout()
    
    save_filename = f'{REC_date}_{REC_time}_Br.png'
    save_path = os.path.join(save_dir, save_filename)
    
    if save_or_not:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close()
    else:
        plt.show()
        db