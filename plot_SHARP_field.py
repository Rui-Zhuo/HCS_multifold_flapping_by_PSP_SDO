import numpy as np
import matplotlib.pyplot as plt
import sunpy.map
import os
import re

NOAA_num = 2796
SHARP_num = 7532
save_or_not = 0

lon_ft = 101.42851491767911
lat_ft = -27.155166512345676

Br_dir = f'E:/Research/Data/SDO/HMI/SHARP/AR{NOAA_num}/'
save_dir = 'E:/Research/Work/HCS_multifold_flapping_by_PSP_SDO/20210117/SHARP/'

# (hmi.sharp_cea_720s.SHARP_num.yyyymmdd_hhMMSS_TAI.xxx.fits)
pattern = re.compile(
    rf'hmi\.sharp_cea_720s\.{SHARP_num}\.(\d{{8}})_(\d{{6}})_TAI\.[A-Za-z]+\.fits'
)

Br_files = [f for f in os.listdir(Br_dir) if f.endswith('.fits')]

for Br_fn in Br_files:
    match = pattern.match(Br_fn)
    if not match:
        continue
    
    REC_date = match.group(1)
    REC_time = match.group(2)
    
    if REC_date != '20210117':
        continue
    
    mag_fn = f'hmi.sharp_cea_720s.{SHARP_num}.{REC_date}_{REC_time}_TAI.magnetogram.fits'
    
    Br_map = sunpy.map.Map(os.path.join(Br_dir, Br_fn))
    
    mag_map = sunpy.map.Map(os.path.join(Br_dir, mag_fn))
    
    Br = Br_map.data
    mag = mag_map.data
    
    fig, ax = plt.subplots(figsize=(8, 3))
    
    im1 = ax.pcolormesh(mag, cmap='gray')
    cbar1 = fig.colorbar(im1, ax=ax, label='magnetogram (G)')
    
    # ax.plot(rect_x, rect_y, 'k-', linewidth=2) 
    
    # levels = [-600, -400, -200, 200, 400, 600]
    # im2 = ax.contour(Br, levels=levels, cmap='bwr')
    # cbar2 = fig.colorbar(im2, ax=ax, label='Br (G)', extend='both')
    
    ax.axis('equal')
    plt.title(f'{REC_date}_{REC_time}', fontsize=14)
    plt.tight_layout()
    
    save_filename = f'{REC_date}_{REC_time}_Br_mag.png'
    save_path = os.path.join(save_dir, save_filename)
    
    if save_or_not:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close()
    else:
        plt.show()
        db