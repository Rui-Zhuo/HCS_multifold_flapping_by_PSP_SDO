import os

import astropy.units as u
from astropy.time import Time
import sunpy.map
from sunpy.net import Fido
from sunpy.net import attrs as a

save_dir = 'E:/Research/Data/SDO/HMI/FullDisk/Magnetogram/4096/20210116-20210117/'
jsoc_email = 'ruizhuo@pku.edu.cn'

# Iterating for time range and parameters
year = '2021'
month = '01'
for day in ['16', '17']:
    for hour in ['00', '01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', \
                    '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23']:
        for minute in ['00', '12', '24', '36', '48']:
            name_record = year + month + day + '_' + hour + minute + '00'
            print('Downloading ' + name_record + ' with Fido')        

            sear_time = Time(f"{year}-{month}-{day}T{hour}:{minute}:00", 
                            scale='utc', format='isot')

            query = Fido.search(
                a.Time(sear_time - 1*u.s, sear_time + 1*u.s),
                a.Sample(360*u.s),
                a.jsoc.Series.hmi_m_720s,
                a.jsoc.Notify(jsoc_email)
            )
            print(query)

            # Checking download    
            if len(query) > 0 and len(query[0]) > 0:
                file_info = query[0][0]
                T_REC = file_info['T_REC']
                
                save_fn = f'hmi.m_720s.{T_REC[0:4]}{T_REC[5:7]}{T_REC[8:10]}_{T_REC[11:13]}{T_REC[14:16]}{T_REC[17:19]}_TAI.1.magnetogram.fits'
                save_path = os.path.join(save_dir, save_fn)
                
                if os.path.exists(save_path):
                    print(f'File {save_fn} already exists, skipping download')
                    continue
            
                files = Fido.fetch(query[0][0], path=save_dir, overwrite=False)
                print(f'Saved to {files}')
            else:
                print(f'No data found')
            