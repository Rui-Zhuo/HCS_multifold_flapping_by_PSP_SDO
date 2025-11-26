import requests
from pathlib import Path

save_dir = Path('E:/Research/Data/SDO/HMI/SHARP/AR2796/')
# Iterating for time range and parameters
for date in ['0118', '0119']:
    for hour in ['00', '01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', \
                    '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23']:
        for minute in ['00', '12', '24', '36', '48']:
            for param in ['Bp', 'Bt', 'Br', 'magnetogram']:
                # Getting time record
                time_record = '2021' + date + '_' + hour + minute + '00'
                print('downloading ' + time_record + ' ' + param + ' with requests')
                # Setting download url (submitted request on JSOC first)
                url = 'https://jsoc1.stanford.edu/SUM38/D1937912030/S00000/hmi.sharp_cea_720s.7532.' \
                    + time_record + '_TAI.' + param + '.fits'
                # Downloading
                r = requests.get(url)
                file_name = url.split('/')[-1]
                file_path = save_dir / file_name
                with open(file_path, 'wb') as file:
                    file.write(r.content)
                # Checking download
                print('saved ' + time_record + ' ' + param)