"""
Configuration for rf tests.
"""
events = None

stations = None

data_path = '~/obspyDMT-data/2010-01-01_2010-01-10_5.5_9.9'
data = '{eventid}/BH_RAW/{net}.{sta}.{loc}.{cha}'

output_path = '~/obspyDMT-data/2010-01-01_2010-01-10_5.5_9.9'
rf = '{eventid}/RF/{net}.{sta}.{loc}.{cha}'
mout = '{eventid}/RF_{motype}MOUT/{net}.{sta}.{loc}.{cha}'
mean = 'RF_SUM/RF_SUM_{net}.{sta}.{loc}.{cha}'
format = 'Q'
