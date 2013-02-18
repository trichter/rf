import rf

# First get data by running this command in the console
# obspyDMT --identity 'TA.Z30A.*.BHZ' --min_mag '5.8' --min_date '2011-01-01' --max_date '2011-01-10' --event_catalog 'IRIS' --arc 'N'

# Set data path of obspyDMT here or in conf.py
rf.set_paths('~/obspyDMT-data/2011-01-01_2011-01-10_5.8_9.9')

# Convert pickled events to catalog file events.xml
rf.convert_dmteventfile()

# Select events for RF
rf.create_rfeventsfile(filters=[])

# Calculate RFs. They will be in the RF directory of the corresponding event
rf.rf('dmt', downsample=10)
