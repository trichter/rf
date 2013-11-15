### Example configuration file for batch usage of rf package

# Method can be 'client' or 'dmt'.
# Method 'client' circles through all events and stations loaded from the
# events and stations files and requests the necessary data from a function
# defined in this config file.
# Method 'dmt' is designed to be used with obspyDMT. It loads files in
# obspyDMT data directory and loads events/stations corresponding
# to the data
method = 'client'

phase = 'P'

# File names of events and station files .
# This files have to exist for method 'client', but not necessarily for 'dmt'.
# If files do not exist, they will be created.
events = 'test_events.xml'
stations = 'test_stations.xml'

### Options for input and output

# data_path is only necessary for method 'dmt'
data_path = 'test_data/'

# definition of function get_waveform and request_window is only necessary for
# method 'client'
# Note: You can easily write a 'get_waveform' function delivering data from
# your local hard drive or you can use the powerful ObsPy clients.

# Example for usage of ObsPy client:
# import obspy.iris
# client = obspy.iris.Client()
# def getwaveform(station, t1, t2):
#     return client.getWaveform('TA', station, '', 'BH?', t1, t2)

# client delivering from local hard disc
def get_waveform(station, t1, t2):
    pass
request_window = (-50, 150)

# Output paths. The output paths also determine if multiple files are
# written or not. For SAC each trace is written to one file, anyway.
# Supported modes are: Files for each station; for each station and year;
# for each station and event.
# Valid expressions: {method}, {eventid}, {net}, {sta}, {loc}, {cha}, {year},
# {month}, {day}, {event_time}, {magnitude}, {motype} (only for moveout)
rf = 'test_output/{method}RF/{net}.{sta}.{loc}.{cha}'
mout = 'test_output/{method}RF_{motype}MOUT/{net}.{sta}.{loc}.{cha}'
mean = 'test_output/{method}RF_SUM/SUM_{net}.{sta}.{loc}.{cha}'

# format should be Q or SAC
format = 'Q'

### Options for receiver function calculation

# general
filter = {'type': 'bandpass', 'freqmin': 0.01, 'freqmax': 2.}
windowP = (-30, 100)  #trim to window around P-onset (P-receiver function)
windowS = (-80, 50)  #trim to window around S-onset (S-receiver function)
downsample = 10

# rotate can be a function accepting a stream for your own
# rotation function, this function has to return the source
# component. E.g.:
# def rotate(stream):
#    code to rotate stream
#    return 'L' # return name of source component
# But it is easier to use built in rotation 'ZNE->LQT' or 'ZNE->LRT'.
rotate = 'ZNE->LQT'
# deconvolve is 'time' or 'freq' for time or frequency domain deconvolution
deconvolve = 'time'

# time domain deconvolution options
spiking = 1.  # spiking factor for noise suppression

# frequency domain deconvolution options
water = 0.05  # water level for noise suppression
gauss = 2.  # low pass Gauss filter with corner frequency in Hz

# Additionally you can set the following options for the deconvolution.
# If not set, reasonable values are used depending on options 'method'
# and 'deconvolve'
# winsrc = (-10, 30, 5) # window for source function, (start, end, tapering)
# next two options only for deconvolve='time':
# winrsp = (-20, 80) # window for response
# winrf= (-20, 80) # desired window for receiver function
# next option only for deconvolve='freq':
# tshift = 10 # shift of rf. Peak on L/Z component will be at this position
