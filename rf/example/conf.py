# Example configuration file for batch usage of rf package
# It is possible to use the variables 'method' (equals 'P' or 'S') and
# 'phase' (equals 'P', 'PP' or 'PS', etc.).


### Options for input and output ###

# File name of events file in QuakeML format and filename for inventory
# of stations in StationXML format.
# A receiver function is possibly created for each event-station combination.
events = 'example_events.xml'
inventory = 'example_inventory.xml'

# Root output path and format of output waveform files
# Supported formats: 'Q', 'SAC', 'H5' (obspyh5)
path = '.'
format = 'Q'

# You have to declare a function init which returns a function.
# This returned function gets called with the optional arguments
# network, station, location, channel, seed_id (all strings),
# starttime, endtime (UTCDateTimes) and event (Obspy event).
# It has to return all 3 components of the requested
# data in a stream object.
#
# Example for usage of ObsPy client:
#
# def init()
#     from obspy.fdsn import Client
#     client = Client()
#     def get_waveform(network, station, location, channel, starttime, endtime,
#                      **kwargs):
#         kws = {'network': network, 'station': station, 'location': location,
#                'channel': channel,
#                'starttime': starttime, 'endtime': endtime}
#         return client.getWaveform(**kws)
#     return get_waveform
#
# For this example configuration get_waveform delivers the data from a
# local waveform file


def init():
    from obspy import read
    data = read('example_data.mseed')

    def get_waveform(network, station, location, channel, seed_id,
                     starttime, endtime, event):
        st = data.select(network=network, station=station, location=location)
        st = st.slice(starttime, endtime)
        if len(st) == 3:
            return st
        print st
    return get_waveform

# Request this time window around P-onset or S-onset
request_window = (-50, 150) if method == 'P' else (-100, 50)


### Options for receiver function calculation ###

# Events with an epicentral distance not in dist_range will be discarded
# for a specific station
if phase == 'P':
    dist_range = (30, 90)
elif phase == 'S':
    dist_range = (60, 85)

# Filter stream.
if method == 'P':
    filter = {'type': 'bandpass', 'freqmin': 0.01, 'freqmax': 2.}
else:
    filter = {'type': 'bandpass', 'freqmin': 0.01, 'freqmax': 0.5}


# Trim window around P-onset (P-receiver function)
# or S-onset (S-receiver function)
window = (-30, 100) if method == 'P' else (-80, 50)

# Downsample stream to this frequency in Hz.
downsample = 10

# Roate stream with this method.
# rotate should be one of 'ZNE->LQT' or 'NE->RT'
# or your own rotation function accepting a stream as arg.
rotate = 'ZNE->LQT'

# Deconvolve stream with this method.
# deconvolve is one of 'time' or 'freq' for time or
# frequency domain deconvolution.
# source_component is the component used as source for the deconvolution.
deconvolve = 'time'
source_component = 'L'

# time domain deconvolution options
spiking = 1.  # spiking factor for noise suppression

# frequency domain deconvolution options
water = 0.05  # water level for noise suppression
gauss = 2.  # low pass Gauss filter with corner frequency in Hz

# window for source function relative to onset, (start, end, tapering)
winsrc = (-10, 30, 5)


### Options for other routines defined as dictionaries ###

kwargs_moveout = {}  # See documentation of RFStream.moveout

kwargs_plot = {'trace_height': 0.2, 'norm': 4., 'fill': True,
               'window': (-5, 22)}  # See documentation of RFStream.plot_rf
