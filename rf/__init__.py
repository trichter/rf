"""
rf: Receiver function calculation in seismology
===============================================

This module heavily depends on ObsPy obspy.org
The main functionality is provided by the class `~rf.rfstream.RFStream` which
is derived from the Stream class of ObsPy.
http://docs.obspy.org/packages/autogen/obspy.core.stream.Stream.html

:copyright:
    Tom Richter
:license:
    MIT

Usage
-----

The canonical way to load a waveform file into an RFStream object is to read
it with ObsPy and pass the Obspy Stream to the generator as a kwarg:

>>> from obspy import read
>>> from rf import RFStream
>>> stream = RFStream(read('infile.SAC')

The stream is again written to disc as usual by its write method:
>>> stream.write('outfile', 'SAC')

The RFStream object inherits a lot of useful methods from its ObsPy ancestor
(e.g. filter, taper, simulate, ...).

The module automatically maps important (for rf calculation) header information
from the stats object attached to every trace to the format specific headers.
At the moment only SAC and SH/Q headers are supported. When initializing an
RFStream the header information in the format specific headers are written to
the stats object and before writing the information stored in the stats object
is written back to the format specific headers. In this way the important
header information is guaranteed to be saved in the waveform files.
The following table reflects the mapping:
stats    station_latitude station_longitude station_elevation
         event_latitude event_longitude event_depth event_magnitude event_time
         onset distance back_azimuth inclination slowness
SH/Q     DCVREG DCVINCI PWDW
         LAT LON DEPTH MAGNITUDE ORIGIN
         P-ONSET DISTANCE AZIMUTH INCI SLOWNESS
SAC      stla stlo stel
         evla evlo evdp mag o
         a gcarc baz user0 user1
Note that Q-file headers DCVREG DCVINCI PWDW are used for the station
information, because the Q format has a shortage of predefined headers.

The first task when calculating receiver functions is calculating some ray
specific values like azimuth and epicentral distance. An appropriate stats
dictionary can be calculated with `~rf.rfstats`:

>>> stats = rfstats(station=station, event=event, phase='P', dist_range=(30,90))
>>> for tr in stream:
>>>     tr.stats.update(stats)

or if the station and event information is alread stored in the stats object:

>>> for tr in stream:
>>>     rfstats(stats=tr.stats)

Now Ps receiver function calculation is as easy as:

>>> stream = RFStream(read('infile.SAC')
>>> stream.filter('bandpass', freqmin=0.05, freqmax=1.)
>>> stream.rf()
>>> stream.write('rf', 'Q')

rf can also calculate Sp receiver functions (not much tested):
>>> stream.rf(method='S')

When calling stream.rf the following operations are performed depending on
the given kwargs:
  * filtering
  * trimming data to window relative to onset
  * downsampling
  * rotation
  * deconvolution
Please see `~RFStream.rf` for a more detailed description.

RFStream provides the possibility to perform moveout correction
and piercing point calculation. Xiahou Yuan providing his Fortran
scripts.

Additionally rf provides a function `~rf.rf_batch` which does nearly all
the work for you.

Please feel free to request features, report bugs or contribute some code on
GitHub.
https://github.com/trichter/rf
"""


#try:
#    EXAMPLE_CONFIG = False
#    import config as conf
#except ImportError:
#    import warnings
#    warnings.warn("Didn't find file config.py. Using example configuration.")
#    EXAMPLE_CONFIG = True
#    from rf import config_example as conf

from rfstream import RFStream, rfstats
#from batch import rf_batch
#from io import convert_dmteventfile, create_rfeventsfile, set_paths


