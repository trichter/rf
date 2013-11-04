"""
rf: Receiver function calculation in seismology
===============================================

This module heavily depends on `Obspy <http://www.obspy.org/>`_.
The main functionality is provided by the class :class:`~rf.rfstream.RFStream`
which is derived from ObsPy's :class:`~obspy.core.stream.Stream` class.

:copyright:
    Tom Richter
:license:
    MIT

Installation
------------

After installing `Obspy <http://www.obspy.org/>`_, its dependencies and
`toeplitz <https://github.com/trichter/toeplitz>`_ the package can be installed
with pip by running
::
    pip install rf

Alternatively download the source code and run
::
    python setup.py install

The tests can be run with the script
::
    rf-runtests

Usage
-----

The canonical way to load a waveform file into an RFStream object is to read
it with ObsPy and pass the Obspy Stream to the generator as a kwarg:

>>> from obspy import read
>>> from rf import RFStream
>>> stream = RFStream(stream=read('infile.SAC'))

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

=================  =========  =====
stats              SH/Q       SAC
=================  =========  =====
station_latitude   DCVREG     stla
station_longitude  DCVINCI    stlo
station_elevation  PWDW       stel
event_latitude     LAT        evla
event_longitude    LON        evlo
event_depth        DEPTH      evdp
event_magnitude    MAGNITUDE  mag
event_time         ORIGIN     o
onset              P-ONSET    a
distance           DISTANCE   gcarc
back_azimuth       AZIMUTH    baz
inclination        INCI       user0
slowness           SLOWNESS   user1
=================  =========  =====

.. note::
    Q-file headers DCVREG DCVINCI PWDW are used for the station
    information, because the Q format has a shortage of predefined headers.

The first task when calculating receiver functions is calculating some ray
specific values like azimuth and epicentral distance. An appropriate stats
dictionary can be calculated with :func:`~rf.rfstream.rfstats`:

>>> from rf import rfstats
>>> stats = rfstats(station=station, event=event, phase='P', dist_range=(30,90))
>>> for tr in stream:
>>>     tr.stats.update(stats)

or if the station and event information is already stored in the stats object:

>>> for tr in stream:
>>>     rfstats(stats=tr.stats)

Now P receiver functions can be calculated by

>>> stream.filter('bandpass', freqmin=0.05, freqmax=1.)
>>> stream.rf()
>>> stream.write('rf', 'Q')

rf can also calculate S receiver functions (not much tested):

>>> stream.rf(method='S')

When calling stream.rf the following operations are performed depending on
the given kwargs:
  * filtering
  * trimming data to window relative to onset
  * downsampling
  * rotation
  * deconvolution

Please see :func:`~rfstream.RFStream.rf` for a more detailed description.
RFStream provides the possibility to perform moveout correction
and piercing point calculation. rf is going to provide a function
:func:`~batch.rf_batch` which will run all the necessary steps.


Please feel free to request features, report bugs or contribute some code on
`GitHub <https://github.com/trichter/rf/>`_.
"""


try:
    EXAMPLE_CONFIG = False
    import rfconf
except ImportError:
    import warnings
    warnings.warn("Didn't find file rfconf.py. Using example configuration.")
    EXAMPLE_CONFIG = True
    from rf import rfconf_example as rfconf

from rfstream import RFStream, rfstats
#from batch import rf_batch
#from io import convert_dmteventfile, create_rfeventsfile, set_paths


