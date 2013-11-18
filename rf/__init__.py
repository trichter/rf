"""
rf: Receiver function calculation in seismology
===============================================

This module heavily depends on ObsPy_.
The main functionality is provided by the class :class:`~rf.rfstream.RFStream`
which is derived from ObsPy's :class:`~obspy.core.stream.Stream` class.

:copyright:
    Tom Richter
:license:
    MIT

Installation
------------

Install ObsPy_, its dependencies and pip_, eg. by ::

    sudo apt-get install python-obspy python-pip

rf can then be installed by ::

    pip install rf

The tests can be run with the script ::

    rf-runtests

Manual installation of dev
--------------------------

Install

    * ObsPy_ and its dependencies,
    * toeplitz_ for time domain deconvolution,
    * geographiclib_ for ppoint calculation.

Then download the source code from GitHub_ eg. by ::

    mkdir rf
    git clone https://github.com/trichter/rf.git rf

Now install rf with ::

    python setup.py install

.. _ObsPy: http://www.obspy.org/
.. _pip: http://www.pip-installer.org/
.. _toeplitz: https://github.com/trichter/toeplitz/
.. _geographiclib: https://pypi.python.org/pypi/geographiclib/
.. _GitHub: https://github.com/trichter/rf/

Basic Usage
-----------

The canonical way to load a waveform file into an RFStream object is to use
the :func:`~obspy.core.stream.read` function of ObsPy and pass
the Obspy Stream to the generator of :class:`~rf.rfstream.RFStream` as a
parameter for the stream kwarg:

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
>>> stats = rfstats(station=station, event=event, phase='P',dist_range=(30,90))
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

Please see :meth:`RFStream.rf() <rf.rfstream.RFStream.rf>`
for a more detailed description.
RFStream provides the possibility to perform moveout correction
and piercing point calculation. rf is going to provide a function
:func:`~rf.batch.rf_batch` which will run all the necessary steps.

Batch Usage
-----------

TODO

Miscellaneous
-------------

Please feel free to request features, report bugs or contribute code on
GitHub_. The code is continiously tested by travis-ci
which reports the latest |build hopefully passing|.

.. |build hopefully passing| image::
    https://api.travis-ci.org/trichter/rf.png?branch={travis_version}
    :target: https://travis-ci.org/trichter/rf
"""
# Suggest people to cite rf.

from _version import __version__
from rfstream import RFStream, rfstats

# get image for this version from travis-ci
if 'dev' in __version__:
    _travis_version = 'master'
else:
    _travis_version = 'v' + __version__
__doc__ = __doc__.format(travis_version=_travis_version)
