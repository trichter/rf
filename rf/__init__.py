"""
rf: Receiver function calculation in seismology
===============================================

This module heavily depends on ObsPy_.
The main functionality is provided by the class :class:`~rf.rfstream.RFStream`
which is derived from ObsPy's :class:`~obspy.core.stream.Stream` class.

:Author: Tom Eulenfeld
:License: MIT
:Documentation: http://rf.readthedocs.org/ (installation instructions, usage,
    examples and module reference)
:Project page: https://github.com/trichter/rf (bug reports and feature requests
    via issue tracker)
:Pypi page: https://pypi.python.org/pypi/rf
:Test status: |buildstatus|

Installation
------------

Install ObsPy_, its dependencies and pip_, e.g. by (Ubuntu) ::

    sudo apt-get install python-obspy python-pip

rf can then be installed by ::

    pip install rf

The tests can be run with the script ::

    rf-runtests

Hdf5 file support can be installed with the obspyh5_ package.

Manual installation of dev
--------------------------

Install

    * ObsPy_ and its dependencies,
    * toeplitz_ for time domain deconvolution,
    * geographiclib_ for ppoint calculation,
    * optionally obspyh5_ for hdf5 file support.

Then download the source code from GitHub_ eg. by ::

    mkdir rf
    git clone https://github.com/trichter/rf.git rf

Now install rf with ::

    python setup.py install

Using the underlying Python module
-----------------------------------

The canonical way to load a waveform file into a RFStream is to use
the :func:`~rf.rfstream.read_rf` function.

>>> from rf import read_rf
>>> stream = read_rf('myfile.SAC')

If you already have an ObsPy Stream and you want to turn it into a RFStream
use the generator of :class:`~rf.rfstream.RFStream`:

>>> from rf import RFStream
>>> stream = RFStream(obspy_stream)

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

.. note::
    Alternatively the hdf5 file format can be used. It is supported via the
    obspyh5 package. In this case all supported stats
    entries are automatically attached to the stored data.

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

Please see :meth:`RFStream.rf() <rf.rfstream.RFStream.rf>`
for a more detailed description.
RFStream provides the possibility to perform moveout correction
and piercing point calculation.

Batch Usage
-----------

rf provides a command line utility 'rf' which runs all the necessary steps
to perform receiver function calculation. All you need is a StationXML
file with your stations as an inventory and a QuakeML file with the events
you want to analyze.

The command ::

    rf create

creates a template configuration file in the current directory. This file is
in JSON format and well documented.
After adapting the file to your needs you can use the various
subcommands of rf to perform different tasks (e.g. receiver function
calculation, plotting).

To create the tutorial with a small included dataset and working configuration
you can use ::

    rf create --tutorial

Now start using rf ..., e.g. ::

    rf calc
    rf moveout
    rf plot Prf_Ps
    rf --moveout Psss moveout
    rf plot Prf_Psss

Miscellaneous
-------------

Please feel free to request features, report bugs or contribute code on
GitHub_. The code is continiously tested by travis-ci.


.. _ObsPy: http://www.obspy.org/
.. _pip: http://www.pip-installer.org/
.. _obspyh5: https://github.com/trichter/obspyh5/
.. _toeplitz: https://github.com/trichter/toeplitz/
.. _geographiclib: https://pypi.python.org/pypi/geographiclib/
.. _GitHub: https://github.com/trichter/rf/
.. |buildstatus| image:: https://api.travis-ci.org/trichter/rf.png?
    branch=master
   :target: https://travis-ci.org/trichter/rf
"""
# Suggest people to cite rf.

from _version import __version__
from rfstream import read_rf, RFStream, rfstats

if 'dev' not in __version__:  # get image for correct version from travis-ci
    _travis_version = 'v' + __version__
    __doc__ = __doc__.replace('branch=master', 'branch=v%s' % __version__)
