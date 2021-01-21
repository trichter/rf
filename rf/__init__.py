# Copyright 2013-2019 Tom Eulenfeld, MIT license
"""
rf Documentation
================

rf is a Python framework for receiver function analysis.
Read and write support of necessary metadata is provided for
SAC, SeismicHandler and HDF5 waveform files.
Data is handled by the ``RFStream`` class which inherits a lot of useful
methods from its ObsPy ancestor ``Stream``,
but also has some unique methods necessary for receiver function calculation.


Method
------

The receiver function method is a popular technique to investigate crustal and
upper mantle velocity discontinuities. Basic concept of the method is that a
small part of incident P-waves from a teleseismic event gets converted to
S-waves at significant discontinuities under the receiver
(for P-receiver functions).
These converted Ps phases arrive at the station after the main P phase.
The response function of the receiver side (receiver function) is constructed
by removing the source and deep mantle propagation effects.
Firstly, the S-wave field is separated from the P-wave field by a rotation
from the station coordinate system (ZNE - vertical, north, east)
to the wave coordinate system (LQT - P-wave polarization,
approx. SV-wave polarization, SH-wave polarization). Alternatively, it is
possible to rotate to the ZRT coordinate system (vertical, radial, transverse).
Secondly, the waveform on the L component is deconvolved from the other
components, which removes source side and propagation effects.
The resulting functions are the Q and T component of the P receiver function.
Multiple reflected waves are also visible in the receiver function.
The conversion points of the rays are called piercing points.

For a more detailed description of the working flow see e.g. chapter 4.1 of
this_ dissertation.

.. image:: _static/2layer_rays.svg
   :height: 250px
   :alt: Ray paths

.. image:: _static/2layer_synrf.svg
   :height: 250px
   :alt: Synthetic receiver function

| *Left*: In a two-layer-model part of the incoming P-wave is converted to a
    S-wave at the layer boundary. Major multiples are Pppp, Ppps and Ppss.
| *Right*: Synthetic receiver function of Q component in a two-layer-model.

Installation
------------

Dependencies of rf are

    * ObsPy_ and some of its dependencies,
    * cartopy, geographiclib, shapely,
    * mtspec_ for multitaper deconvolution,
    * toeplitz_ for faster time domain deconvolution (optional),
    * obspyh5_ for hdf5 file support (optional),
    * tqdm for progress bar in batch processing (optional).

rf can be installed with ::

    pip install rf

The tests can be run with the script ::

    rf-runtests

To install the development version of rf download the source code and run ::

    pip install .

Here are some instructions to install rf into a fresh conda environment::

    conda config --add channels conda-forge
    conda create -n rfenv obspy cartopy geographiclib shapely h5py mtspec tqdm fortran-compiler
    conda activate rfenv
    pip install obspyh5 toeplitz rf
    rf-runtests

The ``fortran-compiler`` metapackage helps finding the right fortran compiler
and its dependencies for any OS. Often, problems with the compilation of Fortran
codes (toeplitz package) can be solved by setting
``export NPY_DISTUTILS_APPEND_FLAGS=1``.

Using the Python module
-----------------------

The main functionality is provided by the class `.RFStream`
which is derived from ObsPy's `~obspy.core.stream.Stream` class.

The canonical way to load a waveform file into a RFStream is to use
the `.read_rf()` function.

>>> from rf import read_rf
>>> stream = read_rf('myfile.SAC')

If you already have an ObsPy Stream and you want to turn it into a RFStream
use the generator of RFStream:

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

=================  =========  ======
stats              SH/Q       SAC
=================  =========  ======
station_latitude   COMMENT    stla
station_longitude  COMMENT    stlo
station_elevation  COMMENT    stel
event_latitude     LAT        evla
event_longitude    LON        evlo
event_depth        DEPTH      evdp
event_magnitude    MAGNITUDE  mag
event_time         ORIGIN     o
onset              P-ONSET    a
type               COMMENT    kuser0
phase              COMMENT    kuser1
moveout            COMMENT    kuser2
distance           DISTANCE   gcarc
back_azimuth       AZIMUTH    baz
inclination        INCI       user0
slowness           SLOWNESS   user1
pp_latitude        COMMENT    user2
pp_longitude       COMMENT    user3
pp_depth           COMMENT    user4
box_pos            COMMENT    user5
box_length         COMMENT    user6
=================  =========  ======

.. note::
    Q-file header COMMENT is used for storing some information, because
    the Q format has a shortage of predefined headers.

.. note::
    Alternatively the hdf5 file format can be used. It is supported via the
    obspyh5 package. In this case all supported stats
    entries are automatically attached to the stored data. The plugin also
    saves entries not listed here (e.g. processing information).

The first task when calculating receiver functions is calculating some ray
specific values like azimuth and epicentral distance. An appropriate stats
dictionary can be calculated with `.rfstats()`:

>>> from rf import rfstats
>>> stats = rfstats(station=station, event=event, phase='P', dist_range=(30,90))
>>> for tr in stream:
>>>     tr.stats.update(stats)

or if the station and event information is already stored in the stats object:

>>> rfstats(stream)

A typical workflow for P receiver function calculation and writing looks like

>>> stream.filter('bandpass', freqmin=0.05, freqmax=1.)
>>> stream.rf()
>>> stream.write('rf', 'Q')

rf can also calculate S receiver functions:

>>> stream.rf(method='S')

When calling stream.rf the following operations are performed depending on
the given kwargs:

    * filtering
    * trimming data to window relative to onset
    * downsampling
    * rotation
    * deconvolution

Please see `.RFStream.rf()`
for a more detailed description.
RFStream provides the possibility to perform moveout correction,
piercing point calculation and profile stacking.

Tutorials
---------

The following Jupyter notebooks can be viewed online or downloaded
to be locally reproduced.

    1. Calculate receiver functions - minimal example (notebook1_)
    2. Calculate receiver functions and stack them by common conversion points
       to create a profile (notebook2_)
    3. Calculate and compare receiver functions calculated with different
       deconvolution methods (notebook3_)

Command line tool for batch processing
--------------------------------------

The rf package provides a command line utility 'rf' which runs all the
necessary steps to perform receiver function calculation.
Use
``rf -h`` and ``rf {your_command} -h`` to discover the interface.
All you need is an inventory file (StationXML) and a file with events
(QuakeML) you want to analyze.

The command ::

    rf create

creates a :ref:`template configuration file <config_label>` in the current
directory. This file is in JSON format and well documented.
After adapting the file to your needs you can use the various
commands of rf to perform different tasks (e.g. receiver function
calculation, plotting).

To create the tutorial with a small included dataset and working configuration
you can use ::

    rf create --tutorial

Now start using rf ..., e.g. ::

    rf data calc myrf
    rf moveout myrf myrfmout
    rf plot myrfmout myrfplot
    rf --moveout-phase Psss moveout myrf myrfPsssmout

.. note::
    Development of the batch module has a lower priority than the Python API.

Miscellaneous
-------------

Please feel free to request features, report bugs or contribute code on
GitHub_. The code is continuously tested by Github Actions.

Citation
--------

If you found this package useful, please consider citing it.

Tom Eulenfeld (2020), rf: Receiver function calculation in seismology, *Journal of Open Source Software*, 5(48), 1808, doi: `10.21105/joss.01808`__.

.. __: https://dx.doi.org/10.21105/joss.01808

.. _this: http://www.diss.fu-berlin.de/diss/servlets/MCRFileNodeServlet/FUDISS_derivate_000000014929/dissertation_richter.pdf
.. _ObsPy: http://www.obspy.org/
.. _pip: http://www.pip-installer.org/
.. _obspyh5: https://github.com/trichter/obspyh5/
.. _toeplitz: https://github.com/trichter/toeplitz/
.. _mtspec: https://github.com/krischer/mtspec

.. _notebook1: http://nbviewer.jupyter.org/github/trichter/notebooks/blob/master/receiver_function_minimal_example.ipynb
.. _notebook2: http://nbviewer.jupyter.org/github/trichter/notebooks/blob/master/receiver_function_profile_chile.ipynb
.. _notebook3: https://nbviewer.jupyter.org/github/hfmark/notebooks/blob/main/rf_comparison.ipynb

.. _GitHub: https://github.com/trichter/rf/
"""

__version__ = '1.0.0'

from rf.profile import get_profile_boxes
from rf.rfstream import read_rf, RFStream, rfstats
from rf.util import iter_event_data, IterMultipleComponents

if 'dev' not in __version__:  # get image for correct version from travis-ci
    __doc__ = __doc__.replace('branch=master', 'branch=v%s' % __version__)
