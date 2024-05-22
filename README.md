# rf
## Receiver function calculation in seismology
[![build status](https://github.com/trichter/rf/workflows/tests/badge.svg)](https://github.com/trichter/rf/actions)
[![docs status](https://readthedocs.org/projects/rf/badge/?version=latest)](https://rf.readthedocs.io)
[![codecov](https://codecov.io/gh/trichter/rf/branch/master/graph/badge.svg)](https://codecov.io/gh/trichter/rf)
[![pypi version](https://img.shields.io/pypi/v/rf.svg)](https://pypi.python.org/pypi/rf)
[![python version](https://img.shields.io/pypi/pyversions/rf.svg)](https://python.org)
[![JOSS](http://joss.theoj.org/papers/10.21105/joss.01808/status.svg)](https://doi.org/10.21105/joss.01808)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4455036.svg)](https://doi.org/10.5281/zenodo.4455036)

##### Documentation: https://rf.readthedocs.io/
##### Tutorials:
  1. Calculate receiver functions - minimal example ([notebook][nb1])
  2. Calculate receiver functions and stack them by common conversion points to create a profile ([notebook][nb2])
  3. Compare different deconvolution methods ([notebook][nb3])
  4. Harmonic deconvolution with synthetics - minimal example ([notebook][nb4])

[nb1]: https://nbviewer.jupyter.org/github/trichter/notebooks/blob/master/receiver_function_minimal_example.ipynb
[nb2]: https://nbviewer.jupyter.org/github/trichter/notebooks/blob/master/receiver_function_profile_chile.ipynb
[nb3]: https://nbviewer.jupyter.org/github/hfmark/notebooks/blob/main/rf_comparison.ipynb
[nb4]: https://nbviewer.jupyter.org/github/hfmark/notebooks/blob/main/rf_harmonics.ipynb

##### Get help and discuss: [ObsPy Related Projects category](https://discourse.obspy.org/c/obspy-related-projects/rf/14) in the ObsPy forum

##### Contribute:

All contributions are welcome ... e.g. report bugs, discuss or add new features.

##### Citation:

If you found this package useful, please consider citing it.

Tom Eulenfeld (2020), rf: Receiver function calculation in seismology, *Journal of Open Source Software*, 5(48), 1808, doi: [10.21105/joss.01808](https://doi.org/10.21105/joss.01808) [[pdf]](https://www.theoj.org/joss-papers/joss.01808/10.21105.joss.01808.pdf)

##### Related receiver function projects

* [seispy](https://github.com/xumi1993/seispy) including hk-stacking
* [RFPy](https://github.com/paudetseis/RfPy) including hk-stacking, harmonic decomposition
* [BayHunter](https://github.com/jenndrei/BayHunter) inversion of receiver functions and surface wave dispersion
* [telewavesim](https://github.com/paudetseis/Telewavesim) synthetics
