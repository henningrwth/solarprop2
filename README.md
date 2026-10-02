SOLARPROP 2.0
=========

Solarprop offers a framework for building models of charge-sign dependent solar modulation of cosmic rays in the heliosphere. The transport equation is solved with the method of stochastic differential equations (SDEs). Models are implemented and the SDE solution is calculated in a backend written in C++, but a Python interface is provided for applying the models to the modulation of cosmic-ray fluxes.

Solarprop is *fast* because the C++ backend is parallelized with OpenMP. It is *accessible* due to its Python bindings, allowing users to set up a model and calculating a modulated flux with a few lines of simple Python code. It is *transparent* and *versatile* because the entire source code is available publicly online and the C++ code is structured such that new models can be added easily.


Installation
------------

It is highly recommended to create and activate a virtual environment first.

After obtaining a copy of the source code, go to the main directory (the one containing this `README.md` file).

Install the lookups of time-dependent data on heliospheric observables:

    pushd data && ./download_data.sh && popd

Then install the `solarprop` python module with `pip`:

    pip install . -v

Useful tools for analyzing cosmic-ray data in conjunction with Solarprop are contained in the `sproptools` package, which is also needed for the model validation scripts. Install it as follows:

    pushd sproptools
    pip install .
    popd

The C++ backend comes with Doxygen documentation. To generate it, run

    doxygen

and to view it, point your browser to `doc/html/index.html`.

Run

    SOLARPROP --config

to test the installation.


Detailed description
--------------------

The physics background of Solarprop is described in detail in an [accompanying preprint](https://arxiv.org/abs/2610.01360). This covers the transport equation governing charge-sign dependent solar modulation, its solution using SDEs, and validation studies with comparisons to previously published models. It also introduces the nomenclature used in the code and the documentation and discusses the efficiency of the parallelization of the code.


Getting started
---------------

Go to the `validation` directory and run a few scripts comparing Solarprop calculations to published model predictions:

    ./test_yamada_fig4.py
    ./test_jokipii_kopriva_fig2.py
    ./test_potgieter_moraal_fig6.py
    ./test_burger_potgieter_fig5.py
    ./test_burger_2012_fig2.py

See the docstrings of these scripts or the [list](#list-of-model-implementations) of model implementations for the references to the original literature.


Using Solarprop
---------------

The following python script illustrates the basic usage of Solarprop's python interface. In the script, we define a local interstellar spectrum and calculate its values for a list of logarithmically equidistant values of kinetic energies. We then choosesa model implementation and set the relevant parameters. Finally, we perform the modulation run to obtain the modulated flux at the energy points used to define the LIS, and print the results:

```python
#! /usr/bin/env python3

import numpy as np

from sproptools.kinematics import Species
import solarprop

# cosmic-ray particle species
proton = Species.from_name('proton')

# analytic form of the local interstellar spectrum (LIS)
def lis(T):
    m = proton.mass
    return np.sqrt(T**2 + 2*T*m) * np.pow(T + m, -3.6) * 10.

# create logarithmically equidistant points in kinetic energy for our spectra
ekin = np.geomspace(0.01, 25.0, num=30)

# calculate the LIS values
lis_T = lis(ekin)

# choose a model implementation and set its options and parameters
parameters = { 'model': 'Burger2012',    # name of the model implementation
               'mass': proton.mass,      # mass of the cosmic-ray particle (in units of GeV)
               'charge': proton.charge,  # charge of the cosmic-ray particle (in units of the proton charge)
               'polarity': -1,           # polarity of the heliospheric magnetic field: A<0
               'particlesPerBin': 1000,  # number of pseudo-particles used for the SDE integration per energy bin
               'angle': 30.,             # tilt angle of the heliospheric current sheet (in degrees)
              }

# create an instance of the Solarprop interface, passing the parameters
s = solarprop.Solarprop(parameters)

# perform the modulation run and obtain the modulated flux at the energy points used to define the LIS
flux_mod = s.modulate(ekin, lis_T)

# print the results
print('Ekin (GeV)   modulated flux (a.u.)   LIS flux (arbitrary units)')
for r in zip(ekin, flux_mod, lis_T):
    print(f'{r[0]:10.4f}   {r[1]:14g}    {r[2]:14g}')
```


List of model implementations
-----------------------------

The following models have already been implemented in Solarprop:

| Model name     | Description                                                                           |
|----------------|---------------------------------------------------------------------------------------|
| `Yamada1998`   | Simple one-dimensional model, published in [Y. Yamada et al., Geophys. Res. Lett. Vol. 25 No. 13 (1998) 2353-2356](https://doi.org/10.1029/98GL51869). (This model was called `ref1` in vanilla Solarprop.) |
| `JK1979`       | Two-dimensional model reproducing the results of [J.R. Jokipii and D.A. Kopriva, ApJ 234 (1979) 384-392](https://adsabs.harvard.edu/pdf/1979ApJ...234..384J) (This model was called `ref2` in vanilla Solarprop.) |
| `PM1985`       | Two-dimensional model reproducing the results of [M.S. Potgieter and H. Moraal, ApJ 294 (1985) 425-440](https://adsabs.harvard.edu/abs/1985ApJ...294..425P) (This model was called `ref3` in vanilla Solarprop.) |
| `BP1989`       | Two-dimensional model reproducing the results of [R.A. Burger and M.S. Potgieter, ApJ 339 (1989) 501-511](https://adsabs.harvard.edu/full/1989apj...339..501b) (This model was called `ref4` in vanilla Solarprop.) To be more precise, this model aims at reproducing Figure 5 in the reference paper. The figure was done with the same parameters as used in [Kota and Jokipii, ApJ 265 (1983) 573](https://doi.org/10.1086/160701), but with the new implementation of current sheet drift by Burger and Potgieter, see also [Burger and Hattingh, Astrophys. Space Sci. 230 (1995) 375-382](https://doi.org/10.1007/BF00658195). |
| `standard2D`   | A simple two-dimensional model introduced in vanilla Solarprop: [R. Kappl, Comp. Phys. Comm, 207 (2016) 386-399](https://doi.org/10.1016/j.cpc.2016.05.025). |
| `Burger2012`   | Three-dimensional model that includes particle drifts along a wavy heliospheric current sheet, for a simple Parker spiral, published in [R.A. Burger, ApJ 760:60 (2012)](https://doi.org/10.1088/0004-637x/760/1/60). |
| `Strauss2012`  | Solarprop implementation of the model described in [R. D. Strauss et al., Astrophys. Space Sci. (2012) 339:223-236](https://doi.org/10.1007/s10509-012-1003-z). The model includes a detailed treatment of particle drifts along the heliospheric current sheet. Note that there are discrepancies between the Solarprop calculations and the results shown in the reference publication. |
| `Helmod2012`   | Solarprop implementation of the model described in [Bobik et al., ApJ 745 (2012) 132](https://doi.org/10.1088/0004-637X/745/2/132). The model uses a Parker spiral magnetic field with polar modifications according to Jokipii and Kóta. Note that there are discrepancies between the Solarprop calculations and the results shown in the reference publication. |
| `Helmod2018`   | Solarprop implementation of the model described in [Boschini et al., Adv. Space Res. 62 (2018) 2859-2879](https://doi.org/10.1016/j.asr.2017.04.017). The model is similar to `Helmod2012` but uses a different spatial dependence for the diffusion coefficient. Note that there are discrepancies between the Solarprop calculations and the results shown in the reference publication. |
| `custom`       | Template for a user-defined model. |

To implement your own model, follow the instructions in `src/models/custom.h` and fill the placeholders in `src/models/custom.cc`.


Model options and parameters
============================

Global options
--------------

The following options and parameters can be used globally, though most of the validation models use fixed values for their parameters:

| option                    | description                                                                                                            |
|---------------------------|------------------------------------------------------------------------------------------------------------------------|
|`model`                    | Model name. (Note that all models have been renamed for clarity. E.g., the `ref1` model is now called `Yamada1998`. Run `SOLARPROP --config` to see the list of available models, and refer to the [list](list-of-model-implementations) of model implementations.)|
|`mass`                     | Mass of the particle (in units of $\mathrm{GeV}/c^2$).         |
|`charge`                   | Charge of the particle (in units of the proton charge).|
|`polarity`                 | Override the polarity of heliospheric magnetic field (+1 or -1).|
|`phase`                    | If used together with the `polarity` option: specify phase of solar activity (ascending: +1, descending: -1, unknown: 0).|
|`timestamp`                | Time stamp (unix time) for solution. (Mutually exclusive with `year`/`month`/`day` options.)|
|`year`                     | Year for solution. (Mutually exclusive with `timestamp` option.)|
|`month`                    | Month in `year` for solution. (Mutually exclusive with `timestamp` option.)|
|`day`                      | Day in `month` for solution. (Mutually exclusive with `timestamp` option.)|
|`dt`                       | (Maximum) time step for SDE integration (in seconds).         |
|`dynamicStep`              | Use dynamic step size? (Implemented for 2D SDEs.)                   |
|`save`                     | Save Green function matrix to disk using the given file name.|
|`green`                    | Read Green function matrix from given file instead of performing the modulation run again.|
|`kappaScaling`             | Scaling factor for the diffusion tensor (applied after all other sources for the diffusion tensor).|
|`BfieldScaling`            | Scaling factor for the magnetic field amplitude.    |
|`particlesPerBin`          | Number of pseudo particles for each energy bin.     |
|`extraBins`                | Number of extra bins for quality improvement. (The default value of 5 is usually sufficient for input spectra with 20 to 30 energy points.)|
|`globalSeed`               | Number to be added to the unique particle number to obtain the random seed for each particle. |
|`runNumber`                | When multiple runs of Solarprop are performed and the results combined, use this option to set a unique run number for each run. This ensures unique random seeds for all pseudo-particles.|
|`angleFile`                | Input file for time-dependent data on the HCS tilt angle.|
|`tiltModel`                | Set to `R` or `L` to choose the model for the determination of the HCS tilt angle. (Default: `R`)|
|`angle`                    | Override value for HCS tilt angle (in degrees).|
|`nmFile`                   | Name of the lookup file for the time-dependent modulation potential as derived from neutron monitor data.|
|`nmValue`                  | Override neutron monitor value (in units of $\mathrm{MV}$).               |
|`kappaFile`                | Filename of a lookup file with values of $\kappa_0$ as a function of time stamp. Format: First column contains unix time stamps, second column contains $\kappa_0$ in units of $\mathrm{cm}^2/\mathrm{s}/\mathrm{GV}$.|
|`kappa0`                   | Override value of the diffusion constant scale $\kappa_0$ (in units of $\mathrm{cm}^2\,\mathrm{s}^{-1}\,\mathrm{GV}^{-1}$). |
|`omniFile`                 | Lookup file containing OmniWEB data with time-dependent values for the magnetic field amplitude and the solar wind speed.|
|`ssnFile`                  | Lookup file with time-dependent values of the smoothed sunspot number.|
|`forceNumericalDerivatives`| Always use numerical calculation of derivatives? Useful for debugging and testing. (Only available for selected models.)|
|`verbosity`                | Amount of information to be printed. Levels `0` to `3`: print increasing amount of information on settings and progress. Level `4`: print final steps for each pseudo-particle. Level `5`: print all intermediate steps.|

Options for individual models
-----------------------------

The `Burger2012` model implementing a wavy heliospheric current sheet takes the following options:

| option                     | description                                                                                     |
|----------------------------|-------------------------------------------------------------------------------------------------|
|`heliosphereBoundary`       | Radial extent of the heliosphere (in units of AU). |
|`kappaExponent`             | Set the exponent $a$ for the scaling of diffusion coefficient with rigidity $R$: $\kappa\propto R^a$. |

The `Strauss2012` model implementing a microscopic treatment of current sheet drifts takes the following options:

| option                     | description                                                                                     |
|----------------------------|-------------------------------------------------------------------------------------------------|
|`heliosphereBoundary`       | Radial extent of the heliosphere (in units of AU). |
|`kappaExponent`             | Set the exponent $a$ for the scaling of diffusion coefficient with rigidity $R$: $\kappa\propto R^a$. |
|`hcsFactor`                 | Fudge factor multiplied to the HCS drift velocity. (The default is 1.)|
|`simpleHcsDrift`            | Use the simplified HCS drift scheme, where an averaged velocity field is used instead of tracing particles along the HCS? The simplified treatment of current sheet drifts is done according to the "BP" model, see section 4.2 of Burger and Hattingh, Astrophys Space Sci 230 (1995) 375-382.|

The `Helmod2012` and `Helmod2018` models using a modified Parker spiral for the magnetic field accept the following options in addition:

| option                     | description                                                                                     |
|----------------------------|-------------------------------------------------------------------------------------------------|
|`highSolarActivityThreshold`| Smoothed sunspot number above which solar activity is considered high.|
|`delta0`                    | Amplitude (numerator) in delta(theta) function for latitudinal component of magnetic field.|
|`vswMin`                    | Override value of solar wind speed at heliospheric equator (in units of $\mathrm{km}/\mathrm{s}$).|
|`Bearth`                    | Override amplitude of magnetic field at Earth (in units of $\mathrm{nT}$).|
|`iotaPolarRegion`           | Enhancement factor for polar perpendicular diffusion coefficient $\kappa_{\perp\theta}$ in polar region.|
|`modelSpecificTimestep`     | Override the time step for the SDE integration by the model prescription of $\mathrm{d}t=r^2/\kappa_{rr}$?|
|`parameterLookup`           | Lookup file specifying parameters as a function of time. See the [section on parameter lookup files](#lookup-files-for-parameters).|
|`driftFactor`               | Drift suppression factor (unity for normal unsuppressed drifts). Set to `lookup` to interpolate the values according to the lookup file specified by the `parameterLookup` option.|
|`rho`                       | Factor between radial perpendicular and parallel diffusion coefficients. Set to `lookup` to interpolate the values according to the lookup file specified by the `parameterLookup` option.|
|`glow`                      | Softening term $g_\mathrm{low}$ for diffusion at low rigidities (`Helmod2018` model only). Set to `lookup` to interpolate the values according to the lookup file specified by the `parameterLookup` option.|
|`kappaScaling`              | Scaling factor for the diffusion tensor, may be set to `lookup` to interpolate the values according to the lookup file specified by the `parameterLookup` option.|


Lookup files for parameters
---------------------------

If the model supports it, time-dependent model parameters can be interpolated for the given point in time based on values contained in a lookup file. The lookup file must have the following structure:

    timestamp  parname1  parname2  parname3
    1306800000 value value value
    1307145600 value value value
    etc

The first column in the header line must always be called `timestamp`. The remaining optons, `parname1` and so on, may be replaced by the parameter names. On the C++ side, the model can then access the
interpolated parameter values by calls like

~~~~c++
double defaultValue = 1.0;
double parname1 = modelInfo.getFloatParameterFromOptionOrLookupOrDefault("parname1", defaultValue, timestamp);
~~~~

provided that the user sets the `parname1` option to "`lookup`" for the modulation run..


Control file structure
----------------------

For the `SOLARPROP` executable and the visualization scripts (for example, `plot_solarprop_trajectories`), model options can be set using a control file, which is then passed using the `-c` option.
Control files with model options have the following structure:

    [solarprop]
    parameter1 = value
    parameter2 = value

See the [list](#model-options-and-parameters) of model options and parameters.


Data files
----------

Solarprop uses a set of data files containing default time-dependent data on
  - the tilt angle of the heliospheric current sheet as computed by [WSO](http://wso.stanford.edu/Tilts.html),
  - the smoothed sunspot number from [SIDC](https://www.sidc.be/SILSO/DATA/SN_ms_tot_V2.0.txt),
  - the amplitude of the heliospheric magnetic field at Earth and the solar wind speed, extracted from [OmniWEB](https://omniweb.gsfc.nasa.gov/form/dx1.html
), and
  - the modulation potential calculated from an analysis of [neutron monitor data](https://cosmicrays.oulu.fi/phi/Phi_mon.txt).

The data files have to be downloaded as the first step in the installation process:

    pushd data
    ./download_data.sh
    popd

The data files will then be installed alongside the python module. To use them for a model that supports time-dependent parameter values, you can use the following
code in your python script:

```python
from importlib import resources
import solarprop

parameters = {
    'angleFile': str(resources.files(solarprop).joinpath('solarprop-data/angle.dat')),
    'ssnFile': str(resources.files(solarprop).joinpath('solarprop-data/ssn.dat')),
    'omniFile': str(resources.files(solarprop).joinpath('solarprop-data/omniweb.dat')),
    'nmFile': str(resources.files(solarprop).joinpath('solarprop-data/nm.dat')),
    }
```

The point in time for the modulated flux then has to be defined using either the `timestamp` or the `year` (and optionally, `month` and `day`) options.


Technical background
====================

The python module for Solarprop is built with the [CMake](https://cmake.org/) build system. The main entry point to the build system is the file `CMakeLists.txt` in the main directory. The python bindings are created using [pybind11](https://pybind11.readthedocs.io), with the main code for this in `src/solarprop.h` and `src/solarprop.cc`. See the Doxygen documentation of the C++ backend for a quick introduction to the structure of the C++ side. Usually, the build will not be invoked by calling `cmake` directly, but by installing Solarprop with `pip` as described above. The installation of the python module uses [scikit-build-core](https://github.com/scikit-build/scikit-build-core) as the build backend, while the `sproptools` module uses `hatchling`. The file `sprop.clang-format` defines a coding style for Solarprop that can be used as input to `clang-format` in order to adapt the coding style automatically, e.g.:

    clang-format --style=file:sprop.clang-format -i src/solarprop.cc

When `cmake` and `make` are invoked manually, a C++ standalone executable `sprop_test` useful for profiling and debugging will be built in addition to the usual python module. See [this section](#configuring-for-qtcreator) for an example.


Directory structure
-------------------

The following is a brief summary of the source code directory structure. Note that there is a `CMakeLists.txt` file for the `CMake` build system in each sub-directory.

    main directory
    |
    |-- calculations: ipython notebooks for SDE derivation and analytic calculations relevant for
    |                 individual models.
    |
    |-- data: default data files, e.g. for smoothed sunspot number vs time, will be installed along
    |         with the python module.
    |
    |-- doc: will be created by the call to doxygen and contains the auto-generated documentation for
    |        the C++ backend, entry point is html/index.html.
    |
    |-- scripts: main SOLARPROP executable (not needed when python interface is used) and
    |            visualization scripts.
    |
    |-- sproptools: a pure-python module providing utilities for the analysis of cosmic-ray data, in
    |               particular in conjunction with Solarprop. Used by the validation scripts.
    |
    |-- src: source code of C++ backend.
    |    |
    |    |-- models: individual model implementation go here.
    |
    |-- validation: python scripts comparing Solarprop calculations to published model results.


Configuring for QtCreator
-------------------------

To configure the project for `QtCreator`, use standalone `cmake` in a dedicated build directory, e.g.:

    cd ..
    mkdir solarprop-build
    cd solarprop-build
    cmake -Dpybind11_DIR=$(python3 -c 'import site; print(site.getsitepackages()[0]);')/pybind11/share/cmake/pybind11 ../solarprop
    cd ../solarprop
    qtcreator CMakeLists.txt

Running `make` in the `solarprop-build` directory will build the usual python module and in addition, the `sprop_test` C++ standalone executable. The coding style as used in Solarprop C++ code can be imported to QtCreator from `Solarprop_coding_style_qtcreator.xml`. For newer versions of QtCreator, the coding style file `sprop.clang-format` has to be used.


System of units
---------------

Solarprop comes with a system of units used consistently throughout the C++ backend. It is defined in `src/units.h`. To apply a unit to a quantity, simply multiply it to the numerical value. When a dimensional quantity is printed, its value should be divided by a meaningful unit and the unit included in the print statement. This is illustrated in the following example:

~~~~~~c++
#include "units.h"
using namespace Units;

double kappaZero = 5.e21 * cm2/s/GV;

std::cout << "Value of the diffusion coefficent kappa0= " << kappaZero / (cm2/s/GV) << " cm2/s/GV" << std::endl;
~~~~~~

Arguments of functions used to pass information from Python code to the C++ backend have a name explicitly stating the required unit, for example:

~~~~~~c++
py::array Solarprop::modulate(const std::vector<double>& ekin_lis_in_GeV, const std::vector<double>& flux_lis)
~~~~~~

Likewise, model options for parameters with units state the required unit, and on extracting the respective value in the C++ backend, the unit is applied immediately, e.g.:

~~~~~~c++
#include "units.h"
using namespace Units;

double kappaZero = modelInfo.getFloatOption("kappa0") * cm2/s/GV;
~~~~~~

Return values from the C++ backend are given in units defined by this unit system and have to be converted manually, though this will often not be necessary due to the choices for the basic units (for example, `GeV = 1.0`), and because the most common return value is a cosmic-ray flux, which will be returned in the same units as the input flux anyway. If a conversion is necessary, the unit system can be used in python code like this:

```python
import solarprop
from solarprop.units import AU

s = solarprop.Solarprop(parameters)
endpoints = s.endpoints(number, ekin)

# histogram the distance of closest approach to the Sun, in astronomical units (AU)
rmin = np.array([t.rmin() for t in endpoints if not t.isInsideHeliosphere()])
plt.hist(rmin / AU)
```

Note that the handling of units is somewhat inconsistent, for reasons of simplicity. For example, `kinematics.Species` from the `sproptools` package uses natural units, i.e. provides particles masses in units of $\mathrm{GeV}$, because this is the unit expected for kinetic energies as Solarprop input anyway, while the Solarprop C++ backend internally uses SI units throughout, i.e. masses are given in units of $\mathrm{GeV}/c^2$.  Particle rigidities are assumed to be in units of $\mathrm{GV}/c$ by `kinematics.Species`.


List of changes in version 2.0
------------------------------

- Add OpenMP support for parallelization of the main particle tracking loop. This results in a speedup close to the number of CPU cores.

- Speedup by avoiding unnecessary calculations for particles in additional energy bins created by rebinning the original binning. A fine binning is needed for the histogram of particle energies at the edge of the heliosphere, to get a numerically accurate Green function matrix. But in vanilla Solarprop, pseudo-particles were also tracked for finer bins created from the TOA binning, even though the results of their tracking were never used. Typical speedup factor of 5.

- Add Python interface.

- Switch to CMake build system.

- Provide a script for downloading up-to-date data files from their respective sources.

- Introduce a system of units, used throughout the C++ side.

- Major code refactoring and restructuring. Rename models. Add factory class for models.

- Generalize 2D SDEs to the case of $\kappa_{r\theta}\neq 0$. Add 3D SDEs for models without polar diffusion, $\kappa_{r\theta}=\kappa_{\theta{}r}=\kappa_{\theta\phi}=\kappa_{\phi\theta}=0$.

- New validation model for 3D SDEs: `Burger2012`.

- Additional model implementations: `Strauss2012`, `Helmod2012`, and `Helmod2018`.

- Option for numerical calculation of derivatives needed for diffusion and drift terms.

- Support for reading local interstellar input spectrum in FITS format was dropped on the C++ side. If needed, users can easily read FITS files on the python side.

- Change SDEs to calculate $\mathrm{d}\cos\theta=-\sin\theta\,\mathrm{d}\theta$ instead of $\mathrm{d}\theta$, to avoid unphysical jumps in $\theta$ when close to zenith/nadir.

- Switch SDEs to use $\mathrm{d}p$ instead of $\mathrm{d}T$, to avoid loss term (see Kopp et al. (2012) and Bobik et al. (2016)).

- Add scripts for visualizing pseudo-particle trajectories as well as the individual contributions to the transport process.

- Add options for saving and loading the Green matrix. Allow subsequent modulation of different local interstellar spectra for the same model setup without re-running the SDE calculations.

- Make sure that pseudo-particle trajectories end right at the heliospheric boundary, avoiding large (unphysical) energy changes for high-energy particles.

- Use dynamic time steps: reduce step size close to the Sun in selected models.

- Switch to completely new interpolation scheme. This avoids numerical problems at intermediate energies:
   - introduce Binning class,
   - convert kinetic energy input points to a binning,
   - create a finer version of the binning in `Spectrum::quality()` and use its logarithmic bin centers as new points,
   - use binning when weighting fluxes,
   - use logarithmic interpolation to calculate intermediate fluxes,
   - smear energies according to a local power law, to solve numerical issues at high energies, where the energy loss is small compared to the distance between flux points.

- Fixes and improvements for existing model implementations:
  - Fix HCS drift term in `ref4` (new name: `BP1989`) model
  - Improve `standard2D` model:
     - Optimize model implementation for speed. Speedup is roughly factor of 2.
     - Improve numerical stability:
       - adjust time step close to Sun to avoid large steps and thus $1/r$ divergence in $\mathrm{d}p$,
       - suppress enormous HCS drift close to Earth, which leads to unphysical results.
  - Yamada model: Calculate $\kappa_0$ from the value of the modulation potential inferred from neutron monitor data. This gives results very similar to the force-field approximation.

- New options for `SOLARPROP` executable: Print configuration information with `SOLARPROP --config`. Add `-p` option ("parameters") to remove need for control file. Takes a comma-separated list of key=value pairs.



Authors
=======

Solarprop was originally developed by Rolf Kappl at the University of Bonn. It was modernized and substantially extended by Henning Gast at RWTH Aachen University.

Homepage of vanilla SOLARPROP:
http://www.th.physik.uni-bonn.de/nilles/people/kappl/
