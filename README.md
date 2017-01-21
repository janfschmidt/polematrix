# polematrix
#### a spin tracking code for electron accelerators using Elegant

*Copyright (C) 2017 Jan Felix Schmidt <janschmidt@mailbox.org>, GNU General Public License*

polematrix is a spin tracking code. It calculates the motion of the spin of ultra
relativistic electrons in particle accelerators. It implements a simple matrix tracking as
multi-thread C++ application. The input can be any lattice file from
[Elegant](http://www.aps.anl.gov/Accelerator_Systems_Division/Accelerator_Operations_Physics/manuals/elegant_latest/elegant.html)
or Mad-X. All required particle trajectories are also imported from the established
particle tracking codes Elegant or Mad-X.

polematrix includes synchrotron radiation via import from Elegant or an own implementation
of longitudinal phasespace. A short description of polematrix is given in [this IPAC'16
contribution](http://jacow.org/ipac2016/papers/mopor046.pdf) and in my phd thesis, which
will be available soon (in German).

Do not hesitate to contact me if you have any questions and please report bugs.
A user guide is coming soon.



## Quick-Info

####Dependencies:

- CMake
- Gnu Scientific Library (GSL)
- Boost program options
- Boost filesystem
- Boost property tree
- Boost random
- [Armadillo](http://arma.sourceforge.net/) library
- [palattice](https://github.com/janfschmidt/palattice) library

So far, polematrix was tested only with GCC compiler under Ubuntu/Debian Linux.

####Installation:

```
cd polematrix/build
cmake ..
make
sudo make install
```



## Usage

```
polematrix [CONFIGURATION FILE] [options]
```
*[CONFIGURATION FILE]* is an easy to read xml file, which holds the tracking parameters.
A template configuration file can be generated with option `--template`

###Allowed options:

####Program modes:

```
-h [ --help ]                  display this help message
-v [ --version ]               display version
-T [ --template ]              create config file template (template.pole) and quit
-R [ --resonance-strengths ]   estimate strengths of depolarizing resonances
```

####Configuration options:

```
-t [ --threads ] arg (=all)    number of threads used for tracking
-o [ --output-path ] arg (=.)  path for output files
-V [ --verbose ]               more output, e.g. each written spin file
-n [ --no-progressbar ]        do not show progress bar during tracking
-a [ --all ]                   write all output (e.g. lattice and orbit)
-s [ --spintune ] arg	       in resonance-strengths mode: calculate for given spin tune only
```
