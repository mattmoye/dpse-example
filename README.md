# dynamical state and parameter estimation example

[![DOI](https://zenodo.org/badge/85078265.svg)](https://zenodo.org/badge/latestdoi/85078265)

This directory is a complete example of an analysis to estimate the parameters
and state variables of a conductance-based neuron model using intracellular
current clamp data. The code and data here were used in "Estimating Parameters
and Predicting Membrane Voltages with Conductance-Based Neuron Models" by C
Daniel Meliza, Mark Kostuk, Hao Huang, Alain Nogaret, Daniel Margoliash, and
Henry D I Abarbanel,
[doi:10.1007/s00422-014-0615-5](http://doi.org/10.1007/s00422-014-0615-5). The
code is free for use under the terms described in the file `LICENSE`. The data
are part of the public domain and free for use.

This file is a summary of the steps you need to run the analysis. For a
more complete description of the algorithm, please consult Bryan Toth's
dissertation, supplied here in the file `tothdissertation0511.pdf`.


## Install prerequisites

### Docker image

The quickest way to get started is to run the software using a pre-built
[Docker](https://www.docker.com/) container. Assuming you have Docker installed,
simply run `docker run -it -v `pwd`:/app --rm dmeliza/dpse` from this directory.
Note that running the model can require a lot of memory, so be sure to allocate
at least 16 GB to the container.

### Manual build

You will need to have g++, IPOPT, Python, and sympy installed. This
example has been tested with g++ 4.3, IPOPT 3.10.2, Python 2.7, and
sympy 0.7.1. running on Debian 7.

All of the prerequisites are available through most OS package managers.
We recommend Debian Linux or OS X. IPOPT requires some facility with the
shell to install, but the instructions at
<http://www.coin-or.org/Ipopt/documentation/node10.html> are reliable if
followed carefully. You need to install BLAS, LAPACK, and the MA57
solver at a minimum. Your mileage may vary with other solvers.

A `Dockerfile` is included for help in building the prerequisites in a
container. You'll need to edit it to set the link to the HSL solvers (which
you'll need to download from somewhere; they are only distributed by request).
Run `docker build -t dpse .` to build from these instructions.

## Specify the model

The model comprises a set of ordinary differential equations, which are
specified in a file called `equations.txt`. This directory contains an
example model with 12 state variables and 73 parameters. The file is
commented. The model distributed here is the 9-current model described
in the manuscript.

## Generate the estimation code

Next, a python program (originally by Bryan Toth) converts the symbolic
description of the model in `equations.txt` into C++ functions.
Essentially this involves computing a lot of derivatives to get the
Jacobian and Hessian. To run this step:

``` {.example}
python model/makecode.py
```

Expect this to take a long time, but you only need to run this when the
`equations.txt` file changes. The makecode script expects to find
`equations.txt` in the current directory.

## Compile the estimation code

Successful completion of the previous step will result in generation of
a number of C++ files and a Makefile. Run `make` to compile the code to
a binary. The Makefile expects to find IPOPT using `pkg-config`. If you
have problems you may need to edit the `COIN_LDFLAGS` and `COIN_CFLAGS`
lines in the Makefile. You will also need a reasonably modern C++
compiler (e.g. gcc 4.3 or later). Older compilers may run out of memory
when used with large models.

## Run the estimation procedure

If the previous step completes successfully, you will have an executable
called `biohh1`. At this point you can adjust the parameter guesses and
bounds by editing `specs.txt`. This file is also commented. Changes to
this file do not require a recompilation of the executable. Make sure it
specifies the location of the input data files correctly, or you'll get
a segmentation fault. You may also want to edit `biohh1.opt` to set
tolerances and maximum number of iterations. Consult the IPOPT
documentation for more information on the parameters in this file. To
run the estimation:

``` {.example}
./biohh1
```

This takes a long time, too. You will see some nice status updates as
the estimates converge. At the end there will be two files you can
inspect. `param.dat` has the parameter estimates; `data.dat` has the
state variables.

## Run forward predictions

The completed model can be used to predict forward in time by
integrating with the estimated state and parameter values. The example
`specs.txt` file uses 1500 ms for the assimilation, leaving quite a bit
of data for prediction. Chris Knowlton's `ipopt_predict` script uses the
equations.txt file and the output from IPOPT to generate forward
predictions. To run the prediction for the rest of the data in the
example:

``` {.example}
python model/ipopt_predict.py -e equations.txt pred.dat 258334
```

The last commandline parameter specifies the number of time points to
integrate forward after the end of the assimilation period.

## Editing the model

To change the equations of motion, edit `equations.txt`. This file is
commented, but the format is very sensitive to bookkeeping errors. It is
critical that the `nY,NP,nU,nI,nF` line match the number of equations,
parameters, control terms, forcing terms, and adjust functions
**exactly**. To change parameter bounds or input file (i.e., the forcing
current and voltage observations), edit `specs.txt`. It is also crucial
that this file match the variable and parameter counts exactly. Any
errors may result in segmentation faults. Sorry.

## Additional data

The `data` subdirectory contains data for all three exemplar neurons
from the manuscript. Each data collection epoch corresponds to two
files, one ending in `v.dat` and the other ending in `i.dat`. The former
file contains a single column of floating point numbers in ASCII format,
which is the recorded voltage in mV, starting at time 0 and sampled at
the rate specified in `specs.txt` (50 kHz). The latter file contains the
injected current (in pA) sampled at the same times.

Three recording epochs are provided for each exemplar neuron. The epoch
is given by the last number in the file name. The files used for data
assimilation and prediction in the manuscript are as follows:

|neuron |   assimilation  |   cross-epoch prediction|
|-------|----------------:|------------------------:|
  N1  |     `20120406_1_3_19` | `20120406_1_3_15`
  N2  |     `20120424_1_1_21` |  `20120424_1_1_15`
  N3  |     `20120116_1_1_6`  |  `20120116_1_1_2`

Currently the `specs.txt` file is set to run the analysis for N1. To
analyze a different neuron, edit `specs.txt` and replace the names of
the data files. Because the file names for the output of the `biohh1`
analysis are fixed, you may wish to copy `specs.txt`, `equations.txt`,
the data files, and the executable to a new directory if you want to run
a different analysis.

To generate cross-epoch predictions, you will need to take the parameter
estimates from the estimation epochs and create a new `specs.txt` with
those values fixed. Then run an estimation procedure on the first 100 ms
or so of data from the new epoch to get updated state estimates for the
new epoch, then predict forward as above.

The full dataset used in the manuscript is too large to include in this
archive, but we are happy to make it available on request. Contact Dan
Meliza (now in the University of Virginia Department of Psychology).
