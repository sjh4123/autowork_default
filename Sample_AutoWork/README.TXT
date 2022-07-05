## DESCRIPTION ##
--------------------------------------------------------------------------------
AUTOWORK is a software designed to integrate all components necessary to perform
automated molecular simulations and trajectory analyses;
PreSQ is an algorithm embedded in AUTOWORK for prediction of quench temperatures
in automated glass formation simulations.


## HOW TO EXECUTE AUTOWORK PROGRAM ##
--------------------------------------------------------------------------------
./AUTOWORK [-i or -inp] <input file name>

Ex. ./AUTOWORK -i initial.inp


## SOFTWARE DOWNLOAD ##
--------------------------------------------------------------------------------
Documentation, examples, and supporting code are available at:
http://blogs.uakron.edu/simmonsgroup


## COMPILING SOURCE CODES ##
--------------------------------------------------------------------------------
Makefile is prepared ready in the same folder; simply type 'make' in a Linux
terminal, navigated to the folder that contains the source code and the makefile


## QUICKSTART EXAMPLE ##
--------------------------------------------------------------------------------
please see initial.inp and sample scripts


## REQUIREMENTS ##
--------------------------------------------------------------------------------
AUTOWORK requires Linux Bash and a C++ compiler that supports C++11

To run parallel simulations, it needs to be run in a batch-queuing scheduler
system. AUTOWORK currently uses SGE (Sun Grid Engine) as the default and only
job scheduling system (will possibly add PBS and others) and LAMMPS (26Jan2014)
as the default MD simulation package

Trajectory analyses are performed by AMDAT


## LICENSE ##
--------------------------------------------------------------------------------
PreSQ & AUTOWORK is available under the terms of the open-source 3-clause BSD
license (See LICENSE.TXT)