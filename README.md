Marine Sampling Simulation
================

<style>
body {text-align: justify}
</style>

<br>

## Introduction

We provide here the code and original dataset to replicate the results
of our study showing the effects of sampling bias on marine
macroecological patterns. As we explain below, we developed a computer
simulation which create a spatial diversity pattern and apply upon it
the empirical spatial distribution of sampling effort. In this
simulation, our geographical domain is the Atlantic Ocean and our target
group, which provide the number of species and the distribution of
sampling effort, is the Ophiuroidea Class. However, such features can be
easily altered once the user changes the input information. The code of
the simulation is available in the ‘Simulation’ folder, the original
Ophiuroidea data with the empirical observations is available in the
‘Data/Basic’ folder, and the R scripts to analyze the simulation
results and compare it with the empirical observations can be found in
the ‘Analyzes’ folder (soon).

<br>

## The simulation

Our simulation (SampSim.exe; built in Delphi) creates a spatial
diversity pattern based on the mid domain effect (contiguous range) or
the area effect (not contiguous range) and apply upon it the empirical
latitudinal distribution of sampling effort of a specific taxonomic
group. The program will ask the user to enter the number of species to
be created and the number of simulations to run. The program also will
ask the directory path to call the input matrices. Two input matrices
are necessary: 1) a matrix of neighborhood among cells (NeigCellMat); 2)
a matrix with the cells and sampling effort of each latitude
(CellLatMat). The sampling effort of each latitude and number of species
within the geographical domain are based on empirical observation. The
code to recreate both matrices are available in the script ‘Input
matrices.R’.

The simulation starts by creating the greographical distribution of all
species. Once the known spatial pattern was created, we apply the
sampling effort based on empirical latitudinal distribution of sampling
events. At each sampling event, the simulation sample two random
species, and generate a database of observed richness. With this
database, which has the goal to emulate the databases in real world, we
can estimate the species richness, inventory completeness, and the
spatial gaps in species distribution. After run all simulations the
program will generate a matrix of konwn species richness by latitude and
cell (MatRichLat/MatRichCell), observed species richness after sampling
(MatObsLat/MatObsCell), estimated species richness after sampling
(MatEstLat/MatEstCell), and inventory completeness (MatCoverageLat) and
spatial gaps (MatGapsLat) by latitude. In these matrices, rows represent
latitudes or cells and columns represent the simulations.

We also provide a second simulation (SampSim\_vES.exe) where the
sampling effort is not based on empirical data, but is equal through the
latitudes. In this alternative simulation the program creates different
scenarios of sampling effort which increases at a given rate until a
maximum number of sampling events provided by the user. The output are
matrices where rows represent latitudes or cells and columns represent
the sampling scenarios. Therefore, the values represent the average of
all simulations (except for MatRichLat/MatRichCell).

<br>

## Data analyses

Once we know the “real” species richness distribution across the
latitudes, we can evaluate the effects of sampling bias on the detection
of the latitudinal diversity pattern. Moreover, we can compare the
simulated sampling results with the empirical ones. All the analyses can
be found in the R script ‘SimAnalyses.R’ (soon).

All the Delphi code and R scripts are commented so that you can follow
and reproduce each step. If something is not clear, please, do not
hesitate in contact me.
