# Improved micro-continuum simulator for capillary-dominated multiphase flow

## Introductions
A diverse range of multiphase flow and transport occurs in multiscale porous media. The multiphase micro-continuum Darcy–Brinkmann–Stokes (DBS) model has been developed to simulate the multiphase flow at both the pore and continuum scales via single-field equations. However, the unacceptable spurious velocities produced by the conventional micro-continuum DBS model present challenges to the modeling of capillary-dominated flow dynamics. This study improves the micro-continuum DBS model to mitigate these spurious velocities at the gas–liquid interface and contact-line regions.

The solver was developed based on the open-source [hybridPorousInterFoam solver](https://github.com/Franjcf/hybridPorousInterFoam.git) which is the first implementation of the multiphase micro-continuum model with the FVM based on OpenFOAM7 and later.

This repository was created by Zhiying Liu and Qianghui XU. 

## Installations 
1. install OpenFOAM 7 as the official [guidance](https://openfoam.org/download/7-source/)

2. go to this repository directory and run `./Allwmake`

## Runing tutorials 
go to each case folder in tutorials and run `./run.sh`

## Key references 
1. F. J. Carrillo, I. C. Bourg, and C. Soulaine, "Multiphase flow modeling in multiscale porous media: An open-source micro-continuum approach," Journal of Computational Physics: X 8, 100073 (2020).
