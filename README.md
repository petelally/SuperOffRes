# SuperOffRes
Super-resolution reconstruction with unbalanced SSFP

## Content
This repository contains a simple multifrequency reconstruction algorithm, as outlined in the paper: Lally PJ, Matthews PM, Bangerter NK. Unbalanced SSFP for super-resolution 
in MRI. Magn Reson Med. 2020;00:1–13. https://doi.org/10.1002/mrm.28593


## Structure
* '/src': source code (written for MATLAB 2017a)
  * '/src/external': external source code, see Acknowledgements
* '/data': example 2D dataset from a Lego brick on a 9.4T Bruker Biospec 94/20 (MATLAB .mat file)
* '/doc': additional documentation 

## Example data
The sequence parameters were as follows: TR/TE=5/2.5ms; N<sub>PE</sub> x N<sub>FE</sub> = 32x128; FOV = 58x40mm<sup>2</sup>; slice thickness = 1.2mm; α=0.4°; 1000 dummy TRs (5s of steady-state preparation); 36 separate images at equidistant phase increments (i.e. 0°,10°,20°,...,350°), with the unbalanced gradient along the PE direction. The complex data array contained in the .mat file therefore has dimensions 32x128x36

## Acknowledgements
The ifft2c and fft2c functions were written by Michael Lustig for ESPIRiT, downloaded from http://people.eecs.berkeley.edu/~mlustig/Software.html on 10/02/2020. All use and distribution rights are as described in the original code.
