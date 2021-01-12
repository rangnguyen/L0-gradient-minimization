INTRODUCTION (MESH DENOISE)
------------
The folder contains two things as follows:
 * exe: contains the executable file for Windows OS (including Win32 only).
 * src: contains our source code implemented in C++. 

HOW TO RUN
------------
Usage is -i <infile> -t <lambda> -o <outfile>
Example: type as follows in command line 
Mesh_Denoising.exe -i pyramid_noise.ply -t 0.5 -o pyramid_denoise.ply

NOTE: can open source code using Visual Studio 2010

REFERENCE
------------
Rang M. H. Nguyen, Michael S. Brown
Fast and Effective L0 Gradient Minimization by Region Fusion
ICCV 2015