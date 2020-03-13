# MonteCarloProteinVolumeApproximation

Nicholas Treynor, Qian Li
ECS 129 

This program is meant to be used to calculate the the volume, surface area, and standard deviation of both surface area and volume of
a protein structure. 

The protein structure in question will be modeled by a union of balls representing atoms.

The program is written in C++ and takes a file name and two integers as command line parameters. 

The input file should store a number denoting the number of balls it will contain on the first line, and the following lines should 
each consist of four floating point numbers separated by whitespace that will denote the X, Y, and Z position of a given ball, 
as well as its radius.

The two integers passed in as input denote the number of points that will be tested for the estimation of volume of the 
protein via the Monte Carlo method, and the number of points that will be utilized in the estimation of the surface area of 
the protein, respectively. 

These command line arguments should be passed in via the form: string ‘fileName,’ int PointsToTest, int PointsToTestOnSurface. 
For example, the input “Proteins.txt 10,000 30,000” would parse the file Proteins.txt for the number of spheres to utilize in our model 
and their values, and then test those spheres with 10,000 points to approximate volume, and then generate 30,000 points on the surface of
those spheres for the sake of approximating the surface area of the protein. 

Proteins.txt has been included as an example file, but can easily be modified to produce different results. For example, changing the
value "2775" on line 1 to "500" will result in the program calculating the volume, surface area, etc of only the first 500 points. 

