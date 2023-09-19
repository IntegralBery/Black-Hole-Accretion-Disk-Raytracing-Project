# Black-Hole-Accretion-Disk-Raytracing-Project
I did this project for a computational methods in physics class in my second year! It raytraces and renders a black hole with its accretion disk which temperature you can specify. Implemented in MATLAB.

The parameters used to generate the image (mass of black hole, radius and temperature of accretion disk, dimensions of image etc.) can be changed at the top of img_gen.m. img_gen.m generates an array (stored in default output.mat) representing the raytraced image before color is added.

img_colorize.m can then be used to render the actual image from output.mat.

Details of the implementation can be found in my project report in proj_report.pdf.
