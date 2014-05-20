THIS DIRECTORY PROVIDES A WORKING ALGORITHM FOR PBC SPIN CHAINS, 
WITH one central CLUSTER and the rest is the environment;
the central bulk cluster rotates clockwise (more stable than anti-clockwise, 
probably because of left-canonical form combined with left matrix-vector 
multiplication in the approximation algorithm)

randenv2 function provides stable results with pn = 2, and ph = 4;

randenv function provides stable results with pn = 2, and ph = 4, but is 
4 times faster than randenv2;


