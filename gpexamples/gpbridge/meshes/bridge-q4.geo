// This variable is a characteristic length and it controles the mesh size around a Point.
// It is possible to specify more than one variable for this purpose.
cl = 0.5;

// These variables control the dimensions
H = 0.5;
L1 = 3.;
L2 = 4.;
P0 = 0.;
P1 = L1;
P2 = L1+L2;
P3 = L1+L2+L1;

// Points contains the x, y and z coordinate and the characteristic length of the Point.
Point(1) = {P0,0,0,cl};
Point(2) = {P1,0,0,cl};
Point(3) = {P2,0,0,cl};
Point(4) = {P3,0,0,cl};
Point(5) = {P3,H,0,cl};
Point(6) = {P2,H,0,cl};
Point(7) = {P1,H,0,cl};
Point(8) = {P0,H,0,cl};

// A Line is basically a connection between two Points. A good practice is to connect the
// Points in a counter-clockwise fashion.
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,1};
Line(9) = {2,7};
Line(10)= {3,6};

// A Line Loop is the connection of Lines that defines an area. Again it is good practice
// to do this in a counter-clockwise fashion.
Line Loop(1) = {1,9,7,8};
Line Loop(2) = {2,10,6,-9};
Line Loop(3) = {3,4,5,-10};

// From the Line Loop it is now possible to create a surface, in this case a Plane Surface.
Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};

// From the Plane Surface a Physical Surface is generated, this makes is possible to only
// save elements which are defined on the area specified by the Line Loop.
Physical Surface(1) = {1,2,3};

// Now, we create Transfinite Surfaces on each Line Loop
Transfinite Surface {1};
Transfinite Surface {2};
Transfinite Surface {3};

// In order to define the density, we use transfinite curves
Transfinite Curve {1,3,5,7}  = L1/cl+1 Using Progression 1;
Transfinite Curve {2,6}      = L2/cl+1 Using Progression 1;
Transfinite Curve {4,8,9,10} = H/cl+1 Using Progression 1;

// Recombine Surface in order to turn the triangular grid into quadrilaterals
Recombine Surface {1,2,3};

