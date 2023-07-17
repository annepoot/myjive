// This variable is a characteristic length and it controles the mesh size around a Point.
// It is possible to specify more than one variable for this purpose.
cl = 2.5;

// These variables control the dimensions
L = 5.;

// Points contains the x, y and z coordinate and the characteristic length of the Point.
Point(1) = {0,0,0,cl};
Point(2) = {L,0,0,cl};
Point(3) = {2*L,0,0,cl};
Point(4) = {2*L,L,0,cl};
Point(5) = {2*L,2*L,0,cl};
Point(6) = {L,2*L,0,cl};
Point(7) = {L,L,0,cl};
Point(8) = {0,L,0,cl};

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
Line(10)= {4,7};

// A Line Loop is the connection of Lines that defines an area. Again it is good practice
// to do this in a counter-clockwise fashion.
Line Loop(1) = {1,9,7,8};
Line Loop(2) = {2,3,10,-9};
Line Loop(3) = {4,5,6,-10};

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

// Recombine Surface in order to turn the triangular grid into quadrilaterals
Recombine Surface {1,2,3};
