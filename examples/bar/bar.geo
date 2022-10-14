// This variable is a characteristic length and it controles the mesh size around a Point.
// It is possible to specify more than one variable for this purpose.
cl = 0.20;

// These variables control the dimensions
L = 10.;

// Points contains the x, y and z coordinate and the characteristic length of the Point.
Point(1) = {0,0,0,cl};
Point(2) = {L,0,0,cl};

// A Line is basically a connection between two Points. 
Line(1) = {1,2};
Physical Curve(1) = {1};

Mesh.ElementOrder = 2;
