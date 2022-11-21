// This variable is a characteristic length and it controles the mesh size around a Point.
// It is possible to specify more than one variable for this purpose.
cl = 0.20;

// These variables control the dimensions
H = 2.;
L = 10.;

// Points contains the x, y and z coordinate and the characteristic length of the Point.
Point(1) = {0,0,0,cl};
Point(2) = {0.5*L,0,0,cl};
Point(3) = {L,0,0,cl};
Point(4) = {L,H,0,cl};
Point(5) = {0.5*L,H,0,cl};
Point(6) = {0,H,0,cl};

// A Line is basically a connection between two Points. A good practice is to connect the
// Points in a counter-clockwise fashion.
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,1};
Line(7) = {2,5};

// A Line Loop is the connection of Lines that defines an area. Again it is good practice
// to do this in a counter-clockwise fashion.
Line Loop(1) = {1,7,5,6};
Line Loop(2) = {2,3,4,-7};

// From the Line Loop it is now possible to create a surface, in this case a Plane Surface.
Plane Surface(1) = {1};
Plane Surface(2) = {2};

// From the Plane Surface a Physical Surface is generated, this makes is possible to only
// save elements which are defined on the area specified by the Line Loop.
Physical Surface('left') = {1};
Physical Surface('right') = {2};
