// This variable is a characteristic length and it controles the mesh size around a Point.
// It is possible to specify more than one variable for this purpose.
cc = 0.5;

// These variables control the dimensions
L = 4.0;
H = 2.0;
R = 0.8;

// Points contains the x, y and z coordinate and the characteristic length of the Point.
Point(1) = {0,0,0,cc};
Point(2) = {L/2,0,0,cc};
Point(3) = {L,0,0,cc};
Point(4) = {L,H/2,0,cc};
Point(5) = {L,H,0,cc};
Point(6) = {L/2,H,0,cc};
Point(7) = {0,H,0,cc};
Point(8) = {0,H/2,0,cc};

Point(10) = {L/2,H/2,0,cc};
Point(11) = {L/2,H/2-R,0,cc};
Point(12) = {L/2+R,H/2,0,cc};
Point(13) = {L/2,H/2+R,0,cc};
Point(14) = {L/2-R,H/2,0,cc};

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

Line(9) = {4,12};
Line(10) = {8,14};

Circle(11) = {11,10,12};
Circle(12) = {12,10,13};
Circle(13) = {13,10,14};
Circle(14) = {14,10,11};

// A Line Loop is the connection of Lines that defines an area. Again it is good practice
// to do this in a counter-clockwise fashion.
Line Loop(1) = {1,2,3,9,-11,-14,-10,8};
Line Loop(2) = {5,6,7,10,-13,-12,-9,4};

// From the Line Loop it is now possible to create a surface, in this case a Plane Surface.
Plane Surface(1) = {1};
Plane Surface(2) = {2};

// From the Plane Surface a Physical Surface is generated, this makes is possible to only
// save elements which are defined on the area specified by the Line Loop.
Physical Surface(1) = {1};
Physical Surface(2) = {2};
