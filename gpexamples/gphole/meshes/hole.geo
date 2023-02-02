// Dimensions of the cantilever itself
H = 2.;
L = 10.;

// Dimensions of the hole
l = 1.;
h = 1.;
left = 8.5;
bot = (H-h)/2;

// Length scale
cl = 0.5;

// Points contain the x, y and z coordinates
Point(1)  = {0, 0, 0};
Point(2)  = {L, 0, 0};
Point(3)  = {L, H, 0};
Point(4)  = {0, H, 0};

Point(11) = {left, bot, 0};
Point(12) = {left+l, bot, 0};
Point(13) = {left+l, bot+h, 0};
Point(14) = {left, bot+h, 0};

Point(21) = {left, 0, 0};
Point(22) = {left+l, 0, 0};
Point(23) = {L, bot, 0};
Point(24) = {L, bot+h, 0};
Point(25) = {left+l, H, 0};
Point(26) = {left, H, 0};
Point(27) = {0, bot+h, 0};
Point(28) = {0, bot, 0};

// A Line is basically a connection between two Points
// A good practice is to connect the Points in a counter-clockwise fasion
Line(1) = {1, 21};
Line(2) = {21, 22};
Line(3) = {22, 2};
Line(4) = {2, 23};
Line(5) = {23, 24};
Line(6) = {24, 3};
Line(7) = {3, 25};
Line(8) = {25, 26};
Line(9) = {26, 4};
Line(10) = {4, 27};
Line(11) = {27, 28};
Line(12) = {28, 1};

Line(21) = {11, 12};
Line(22) = {12, 13};
Line(23) = {13, 14};
Line(24) = {14, 11};

Line(31) = {21, 11};
Line(32) = {22, 12};
Line(33) = {23, 12};
Line(34) = {24, 13};
Line(35) = {25, 13};
Line(36) = {26, 14};
Line(37) = {27, 14};
Line(38) = {28, 11};

// A Line Loop is the connection of Lines that defines an area.
// Again, it is good practice to do this in a counter-clockwise fashion
Line Loop(1) = {1, 31, -38, 12};
Line Loop(2) = {2, 32, -21, -31};
Line Loop(3) = {3, 4, 33, -32};
Line Loop(4) = {5, 34, -22, -33};
Line Loop(5) = {6, 7, 35, -34};
Line Loop(6) = {8, 36, -23, -35};
Line Loop(7) = {9, 10, 37, -36};
Line Loop(8) = {11, 38, -24, -37};

// From the Line Loop it is now possible to cretae a surface, in this case a Plane Surface.
Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};
Plane Surface(6) = {6};
Plane Surface(7) = {7};
Plane Surface(8) = {8};

// From the Plane Surface a Physical Surface is generated.
// This makes it possible to only save elements which are defined on the area specified by the Line Loop.
Physical Surface(1) = {1, 2, 3, 4, 5, 6, 7, 8};

// To create a structured mesh, we need to define a Transfinite Curve on each Line
Transfinite Curve {1, 38, 37, 9} = left/cl+1 Using Progression 1;
Transfinite Curve {2, 21, 23, 8} = l/cl+1 Using Progression 1;
Transfinite Curve {3, 33, 34, 7} = (L-left-l)/cl+1 Using Progression 1;
Transfinite Curve {4, 32, 31, 12} = bot/cl+1 Using Progression 1;
Transfinite Curve {5, 22, 24, 11} = h/cl+1 Using Progression 1;
Transfinite Curve {6, 35, 36, 10} = (H-bot-h)/cl+1 Using Progression 1;

// Now, we create Transfinite Surfaces on each Line Loop
Transfinite Surface {1};
Transfinite Surface {2};
Transfinite Surface {3};
Transfinite Surface {4};
Transfinite Surface {5};
Transfinite Surface {6};
Transfinite Surface {7};
Transfinite Surface {8};

// Recombine Surface in order to turn the triangular grid into quadrilaterals
Recombine Surface {1, 2, 3, 4, 5, 6, 7, 8};
