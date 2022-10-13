
SetFactory("OpenCASCADE");

l = 2;
mid = l/2;
h = 1;
r1 = 0.5;
r2 = 0.15;
hoff = 0.6;
voff = 0.6;
c1 = mid-hoff;
c2 = mid+hoff;

// Main rectangle
Point(1) = {mid, 0, 0};
Point(2) = {mid-r1, 0, 0};
Point(3) = {mid, r1, 0};
Point(4) = {mid+r1, 0, 0};
Point(5) = {l, 0, 0};
Point(6) = {l, h, 0};
Point(7) = {0, h, 0};
Point(8) = {0, 0, 0};

Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Line(3) = {4, 5};
Line(4) = {5, 6};
Line(5) = {6, 7};
Line(6) = {7, 8};
Line(7) = {8, 2};

Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7};

// Left side hole
Point(9) = {c1, voff, 0};
Point(10) = {c1+r2, voff, 0};
Point(11) = {c1, voff+r2, 0};
Point(12) = {c1-r2, voff, 0};
Point(13) = {c1, voff-r2, 0};

Circle(8) = {10, 9, 11};
Circle(9) = {11, 9, 12};
Circle(10) = {12, 9, 13};
Circle(11) = {13, 9, 10};

Curve Loop(2) = {8, 9, 10, 11};

// Right side hole
Point(14) = {c2, voff, 0};
Point(15) = {c2+r2, voff, 0};
Point(16) = {c2, voff+r2, 0};
Point(17) = {c2-r2, voff, 0};
Point(18) = {c2, voff-r2, 0};

Circle(12) = {15, 14, 16};
Circle(13) = {16, 14, 17};
Circle(14) = {17, 14, 18};
Circle(15) = {18, 14, 15};

Curve Loop(3) = {12, 13, 14, 15};

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};

BooleanDifference{ Surface{1}; Delete; }{ Surface{2,3}; Delete; }

Physical Surface(1) = {1};

Mesh.SecondOrderLinear = 1;
