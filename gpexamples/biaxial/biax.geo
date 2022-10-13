
SetFactory("OpenCASCADE");

alpha = 0.2;
len = 2.;
rad = 0.2;
trans = 0.8;
thick2 = 0.2;
len2 = 3.;
thick = 1.;
h = 0.2;

Mesh.CharacteristicLengthMin = h;
Mesh.CharacteristicLengthMax = h;

Point(1) = {0, 0, 0};
Point(2) = {rad*Sin(alpha), -rad*Cos(alpha), 0};
Point(3) = {rad*Sin(alpha), rad*Cos(alpha), 0};
Point(4) = {-len, rad + rad*Sin(alpha) + len*Sin(alpha), 0};
Point(5) = {-len, -rad - rad*Sin(alpha) - len*Sin(alpha) , 0};

Circle(1) = {2, 1, 3};

Line(2) = {3, 4};
Line(3) = {4, 5};
Line(4) = {5, 2};

Curve Loop(1) = {1, 2, 3, 4};

Plane Surface(1) = {1};

Rotate {{0, 0, 1}, {0, 0, 0}, -Pi/4} {
  Surface{1};
}

Translate {-trans, trans, 0} {
  Surface{1};
}

Point(6) = {0, 0, 0};
Point(7) = {0, len2, 0, 0};
Point(8) = {-thick, len2, 0};
Point(9) = {-thick, thick, 0};
Point(10) = {-len2, thick, 0};
Point(11) = {-len2, 0, 0};

Line(5) = {6, 7};
Line(6) = {7, 8};
Line(7) = {8, 9};
Line(8) = {9, 10};
Line(9) = {10, 11};
Line(10) = {11, 6};

Curve Loop(2) = {5, 6, 7, 8, 9, 10};

Plane Surface(2) = {2};

BooleanDifference{ Surface{2}; Delete; }{ Surface{1}; Delete; }

Physical Surface(2) = {2};
