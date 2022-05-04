//+
SetFactory("OpenCASCADE");
Rectangle(1) = {0., 0., 0, 2., 1., 0};
//+
Disk(2) = {1., 0., 0, 0.5, 0.5};
//+
Disk(3) = {0.4, 0.6, 0, 0.15, 0.15};
//+
Disk(4) = {1.6, 0.6, 0, 0.15, 0.15};
//+
BooleanDifference{ Surface{1}; Delete; }{ Surface{2,3,4}; Delete; }
//+
BooleanDifference{ Surface{1}; Delete; }{ Surface{2}; Delete; }