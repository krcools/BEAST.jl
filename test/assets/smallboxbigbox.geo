h = 0.1;
r = 0.25;

Point(1) = {0, 0, 0, h};
Point(2) = {1, 0, 0, h};
Point(3) = {1, 1, 0, h};
Point(4) = {0, 1, 0, h};
Point(5) = {0, 0, 1, h};
Point(6) = {1, 0, 1, h};
Point(7) = {1, 1, 1, h};
Point(8) = {0, 1, 1, h};

Point(9) = {-r,    0, 0, h};
Point(10) = {    0, r, 0, h};
Point(11) = {-r, r, 0, h};


Point(12) = { 0, 0, r, h};
Point(13) = {-r, 0, r, h};
Point(14) = { 0, r, r, h};
Point(15) = {-r, r, r, h};


Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 10};
Line(5) = {10, 1};

Line(6) = {2, 6};
Line(7) = {3, 7};
Line(8) = {4, 8};
Line(9) = {1, 12};
Line(10) = {12, 5};

Line(11) = {5, 6};
Line(12) = {6, 7};
Line(13) = {7, 8};
Line(14) = {8, 5};

Line(15) = {10, 14};
Line(16) = {14, 12};

Line(17) = {1, 9};
Line(18) = {9, 11};
Line(19) = {11, 10};
Line(20) = {9, 13};
Line(21) = {11, 15};

Line(22) = {12, 13};
Line(23) = {13, 15};
Line(24) = {15, 14};


Curve Loop(1) = {1, 2, 3, 4, 5};
Plane Surface(1) = {1};
Curve Loop(2) = {1, 6, -11, -10, -9};
Plane Surface(2) = {2};
Curve Loop(3) = {-6, 2, 7, -12};
Plane Surface(3) = {3};
Curve Loop(4) = {-3, 7, 13, -8};
Plane Surface(4) = {4};
Curve Loop(5) = {-4, 8, 14, -10, -16, -15};
Plane Surface(5) = {5};
Curve Loop(6) = {11, 12, 13, 14};
Plane Surface(6) = {6};

Curve Loop(7) = {-5, 15, 16, -9};
Plane Surface(7) = {7};

Curve Loop(8) = {5, 19, 18, 17};
Plane Surface(8) = {8};
Curve Loop(9) = {-17, 9, 22, -20};
Plane Surface(9) = {9};
Curve Loop(10) = {18, 21, -23, -20};
Plane Surface(10) = {10};
Curve Loop(11) = {19, 15, -24, -21};
Plane Surface(11) = {11};
Curve Loop(12) = {-22, -16, -24, -23};
Plane Surface(12) = {12};

Physical Surface("Gamma1") = {1, 2,  3, 4, 5, 6};
Physical Surface("Gamma2") = {8, 9, 10, 11, 12};
Physical Surface("Gamma3") = {7};
