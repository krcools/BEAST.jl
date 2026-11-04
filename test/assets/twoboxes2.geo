h = 0.1;
//+
Point(1) = {0, 0, 0, h};
//+
Point(2) = {0.5, 0, 0, h};
//+
Point(3) = {1, 0, 0, h};
//+
Point(4) = {0, 1, 0, h};
//+
Point(5) = {0.5, 1, 0, h};
//+
Point(6) = {1, 1, 0, h};
//+
Point(7) = {0, 0, 1, h};
//+
Point(8) = {0.5, 0, 1, h};
//+
Point(9) = {1, 0, 1, h};
//+
Point(10) = {0, 1, 1, h};
//+
Point(11) = {0.5, 1, 1, h};
//+
Point(12) = {1, 1, 1, h};
//+
Line(1) = {1, 4};
//+
Line(2) = {2, 5};
//+
Line(3) = {3, 6};
//+
Line(4) = {7, 10};
//+
Line(5) = {8, 11};
//+
Line(6) = {9, 12};
//+
Line(7) = {7, 8};
//+
Line(8) = {8, 9};
//+
Line(9) = {10, 11};
//+
Line(10) = {11, 12};
//+
Line(11) = {1, 2};
//+
Line(12) = {2, 3};
//+
Line(13) = {4, 5};
//+
Line(14) = {5, 6};
//+
Line(15) = {9, 3};
//+
Line(16) = {8, 2};
//+
Line(17) = {7, 1};
//+
Line(18) = {10, 4};
//+
Line(19) = {11, 5};
//+
Line(20) = {12, 6};
//+
Curve Loop(1) = {7, 5, -9, -4};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {8, 6, -10, -5};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {11, 2, -13, -1};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {12, 3, -14, -2};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {7, 16, -11, -17};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {8, 15, -12, -16};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {10, 20, -14, -19};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {9, 19, -13, -18};
//+
Plane Surface(8) = {8};
//+
Curve Loop(9) = {4, 18, -1, -17};
//+
Plane Surface(9) = {9};
//+
Curve Loop(10) = {5, 19, -2, -16};
//+
Plane Surface(10) = {10};
//+
Curve Loop(11) = {6, 20, -3, -15};
//+
Plane Surface(11) = {11};
//+
Physical Surface("Gamma1") = {5, 1, 9, 8, 3};
//+
Physical Surface("Gamma2") = {6, 2, 7, 4, 11};
//+
Physical Surface("Gamma3") = {10};
