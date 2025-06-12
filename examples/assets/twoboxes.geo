// Gmsh project created on Fri Dec 13 11:47:12 2024
h = 0.2;
//+
Point(1) = {0, 0, 0, h};
//+
Point(3) = {1, 0, 0, h};
//+
Point(4) = {-1, 0, 0, h};
//+
Point(5) = {0, 0, 1, h};
//+
Point(6) = {1, 0, 1, h};
//+
Point(7) = {-1, 0, 1, h};
//+
Point(8) = {0, 1, 0, h};
//+
Point(9) = {1, 1, 0, h};
//+
Point(10) = {-1, 1, 0, h};
//+
Point(11) = {0, 1, 1, h};
//+
Point(12) = {1, 1, 1, h};
//+
Point(13) = {-1, 1, 1, h};
//+
Line(1) = {3, 9};
//+
Line(2) = {1, 8};
//+
Line(3) = {4, 10};
//+
Line(4) = {6, 12};
//+
Line(5) = {5, 11};
//+
Line(6) = {7, 13};
//+
Line(7) = {6, 5};
//+
Line(8) = {5, 7};
//+
Line(9) = {12, 11};
//+
Line(10) = {11, 13};
//+
Line(11) = {3, 1};
//+
Line(12) = {1, 4};
//+
Line(13) = {9, 8};
//+
Line(14) = {8, 10};
//+
Line(15) = {3, 6};
//+
Line(16) = {9, 12};
//+
Line(17) = {1, 5};
//+
Line(18) = {8, 11};
//+
Line(19) = {4, 7};
//+
Line(20) = {10, 13};
//+
Curve Loop(1) = {16, -4, -15, 1};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {18, -9, -16, 13};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {20, -10, -18, 14};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {6, -20, -3, 19};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {11, 2, -13, -1};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {12, 3, -14, -2};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {10, -6, -8, 5};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {9, -5, -7, 4};
//+
Plane Surface(8) = {8};
//+
Curve Loop(9) = {7, -17, -11, 15};
//+
Plane Surface(9) = {9};
//+
Curve Loop(10) = {8, -19, -12, 17};
//+
Plane Surface(10) = {10};
//+
Curve Loop(11) = {2, 18, -5, -17};
//+
Plane Surface(11) = {11};


//classify parts of geometry
Physical Surface("G12") = {11};
Physical Surface("G10") = {3,4,6,7,10};
Physical Surface("G20") = {1,2,5,8,9};
