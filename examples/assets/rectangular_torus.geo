h = 1.0;
H = 1.0;
// Gmsh project created on Sun Nov 06 14:23:16 2022
SetFactory("Built-in");
//+
Point(1) = {2, -2, 0, h};
//+
Point(2) = {2, 2, 0, h};
//+
Point(3) = {-2, 2, 0, h};
//+
Point(4) = {-2, -2, 0, h};
//+
Point(5) = {1, -1, 0, h};
//+
Point(6) = {1, 1, 0, h};
//+
Point(7) = {-1, 1, 0, h};
//+
Point(8) = {-1, -1, 0, h};
//+
Line(1) = {2, 3};
//+
Line(2) = {3, 4};
//+
Line(3) = {4, 1};
//+
Line(4) = {1, 2};
//+
Line(5) = {6, 7};
//+
Line(6) = {7, 8};
//+
Line(7) = {8, 5};
//+
Line(8) = {5, 6};
//+
Point(9) = {2, -2, H, h};
//+
Point(10) = {2, 2, H, h};
//+
Point(11) = {-2, 2, H, h};
//+
Point(12) = {-2, -2, H, h};
//+
Point(13) = {1, -1, H, h};
//+
Point(14) = {1, 1, H, h};
//+
Point(15) = {-1, 1, H, h};
//+
Point(16) = {-1, -1, H, h};
//+
Line(9) = {10, 11};
//+
Line(10) = {11, 12};
//+
Line(11) = {12, 9};
//+
Line(12) = {9, 10};
//+
Line(13) = {14, 15};
//+
Line(14) = {15, 16};
//+
Line(15) = {16, 13};
//+
Line(16) = {13, 14};
//+
Line(17) = {2, 10};
//+
Line(18) = {3, 11};
//+
Line(19) = {4, 12};
//+
Line(20) = {1, 9};
//+
Line(21) = {6, 14};
//+
Line(22) = {7, 15};
//+
Line(23) = {8, 16};
//+
Line(24) = {5, 13};
//+
Curve Loop(1) = {-12, 17, 4, -20};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {1, 18, -9, -17};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {-18, -10, 19, 2};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {-19, -11, 20, 3};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {24, 16, -21, -8};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {-5, -22, 13, 21};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {-6, -23, 14, 22};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {23, 15, -24, -7};
//+
Plane Surface(8) = {8};
//+
Curve Loop(9) = {10, 11, 12, 9};
//+
Curve Loop(10) = {14, 15, 16, 13};
//+
Plane Surface(9) = {9, 10};
//+
Curve Loop(11) = {-3, -4, -1, -2};
//+
Curve Loop(12) = {-8, -5, -6, -7};
//+
Plane Surface(10) = {11, 12};
