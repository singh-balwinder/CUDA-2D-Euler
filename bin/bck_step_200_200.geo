//+
Point(1) = {0.05, 0, 0, 1};
//+
Point(2) = {1, 0, 0, 1};
//+
Point(3) = {1, 0.6, 0, 1};
//+
Point(4) = {1, 1, 0, 1};
//+
Point(5) = {.05, 1, 0, 1};
//+
Point(6) = {0, 1, 0, 1};
//+
Point(7) = {0, 0.6, 0, 1};
//+
Point(8) = {0.05, 0.6, 0,1};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 5};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 7};
//+
Line(7) = {7, 8};
//+
Line(8) = {8, 1};
//+
Line(9) = {8, 3};
//+
Line(10) = {8, 5};
//+
Transfinite Curve {1, 9, 4} = 191 Using Progression 1;
//+
Transfinite Curve {2, 8} = 121 Using Progression 1;
//+
Transfinite Curve {3, 10, 6} = 81 Using Progression 1;
//+
Transfinite Curve {5, 7} = 11 Using Progression 1;
//+
Curve Loop(1) = {1, 2, -9, 8};
//+
Plane Surface(1) = {1};
//+
Recombine Surface {1};
//+
Transfinite Surface {1};
//+
Curve Loop(2) = {9, 3, 4, -10};
//+
Plane Surface(2) = {2};
//+
Transfinite Surface {2};
//+
Curve Loop(3) = {6, 7, 10, 5};
//+
Plane Surface(3) = {3};
//+
Recombine Surface {3};
//+
Transfinite Surface {3};
//+
Recombine Surface {1,2,3};
//+
Physical Curve("INLET", 7) = {6};
//+
Physical Curve("WALL", 5) = {7, 8, 1, 2, 3, 5, 4};
