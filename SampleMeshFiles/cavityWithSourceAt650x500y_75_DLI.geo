// Gmsh project created on Mon Dec 12 19:12:37 2016
//+

lc = 0.0;

Point(1) = {0, 0, 0};
//+
Point(2) = {0.25, 0, 0};
//+
Point(3) = {0.25, 0.2, 0};
//+
Point(4) = {0, 0.2, 0};
//+
Point(5) = {0.0625+.01, .05+.01, 0, lc};
//+
Point(6) = {0.0625+.01, .05-.01, 0, lc};
//+
Point(7) = {0.0625-.01, .05-.01, 0, lc};
//+
Point(8) = {0.0625-.01, .05+.01, 0, lc};
//+
Line(1) = {7, 8};
//+
Line(2) = {8, 5};
//+
Line(3) = {5, 6};
//+
Line(4) = {6, 7};
//+
Line Loop(5) = {2, 3, 4, 1};
//+
Plane Surface(6) = {5};
//+
Extrude {0, 0, .0075} {
  Surface{6};
}
//+
Line(29) = {2, 1};
//+
Line(30) = {1, 4};
//+
Line(31) = {4, 3};
//+
Line(32) = {3, 2};
//+
Line Loop(33) = {32, 29, 30, 31};
//+
Plane Surface(34) = {5, 33};
//+
Extrude {0, 0, .0075} {
  Surface{34};
}
//+
//Physical Volume("TryThis") = {2, 1};
