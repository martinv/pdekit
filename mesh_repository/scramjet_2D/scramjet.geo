lc = 0.1;

Point(1) = {0.0, 0.0, 0.0, lc};
Point(2) = {8.0, 0.0, 0.0, lc};
Point(3) = {8.0, 0.8, 0.0, lc};
Point(4) = {0.0, 2.0, 0.0, lc};

Point(5) = {2.0, 0.7, 0.0, lc};
Point(6) = {4.0, 0.2, 0.0, lc};
Point(7) = {7.0, 0.6, 0.0, lc};
Point(8) = {6.0, 0.7, 0.0, lc};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,5};

Line Loop(1) = {1:4};
Line Loop(2) = {5:8};

Plane Surface(1) = {1,2};

Physical Line("BottomWall") = {1};
Physical Line("TopWall") = {3};
Physical Line("InteriorWall") = {5,6,7,8};
Physical Line("Inlet") = {4};
Physical Line("Outlet") = {2};
Physical Surface("Interior") = {1};
