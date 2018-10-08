lc = 0.1;

Ly = 5.0;

Point(1) = {  0.0,  0.0, 0.0, lc };
Point(2) = {  0.0,  1.0, 0.0, lc };
Point(3) = { -1.0,  0.0, 0.0, lc };
Point(4) = {  0.0, -1.0, 0.0, lc };
Point(5) = {  0.0,   Ly, 0.0, lc };
Point(6) = { -2.0,   Ly, 0.0, lc };
Point(7) = { -2.0,  0.0, 0.0, lc };
Point(8) = { -2.0,  -Ly, 0.0, lc };
Point(9) = {  0.0,  -Ly, 0.0, lc };



Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Line(3) = {2,5};
Line(4) = {5,6};
Line(5) = {6,7};
Line(6) = {7,8};
Line(7) = {8,9};
Line(8) = {9,4};

Line Loop(1) = { -1, 3:8, -2 };

Plane Surface(1) = { 1 };


Physical Line("Wall") = {1,2};
Physical Line("Inlet") = {5,6};

Physical Line("OutletTop") = {3,4};
Physical Line("OutletBottom") = {7,8};

Physical Surface("InnerCells") = { 1 };
