Mesh.Algorithm=7;

lc = 1.0;

x_min = 0.0;
x_max =  1.0;
y_min = 0.0;
y_max = 1.0;

Point(1) = { x_min, y_min, 0.0, lc };
Point(2) = { x_max, y_min, 0.0, lc };
Point(3) = { x_max, y_max, 0.0, lc };
Point(4) = { x_min, y_max, 0.0, lc };

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line Loop(1) = {1,2,3,4};

Plane Surface(1) = {1};

Physical Surface("InnerCells") = {1};
Physical Line("Bottom") = {1};
Physical Line("Right") = {2};
Physical Line("Top") = {3};
Physical Line("Left") = {4};

Transfinite Line(1) = 8;
Transfinite Line(2) = 8;
Transfinite Line(3) = 8;
Transfinite Line(4) = 8;

/*
Transfinite Surface(1) = {1,2,3,4};
Recombine Surface { 1 };
*/

Field[1] = BoundaryLayer;
Field[1].EdgesList = { 1,2,3,4 };
Field[1].FanNodesList = { 1,2,3,4 };
Field[1].hfar = 0.1;
Field[1].hwall_n = 0.001;
Field[1].ratio = 1.3;
Field[1].thickness = 0.01;
Field[1].IntersectMetrics = 1;

Background Field = 1;
