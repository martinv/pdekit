nb_pts_on_edge = 16;
lc = 1.0/(nb_pts_on_edge-1);

Point(1) = {  0.0,  0.0, 0.0, lc };
Point(2) = {  1.0,  0.0, 0.0, lc };
Point(3) = {  0.0,  1.0, 0.0, lc };
Point(4) = { -1.0,  0.0, 0.0, lc };
Point(5) = {  0.0, -1.0, 0.0, lc };

Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};

Line Loop(1)= {1,2,3,4};

Plane Surface(1) = {1};

Physical Line("boundary") = {1,2,3,4};

Physical Surface("interior") = {1};

Transfinite Line(1) = nb_pts_on_edge;
Transfinite Line(2) = nb_pts_on_edge;
Transfinite Line(3) = nb_pts_on_edge;
Transfinite Line(4) = nb_pts_on_edge;

