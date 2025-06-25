nb_pts_on_edge = 16;
lc = 1.0/(nb_pts_on_edge-1);

Point(1) = { 0.0, 0.0, 0.0, lc };
Point(2) = { 1.0, 0.0, 0.0, lc };
Point(3) = { 1.0, 1.0, 0.0, lc };
Point(4) = { 0.0, 1.0, 0.0, lc };

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line Loop(1)= {1,2,3,4};

Plane Surface(1) = {1};

Physical Line("left") = {4};
Physical Line("bottom") = {1};
Physical Line("right") = {2};
Physical Line("top") = {3};

Physical Surface("interior") = {1};

Transfinite Line(1) = nb_pts_on_edge;
Transfinite Line(2) = nb_pts_on_edge;
Transfinite Line(3) = nb_pts_on_edge;
Transfinite Line(4) = nb_pts_on_edge;

/*
Transfinite Surface(1) = {1,2,3,4};
Recombine Surface{1};
*/

