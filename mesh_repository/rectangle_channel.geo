lc = 1.0;
nb_pts_horizontal = 20;
nb_pts_vertical   = 8; 


Point(1) = { 0.0, 0.0, 0.0, lc };
Point(2) = { 4.0, 0.0, 0.0, lc };
Point(3) = { 4.0, 1.0, 0.0, lc };
Point(4) = { 0.0, 1.0, 0.0, lc };

Line(1) = { 1,2 };
Line(2) = { 2,3 };
Line(3) = { 3,4 };
Line(4) = { 4,1 };

Line Loop(1) = { 1,2,3,4 };

Plane Surface(1) = {1};

Transfinite Line(1) = nb_pts_horizontal;
Transfinite Line(3) = nb_pts_horizontal;
Transfinite Line(2) = nb_pts_vertical;
Transfinite Line(4) = nb_pts_vertical;

Physical Line("bottom_wall") = {1};
Physical Line("outlet")      = {2};
Physical Line("top_wall")    = {3};
Physical Line("inlet")       = {4};

Physical Surface("domain")   = {1};

