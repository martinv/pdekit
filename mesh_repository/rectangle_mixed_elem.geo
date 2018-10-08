nb_pts_per_edge = 5;

L = 1.0;

lc = L/nb_pts_per_edge;

Point(1) = { 0.0 ,0.0,0.0,lc};
Point(2) = { L   ,0.0,0.0,lc};
Point(3) = {2.0*L,0.0,0.0,lc};
Point(4) = {2.0*L, L ,0.0,lc};
Point(5) = {  L  , L ,0.0,lc};
Point(6) = { 0.0 , L ,0.0,lc};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,1};
Line(7) = {2,5};

Line Loop(1) = {1,7,5,6};
Line Loop(2) = {2,3,4,-7};

Plane Surface(1) = {1};
Plane Surface(2) = {2};

Transfinite Line(2) = nb_pts_per_edge;
Transfinite Line(3) = nb_pts_per_edge;
Transfinite Line(4) = nb_pts_per_edge;
Transfinite Line(7) = nb_pts_per_edge;

Transfinite Surface(2) = {2,3,4,5};
Recombine Surface(2);

Physical Line("Bottom") = {1,2};
Physical Line("Right")  = {3};
Physical Line("Top")    = {4,5};
Physical Line("Left")   = {6};

Physical Surface("Tri")  = {1};
Physical Surface("Quad") = {2};

