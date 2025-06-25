lc = 0.050;
L = 2.0;
H = 2.0;

nb_pts_hor = 11;
nb_pts_vert = 12;


Point(0) = {0.0,0.0,0.0,lc};
Point(1) = { L, 0.0,0.0,lc};
Point(2) = { L ,H,  0.0,lc};
Point(3) = { 0 ,H,  0.0,lc};
Point(4) = {-L ,H  ,0.0,lc};
Point(5) = {-L ,0.0,0.0,lc};


Line(1) = {0,1};
Line(2) = {1,2};
Line(3) = {2,3};
Line(4) = {3,4};
Line(5) = {4,5};
Line(6) = {5,0};

Line Loop(1) = {1,2,3,4,5,6};

Plane Surface(1) = {1};

Transfinite Line(1) = nb_pts_hor;
Transfinite Line(3) = nb_pts_hor;
Transfinite Line(4) = nb_pts_hor;
Transfinite Line(6) = nb_pts_hor;
Transfinite Line(2) = nb_pts_vert;
Transfinite Line(5) = nb_pts_vert;

Physical Line("inlet") = {6}; 
Physical Line("outlet") = {1}; 
Physical Line("farfield") = {3,5}; 
Physical Surface("fluid") = {1};
