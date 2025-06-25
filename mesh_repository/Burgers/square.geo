// 2D mesh algorithm (1=MeshAdapt, 2=Automatic, 5=Delaunay, 6=Frontal, 7=bamg, 8=delquad)
// Mesh.Algorithm = 8;

lc = 0.05;
nb_pts_on_line_x = 25;
nb_pts_on_line_y = 25;

Point(1) = { 0.0, 0.0, 0.0, lc };
Point(2) = { 1.0, 0.0, 0.0, lc };
Point(3) = { 1.0, 1.0, 0.0, lc };
Point(4) = { 0.0, 1.0, 0.0, lc };

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line Loop(1) = {1:4};

Plane Surface(1) = {1};

Physical Line("bottom") = {1};
Physical Line("right")  = {2};
Physical Line("top")    = {3};
Physical Line("left")   = {4};

Physical Surface("interior") = {1};

Transfinite Line(1) = nb_pts_on_line_x;
Transfinite Line(2) = nb_pts_on_line_y;
Transfinite Line(3) = nb_pts_on_line_x;
Transfinite Line(4) = nb_pts_on_line_y;

/*
Transfinite Surface{1} = {1,2,3,4};
Recombine Surface{1};
*/
