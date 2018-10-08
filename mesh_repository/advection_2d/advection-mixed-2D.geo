// 1=MeshAdapt, 2=Automatic, 5=Delaunay, 6=Frontal, 7=bamg, 8=delquad
Mesh.Algorithm = 2;

nb_pts_per_edge = 10;

L = 2.0;

lc = L/nb_pts_per_edge;

Point(1) = { -L,  0.0, 0.0,lc};
Point(2) = { 0.0, 0.0, 0.0,lc};
Point(3) = {  L,  0.0, 0.0,lc};
Point(4) = {  L,   L , 0.0,lc};
Point(5) = { 0.0,  L , 0.0,lc};
Point(6) = { -L,   L , 0.0,lc};

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

//Transfinite Surface(2) = {2,3,4,5};
Recombine Surface(2);

Physical Line("inlet") = {1};
Physical Line("outlet") = {2};
Physical Line("farfield")  = {4,6};

Physical Surface("fluid")  = {1,2};
