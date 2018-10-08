Mesh.Algorithm = 5;    // 1=MeshAdapt, 2=Automatic, 5=Delaunay, 6=Frontal, 7=bamg, 8=delquad
Mesh.Algorithm3D = 1;  // 1=Delaunay, 4=Frontal, 5=Frontal Delaunay, 6=Frontal Hex, 7=MMG3D, 9=R-tree

L = 1.0;

lc = L/10;

nb_pts_x = 9;
nb_pts_y = 9;
nb_pts_z = 9;

// ----------------------------------------------------------------------------
// Points
// ----------------------------------------------------------------------------

Point(1)  = { 0.0 , 0.0, -L/2, lc};
Point(2)  = { L   , 0.0, -L/2, lc};
Point(3)  = { 2*L , 0.0, -L/2, lc};
Point(4)  = { 3*L , 0.0, -L/2, lc};

Point(5)  = { 0.0 , 0.0,  L/2, lc};
Point(6)  = { L   , 0.0,  L/2, lc};
Point(7)  = { 2*L , 0.0,  L/2, lc};
Point(8)  = { 3*L , 0.0,  L/2, lc};

Point(9)  = { 0.0 ,  L , -L/2, lc};
Point(10) = { 3*L ,  L , -L/2, lc};
Point(11) = { 0.0 ,  L ,  L/2, lc};
Point(12) = { 3*L ,  L ,  L/2, lc};

Point(100) = { 1.5*L, -1.2*L, -L/2, lc };
Point(101) = { 1.5*L, -1.2*L,  L/2, lc };

// ----------------------------------------------------------------------------
// Lines
// ----------------------------------------------------------------------------

Line(1)   = {1,2};
Circle(2) = {2, 100, 3};
Line(3)   = {3,4};

Line(4)   = {5,6};
Circle(5) = {6, 101, 7};
Line(6)   = {7,8};

Line(7)   = {1,5};
Line(8)   = {2,6};
Line(9)   = {3,7};
Line(10)  = {4,8};

Line(11)  = {9,11};
Line(12)  = {11,12};
Line(13)  = {10,12};
Line(14)  = {9,10};

Line(15)  = {5,11};
Line(16)  = {1,9};
Line(17)  = {8,12};
Line(18)  = {4,10};

// ----------------------------------------------------------------------------
// Line loops
// ----------------------------------------------------------------------------

// Bottom wall
Line Loop(1) = {7,4,-8,-1};
Line Loop(2) = {8,5,-9,-2};
Line Loop(3) = {9,6,-10,-3};

// Left side
Line Loop(4) = {1,2,3,18,-14,-16};
// Right side
Line Loop(5) = {4,5,6,17,-12,-15};
// Inlet
Line Loop(6) = {7,15,-11,-16};
// Outlet
Line Loop(7) = {-10, 18, 13, -17};
// Top wall
Line Loop(8) = {11,12,-13,-14};

// ----------------------------------------------------------------------------
// Surfaces
// ----------------------------------------------------------------------------

Plane Surface(1) = {1};
Ruled Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};
Plane Surface(6) = {6};
Plane Surface(7) = {7};
Plane Surface(8) = {8};

// ----------------------------------------------------------------------------
// Volumes
// ----------------------------------------------------------------------------

Surface Loop(1) = {1:8};
Volume(1) = {1};

// ----------------------------------------------------------------------------
// Physical entities
// ----------------------------------------------------------------------------
Physical Surface("Inlet") = {6};
Physical Surface("Outlet") = {7};
Physical Surface("SymmetryLeft") = {4};
Physical Surface("SymmetryRight") = {5};
Physical Surface("Bottom") = {1,2,3};
Physical Surface("Top") = {8};

Physical Volume("Fluid") = {1};

// ----------------------------------------------------------------------------
// Force placing of points on geometry
// ----------------------------------------------------------------------------
Transfinite Line(1) = nb_pts_x;
Transfinite Line(2) = nb_pts_x;
Transfinite Line(3) = nb_pts_x;
Transfinite Line(4) = nb_pts_x;
Transfinite Line(5) = nb_pts_x;
Transfinite Line(6) = nb_pts_x;

Transfinite Line{15} = nb_pts_y;
Transfinite Line{16} = nb_pts_y;
Transfinite Line{17} = nb_pts_y;
Transfinite Line{18} = nb_pts_y;

Transfinite Line(7)  = nb_pts_z;
Transfinite Line(10) = nb_pts_z;
Transfinite Line(11) = nb_pts_z;
Transfinite Line(13) = nb_pts_z;


