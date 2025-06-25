lc = 0.050;
Lx = 2.0;
Ly = 1.0;
Lz = 1.0;

factor = 1;
nb_pts_x = factor*10;
nb_pts_y = factor*10;
nb_pts_z = factor*10;


Point(1)  = { -Lx, 0.0,   -Lz, lc};
Point(2)  = { 0.0, 0.0,   -Lz, lc};
Point(3)  = {  Lx, 0.0,   -Lz, lc};
Point(4)  = {  Lx, 2.*Ly, -Lz, lc};
Point(5)  = { 0.0, 2.*Ly, -Lz, lc};
Point(6)  = { -Lx, 2.*Ly, -Lz, lc};
Point(7)  = { -Lx, 0.0,    Lz, lc};
Point(8)  = { 0.0, 0.0,    Lz, lc};
Point(9)  = {  Lx, 0.0,    Lz, lc};
Point(10) = {  Lx, 2.*Ly,  Lz, lc};
Point(11) = { 0.0, 2.*Ly,  Lz, lc};
Point(12) = { -Lx, 2.*Ly,  Lz, lc};

Line(1)  = {1,2};
Line(2)  = {2,5};
Line(3)  = {5,6};
Line(4)  = {6,1};
Line(5)  = {2,3};
Line(6)  = {3,4};
Line(7)  = {4,5};
Line(8)  = {7,8};
Line(9)  = {8,11};
Line(10) = {11,12};
Line(11) = {12,7};
Line(12) = {8,9};
Line(13) = {9,10};
Line(14) = {10,11};
Line(15) = {1,7};
Line(16) = {2,8};
Line(17) = {3,9};
Line(18) = {4,10};
Line(19) = {5,11};
Line(20) = {6,12};


Line Loop(1)  = {1,2,3,4};
Line Loop(2)  = {5,6,7,-2};
Line Loop(3)  = {8,9,10,11};
Line Loop(4)  = {12,13,14,-9};
Line Loop(5)  = {1,16,-8,-15};
Line Loop(6)  = {5,17,-12,-16};
Line Loop(7)  = {6,18,-13,-17};
Line Loop(8)  = {7,19,-14,-18};
Line Loop(9)  = {3,20,-10,-19};
Line Loop(10) = {4,15,-11,-20};
Line Loop(11) = {2,19,-9,-16};

For i In {1:11}
  Plane Surface(i) = {i};
EndFor

Surface Loop(1) = {1,5,3,9,10,11};
Surface Loop(2) = {2,8,4,6,11,7};

Volume(1) = {1};
Volume(2) = {2};

Transfinite Line(1)  = nb_pts_x + 1;
Transfinite Line(3)  = nb_pts_x + 1;
Transfinite Line(5)  = nb_pts_x + 1;
Transfinite Line(7)  = nb_pts_x + 1;
Transfinite Line(8)  = nb_pts_x + 1;
Transfinite Line(10) = nb_pts_x + 1;
Transfinite Line(12) = nb_pts_x + 1;
Transfinite Line(14) = nb_pts_x + 1;


Transfinite Line(2)  = nb_pts_y + 1;
Transfinite Line(4)  = nb_pts_y + 1;
Transfinite Line(6)  = nb_pts_y + 1;
Transfinite Line(9)  = nb_pts_y + 1;
Transfinite Line(11) = nb_pts_y + 1;
Transfinite Line(13) = nb_pts_y + 1;

Transfinite Line(15) = nb_pts_z + 1;
Transfinite Line(16) = nb_pts_z + 1;
Transfinite Line(17) = nb_pts_z + 1;
Transfinite Line(18) = nb_pts_z + 1;
Transfinite Line(19) = nb_pts_z + 1;
Transfinite Line(20) = nb_pts_z + 1;

Physical Surface("inlet") = {5}; 
Physical Surface("outlet") = {6}; 
//Physical Surface("farfield") = {1,2,10,9,8,7,3,4};
Physical Surface("farfield") = {10,8};

Physical Volume("fluid") = {1,2};

