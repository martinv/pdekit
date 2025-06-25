L = 15.0; // Half of the length of the channel
W = 5.0;  // Half of the width of the channel
H = 5.0;  // Half of the height of the channel

nb_points_length = 30;
nb_points_width  = 12;
nb_points_height = 12;

lc = 2.0*L/nb_points_length;

Point(1)  = {  W, -L, -H };
Point(2)  = {  W,  L, -H };
Point(3)  = { -W,  L, -H };
Point(4)  = { -W, -L, -H };
Point(5)  = {  W, -L,  H };
Point(6)  = {  W,  L,  H };
Point(7)  = { -W,  L,  H };
Point(8)  = { -W, -L,  H };

Line(1)  = { 1, 2 };
Line(2)  = { 2, 3 };
Line(3)  = { 3, 4 };
Line(4)  = { 4, 1 };
Line(5)  = { 5, 6 };
Line(6)  = { 6, 7 };
Line(7)  = { 7, 8 };
Line(8)  = { 8, 5 };
Line(9)  = { 1, 5 };
Line(10) = { 2, 6 };
Line(11) = { 3, 7 };
Line(12) = { 4, 8 };

Line Loop(1) = { 1, 2, 3, 4 };
Line Loop(2) = { 5, 6, 7, 8 };
Line Loop(3) = { 4, 9, -8, -12 };
Line Loop(4) = { 1, 10, -5, -9 };
Line Loop(5) = { 2, 11, -6, -10 };
Line Loop(6) = { 3, 12, -7, -11 };

Plane Surface(1) = { 1 };
Plane Surface(2) = { 2 };
Plane Surface(3) = { 3 };
Plane Surface(4) = { 4 };
Plane Surface(5) = { 5 };
Plane Surface(6) = { 6 };

Surface Loop(1) = { 1:6 };
Volume(1) = { 1 };

Physical Surface("Bottom") = { 1 };
Physical Surface("Top")    = { 2 };
Physical Surface("Inlet")  = { 3 };
Physical Surface("Outlet") = { 5 };
Physical Surface("SymmetryLeft") = { 6 };
Physical Surface("SymmetryRight") = { 4 };

Physical Volume("InnerCells") = { 1 };

Transfinite Line(1) = nb_points_length;
Transfinite Line(3) = nb_points_length;
Transfinite Line(5) = nb_points_length;
Transfinite Line(7) = nb_points_length;

Transfinite Line(4) = nb_points_width;
Transfinite Line(2) = nb_points_width;
Transfinite Line(6) = nb_points_width;
Transfinite Line(8) = nb_points_width;

Transfinite Line(9)  = nb_points_height;
Transfinite Line(10) = nb_points_height;
Transfinite Line(11) = nb_points_height;
Transfinite Line(12) = nb_points_height;
