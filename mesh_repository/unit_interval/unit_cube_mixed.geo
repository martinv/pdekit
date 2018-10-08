Mesh.CharacteristicLengthFromPoints = 1;

nb_pts_on_edge = 4;
lc = 1.0/(nb_pts_on_edge-1);

Point(1) = { 0.0, 0.0, 0.0, lc };
Point(2) = { 1.0, 0.0, 0.0, lc };
Point(3) = { 1.0, 1.0, 0.0, lc };
Point(4) = { 0.0, 1.0, 0.0, lc };
Point(5) = { 0.0, 0.0, 1.0, lc };
Point(6) = { 1.0, 0.0, 1.0, lc };
Point(7) = { 1.0, 1.0, 1.0, lc };
Point(8) = { 0.0, 1.0, 1.0, lc };

Point(9) = { 0.25, 0.25, 1.0, lc };
Point(10) = { 0.75, 0.25, 1.0, lc };
Point(11) = { 0.75, 0.75, 1.0, lc };
Point(12) = { 0.25, 0.75, 1.0, lc };

Line(1)  = {1,2};
Line(2)  = {2,3};
Line(3)  = {3,4};
Line(4)  = {4,1};
Line(5)  = {5,6};
Line(6)  = {6,7};
Line(7)  = {7,8};
Line(8)  = {8,5};
Line(9)  = {1,5};
Line(10) = {2,6};
Line(11) = {3,7};
Line(12) = {4,8};

Line(13) = {9, 10};
Line(14) = {10, 11};
Line(15) = {11, 12};
Line(16) = {12, 9};

Line Loop(1)= {-4,-3,-2,-1};
Line Loop(2)= {5,6,7,8};
Line Loop(3)= {1,10,-5,-9};
Line Loop(4)= {2,11,-6,-10};
Line Loop(5)= {3,12,-7,-11};
Line Loop(6)= {4,9,-8,-12};

Line Loop(7) = {13,14,15,16};


Plane Surface(1) = {1};
Plane Surface(2) = {2,7};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};
Plane Surface(6) = {6};

Plane Surface(7) = {7};


For iLine In {13:16}
  Transfinite Line(iLine) = nb_pts_on_edge;
EndFor

Transfinite Surface(7) = {9, 10, 11, 12};
Recombine Surface{7};

Surface Loop(1) = {1:7};
Volume(1) = {1};

Physical Surface("bottom") = {1};
Physical Surface("top") = {2,7};
Physical Surface("wall") = {3:6};

Physical Volume("interior") = {1};

