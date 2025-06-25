Mesh.CharacteristicLengthFromPoints = 1;

nb_pts_on_edge = 16;
lc = 1.0/(nb_pts_on_edge-1);

Point(1) = { 0.0, 0.0, 0.0, lc };
Point(2) = { 1.0, 0.0, 0.0, lc };
Point(3) = { 1.0, 1.0, 0.0, lc };
Point(4) = { 0.0, 1.0, 0.0, lc };
Point(5) = { 0.0, 0.0, 1.0, lc };
Point(6) = { 1.0, 0.0, 1.0, lc };
Point(7) = { 1.0, 1.0, 1.0, lc };
Point(8) = { 0.0, 1.0, 1.0, lc };

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

Line Loop(1)= {-4,-3,-2,-1};
Line Loop(2)= {5,6,7,8};
Line Loop(3)= {1,10,-5,-9};
Line Loop(4)= {2,11,-6,-10};
Line Loop(5)= {3,12,-7,-11};
Line Loop(6)= {4,9,-8,-12};

For iSurf In {1:6}
  Plane Surface(iSurf) = {iSurf};
EndFor

Surface Loop(1) = {1:6};
Volume(1) = {1};

Physical Surface("bottom") = {1};
Physical Surface("top") = {2};
Physical Surface("wall") = {3:6};

Physical Volume("interior") = {1};

For iLine In {1:12}
  Transfinite Line(iLine) = nb_pts_on_edge;
EndFor

/*
// Structured mesh

Transfinite Surface(1) = {1,2,3,4};
Transfinite Surface(2) = {5,6,7,8};
Transfinite Surface(3) = {1,2,6,5};
Transfinite Surface(4) = {2,3,7,6};
Transfinite Surface(5) = {3,4,8,7};
Transfinite Surface(6) = {4,1,5,8};

For iSurf In {1:6}
  Recombine Surface{iSurf};
EndFor

Transfinite Volume(1) = {1:8};
*/
