lc = 0.050;
L = 1.0;

factor = 1;
nb_pts = factor*30;


Point(1) = { 0.0  , 0.0, 0.0, lc };
Point(2) = {  L   , 0.0, 0.0, lc };
Point(3) = { 1.1*L,  L , 0.0, lc };
Point(4) = {-0.1*L,  L , 0.0, lc };
Point(5) = { 0.0  , 0.0,  L , lc };
Point(6) = {  L   , 0.0,  L , lc };
Point(7) = { 1.1*L,  L ,  L , lc };
Point(8) = {-0.1*L,  L ,  L , lc };


Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,5};

Line(9)  = {1,5};
Line(10) = {2,6};
Line(11) = {3,7};
Line(12) = {4,8};

Line Loop(1) = {1,2,3,4};
Line Loop(2) = {5,6,7,8};
Line Loop(3) = {1,10,-5,-9};
Line Loop(4) = {2,11,-6,-10};
Line Loop(5) = {3,12,-7,-11};
Line Loop(6) = {4,9,-8,-12};

For i In {1:6}
  Plane Surface(i) = {i};
EndFor

For i In {1:12}
  Transfinite Line(i) = nb_pts + 1;
EndFor

Surface Loop(1) = {1:6};
Volume(1) = {1};

Physical Surface("bottom") = {1};
Physical Surface("top") = {2};
Physical Surface("inlet") = {3};
Physical Surface("right") = {4};
Physical Surface("outlet") = {5};
Physical Surface("letf") = {6};

Physical Volume("domain") = {1};



