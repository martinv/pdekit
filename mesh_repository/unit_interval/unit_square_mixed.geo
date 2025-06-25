nb_pts_on_edge = 16;
lc = 1.0/(nb_pts_on_edge-1);

Point(1) = { 0.0, 0.0, 0.0, lc };
Point(2) = { 1.0, 0.0, 0.0, lc };
Point(3) = { 1.0, 1.0, 0.0, lc };
Point(4) = { 0.0, 1.0, 0.0, lc };
Point(5) = { 0.5, 0.0, 0.0, lc };
Point(6) = { 1.0, 0.5, 0.0, lc };
Point(7) = { 0.5, 1.0, 0.0, lc };
Point(8) = { 0.0, 0.5, 0.0, lc };

Line(1)  = {1,5};
Line(2)  = {5,2};
Line(3)  = {2,6};
Line(4)  = {6,3};
Line(5)  = {3,7};
Line(6)  = {7,4};
Line(7)  = {4,8};
Line(8)  = {8,1};
Line(9)  = {8,5};
Line(10) = {5,6};
Line(11) = {6,7};
Line(12) = {7,8};


Line Loop(1)= {1,-9,8};
Line Loop(2)= {2,3,-10};
Line Loop(3)= {4,5,-11};
Line Loop(4)= {6,7,-12};
Line Loop(5)= {9,10,11,12};

For isurf In {1:5}
 Plane Surface(isurf) = {isurf};
EndFor

For iline In {1:12}
 Transfinite Line(iline) = nb_pts_on_edge;
EndFor

Transfinite Surface(5) = {5,6,7,8};
Recombine Surface{5};

Physical Line("left") = {7,8};
Physical Line("bottom") = {1,2};
Physical Line("right") = {3,4};
Physical Line("top") = {5,6};

Physical Surface("interior") = {1:6};

