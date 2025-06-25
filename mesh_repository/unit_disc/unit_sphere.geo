nb_pts_on_edge = 16;
lc = 1.0/(nb_pts_on_edge-1);

Point(1) = {  0.0,  0.0,  0.0, lc };
Point(2) = {  1.0,  0.0,  0.0, lc };
Point(3) = {  0.0,  1.0,  0.0, lc };
Point(4) = { -1.0,  0.0,  0.0, lc };
Point(5) = {  0.0, -1.0,  0.0, lc };
Point(6) = {  0.0,  0.0, -1.0, lc };
Point(7) = {  0.0,  0.0,  1.0, lc };

Circle(1)  = {2,1,3};
Circle(2)  = {3,1,4};
Circle(3)  = {4,1,5};
Circle(4)  = {5,1,2};
Circle(5)  = {2,1,7};
Circle(6)  = {3,1,7};
Circle(7)  = {4,1,7};
Circle(8)  = {5,1,7};
Circle(9)  = {2,1,6};
Circle(10) = {3,1,6};
Circle(11) = {4,1,6};
Circle(12) = {5,1,6};

Line Loop(1) = {1,6,-5};
Line Loop(2) = {2,7,-6};
Line Loop(3) = {3,8,-7};
Line Loop(4) = {4,5,-8};
Line Loop(5) = {1,10,-9};
Line Loop(6) = {2,11,-10};
Line Loop(7) = {3,12,-11};
Line Loop(8) = {4,9,-12};

For iSurf In {1:8}
 Ruled Surface(iSurf) = {iSurf};
EndFor

Surface Loop(1) = {1:8};
Volume(1) = {1};

Physical Surface("boundary") = {1:8};

Physical Volume("interior") = {1};

For iLine In {1:12}
 Transfinite Line(iLine) = nb_pts_on_edge;
EndFor

