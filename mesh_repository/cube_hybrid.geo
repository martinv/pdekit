lc = 1.2;

Point(1) = {-1.,-1.,-1., lc };
Point(2) = { 1.,-1.,-1., lc };
Point(3) = { 1., 1.,-1., lc };
Point(4) = {-1., 1.,-1., lc };
Point(5) = {-1.,-1., 1., lc };
Point(6) = { 1.,-1., 1., lc };
Point(7) = { 1., 1., 1., lc };
Point(8) = {-1., 1., 1., lc };

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,5};
Line(9) = {1,5};
Line(10) = {2,6};
Line(11) = {3,7};
Line(12) = {4,8};

Line Loop(1) = {1,2,3,4};
Line Loop(2) = {5,6,7,8};
Line Loop(3) = {1,10,-5,-9};
Line Loop(4) = {2,11,-6,-10};
Line Loop(5) = {3,12,-7,-11};
Line Loop(6) = {4,9,-8,-12};

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};
Plane Surface(6) = {6};

Surface Loop(1) = {1:6};
Volume(1) = {1};

Physical Surface("Wall") = {1:6};
Physical Volume("Volume") = {1};

Transfinite Line(1) = 4;
Transfinite Line(2) = 4;
Transfinite Line(3) = 4;
Transfinite Line(4) = 4;

Transfinite Surface(1) = {1,2,3,4};
Recombine Surface(1);


