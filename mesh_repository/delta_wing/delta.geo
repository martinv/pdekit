lc = 0.07;
lc_ff = 1;

// Mesh.Algorithm = 6;

x_min = -6.5;
x_max =  6.5;
y_min =  0.0;
y_max =  5.3;
z_min = -5.3;
z_max =  5.3;

// Wing geometry
Point(1) = {1.0, 0.0, 0.0, lc};
Point(2) = {1.0, 0.0,  -0.024416137, lc};
Point(3) = {1.0, 0.22565924, -0.024416137, lc};
Point(4) = {1.0, 0.26794922,  0.0, lc};
Point(5) = {0.15782839, 0.0, -0.024416137, lc};
Point(6) = {0.0, 0.0, 0.0, lc};

// Bounding box
Point(10) = {x_max, y_min, z_min, lc_ff};
Point(11) = {x_max, y_max, z_min, lc_ff};
Point(12) = {x_max, y_max, z_max, lc_ff};
Point(13) = {x_max, y_min, z_max, lc_ff};

Point(14) = {x_min, y_min, z_min, lc_ff};
Point(15) = {x_min, y_max, z_min, lc_ff};
Point(16) = {x_min, y_max, z_max, lc_ff};
Point(17) = {x_min, y_min, z_max, lc_ff};

// Wing geometry - lines
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line(5) = {4,6};
Line(6) = {3,5};
Line(7) = {2,5};
Line(8) = {1,6};
Line(9) = {5,6};

// Bounding box - lines
Line(20) = {10,11};
Line(21) = {11,12};
Line(22) = {12,13};
Line(23) = {13,10};

Line(24) = {14,15};
Line(25) = {15,16};
Line(26) = {16,17};
Line(27) = {17,14};

Line(28) = {11,15};
Line(29) = {14,10};
Line(30) = {12,16};
Line(31) = {17,13};

// Line loops - wing

Line Loop(1) = {1,2,3,4};
Line Loop(2) = {-4,5,-8};
Line Loop(3) = {8,-9,-7,-1};
Line Loop(4) = {-6,-2,7};
Line Loop(5) = {6,9,-5,-3};

// Line loops - bounding box
Line Loop(10) = {20,21, 22, 23};
Line Loop(11) = {28,25,-30,-21};
Line Loop(12) = {-24, -27, -26, -25};
Line Loop(13) = {29, -23, -31, 27};
Line Loop(14) = {24,-28, -20, -29};
Line Loop(15) = {30, 26, 31, -22};

// Wing surfaces
Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(4) = {4};
Plane Surface(5) = {5};

// Farfield surfaces
Plane Surface(10) = {10};
Plane Surface(11) = {11};
Plane Surface(12) = {12};
Plane Surface(13) = {13,3};
Plane Surface(14) = {14};
Plane Surface(15) = {15};

// Volume
Surface Loop(1) = {1,2,4,5,10:15};
Volume(1) = {1};

Physical Surface("Wall") = {1,2,4,5};
Physical Surface("Symmetry") = {13};
Physical Surface("Farfield") = {10,11,12,13,14,15};
Physical Volume("Fluid") = {1};

