Include "naca_contour.geo";

// =============================================================
// USER-DEFINED PARAMETERS
// =============================================================
// Length of the domain in front of the cylinder
L1 = 15;

// Length of the domain in behind the cylinder
L2 = 20;

// Half-height of the domain
H = 10;

// Extrusion of the airfoil contour
out[] = Extrude{0,0,2*W}{ Line{1:6}; };

// Points of the left wall (which contains the original airfoil contour)
Point(10001) = { -L1 , -H, -W, lc };
Point(10002) = {  L2 , -H, -W, lc };
Point(10003) = {  L2 ,  H, -W, lc };
Point(10004) = { -L1 ,  H, -W, lc };


// Points of the right wall (which contains the extruded airfoil contour)
Point(10005) = { -L1 , -H, W, lc };
Point(10006) = {  L2 , -H, W, lc };
Point(10007) = {  L2 ,  H, W, lc };
Point(10008) = { -L1 ,  H, W, lc };

// =============================================================
// LINES
// =============================================================

// Lines of the left wall (contains the original airfoil contour)
Line(101) = {10001,10002};
Line(102) = {10002,10003};
Line(103) = {10003,10004};
Line(104) = {10004,10001};

// Lines of the right wall (contains the extruded airfoil contour)
Line(105) = {10005,10006};
Line(106) = {10006,10007};
Line(107) = {10007,10008};
Line(108) = {10008,10005};

// Lines connecting the left and right walls
Line(109) = { 10001, 10005 };
Line(110) = { 10002, 10006 };
Line(111) = { 10003, 10007 };
Line(112) = { 10004, 10008 };

// =============================================================
// SURFACES
// =============================================================

// Airfoil on the left
Line Loop(1) = {1:6};
// Left wall
Line Loop(2) = {101,102,103,104};
Plane Surface(1) = {1,2};

// Airfoil on the right
Line Loop(3) = { 7,11,15,19,23,27};
// Right wall
Line Loop(4) = {105,106,107,108};
Plane Surface(2) = {3,4};

// Rear wall (outlet)
Line Loop(5) = {102,111,-106,-110};
Plane Surface(3) = {5};

// Front wall (inlet)
Line Loop(6) = {-108,-112,104,109};
Plane Surface(4) = {6};

// Bottom wall
Line Loop(7) = {105,-110,-101,109};
Plane Surface(5) = {7};

// Top wall
Line Loop(8) = {-107,-111,103,112};
Plane Surface(6) = {8};


// =============================================================
// VOLUMES
// =============================================================
// Walls
Surface Loop(1) = {1:6};
// Surface of the airfoil
Surface Loop(2) = {10,14,18,22,26,30};

Volume(1) = {1,2};


// =============================================================
// MESH SIZES
// =============================================================

nb_pts_le_segment = 6;
nb_pts_te_segment = 15;
nb_pts_middle     = 20;

Transfinite Line{1} = nb_pts_le_segment;
Transfinite Line{7} = nb_pts_le_segment;

Transfinite Line{2}  = nb_pts_middle Using Bump 0.90;
Transfinite Line{11} = nb_pts_middle Using Bump 0.90;

Transfinite Line{3}  = nb_pts_te_segment;
Transfinite Line{15} = nb_pts_te_segment;

Transfinite Line{4}  = nb_pts_te_segment;
Transfinite Line{19} = nb_pts_te_segment;

Transfinite Line{5}  = nb_pts_middle Using Bump 0.90;
Transfinite Line{23} = nb_pts_middle Using Bump 0.90;

Transfinite Line{6}  = nb_pts_le_segment;
Transfinite Line{27} = nb_pts_le_segment;

nb_pts_on_span = 60;

Transfinite Line{17} = nb_pts_on_span;
Transfinite Line{13} = nb_pts_on_span;
Transfinite Line{9}  = nb_pts_on_span;
Transfinite Line{8}  = nb_pts_on_span;
Transfinite Line{25} = nb_pts_on_span;
Transfinite Line{21} = nb_pts_on_span;

Transfinite Surface{18};
Transfinite Surface{14};
Transfinite Surface{10};
Transfinite Surface{30};
Transfinite Surface{26};
Transfinite Surface{22};

R = 0.5;

Field[1] = Attractor;
Field[1].FacesList = {10,14,18,22,26,30};

Field[2] = Threshold;
Field[2].DistMin = 1.5*R;
Field[2].DistMax = 5.0*R;
Field[2].IField = 1;
Field[2].LcMin = 0.05*lc;
Field[2].LcMax = lc;
Field[2].Sigmoid = 0;
Field[2].StopAtDistMax = 0;


Field[3] = Box;
Field[3].XMin = -2.0 * R;
Field[3].XMax =  8.0 * R;
Field[3].YMin = -2.0 * R;
Field[3].YMax =  2.0 * R;
Field[3].ZMin = -W;
Field[3].ZMax =  W;
Field[3].VIn =   0.1*lc;
Field[3].VOut =  1.0*lc;

Field[4] = Min;
Field[4].FieldsList = {2,3};

Background Field = 4;

// =============================================================
// PHYSICAL ENTITIES
// =============================================================
Physical Surface("Airfoil") = {10,14,18,22,26,30};
Physical Surface("Inlet") = {4};
Physical Surface("Outlet") = {3};
Physical Surface("SymmetryLeft") = {1};
Physical Surface("SymmetryRight") = {2};
Physical Surface("Bottom") = {5};
Physical Surface("Top") = {6};

Physical Volume("Fluid") = {1};


