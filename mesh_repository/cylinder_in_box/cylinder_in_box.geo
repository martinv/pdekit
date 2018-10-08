lc = 10.0;

// =============================================================
// USER-DEFINED PARAMETERS
// =============================================================

// Radius of the cylinder
R = 1.0;

// Half-width of the domain
W = 5.0;

// Length of the domain in front of the cylinder
L1 = 15;

// Length of the domain in behind the cylinder
L2 = 20;

// Half-height of the domain
H = 10;

// =============================================================
// POINTS
// =============================================================


// Points of the first circle
Point(1) = { 0.0, 0.0, W, lc };
Point(2) = { 0.0,  R,  W, lc };
Point(3) = { -R , 0.0, W, lc };
Point(4) = { 0.0, -R , W, lc };
Point(5) = {  R , 0.0, W, lc };

// Points of the left wall (which contains the first circle)
Point(6) = { -L1 , -H, W, lc };
Point(7) = {  L2 , -H, W, lc };
Point(8) = {  L2 ,  H, W, lc };
Point(9) = { -L1 ,  H, W, lc };


// Points of the second circle
Point(10) = { 0.0, 0.0, -W, lc };
Point(11) = { 0.0,  R,  -W, lc };
Point(12) = { -R , 0.0, -W, lc };
Point(13) = { 0.0, -R , -W, lc };
Point(14) = {  R , 0.0, -W, lc };

// Points of the left wall (which contains the first circle)
Point(15) = { -L1 , -H, -W, lc };
Point(16) = {  L2 , -H, -W, lc };
Point(17) = {  L2 ,  H, -W, lc };
Point(18) = { -L1 ,  H, -W, lc };

// =============================================================
// LINES
// =============================================================

// Arcs of the first circle
Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};

// Lines of the left wall (contains the first circle)
Line(5) = {6,7};
Line(6) = {7,8};
Line(7) = {8,9};
Line(8) = {9,6};

// Arcs of the second circle
Circle(9)  = {11,10,12};
Circle(10) = {12,10,13};
Circle(11) = {13,10,14};
Circle(12) = {14,10,11};


// Lines of the right wall (contains the second circle)
Line(13) = {15,16};
Line(14) = {16,17};
Line(15) = {17,18};
Line(16) = {18,15};

// Lines connecting the left and right boundaries
Line(17) = {2,11};
Line(18) = {3,12};
Line(19) = {4,13};
Line(20) = {5,14};

Line(21) = {6,15};
Line(22) = {7,16};
Line(23) = {8,17};
Line(24) = {9,18};


// =============================================================
// SURFACES
// =============================================================

// Surface of the cylinder
Line Loop(1) = {1,18,-9,-17};
Line Loop(2) = {2,19,-10,-18};
Line Loop(3) = {3,20,-11,-19};
Line Loop(4) = {4,17,-12,-20};

Line Loop(5) = {1,2,3,4};
Line Loop(6) = {9,10,11,12};

Ruled Surface(1) = {1};
Ruled Surface(2) = {2};
Ruled Surface(3) = {3};
Ruled Surface(4) = {4};

// Left wall
Line Loop(7) = {5,6,7,8};
Plane Surface(7) = {7,5};

// Rear wall (outlet)
Line Loop(8) = {22,14,-23,-6};
Plane Surface(8) = {8};

// Right wall
Line Loop(9) = {13,14,15,16};
Plane Surface(9) = {9,6};

// Front wall (inlet)
Line Loop(10) = {8,21,-16,-24};
Plane Surface(10) = {10};

// Bottom wall
Line Loop(11) = {5,22,-13,-21};
Plane Surface(11) = {11};

// Top wall
Line Loop(12) = {-7,23,15,-24};
Plane Surface(12) = {12};


// =============================================================
// VOLUMES
// =============================================================
// Walls
Surface Loop(1) = {7,8,9,10,11,12};
// Surface of the cylinder
Surface Loop(2) = {1,2,3,4};

Volume(1) = {1,2};


// =============================================================
// MESH SIZES
// =============================================================
Field[1] = Attractor;
Field[1].FacesList = {1,2,3,4};

Field[2] = Threshold;
Field[2].DistMin = 1.5*R;
Field[2].DistMax = 5.0*R;
Field[2].IField = 1;
Field[2].LcMin = 0.25*R;
Field[2].LcMax = 1.3*R;
Field[2].Sigmoid = 0;
Field[2].StopAtDistMax = 0;

Field[3] = Box;
Field[3].XMin = -2.0 * R;
Field[3].XMax =  8.0 * R;
Field[3].YMin = -2.0 * R;
Field[3].YMax =  2.0 * R;
Field[3].ZMin = -W;
Field[3].ZMax =  W;
Field[3].VIn =   0.3*R;
Field[3].VOut =  1.2*R;

Field[4] = Min;
Field[4].FieldsList = {2,3};

Background Field = 4;

// =============================================================
// PHYSICAL ENTITIES
// =============================================================
Physical Surface("CylinderWall") = {1,2,3,4};
Physical Surface("Inlet") = {10};
Physical Surface("Outlet") = {8};
Physical Surface("SymmetryLeft") = {7};
Physical Surface("SymmetryRight") = {9};
Physical Surface("Bottom") = {11};
Physical Surface("Top") = {12};

Physical Volume("Fluid") = {1};


