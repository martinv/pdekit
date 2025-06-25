// ----------------------------------------------------------------------------
// Global variables: characteristic length on body and farfield
// ----------------------------------------------------------------------------
lc_body     = 0.05;
lc_farfield = 4.0;
delta_y = 0.01;

// ----------------------------------------------------------------------------
// FarField: a sphere with radius = 10 and center at {1,0,0} to represent the farfield boundary
// Helper variables
// ----------------------------------------------------------------------------
ff_radius = 10.0;
pi_half=Pi/2;
xc = 1.0;
yc = 0.0;
zc = 0.0;
ff_center[] = {xc,yc,zc};

// Minimum and maximum char. lengths
Mesh.CharacteristicLengthMin = 0.1*lc_body;
Mesh.CharacteristicLengthMax = lc_farfield;

// Control point distribution on the elliptic part of the body
nb_pts_ellipse = 10;
// Progression to cluster the points towards the tip of the body
ellipse_progr = 1.2;
// Number of points on the middle stripe
nb_pts_mid = 2; 

//   Mesh.Algorithm3D=1; // for Tetgen+Delaunay (default)
//   Mesh.Algorithm3D=4; // for Netgen
Mesh.Algorithm = 2; //6; // 1=MeshAdapt, 2=Automatic, 5=Delaunay, 6=Frontal, 7=bamg, 8=delquad
Mesh.Algorithm3D = 4;

// Optimize 2D mesh?
// Mesh.Optimize=1;

// Optimize 3D mesh?
// Mesh.OptimizeNetgen=1;

// Optimize high order meshes? Default: 0
Mesh.HighOrderOptimize = 1;

// Minimum threshold for high order element optimization, default: 0.1
Mesh.HighOrderThresholdMin = 0.5;

// Maximum threshold for high order element optimization, default: 2.0
Mesh.HighOrderThresholdMax = 2.0;

// Try to fix flipped surface mesh elements in high-order optimizer, default: 0
Mesh.HighOrderOptPrimSurfMesh = 1;

Point(1) =  {1/4, -delta_y, 0, lc_body};
Point(2) =  {1/4,  delta_y, 0, lc_body};

Point(10) = {0.0, -delta_y, 0.0, lc_body};
Point(11) = {1/4, -delta_y, 0.05, lc_body};
Point(12) = {1/3, -delta_y, Sqrt(2.)/30, lc_body};
Point(13) = {1.0, -delta_y, 0.0, lc_body};

Point(14) = {1.0,  delta_y, 0.0, lc_body};
Point(15) = {1/3,  delta_y, Sqrt(2.)/30, lc_body};
Point(16) = {1/4,  delta_y, 0.05, lc_body};
Point(17) = {0.0,  delta_y, 0.0, lc_body};

// Ellipse(Nb. of ellipse) = {Start, Center, Main axis, Stop};
Ellipse(1) = {10, 1, 11, 12};
// Line {Start vertex, End vertex}
Line(2) = {12, 13};

Line(3) = {13, 14};

Line(4) = {14, 15};
// Ellipse(Nb. of ellipse) = {Start, Center, Main axis, Stop};
Ellipse(5) = {17, 2, 16, 15};
Line(6) = {17, 10};

Line(7) = {12, 15};

Line Loop(1) = {1,7, -5, 6};
Line Loop(2) = {2,3,4,-7};

Ruled Surface(1) = {1};
Ruled Surface(2) = {2};

mid_top = {1,2};
// Extrusion by rotation
// Extrude {{rotation axis}, {center of rotation}, angle} {lines to extrude};
sides_top_right[] = Extrude {{1, 0, 0}, {1, -0.01, 0},  Pi/2}{ Line{ 1, 2}; };
sides_top_left[]    = Extrude {{1, 0, 0}, {1,  0.01, 0}, -Pi/2}{ Line{ 4, 5}; };

// Symmetry {a,b,c,d} {Duplicata{Surface numbers}}; 
// Symmetry with respect to the plane given by ax+by+cz+d=0
// Duplicata indicates that we mean to preserve the source surface
mid_bottom[]   = Symmetry {0, 0, 1, 0} {Duplicata{Surface{1,2};}};
sides_bottom_right[] = Symmetry {0, 0, 1, 0} {Duplicata{Surface{sides_top_right[{1,4}]};}};
sides_bottom_left[]  = Symmetry {0, 0, 1, 0} {Duplicata{Surface{sides_top_left[{1,4}]};}};

// For i In {0:(#sides_bottom_right[]-1)}
//   Printf("i = %g, entity idx = %g",i, sides_bottom_right[i]);
// EndFor


// ---------------------------------------------------------------------------- 
// Creation of farfield sphere
// ---------------------------------------------------------------------------- 

i=0;
// Creation of points necessary for the sphere. The first point is the center.
i++; pSphere[i] = newp; Point(pSphere[i]) = {xc          , yc          , zc          , lc_farfield};
i++; pSphere[i] = newp; Point(pSphere[i]) = {xc-ff_radius, yc          , zc          , lc_farfield};
i++; pSphere[i] = newp; Point(pSphere[i]) = {xc       ,    yc+ff_radius, zc          , lc_farfield};
i++; pSphere[i] = newp; Point(pSphere[i]) = {xc       ,    yc          , zc+ff_radius, lc_farfield};

j=0;
// Create the shell.
// Circle(Circle number) = {start, center, stop}; 
// The circle angle cannot be larger than pi/2
j++; lSphere[j] = newl; Circle(lSphere[j]) = {pSphere[2], pSphere[1], pSphere[3]};
j++; lSphere[j] = newl; Circle(lSphere[j]) = {pSphere[3], pSphere[1], pSphere[4]};
j++; lSphere[j] = newl; Circle(lSphere[j]) = {pSphere[4], pSphere[1], pSphere[2]};

// Line Loop of circles...
lShell = newll; Line Loop(lShell)     = {lSphere[{1:3}]};
// This defines 1/8 of the spherical surface.
sShell = news;  Ruled Surface(sShell) = {lShell};


// create remaining 7/8 inner shells
// Rotate {{Axe de rotation},{coord d'un point de l'axe}, angle} {N° de l'objet à rotater};
// Encore une fois, Duplicata pour conserver l'original et une liste devant pour recevoir le 
// numéro de l'objet obtenu.
t1[] = Rotate {{0,0,1},{xc,yc,zc},pi_half}   {Duplicata{Surface{sShell};}};
t2[] = Rotate {{0,0,1},{xc,yc,zc},pi_half*2} {Duplicata{Surface{sShell};}};
t3[] = Rotate {{0,0,1},{xc,yc,zc},pi_half*3} {Duplicata{Surface{sShell};}};
t4[] = Rotate {{0,1,0},{xc,yc,zc},-pi_half}  {Duplicata{Surface{sShell};}};
t5[] = Rotate {{0,0,1},{xc,yc,zc},pi_half}   {Duplicata{Surface{t4[0]};}};
t6[] = Rotate {{0,0,1},{xc,yc,zc},pi_half*2} {Duplicata{Surface{t4[0]};}};
t7[] = Rotate {{0,0,1},{xc,yc,zc},pi_half*3} {Duplicata{Surface{t4[0]};}};

// ---------------------------------------------------------------------------- 
// Surface loops and volume
// ---------------------------------------------------------------------------- 

// Wall
Surface Loop(1) = {mid_top[{0,1}],   sides_top_right[{1,4}],    sides_top_left[{1,4}],
                  mid_bottom[{0,1}], sides_bottom_right[{0,1}], sides_bottom_left[{0,1}]};

// Farfield
Surface Loop(2) = {sShell, t1[0], t2[0], t3[0], t4[0], t5[0], t6[0], t7[0]};

// Volume
Volume(1) = {2,1};

// ---------------------------------------------------------------------------- 
// Control point distribution on the body with transfinite commands
// ---------------------------------------------------------------------------- 

// Number of points close to the front of the body
/*
Transfinite Line(1) = nb_pts_ellipse Using Progression ellipse_progr;
Transfinite Line(5) = nb_pts_ellipse Using Progression ellipse_progr;
Transfinite Line(8) = nb_pts_ellipse Using Progression ellipse_progr;
Transfinite Line(17) = nb_pts_ellipse Using Progression ellipse_progr;
Transfinite Line(21) = nb_pts_ellipse Using Progression ellipse_progr;
Transfinite Line(-23) = nb_pts_ellipse Using Progression ellipse_progr;

// Number of points across the middle stripe
Transfinite Line(3) = nb_pts_mid;
Transfinite Line(6) = nb_pts_mid;
Transfinite Line(7) = nb_pts_mid;
Transfinite Line(22) = nb_pts_mid;
*/

// Front stagnation point:
Field[1] = Attractor;
Field[1].EdgesList = { 6 };

// Trailing edge
Field[2] = Attractor;
Field[2].EdgesList = { 3 };

// Local refinement around the front stagnation point:
Field[5] = Threshold;
Field[5].IField = 1;
Field[5].DistMin = 0.1;
Field[5].DistMax = 0.3;
Field[5].LcMin = 0.01;
Field[5].LcMax = 2.0;
Field[5].Sigmoid = 1;
//Don't impose anything outside DistMax:
Field[5].StopAtDistMax = 1;

// Local refinement around the trailing edge
Field[7] = Threshold;
Field[7].IField = 2;
Field[7].DistMin = 0.1;
Field[7].DistMax = 0.3;
Field[7].LcMin = 0.01;
Field[7].LcMax = 2.0;
Field[7].Sigmoid = 1;
//Don't impose anything outside DistMax:
Field[7].StopAtDistMax = 1;

Field[10] = Min;
Field[10].FieldsList = { 5, 7 };
Background Field = 10;


// ---------------------------------------------------------------------------- 
// Physical entities
// ---------------------------------------------------------------------------- 

Physical Surface("WallTop") = {mid_top[{0,1}], sides_top_right[{1,4}], sides_top_left[{1,4}]};
Physical Surface("WallBottom") = {mid_bottom[{0,1}], sides_bottom_right[{0,1}], sides_bottom_left[{0,1}]};
Physical Surface("Farfield") = {sShell, t1[0], t2[0], t3[0], t4[0], t5[0], t6[0], t7[0]};
Physical Volume("Fluid") = {1};

