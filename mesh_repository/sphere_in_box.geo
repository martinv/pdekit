// Mesh of a sphere in a box
// Run from command line as:
// gmsh -3 -order 2 sphere_in_box.geo -o sphere_p2.msh
// If you have a recent version of gmsh (April 2012), you can run as
// gmsh -3 -order 2 -hoOptimize sphere_in_box.geo -o sphere_p2.msh

// *******************************************************************************
// CHANGE ANY OF THE FOLLOWING PARAMETERS TO CHANGE THE PROPORTIONS OF THE DOMAIN
// *******************************************************************************
// NOTE: The center of the sphere is located in origin (0,0,0)
// Radius of the sphere
   RSphere = 1.;
// Length of the domain in front of the sphere (negative x-axis):
   Lx1 = 10.;
// Length of the domain behind the sphere (positive x-axis)
   Lx2 = 20.0;
// Half the width of the box (dimension in the y-direction)
   Ly = 5.;
// Half the height of the box (dimension in the z-direction)
   Lz = 5.;

// *******************************************************************************
// PART a) - DEFINE GEOMETRY OF A SPHERE
// *******************************************************************************

// Characteristic mesh size (lc) of points on the sphere
// The smaller is the lc, the finer is the mesh
lcSphere = .1 * RSphere;

// Characteristic mesh size in the domain around the sphere
lcDom = .7*RSphere;



// ============================================
// POINTS
// ============================================
// Each point of the geometry is defined by 3 coordinates (even for 2d mesh)
// and characteristic mesh size
// Point 1 is the center of the sphere
Point(1) = {    0,      0,     0,    lcSphere };

Point(2) = { RSphere,   0,     0,    lcSphere };
Point(3) = {    0,   RSphere,  0,    lcSphere };
Point(4) = { -RSphere,  0,     0,    lcSphere };
Point(5) = {    0,  -RSphere,  0,    lcSphere };
Point(6) = {    0,      0, -RSphere, lcSphere};
Point(7) = {    0,      0,  RSphere, lcSphere};



// ============================================
// LINES (straight lines, circles, splines ...)
// ============================================
// 'Circle' is any circular arc whose angle is smaller
// than pi. If you need an arc with angle bigger than
// pi, you have to compose it of several 'smaller' arcs
// Syntax for circle: Circle( CIRCLE_NR ) = { STARTING_PT, CENTER_PT, END_PT };
Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};

Circle(5) = {3,1,6};
Circle(6) = {6,1,5};
Circle(7) = {5,1,7};
Circle(8) = {7,1,3};
Circle(9) = {2,1,7};
Circle(10) = {7,1,4};
Circle(11) = {4,1,6};
Circle(12) = {6,1,2};



// ============================================
// LINE LOOPS:
// ============================================
// A line loop is a sequence of lines which is closed: 
// The starting point of a line in line loop has to be 
// the ending point of the preceding line in line loop
// If line has inverse orientation compared to what is needed, we list
// it with a '-' sign.
// Example: suppose we have points 1,2,3 and three arcs ('Circles')
// connecting these points such that they form a closed curve:
// Circle(10) goes from point 1 to point 2
// Circle(11) goes from point 3 to point 2
// Circle(12) goes from point 3 to point 1
// Then the line loop can be defined as 
// Line Loop(100) = {10,-11,12}

// Note that line loops are also types of lines, and therefore their numbers 
// must be distinct from the numbers of all other lines, arcs, splines ...
// Hence line loop numbering cannot start with 1, because there is already
// a circle with number 1
// On the other hand, point numbering is independent of line numbering
// Line numbering is independent of surface numbering, which does not 
// interfere in any way with volume numbering ...
// Numbering need not be continuous and need not start with 1
Line Loop(21) = {2,8,-10};
Line Loop(22) = {10,3,7};
Line Loop(23) = {-8,-9,1};
Line Loop(24) = {-11,-2,5};
Line Loop(25) = {-5,-12,-1};
Line Loop(26) = {-3,11,6};
Line Loop(27) = {-7,4,9};
Line Loop(28) = {-4,12,-6};


// Surface patches which create the surface
// of the sphere. 'Ruled' surface is a surface 
// such that there's a 1-to-1 mapping between
// this surface and a 2D plane
// There also exists a 'Plane Surface', which
// is a 'flat' 2D surface (possibly embedded in 3D space)


Ruled Surface(1) = {21};
Ruled Surface(2) = {22};
Ruled Surface(3) = {23};
Ruled Surface(4) = {24};
Ruled Surface(5) = {25};
Ruled Surface(6) = {26};
Ruled Surface(7) = {27};
Ruled Surface(8) = {28};

// Surface loop representing the sphere
Surface Loop(10) = {1,2,3,4,5,6,7,8};


// *******************************************************************************
// PART b) - BOX ENCLOSING THE SPHERE
// *******************************************************************************
// Numbering of all entities in this section starts with 101 so that it does not
// interfere with numbers of entities used to define the sphere in PART a)

Point(101) = { -Lx1, -Ly, -Lz, lcDom };
Point(102) = {  Lx2, -Ly, -Lz, lcDom };
Point(103) = {  Lx2,  Ly, -Lz, lcDom };
Point(104) = { -Lx1,  Ly, -Lz, lcDom };
Point(105) = { -Lx1, -Ly,  Lz, lcDom };
Point(106) = {  Lx2, -Ly,  Lz, lcDom };
Point(107) = {  Lx2,  Ly,  Lz, lcDom };
Point(108) = { -Lx1,  Ly,  Lz, lcDom };

// Lines forming the bottom surface of the box
Line(101) = {101,102};
Line(102) = {102,103};
Line(103) = {103,104};
Line(104) = {104,101};

// Lines forming the top surface of the box
Line(105) = {105,106};
Line(106) = {106,107};
Line(107) = {107,108};
Line(108) = {108,105};

// Vertical lines connecting bottom and top
Line(109) = {101,105};
Line(110) = {102,106};
Line(111) = {103,107};
Line(112) = {104,108};

// Line loops:
// Front face (inlet)
Line Loop(121) = {104,109,-108,-112};
// End of the domain (outlet)
Line Loop(122) = {102,111,-106,-110};
// Bottom of the domain
Line Loop(123) = {101,102,103,104};
// Top of the domain
Line Loop(124) = {105,106,107,108};
// Right face:
Line Loop(125) = {101,110,-105,-109};
// Left face:
Line Loop(126) = {103,112,-107,-111};


// Surfaces based on the previously defined line loops:
Plane Surface(101) = {121}; // Inlet
Plane Surface(102) = {122}; // Outlet
Plane Surface(103) = {123}; // Wall
Plane Surface(104) = {124}; // Wall
Plane Surface(105) = {125}; // Wall
Plane Surface(106) = {126}; // Wall


Surface Loop(110) = {101,102,103,104,105,106};

// Volume which represents the 'fluid' mesh
// The rhs contains the list of surface loops
// that form the volume: if only one number is given,
// then the volume is equal to the surface loop
// If multiple numbers are given, then the first 
// surface loop represents the volume and the following 
// surface loops represent holes in this volume. In our
// case, there is one hole in the domain: the sphere (surface loop nr. 10)
Volume(1) = {110,10};

// ******************************************************************
// Physical names - useful to mark different parts of the domain
// with different names so that the solver can apply different boundary 
// conditions to them

Physical Surface("Sphere") = {1,2,3,4,5,6,7,8}; // 1,...,8 are the numbers
                                                // of ruled surfaces which form
                                                // the wall of the sphere
                                                
Physical Surface("Inlet") = {101};
Physical Surface("Outlet") = {102};
Physical Surface("Wall") = {103,104,105,106};

Physical Volume("Fluid") = {1};

// *******************************************************************************
// PART c) - Mesh size fields for refinement around and behind the sphere
// *******************************************************************************

// First, we define an attractor (a point). In our case, the attractor will be the
// center of the sphere
Field[1] = Attractor;
Field[1].NodesList = {1}; // Node 1 is the center of the sphere

// Math eval - give a mathematical expression 
// which defines the mesh size
//Field[2] = MathEval;
//Field[2].F = "0.7";
////Field[2].F = "0.3*x";

// Based on attractor, we can define a threshold: when a distance
// from the attractor is smaller than DistMin, prescribe LcMin as 
// mesh element size
// When a distance is bigger that DistMax from the attractor, 
// prescribe LcMax as mesh element size
// Between DistMin and DistMax, interpolate the mesh size linearly
Field[3] = Threshold;
Field[3].IField  = 1; // Threshold is measured from field 1 - the attractor
Field[3].DistMin = 2.*RSphere;
Field[3].DistMax = 5.*RSphere;
Field[3].LcMin   = 2.*lcSphere;
Field[3].LcMax   = lcDom;
//Field[3].Sigmoid = 1;


// A box (given by 3 coordinates of its 2 diagonal corners)
// We prescribe the desired mesh size inside and outside the box
Field[4] = Box;
Field[4].XMin = -2.*RSphere;
Field[4].YMin = -2.*RSphere;
Field[4].ZMin = -2.*RSphere;
Field[4].XMax =  7.*RSphere;
Field[4].YMax =  2.*RSphere;
Field[4].ZMax =  2.*RSphere;
Field[4].VIn  =  2.*lcSphere;
Field[4].VOut =  2.*lcDom;

// A field which computes the mesh size as the minimum 
// of all mesh sizes prescribed by the previous field
Field[10] = Min;
Field[10].FieldsList = {3,4};

// Define which field should be used as 'background field' 
// Background field is the one on whose basis gmsh computes
// element sizes

// Only threshold:
//Background Field = 3;

// Only box:
//Background Field = 4;

// Minimum of mesh size prescribed by threshold
// and box at every point of the mesh:
Background Field = 10;



