/********************************************************************* 
 *
 *  Geometry for mesh of NACA0012 airfoil
 * 
 *  Points of airfoil surface were defined in included file
 *  'points.geo'
 *  Here we define farfield and characteristic lenghts
 *  in the domain
 *  Generate with clscale = 2.1 and del2d algorithm
 *
 *********************************************************************/

//Original value for Adigma meshes: lc = 0.023
lc = 0.023;

Include "points.geo";

PI = 3.14159265358979323846;

///Set your number of farfield points. This number should be EVEN
nb_farfield_pts = 15;


R = 50; //Radius of farfield
dphi = 2*PI/nb_farfield_pts;

//characteristic length of points close to the airfoil
lc_af = 0.030;

//characteristic length of farfield points
lc_ff = 2*PI*R/nb_farfield_pts;

//Points defining farfield:
For iFPt In {1:nb_farfield_pts}
      Point(nb_airfoil_pts+iFPt) = {R*Cos((iFPt-1)*dphi), R*Sin((iFPt-1)*dphi), 0.0, lc_ff};
EndFor

//Airfoil surface:
Spline(1) = {le_idx:te_idx:-1};
Spline(2) = {1,nb_airfoil_pts:le_idx:-1};

//Physical Line( "UpperAirfoil" ) = {1};
//Physical Line( "LowerAirfoil" ) = {2};
Physical Line( "Wall" ) = {1,2};

//Farfield:
Spline(3) = {nb_airfoil_pts+1:nb_airfoil_pts+nb_farfield_pts, nb_airfoil_pts+1};
//Physical Line( "FField" ) = {3};
Physical Line( "Farfield" ) = {3};

Line Loop(1) = {3};
Line Loop(2) = {1,2};

Plane Surface(1) = {1,2};
Physical Surface( "InnerCells" ) = {1};

//Refinement around the airfoil:
//Point(5000) = {0.0, 0.0, 0.0, lc_af};
//Point(5001) = {1.0, 0.0, 0.0, lc_af};
// Point(6000) = {0.5, 0.0, 0.0, 0.05};

//Leading edge:
Field[0] = Attractor;
Field[0].NodesList = { le_idx };

//Trailing edge:
Field[1] = Attractor;
Field[1].NodesList = { te_idx };

//Suction and pressure and sides of the airfoil:
Field[2] = Attractor;
Field[2].EdgesList = { 1, 2 };

//Linearly increasing element size through the whole domain:
Field[3] = MathEval;
Field[3].F = "0.01 + 0.50*((x^2+y^2)^0.65)*(20.0-0.01)/30";
Field[4] = MathEval;
Field[4].F = "0.01 + 0.50*(((x-1)^2+y^2))^0.65*(20.0-0.01)/30";

//Refinement around leading edge:
Field[5] = MathEval;
Field[5].F = "0.007 + 1.9*(x^2+y^2)^1.1";

//Refinement around trailing edge:
Field[6] = MathEval;
Field[6].F = "0.010 + 1.9*((x-1)^2+y^2)^1.1";

// Local refinement near the airfoil:
Field[7] = Threshold;
Field[7].IField = 2;
//Field[7].DistMin = 0.25;
Field[7].DistMin = 0.2;
Field[7].DistMax = 5.0;
Field[7].LcMin = lc_af;
Field[7].LcMax = lc_ff;
Field[7].Sigmoid = 1;
//Don't impose anything outside DistMax:
Field[7].StopAtDistMax = 1;

Field[8] = Box;
Field[8].VIn = 3.0*lc_af;
Field[8].VOut = lc_ff;
Field[8].XMin = 0.0;
Field[8].XMax = 1.2;
Field[8].YMin = -0.7;
Field[8].YMax = 0.7;

Field[10] = Min;
Field[10].FieldsList = { 3,4,5,6,7,8 };
Background Field = 10;

// Transfinite Line{1,2} = 30 Using Bump 0.1;
Transfinite Line{3} = nb_farfield_pts;
