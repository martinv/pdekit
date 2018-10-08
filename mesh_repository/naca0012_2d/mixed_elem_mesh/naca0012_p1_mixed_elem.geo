/********************************************************************* 
 *
 *  Geometry for mesh of NACA0012 airfoil
 * 
 *  Points of airfoil surface were defined in included file
 *  'points.geo'
 *  Here we define farfield and characteristic lenghts
 *  in the domain
 *
 *********************************************************************/

Geometry.Tolerance = 1e-10;

Include "points.geo";
Include "offset1.geo";
Include "offset2.geo";


PI = 3.14159265358979323846;

///Set your number of farfield points. This number should be EVEN
nb_farfield_pts = 20;

//characteristic length of farfield points
lc_ff = 50.0;

R = 50; //Radius of farfield
dphi = 2*PI/nb_farfield_pts;

//Points defining farfield:
For iFPt In {1:nb_farfield_pts}
      Point(nb_airfoil_pts+iFPt) = {R*Cos((iFPt-1)*dphi), R*Sin((iFPt-1)*dphi), 0.0, lc_ff};
EndFor

//Airfoil surface:
//Spline(1) = {le_idx:te_idx:-1};
//Spline(2) = {1,nb_airfoil_pts:le_idx:-1};
idx_patch1 = 120;
idx_patch2 = 200;
idx_patch3 = 381;
idx_patch4 = 562;
idx_patch5 = 642;

Spline(1) = {1:idx_patch1};
Spline(2) = {idx_patch1:idx_patch2};
Spline(3) = {idx_patch2:idx_patch3};
Spline(4) = {idx_patch3:idx_patch4};
Spline(5) = {idx_patch4:idx_patch5};
Spline(6) = {idx_patch5:nb_airfoil_pts,1};


idx_offset1 = 999;
idx_offset2 = 1999;

Line(7)  = {te_idx,te_idx+idx_offset1};
Line(8)  = {idx_patch1,idx_patch1+idx_offset1};
Line(9)  = {idx_patch2,idx_patch2+idx_offset1};
Line(10) = {idx_patch3,idx_patch3+idx_offset1};
Line(11) = {idx_patch4,idx_patch4+idx_offset1};
Line(12) = {idx_patch5,idx_patch5+idx_offset1};
Line(13) = {te_idx,te_idx+nb_airfoil_pts+idx_offset1};

Spline(14) = {idx_patch1+idx_offset1:te_idx+idx_offset1:-1};
Spline(15) = {idx_patch2+idx_offset1:idx_patch1+idx_offset1:-1};
Spline(16) = {idx_patch3+idx_offset1:idx_patch2+idx_offset1:-1};
Spline(17) = {idx_patch4+idx_offset1:idx_patch3+idx_offset1:-1};
Spline(18) = {idx_patch5+idx_offset1:idx_patch4+idx_offset1:-1};
Spline(19) = {te_idx+nb_airfoil_pts+idx_offset1:idx_patch5+idx_offset1:-1};


Line(20) = {te_idx+idx_offset1,te_idx+idx_offset2};
Line(21) = {idx_patch1+idx_offset1,idx_patch1+idx_offset2};
Line(22) = {idx_patch2+idx_offset1,idx_patch2+idx_offset2};
Line(23) = {idx_patch3+idx_offset1,idx_patch3+idx_offset2};
Line(24) = {idx_patch4+idx_offset1,idx_patch4+idx_offset2};
Line(25) = {idx_patch5+idx_offset1,idx_patch5+idx_offset2};
Line(26) = {te_idx+nb_airfoil_pts+idx_offset1,te_idx+nb_airfoil_pts+idx_offset2};


Spline(27) = {te_idx+idx_offset2:idx_patch1+idx_offset2};
Spline(28) = {idx_patch1+idx_offset2: idx_patch2+idx_offset2};
Spline(29) = {idx_patch2+idx_offset2:idx_patch3+idx_offset2};
Spline(30) = {idx_patch3+idx_offset2:idx_patch4+idx_offset2};
Spline(31) = {idx_patch4+idx_offset2:idx_patch5+idx_offset2};
Spline(32) = {idx_patch5+idx_offset2:te_idx+nb_airfoil_pts+idx_offset2};


//Farfield:
Spline(50) = {nb_airfoil_pts+1:nb_airfoil_pts+nb_farfield_pts, nb_airfoil_pts+1};

//Physical Line( "UpperAirfoil" ) = {1,2,3};
//Physical Line( "LowerAirfoil" ) = {4,5,6};
Physical Line( "Wall" ) = {1:6};
Physical Line( "Farfield" ) = {50};

//Wall of the airfoil:
//Line Loop(0) = {1:6};

//Viscous layer:
//Line Loop(1) = {1,8,14,-7};
Line Loop(1) = {-1,7,-14,-8};
Plane Surface(1) = {1};

//Line Loop(2) = {2,9,15,-8};
Line Loop(2) = {-2,8,-15,-9};
Plane Surface(2) = {2};

//Line Loop(3) = {3,10,16,-9};
Line Loop(3) = {-3,9,-16,-10};
Plane Surface(3) = {3};

//Line Loop(4) = {4,11,17,-10};
Line Loop(4) = {-4,10,-17,-11};
Plane Surface(4) = {4};

//Line Loop(5) = {5,12,18,-11};
Line Loop(5) = {-5,11,-18,-12};
Plane Surface(5) = {5};

//Line Loop(6) = {6,13,19,-12};
Line Loop(6) = {-6,12,-19,-13};
Plane Surface(6) = {6};

Line Loop(7) = {14,20,27,-21};
Plane Surface(7) = {7};

Line Loop(8) = {15,21,28,-22};
Plane Surface(8) = {8};

Line Loop(9) = {16,22,29,-23};
Plane Surface(9) = {9};

Line Loop(10) = {17,23,30,-24};
Plane Surface(10) = {10};

Line Loop(11) = {18,24,31,-25};
Plane Surface(11) = {11};

Line Loop(12) = {19,25,32,-26};
Plane Surface(12) = {12};

//The rest of the domain is meshed by unstructured mesh:

	//Boundary layer envelope:
	Line Loop(13) = {27:32,-26,-13,7,20};

	//Farfield:
	Line Loop(14) = {50};

	Plane Surface(14) = {14,13};

Physical Surface( "InnerCells" ) = {1:12,14};



nb_pts_patch_te = 14;
nb_pts_patch_mid = 40;
nb_pts_patch_le = 10;

ratio_patch_te  = 1.2;
ratio_patch_mid = 1.05;
ratio_patch_le  = 1.2;


//SUCTION SIDE:
//Trailing edge
Transfinite Line{1,27} = nb_pts_patch_te Using Progression 1.2;
Transfinite Line{14}   = nb_pts_patch_te Using Progression 1.0/1.2;
//Mid airfoil
Transfinite Line{2,28} = nb_pts_patch_mid Using Progression 1.0/1.05;
Transfinite Line{15}   = nb_pts_patch_mid Using Progression 1.05;
//Leading edge
Transfinite Line{3} = nb_pts_patch_le Using Progression 1.0/ratio_patch_le;
Transfinite Line{16}   = nb_pts_patch_le Using Progression 0.95*ratio_patch_le;
Transfinite Line{29} = nb_pts_patch_le Using Progression 1.0/(0.90*ratio_patch_le);

//PRESSURE SIDE
//Leading edge
Transfinite Line{4} = nb_pts_patch_le Using Progression ratio_patch_le;
Transfinite Line{17}   = nb_pts_patch_le Using Progression 1.0/(0.95*ratio_patch_le);
Transfinite Line{30} = nb_pts_patch_le Using Progression 0.9*ratio_patch_le;
//Mid airfoil
Transfinite Line{5,31} = nb_pts_patch_mid Using Progression 1.05;
Transfinite Line{18}   = nb_pts_patch_mid Using Progression 1.0/1.05;
//Trailing edge
Transfinite Line{6,32} = nb_pts_patch_te Using Progression 1.0/1.2;
Transfinite Line{19}   = nb_pts_patch_te Using Progression 1.2;

nb_layers1 = 3;
nb_layers2 = 4;

Transfinite Line {7:13} = nb_layers1+1;
Transfinite Line {20:26} = nb_layers2+1;

//Farfield
Transfinite Line{50} = 20;


Transfinite Surface {1}  = { } Left;
Transfinite Surface {2}  = { } Left;
Transfinite Surface {3}  = { } Right;
Transfinite Surface {4}  = { } Left;
Transfinite Surface {5}  = { } Right;
Transfinite Surface {6}  = { } Right;
Transfinite Surface {7}  = { } Left;
Transfinite Surface {8}  = { } Left;
Transfinite Surface {9}  = { } Right;
Transfinite Surface {10} = { } Left;
Transfinite Surface {11} = { } Right;
Transfinite Surface {12} = { } Right;

For isurface In {1:12}
  //Transfinite Surface{isurface};
   Recombine Surface(isurface);
EndFor



//Leading edge:
Field[0] = Attractor;
Field[0].NodesList = { le_idx+idx_offset2 };

//Leading edge:
Field[1] = Attractor;
Field[1].NodesList = { te_idx };

//Suction and pressure and sides of the airfoil:
Field[2] = Attractor;
Field[2].EdgesList = { 1, 2 };
// Field[2].NodesList = { 6000 };

//Linearly increasing element size through the whole domain:
Field[3] = MathEval;
Field[3].F = "0.01 + 0.55*(x^2+y^2)^(0.5)*(20-0.01)/30";
Field[4] = MathEval;
Field[4].F = "0.01 + 0.55*((x-1)^2+y^2)^(0.5)*(20-0.01)/30";

//Refinement around leading edge:
Field[5] = MathEval;
Field[5].F = "0.005 + 1.9*(x^2+y^2)^1.1";
//Refinement around trailing edge:
Field[6] = MathEval;
Field[6].F = "0.007 + 1.9*((x-1.0)^2+y^2)^1.1";

// Local refinement near the airfoil:
Field[7] = Threshold;
Field[7].IField = 2;
Field[7].DistMin = 0.2;
Field[7].DistMax = 0.5;
Field[7].LcMin = 0.06;
Field[7].LcMax = 30.0;
//Field[7].Sigmoid = 1;

/*
Field[8] = Box;
Field[8].XMin = 1.0;
Field[8].XMax = 1.4;
Field[8].YMin = -0.03;
Field[8].YMax = 0.03;
Field[8].VIn = 0.008;
Field[8].VOut = 10.0;
*/

//Refinement in wake area
Point(10000) = {1.4,0.0,0.0,lc};
Line(10000) = {te_idx,10000};

Field[8] = Attractor;
Field[8].EdgesList = { 10000 };

Field[9] = Threshold;
Field[9].IField = 8;
Field[9].DistMin = 0.03;
Field[9].DistMax = 0.1;
Field[9].LcMin = 0.006;
Field[9].LcMax = 15.0;
//Field[9].Sigmoid = 1;





Field[10] = Min;
Field[10].FieldsList = { 3,4,5,6,7 };
Background Field = 10;



//Transfinite Line{1,2} = 80 Using Bump 0.01;




