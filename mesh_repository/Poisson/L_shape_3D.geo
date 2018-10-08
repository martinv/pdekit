/******************************
   L-shape in 3D
******************************/
lc = 0.050;
L1 = 2.0; // x
L2 = 2.0; // y
H  = 1.0; // z
nb_pts_on_side = 10;

Point(1) = {0, 0, 0, lc};
Point(2) = {L1, 0, 0, lc};
Point(3) = {L1, L2, 0, lc};
Point(4) = {0.5*L1, L2, 0, lc};
Point(5) = {0.5*L1, 0.5*L2, 0, lc};
Point(6) = {0, 0.5*L2, 0, lc};
Point(7) = {0, 0, H, lc};
Point(8) = {L1, 0, H, lc};
Point(9) = {L1, L2, H, lc};
Point(10) = {0.5*L1, L2, H, lc};
Point(11) = {0.5*L1, 0.5*L2, H, lc};
Point(12) = {0, 0.5*L2, H, lc};

// Bottom surface
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,1};

// Top surface
Line(7) = {7,8};
Line(8) = {8,9};
Line(9) = {9,10};
Line(10) = {10,11};
Line(11) = {11,12};
Line(12) = {12,7};

// Vertical lines connecting top and bottom
Line(13) = {1,7};
Line(14) = {2,8};
Line(15) = {3,9};
Line(16) = {4,10};
Line(17) = {5,11};
Line(18) = {6,12};

// Line loop for bottom surface
Line Loop(1) = {1:6};

// Line loop for top surface
Line Loop(2) = {7:12};

// Line loops on the sides
Line Loop(3) = {1,14,-7,-13};
Line Loop(4) = {2,15,-8,-14};
Line Loop(5) = {3,16,-9,-15};
Line Loop(6) = {4,17,-10,-16};
Line Loop(7) = {5,18,-11,-17};
Line Loop(8) = {6,13,-12,-18};

// Make surfaces from loops
For iSurf In {1:8}
  Plane Surface(iSurf) = {iSurf};
EndFor

// Define volume
Surface Loop(1) = {1:8};
Volume(1) = {1};

// Define physical entities
Physical Surface("Bottom") = {1};
Physical Surface("Top") = {2};
Physical Surface("Wall") = {3:8};

Physical Volume("InnerCells") = {1};

// Transfinite mesh - lines
For iLine In {1:2}
       Transfinite Line{iLine} = 2*nb_pts_on_side;
EndFor
For iLine In {7:8}
       Transfinite Line{iLine} = 2*nb_pts_on_side;
EndFor

For iLine In {3:6}
       Transfinite Line{iLine} = nb_pts_on_side;
EndFor
For iLine In {9:12}
       Transfinite Line{iLine} = nb_pts_on_side;
EndFor

For iLine In {13:18}
       Transfinite Line{iLine} = nb_pts_on_side;
EndFor


//Transfinite Surface(1) = {1,2,3,4};

