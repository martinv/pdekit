/******************************
          Rectangle
******************************/
lc = 0.050;
L = 2.0; // x
H = 2.0; // y
nb_pts_on_side = 10;

Point(1) = {0, 0, 0, lc};
Point(2) = {L, 0, 0, lc};
Point(3) = {L, H, 0, lc};
Point(4) = {0.5*L, H, 0, lc};
Point(5) = {0.5*L, 0.5*H, 0, lc};
Point(6) = {0, 0.5*H, 0, lc};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,1};

Line Loop(1) = {1:6};
Plane Surface(1) = {1};

Physical Line("Wall1") = {1};
Physical Line("Wall2") = {2};
Physical Line("Wall3") = {3};
Physical Line("Wall4") = {4};
Physical Line("Wall5") = {5};
Physical Line("Wall6") = {6};
Physical Surface("InnerCells") = {1};


For iLine In {1:2}
       Transfinite Line{iLine} = 2*nb_pts_on_side;
EndFor

For iLine In {3:6}
       Transfinite Line{iLine} = nb_pts_on_side;
EndFor

//Transfinite Surface(1) = {1,2,3,4};

