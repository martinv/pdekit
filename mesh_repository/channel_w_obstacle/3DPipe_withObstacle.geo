lc1 = 0.4;
lc2 = 0.2;

R  = 2;   // Radius of the pipe
D  = 2*R;   // Diameter of the pipe
L  = 5*D; // Length of the pipe
h  = R + 0.2*R;  // Height of the obstacle
t = R/20; //Half tichkness of the obstacle


//****************************
// Nodes & Lines creation
//****************************

//Points on Inlet surface
Point(1) = {0, -R, -R, lc1} ;
Point(2) = {0, R, -R, lc1} ;
Point(3) = {0, R,  R, lc1} ;
Point(4) = {0,-R,  R, lc1} ;

//Points on Outlet surface
Point(10) = {L, -R, -R, lc1} ;
Point(20) = {L, R, -R, lc1} ;
Point(30) = {L, R,  R, lc1} ;
Point(40) = {L,-R,  R, lc1} ;

//Points on front face of the obstcale
Point(5) = {(L/2 - t),-t, -R, lc2} ;
Point(6) = {(L/2 - t), t, -R, lc2} ;
Point(7) = {(L/2 - t), t, -R+h, lc2} ;
Point(8) = {(L/2 - t),-t, -R+h, lc2} ;

//Points on back face of the obstcale
Point(50) = {(L/2 + t),-t, -R, lc2} ;
Point(60) = {(L/2 + t), t, -R, lc2} ;
Point(70) = {(L/2 + t), t, -R+h, lc2} ;
Point(80) = {(L/2 + t),-t, -R+h, lc2} ;

//Horizontale lines joining the Inlet and Outlet surfaces of the channel
Line(1) = {1,10} ;
Line(2) = {2,20} ;
Line(3) = {3,30} ;
Line(4) = {4,40} ;

//Horizontale lines joining the front and back surfaces of the obstacle
Line(5) = {5,50} ;
Line(6) = {6,60} ;
Line(7) = {7,70} ;
Line(8) = {8,80} ;

//Lines around the Inlet of the channel
Line(11) = {1,2} ;
Line(12) = {2,3} ;
Line(13) = {3,4} ;
Line(14) = {4,1} ;

//Lines around the Outlet of the channel
Line(21) = {10,20} ;
Line(22) = {20,30} ;
Line(23) = {30,40} ;
Line(24) = {40,10} ;

//Lines around the front and back faces of the obstacle
Line(15) = {5,6} ;
Line(16) = {6,7} ;
Line(17) = {7,8} ;
Line(18) = {8,5} ;
Line(25) = {50,60} ;
Line(26) = {60,70} ;
Line(27) = {70,80} ;
Line(28) = {80,50} ;
//Line(35) = {1,2} ;
//Line(36) = {2,3} ;
//Line(37) = {3,4} ;
//Line(38) = {4,1} ;


//****************************
// Surfaces creation
//****************************

//Around the channel
Line Loop(1) = {11,12,13,14}; //Inlet channel
Line Loop(2) = {21,22,23,24}; //Outlet channel

Line Loop(3) = {1, -24,-4, 14};
Line Loop(4) = {4, -23,-3, 13};
Line Loop(5) = {3, -22,-2, 12};
Line Loop(6) = {2, -21,-1, 11};

//Around the obstacle
Line loop(7) =  {15,16,17,18};
Line loop(8) =  {25,26,27,28};

Line Loop(9) = {5, -28,-8, 18};
Line Loop(10) = {8, -27,-7, 17};
Line Loop(11) = {7, -26,-6, 16};
Line Loop(12) = {6, -25,-5, 15};



//Inlet and Outlet surfaces of the channel
Plane Surface(1) = {1};
Plane Surface(2) = {2};

//Lateral (side) surfaces of the channel
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};
//Plane Surface(6) = {6};

//Inlet and Outlet surfaces of the obstacle
Plane Surface(7) = {7};
Plane Surface(8) = {8};

//Lateral (side) surfaces of the obstacle
Plane Surface(9) = {9};
Plane Surface(10) = {10};
Plane Surface(11) = {11};
Plane Surface(12) = {6,12};

//****************************
// Volume creation
//****************************

Surface Loop(13) = {1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12}; //Missing surf 6 which is the base of the obstacle
Volume(1) = {13};


//****************************
// Physical entites creation
//****************************

Physical Surface("inlet") = {1};  //Inlet
Physical Surface("outlet") = {2}; //Outlet
Physical Surface("obstacle") = {7,8,9,10,11}; //Walls of the obstacle
Physical Surface("walls") = {3,4,5,12}; //Walls of the channel

Physical Volume("interior") = {1}; //Volume



