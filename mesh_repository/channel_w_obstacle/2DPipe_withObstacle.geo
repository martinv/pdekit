lc = 0.4;
R  = 2;   // Radius of the pipe
D  = 2*R;   // Diameter of the pipe
L  = 5*D; // Length of the pipe
h  = R + 0.2*R;  // Height of the obstacle
t = R/20; //Half tichkness of the obstacle

//****************************
// Nodes & Lines creation
//****************************
//Points of Channel geometry
Point(1) = {0, -R, 0, lc} ;
Point(2) = { L, -R, 0, lc} ;
Point(3) = { L,  R, 0, lc} ;
Point(4) = {0,  R, 0, lc} ;

//Points of Obstacle geometry
Point(5) = {(L/2 - t), -R, 0, lc} ;
Point(6) = {(L/2 + t), -R, 0, lc} ;
Point(7) = {(L/2 + t), (-R + h), 0, lc} ;
Point(8) = {(L/2 - t), (-R + h), 0, lc} ;

//Lines joining the points of the domain =  channel - the obstacle
Line(1) = {1,5} ;
Line(2) = {5,8} ;
Line(3) = {8,7} ;
Line(4) = {7,6} ;
Line(5) = {6,2} ;
Line(6) = {2,3} ;
Line(7) = {3,4} ;
Line(8) = {4,1} ;



//****************************
// Surfaces creation
//****************************

//Surface corresponding to the domain = channel- the obstacle
Line Loop(1) = {1,2,3,4,5,6,7,8};
Plane Surface(1) = {1};



//****************************
// Physical entites creation
//****************************


Physical Line("inlet")    = {8}; //Inlet
Physical Line("outlet")   = {6}; //Outlet
Physical Line("obstacle") = {2,3,4}; //Obstacle
Physical Line("bottom")   = {1,5}; //BottomWall
Physical Line("top")      = {7}; //TopWall

Physical Surface("interior") = {1}; //Domain = Channel - Obstacle



