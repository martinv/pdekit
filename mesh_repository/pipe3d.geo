//-*- C++ -*-

lc = 0.0001; //0.001
h = 0.001;// 0.004
R= 0.02;
L = 0.2;


// ############# POINTS ##############
//Base
Point(0) = {0, 0, 0, lc};
Point(1) = {R, 0, 0, lc};
Point(2) = {0, R, 0, lc};
Point(3) = {-R, 0, 0, lc};
Point(4) = {0, -R, 0, lc};
Point(5) = {R+h, 0, 0, lc};
Point(6) = {0, R+h, 0, lc};
Point(7) = {-R-h, 0, 0, lc};
Point(8) = {0, -R-h, 0, lc};


//############### LINES ##############
//Base
Circle(1) = {1,0,2};
Circle(2) = {2,0,3};
Circle(3) = {3,0,4};
Circle(4) = {4,0,1};
Circle(5) = {5,0,6};
Circle(6) = {6,0,7};
Circle(7) = {7,0,8};
Circle(8) = {8,0,5};

Line (9) = {1,5};
Line (10) = {2,6};
Line (11) = {3,7};
Line (12) = {4,8};

// ############# LINE LOOPS ############

Line Loop(1) = {-1,9,5,-10}; //on lines
Line Loop(2) = {-2,10,6,-11};
Line Loop(3) = {-3,11,7,-12};
Line Loop(4) = {-4,12,8,-9};
Line Loop(5) = {1,2,3,4};

// ############ SURFACES ############

Plane Surface(1) = {1}; //on line loops
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};

// ############ TRANSFINITE LINES ############

Transfinite Line{1,2,3,4,5,6,7,8}=40;//20 //on lines or circle
Transfinite Line{9,10,11,12}=10; //5
//Transfinite Line{2,4}=5;

// ############ TRANSFINITE SURFACES ############

Transfinite Surface { 1 } = { 1,5,6,2} ; //on points
Recombine Surface { 1 } ;
Transfinite Surface { 2 } = {2,6,7,3} ;
Recombine Surface { 2 } ;
Transfinite Surface { 3 } = {3,7,8,4} ;
Recombine Surface { 3 } ;
Transfinite Surface { 4 } = {4,8,5,1} ;
Recombine Surface { 4 } ;

//############# EXTRUSION TO 3D ###############

z = 10*R;
  nb_layers = 60; //30

///This loop automatically generates 4 volumes whose indexes should be 1,2,3,4
///These volumes will create physical volume 'PipeWall'
For i In {1:4}
 out[] = Extrude{0, 0, z}{Surface{i}; Layers{nb_layers}; Recombine; }; 
 //Create temporary arrays to store indexes of extruded surfaces
 IndexWallFaceOutlet[i-1] = out[0];
 IndexInnerWall[i-1] = out[2];
 IndexOuterWall[i-1] = out[4];
EndFor

///This loop generates one volume whose number is 5 and which is the inside of the pipe
//This command gives error
//'Cannot tetrahedralize volume with quadrangles' :
// out[] = Extrude{0, 0, z}{Surface{5};}; 

//Prismatic elements inside the pipe
out[] = Extrude{0, 0, z}{Surface{5}; Layers{nb_layers}; Recombine; }; 

IndexOutlet = out[0];

//####### PHYSICAL SURFACES AND VOLUMES #########
Physical Surface("WallOnInlet") = {1,2,3,4};
Physical Surface("WallOnOutlet")= {IndexWallFaceOutlet[0],IndexWallFaceOutlet[1],
                                   IndexWallFaceOutlet[2],IndexWallFaceOutlet[3]};
Physical Surface("Inlet") = {5};
Physical Surface("Outlet") = {IndexOutlet};
Physical Surface("InnerWall") = {IndexInnerWall[0],IndexInnerWall[1],
                                 IndexInnerWall[2],IndexInnerWall[3]};
Physical Surface("OuterWall") = {IndexOuterWall[0],IndexOuterWall[1],
                                 IndexOuterWall[2],IndexOuterWall[3]};
Physical Volume("PipeWall") = {1,2,3,4};
Physical Volume("PipeInside") = {5};
