Mesh.RecombineAll=1;
Mesh.Algorithm=8; // delquad

lc = 0.10;
L = 2.0;
H = 2.0;

nb_pts_hor = 20;
nb_pts_vert =20;


Point(1) = {0.0,0.0,0.0,lc};
Point(2) = { L ,0.0,0.0,lc};
Point(3) = { L,  H ,0.0,lc};
Point(4) = {0.0, H ,0.0,lc};
Point(5) = {-L , H, 0.0,lc};
Point(6) = {-L ,0.0,0.0,lc};


Line(1) = {6,1};
Line(2) = {1,2};
Line(3) = {2,3};
Line(4) = {3,4};
Line(5) = {4,5};
Line(6) = {5,6};

Line Loop(1) = {1,2,3,4,5,6};

Plane Surface(1) = {1};


Step = 1.008;
InvStep = 1.0/Step;


Transfinite Line(1) = nb_pts_hor + 1 Using Progression Step;
Transfinite Line(2) = nb_pts_hor + 1 Using Progression InvStep;
Transfinite Line(4) = nb_pts_hor + 1 Using Progression InvStep;
Transfinite Line(5) = nb_pts_hor + 1 Using Progression Step;
Transfinite Line(3) = nb_pts_vert Using Progression Step;
Transfinite Line(6) = nb_pts_vert Using Progression InvStep;


//Transfinite Surface(1) = {6,1,4,5};
//Transfinite Surface(2) = {1,2,3,4};

//Recombine Surface(1);
//Recombine Surface(2);

Physical Line("inlet") = {1}; 
Physical Line("outlet") = {2}; 
Physical Line("farfield") = {4,6}; 
Physical Surface("fluid") = {1};
