Include "bump_contour.geo";

Mesh.Algorithm = 2; //1=MeshAdapt, 2=Automatic, 5=Delaunay, 6=Frontal, 7=bamg, 8=delquad
Mesh.Algorithm3D = 1; //1=Delaunay, 4=Frontal, 5=Frontal Delaunay, 6=Frontal Hex, 7=MMG3D
Mesh.CharacteristicLengthFromPoints = 1;
Mesh.CharacteristicLengthExtendFromBoundary = 1;

// =============================================================
// SURFACES
// =============================================================

// Left wall
Line Loop(1) = {1:6};
Plane Surface(1) = {1};

// Extrusion of the domain contour
//out[] = Extrude{0,0,2*W}{ Line{1:4}; };
out[] = Extrude{0,0,2*W}{ Surface{1}; };

// Remark: the extrusion directly creates a volume, because we
// extruded a surface !!!

dist = 1.5;

Field[1] = Box;
Field[1].XMin = 2.0 - dist;
Field[1].XMax = 2.0 + dist;
Field[1].YMin = 0.0;
Field[1].YMax = 0.4;
Field[1].ZMin = -W;
Field[1].ZMax =  W;
Field[1].VIn  = 0.5*lc;
Field[1].VOut = 0.8*lc;

Field[2] = Min;
Field[2].FieldsList = {1};

Background Field = 2;

// =============================================================
// PHYSICAL ENTITIES
// =============================================================
Physical Surface("Inlet") = {37};
Physical Surface("Outlet") = {29};
Physical Surface("SymmetryLeft") = {1};
Physical Surface("SymmetryRight") = {38};
Physical Surface("Bottom") = {17,21,25};
Physical Surface("Top") = {33};

Physical Volume("Fluid") = {1};

