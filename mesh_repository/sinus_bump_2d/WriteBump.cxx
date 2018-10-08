#include<iostream>
#include<iomanip>
#include<fstream>
#include<cmath>

// #include<cstdlib>



int main ()
{
const double PI = 3.14159265358979323846;
std::ofstream outFile;
const int nPts = 100;

const double dx = 2.0/(nPts-1);


outFile.precision(10);
outFile.setf(std::ios_base::scientific);
outFile.open("bump.geo");
outFile << "lc = 20;\n\n";
outFile << "Point(1) = {0.0, 0.0, 0.0, lc};\n";
outFile << "Point(2) = {4.0, 0.0, 0.0, lc};\n";
outFile << "Point(3) = {4.0, 1.0, 0.0, lc};\n";
outFile << "Point(4) = {0.0, 1.0, 0.0, lc};\n";

for(int i=0;i<nPts;++i) {

const double x = 1.0 + i*dx;
const double y = 0.1*std::sin(PI*(x-1.5))+0.1;

outFile << "Point(" << 5+i << ") = {" << x << ", " << y << ", 0.0, lc};\n";
}

outFile << "\n";

outFile << "Line(1) = {1,5};\n";
outFile << "Spline(2) = {5:" << 4+nPts << "};\n";
outFile << "Line(3) = {" << 4+nPts << ",2};\n";
outFile << "Line(4) = {2,3};\n";
outFile << "Line(5) = {3,4};\n";
outFile << "Line(6) = {4,1};\n";

outFile << "\n";

outFile << "Physical Line( \"LowerWall\" ) = { 1, 2, 3 };\n";
outFile << "Physical Line( \"SubOutlet\" ) = { 4 };\n";
outFile << "Physical Line( \"UpperWall\" ) = { 5 };\n";
outFile << "Physical Line( \"SubInlet\" ) = { 6 };\n";

outFile << "\n";

outFile << "Line Loop(1) = {1,2,3,4,5,6};\n";

outFile << "\n";

outFile << "Plane Surface(1) = {1};\n";
outFile << "Physical Surface( \"InnerCells\" ) = {1};\n";

outFile << "\n";

outFile << "pts_hor1 = 4; //straight line segments on bottom wall\n";
outFile << "pts_hor2 = 7; //bump\n";
outFile << "pts_vert = 4;\n\n";


outFile << "Transfinite Line{1,3} = pts_hor1;\n";
outFile << "Transfinite Line{2} = pts_hor2;\n";
outFile << "Transfinite Line{5} = 2*pts_hor1 + pts_hor2 Using Bump 0.8;\n\n";

outFile << "Transfinite Line{4} = pts_vert;\n";
outFile << "Transfinite Line{6} = pts_vert;\n";

outFile.close();

return 0;
}
