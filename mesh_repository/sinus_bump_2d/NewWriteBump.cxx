#include<iostream>
#include<iomanip>
#include<fstream>
#include<cmath>

// #include<cstdlib>



int main ()
{
const double PI = 3.14159265358979323846;
std::ofstream outFile;
const int nPts = 200;

const double dx = 2.0/(nPts-1);


outFile.precision(10);
outFile.setf(std::ios_base::scientific);
outFile.open("bump.geo");
outFile << "lc = 20;\n\n";
outFile << "Point(1) = {0.0, 0.0, 0.0, lc};\n";
outFile << "Point(2) = {0.25, 0.0, 0.0, lc};\n";
outFile << "Point(3) = {0.50, 0.0, 0.0, lc};\n";
outFile << "Point(4) = {0.75, 0.0, 0.0, lc};\n";
outFile << "Point(5) = {0.98, 0.0, 0.0, lc};\n";
outFile << "Point(6) = {0.99, 0.0, 0.0, lc};\n";

for(int i=0;i<nPts;++i) {

const double x = 1.0 + i*dx;
const double y = 0.1*std::sin(PI*(x-1.5))+0.1;

outFile << "Point(" << 7+i << ") = {" << x << ", " << y << ", 0.0, lc};\n";
}

outFile << "Point(" << nPts + 7 << ") = {3.01, 0.0, 0.0, lc};\n";
outFile << "Point(" << nPts + 8 << ") = {3.02, 0.0, 0.0, lc};\n";

outFile << "Point(" << nPts + 9 << ") = {3.25, 0.0, 0.0, lc};\n";
outFile << "Point(" << nPts + 10 << ") = {3.50, 0.0, 0.0, lc};\n";
outFile << "Point(" << nPts + 11 << ") = {3.75, 0.0, 0.0, lc};\n";
outFile << "Point(" << nPts + 12 << ") = {4.0, 0.0, 0.0, lc};\n";


outFile << "Point(" << nPts + 13 << ") = {4.0, 1.0, 0.0, lc};\n";
outFile << "Point(" << nPts + 14 << ") = {0.0, 1.0, 0.0, lc};\n";

outFile << "\n";

outFile << "Spline(1) = {1:" << 12+nPts << "};\n";
outFile << "Line(2) = {" << 12+nPts << "," << 13+nPts << "};\n";
outFile << "Line(3) = {" << 13+nPts << "," << 14+nPts << "};\n";
outFile << "Line(4) = {" << 14+nPts << ",1};\n";

outFile << "\n";

outFile << "Physical Line( \"LowerWall\" ) = { 1 };\n";
outFile << "Physical Line( \"SubOutlet\" ) = { 2 };\n";
outFile << "Physical Line( \"UpperWall\" ) = { 3 };\n";
outFile << "Physical Line( \"SubInlet\" ) = { 4 };\n";

outFile << "\n";

outFile << "Line Loop(1) = {1,2,3,4};\n";

outFile << "\n";

outFile << "Plane Surface(1) = {1};\n";
outFile << "Physical Surface( \"InnerCells\" ) = {1};\n";

outFile << "\n";

outFile << "pts_hor = 15; //number of points on horizontal wall\n";
outFile << "pts_vert = 4; //number of points on vertical wall\n\n";


outFile << "Transfinite Line{1} = pts_hor;\n";
outFile << "Transfinite Line{2} = pts_vert;\n";
outFile << "Transfinite Line{3} = pts_hor Using Bump 0.8;\n\n";
outFile << "Transfinite Line{4} = pts_vert;\n";

outFile << "Transfinite Surface { 1 } = {  } Right;\n";

outFile.close();

return 0;
}
