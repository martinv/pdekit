#include<iostream>
#include<cstdlib>
#include<string>
#include<fstream>
#include<iomanip>
#include<cmath>
#include<vector>

int main(int argc, char *argv[]) {


  struct Node { double x; double y; };

  Node* Coords;
  Node* Normals;

  std::string inputFileName;
  std::string outputFileName;
  std::ifstream inFile;
  std::ofstream outFile;
  double nx,ny;
 

 const int idx_offset = 2000;
 const double dist = 0.03;
 unsigned int NPts;

 if(argc != 3) {
    std::cerr << "Wrong number of arguments.";
    std::cerr << "Usage: ./extrude input_file output_file\n";
    std::cerr << "Exiting.\n";
    exit(1);
  }
  else
  {
    inputFileName = std::string(argv[1]);
    outputFileName = std::string(argv[2]);
  }

 inFile.open(inputFileName.c_str());

 inFile >> NPts;

 Coords = new Node[NPts];
 Normals = new Node[2*NPts+1];

 for(uint i=0;i<NPts;++i) {
	inFile >> Coords[i].x;
	inFile >> Coords[i].y;

 }


int idx = 0;

while(Coords[idx].x > 1.0) {

 const double dx = Coords[idx+1].x - Coords[idx].x;
 const double dy = Coords[idx+1].y - Coords[idx].y;

 const double invnorm = 1.0/std::sqrt(dx*dx+dy*dy);

 //Compute normal:
 Normals[idx].x =  invnorm*dist*dy;
 Normals[idx].y = -invnorm*dist*dx;

// std::cout << "N = [" << Normals[idx].x << "," << Normals[idx].y << "]";
// std::cin.get();

 //Store the coordinate of the point moved in the normal direction:
 Normals[idx].x += Coords[idx].x;
 Normals[idx].y += Coords[idx].y;

 idx++;
}

while(idx < NPts) {

 const double x = Coords[idx].x;

 Normals[idx].x = -0.6*(0.14845/std::sqrt(x) - 0.126 - 0.7032*x + 0.8529*x*x - 0.406*x*x*x); 
 Normals[idx].y = -1.0;

 const double invnorm = 1.0/std::sqrt(Normals[idx].x*Normals[idx].x+Normals[idx].y*Normals[idx].y);

 //Compute normal:
 Normals[idx].x =  invnorm*dist*Normals[idx].x;
 Normals[idx].y = -invnorm*dist*Normals[idx].y;

 //Store the coordinate of the point moved in the normal direction:
 Normals[idx].x += Coords[idx].x;
 Normals[idx].y += Coords[idx].y;

 idx++;
}










 //Leading edge normal:
 Normals[idx].x = -dist;
 Normals[idx].y = 0.0;

 
 for(int i=0;i<NPts;++i) {
 Normals[2*NPts-i].x =  Normals[i].x;
 Normals[2*NPts-i].y = -Normals[i].y;
 }


 inFile.close(); 


 outFile.open(outputFileName.c_str()) ;
 outFile.precision(15);

 for(int idx = 0; idx < 2*NPts+1; idx++){
  
  outFile << "Point(" << idx_offset + idx << ") = {" << Normals[idx].x << ", " << Normals[idx].y;
  outFile << ", 0.0, lc};\n";



 }


 delete [] Coords;
 delete [] Normals;

return 0;
}
