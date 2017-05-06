#include "stlcut.h"
using namespace std;
int main(int argc, char **argv) 
{
  if ((argc != 2 && argc!=6 ) ||argv[1]==string("help")) {
    
    //if(argc<2)
     // {
        cerr << "Usage: " << argv[0] << " file.stl a b c d (plane-optional or 0 0 1 1 will be used) " << endl<<"or"<<endl<< argv[0]<<" tests"<<endl;
        cerr<< "Plane is in a form of a plane normal (a, b, c) and d = distance from origin."<<endl;
      return 1;
  //}

}

   stl_plane plane=stl_plane(0,0,1,-1);
  if(argc==6)
  {
    plane=stl_plane(atof( argv[2]),atof(argv[3]),atof(argv[4]),atof(argv[5]));
  }
 
  Mesh mesh;
  mesh.open(argv[1]);
  mesh.cut(plane); 
                              
  if(mesh.createBorderPolylines())
    {
      mesh.findHoles();
      mesh.triangulateCut();
      mesh.save();
    }

  mesh.close();

  

return 0;
}