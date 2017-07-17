#include "stlcut.h"
#include <signal.h>
#include <setjmp.h>
using namespace std;
jmp_buf buf;
sigset_t signal_set;
int numOfSegv=0;

void segv_handler(int s)
{
    switch(s)
    {

        case SIGSEGV:
        printf("\nSegmentation fault signal caught! Attempting recovery..");
        numOfSegv++;
        longjmp(buf, numOfSegv);
        break;
    }
}
void meshUnitTests()
{
  Mesh mesh;
    bool result = mesh.runUnitTests();
    if(result)
      cout<<endl<<"All tests completed without a problem."<<endl;
    else
      cout<<endl<<"Some tests failed."<<endl;
}

void acquireSaveName(string& name)
{
  name = "";
  cout<<"Please enter the new name of meshes ( \"_1.stl\" and \"_2.stl\" will be added)"<<endl;
  cout<<"Or pres Enter for default name (Cut_Mesh_1/2.stl)"<<endl;
  while(true)
  {
    getline(cin,name);//cin>>name;
    if(! (find_if(name.begin(), name.end(), 
        [](char c) { return !(isalnum(c) || (c == ' ') || (c == '_' )); }) == name.end() ) )
    {
      cout<<"Plese use only alphanumeric characters, space and underscore"<<endl;
      name="";
      continue;
    }
    else
      return;
  }
}



int main(int argc, char **argv) 
{
  double error_correction = 0.00015;
  string name = "";
  if ((argc != 2 && argc != 6 ) || argv[1]==string(" --help ")) 
  {
    
    //if(argc<2)
     // {
        cout << "Usage: " << argv[0] << " file.stl a b c d (plane-optional or 0 0 1 1 will be used) " << endl<<"or"<<endl<< argv[0]<<" tests"<<endl;
        cout<< "Plane is in a form of a plane normal (a, b, c) and d = distance from origin."<<endl;
        cout<< "Cutting through non-manifold object might result in error or wrong triangulation of cutted area."<<endl;
      return 1;
  //}

  }
  if(argv[1] == string("tests"))
  {
    meshUnitTests();
    //meshIntegrationTests();
    return 0;
  }

  stl_plane plane = stl_plane(0, 0, 1, -1);

  if(argc == 6)
  {
    plane = stl_plane(atof( argv[2]), atof(argv[3]), atof(argv[4]), atof(argv[5]));
  }

  
  //Setting jump in case of poly2tri segfault
  setjmp(buf);
  sigemptyset(&signal_set);
  sigaddset(&signal_set, SIGSEGV); 
  sigprocmask(SIG_UNBLOCK, &signal_set, NULL); //clearnig segfault signal
  signal(SIGSEGV, segv_handler);

  if(numOfSegv > 6)
  {
     cerr<<"STLCUT wasnt able to made this cut. Try changing the plane position slightly and make sure that your model is 2-manifold."<<endl;
    return 1;
  }
  else
  {
    Mesh mesh;
    mesh.openStl(argv[1]);

    if(numOfSegv > 0) 
    {
      plane = stl_plane(atof( argv[2]),atof(argv[3]),atof(argv[4]),(-1)*plane.d + error_correction);
      error_correction=error_correction>0?(-1)*error_correction:error_correction*10;//(-1)*error_correction * 10;
      cout<<endl<<"Recovered from segmentation fault."<<endl;
    }

    if (mesh.cut(plane))
    {
      acquireSaveName(name);
      mesh.save(name);
    }
  }
  return 0;
}
