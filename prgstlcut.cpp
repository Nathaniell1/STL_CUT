#include "stlcut.h"
#include <signal.h>
#include <setjmp.h>
using namespace std;

jmp_buf buf;
sigset_t signal_set;
int num_of_segv=0;

void segv_handler(int s)
{

    switch(s)
    {

        case SIGSEGV:
        printf("\nSegmentation fault signal caught! Attempting recovery..");
        num_of_segv++;
        longjmp(buf, num_of_segv);
        break;
    }


}
void meshUnitTests()
{
  Mesh mesh;
    bool result=mesh.runTests();
    if(result)
      cout<<endl<<"All tests completed without a problem."<<endl;
    else
      cout<<endl<<"Some tests failed."<<endl;
}

int main(int argc, char **argv) 
{
  double error_correction=0.00015;
  if ((argc != 2 && argc!=6 ) ||argv[1]==string("--help")) 
  {
    
    //if(argc<2)
     // {
        cout << "Usage: " << argv[0] << " file.stl a b c d (plane-optional or 0 0 1 1 will be used) " << endl<<"or"<<endl<< argv[0]<<" tests"<<endl;
        cout<< "Plane is in a form of a plane normal (a, b, c) and d = distance from origin."<<endl;
        cout<< "Cutting through non-manifold object might result in error or wrong triangulation of cutted area."<<endl;
      return 1;
  //}

  }
  if(argv[1]==string("tests"))
  {
    meshUnitTests();
    //meshIntegrationTests();
    return 0;
  }

  stl_plane plane=stl_plane(0,0,1,-1);
  if(argc==6)
  {
    plane=stl_plane(atof( argv[2]),atof(argv[3]),atof(argv[4]),atof(argv[5]));
  }

  //Setting jump in case of poly2tri segfault
  setjmp(buf);
  sigemptyset(&signal_set);
  sigaddset(&signal_set, SIGSEGV); 
  sigprocmask(SIG_UNBLOCK, &signal_set, NULL); //clearnig segfault signal
  signal(SIGSEGV, segv_handler);
  if(num_of_segv==0)//!setjmp(buf))
  {
    Mesh mesh;
    mesh.open(argv[1]);
    /*  
    mesh.cut(plane); 
    mesh.poly2triTest();
    mesh.save();
    */
    mesh.cut(plane); 
    cout<<"1"<<endl;                         
    if(mesh.createBorderPolylines())
      {
        cout<<"2"<<endl;
        mesh.findHoles();
        cout<<"3"<<endl;
        mesh.triangulateCut();
        cout<<"4"<<endl;
        mesh.save();
      }
  
    mesh.close();
  }
  else
  {
    if(num_of_segv>6)
    {
      cerr<<"STLCUT wasnt able to made this cut. Try changing the plane position slightly and make sure that your model is manifold."<<endl;
      return 1;
    }
    cout<<endl<<"recovered from segfault"<<endl;
    int x;
    cin>>x;
    Mesh mesh2;
    cout<<argv[1]<<endl;
    mesh2.open(argv[1]);
    /*  
    mesh.cut(plane); 
    mesh.poly2triTest();
    mesh.save();
    */
    //plane.d=-0.5;//(plane.d)-0.099;
    plane=stl_plane(atof( argv[2]),atof(argv[3]),atof(argv[4]),(-1)*plane.d + error_correction);
    error_correction=error_correction>0?(-1)*error_correction:error_correction*10;//(-1)*error_correction * 10;
    cout<<"planeD je: "<<plane.d<<endl;
    mesh2.cut(plane); 
    cout<<"1"<<endl;                         
    if(mesh2.createBorderPolylines())
    {
      cout<<"2"<<endl;
      mesh2.findHoles();
      cout<<"3"<<endl;
      mesh2.triangulateCut();
      cout<<"4"<<endl;
      mesh2.save();
    }
    mesh2.close();
  }
return 0;
}