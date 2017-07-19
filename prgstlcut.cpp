#include "stlcut.h"
#include <getopt.h>
#include <unistd.h>
#include <string.h>
using namespace std;
/*
jmp_buf buf;
sigset_t signal_set;
int numOfSegv=0;
*//*
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
}*/
void meshUnitTests()
{
  Mesh mesh;
    bool result = mesh.runUnitTests();
    if(result)
      cout<<endl<<"All tests completed without a problem."<<endl;
    else
      cout<<endl<<"Some tests failed."<<endl;
}

void acquireSaveName(string& name, bool silent)
{
  name = "";
  if(!silent)
    {
      cout<<"Please enter the new name of meshes ( \"_1.stl\" and \"_2.stl\" will be added)"<<endl;
      cout<<"Or pres Enter for default name (Cut_Mesh_1/2.stl)"<<endl;
    }
  while(true)
  {
    getline(cin,name);//cin>>name;
    if(! (find_if(name.begin(), name.end(), 
        [](char c) { return !(isalnum(c) || (c == ' ') || (c == '_' )); }) == name.end() ) )
    {
      if(!silent)cout<<"Plese use only alphanumeric characters, space and underscore"<<endl;
      name="";
      continue;
    }
    else
      return;
  }
}

void printHelp(const char * name)
{
  cout << "StlCut is program for cuting STL meshes using plane. If succesful, two new STL files will be created."<<endl<<endl;
  cout << "\e[1m-Usage:\e[0m- " << name << " file.stl a b c d (optional) -- if not provided 0 0 1 1 will be used. " << endl<<endl;
  cout << "Plane is in a form of a plane normal vector (a, b, c) and d = distance from origin."<<endl;
  cout << "Cutting through non-manifold object might result in error or wrong triangulation of cut area."<<endl<<endl;
  cout << "To run tests use: " << name << " tests"<<endl<<endl;
  cout << "Options:"<<endl;
  cout << "\e[1m--help (-h)\e[0m " << "Displays this help."<<endl<<endl;
  cout << "\e[1m--silent (-s)\e[0m "<< "StlCut won't write anything to command line. Beware, that it won't ask you for new name of meshes, but will still require it."<<endl<<endl;
  cout << "\e[1m--error-recovery=<true/false> (-e)\e[0m "<< "If this option is false, StlCut will not try to recover from possible errors (and deliver meshes cut with slightly changed plane position). Short version -e sets this option to false. Error-recovery is true in defult."<<endl;
}



 struct options {
        bool silent;
        bool help;
        bool error_recovery;
};

int main(int argc, char ** argv)
{
  string name = "";
  stl_plane plane = stl_plane(0, 0, 1, 1);
  int op;
  struct options Option = { false, false, true };
  const struct option longopts[] = 
  {
    { "silent",         no_argument,       NULL, 'S' },
    { "help",           no_argument,       NULL, 'H' },
    { "error-recovery", optional_argument, NULL, 'E' },
    { NULL,         no_argument,       NULL,  0  }
  };

  if (argc >= 6 && argv[argc-1] != string("tests"))
  {
      plane = stl_plane(atof( argv[argc-4]), atof(argv[argc-3]), atof(argv[argc-2]), atof(argv[argc-1]));
      argc -=4;
  }
  do 
  {   
      op = getopt_long(argc, argv, "she::", longopts, NULL);
      
      switch(op) 
      {
        case -1:  break;         
        case ':': cout << "Missing argument for :"<<endl; return 2;
        case 'H':
        case 'h': Option.help = true; break;
        case 'S':
        case 's': Option.silent = true; break;
        case 'e': Option.error_recovery = false; break;
        case 'E': Option.error_recovery = false;
                  if(optarg != NULL)
                  {
                    if(strcmp(optarg,"true") == 0)
                    {
                      Option.error_recovery = true;
                      break;
                    }
                    if(strcmp(optarg,"false") != 0) 
                      cout << "Invalid argument for error-recovery, use true or false"<<endl;
                    return 1;
                  }
      }    
  } while(op != -1);

//////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
  if(Option.help || (argc-1) !=  optind )
  {
    printHelp(argv[0]);
    return 0;
  }

  if(argv[optind] == string("tests"))
  {
    meshUnitTests();
    //meshIntegrationTests();
    return 0;
  }

  Mesh mesh;
  mesh.setOptions(Option.silent, Option.error_recovery);
  mesh.openStl(argv[optind]);
  if (mesh.cut(plane))
    {
      acquireSaveName(name, Option.silent);
      mesh.save(name);
    }

  return 0;
}













/*
int main(int argc, char **argv) 
{
  //double error_correction = 0.00015;
  string name = "";
  /*
  if ((argc != 2 && argc != 6 ) || argv[1]==string(" --help ")) 
  {
    
    //if(argc<2)
     // {
        cout << "Usage: " << argv[0] << " file.stl a b c d (plane-optional or 0 0 1 1 will be used) " << endl<<"or"<<endl<< argv[0]<<" tests"<<endl;
        cout<< "Plane is in a form of a plane normal (a, b, c) and d = distance from origin."<<endl;
        cout<< "Cutting through non-manifold object might result in error or wrong triangulation of cutted area."<<endl;
      return 1;
  //}
  }*/
//-------------------------------
/*
  if(argv[optind] == string("tests"))
  {
    meshUnitTests();
    //meshIntegrationTests();
    return 0;
  }
stl_plane plane = stl_plane(0, 0, 1, -1);




 //-----------------
/*
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
*/
//----------------------
  /*
  //Setting jump in case of poly2tri segfault
  setjmp(buf);
  sigemptyset(&signal_set);
  sigaddset(&signal_set, SIGSEGV); 
  sigprocmask(SIG_UNBLOCK, &signal_set, NULL); //clearnig segfault signal
  signal(SIGSEGV, segv_handler);
*/
  /*if(numOfSegv > 6)
  {
     cerr<<"STLCUT wasnt able to made this cut. Try changing the plane position slightly and make sure that your model is 2-manifold."<<endl;
    return 1;
  }
  else
  {*//*
    Mesh mesh;
    mesh.openStl(argv[optind]);//argv[1]);

    /*if(numOfSegv > 0) 
    {
      plane = stl_plane(atof( argv[2]),atof(argv[3]),atof(argv[4]),(-1)*plane.d + error_correction);
      error_correction=error_correction>0?(-1)*error_correction:error_correction*10;//(-1)*error_correction * 10;
      cout<<endl<<"Recovered from segmentation fault."<<endl;
    }
*/
/*
    if (mesh.cut(plane))
    {
      acquireSaveName(name);
      mesh.save(name);
    }
  //}
  return 0;
}
*/