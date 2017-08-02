#include "stlcut.h"
#include <getopt.h>
#include <unistd.h>
#include <string.h>
using namespace std;



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
      name = "";
      continue;
    }
    else
      return;
  }
}

void printHelp(const char * name)
{
  cout << "StlCut is program for cuting STL meshes using plane. If succesful, two new STL files will be created."<<endl<<endl;
  cout << "\e[1m-Usage:\e[0m- " << name << " file.stl a b c d (optional) -- if not provided 0 0 1 0 will be used. " << endl<<endl;
  cout << "Plane is in a form of a plane normal vector (a, b, c) and d = distance from origin."<<endl;
  cout << "Cutting through non-manifold object might result in error or wrong triangulation of cut area."<<endl<<endl;
  cout << "Options:"<<endl;
  cout << "\e[1m--help (-h)\e[0m " << "Displays this help."<<endl<<endl;
  cout << "\e[1m--silent (-s)\e[0m "<< "StlCut won't write anything to command line. Beware, that it won't ask you for new name of meshes, but will still require it."<<endl<<endl;
  cout << "\e[1m--error-recovery=<true/false> (-e)\e[0m "<< "If this option is false, StlCut will not try to recover from possible errors (recovery leads to meshes cut with slightly changed plane position). Short version -e sets this option to false. Error-recovery is set true in default."<<endl;
}



 struct options {
        bool silent;
        bool help;
        bool error_recovery;
};

int main(int argc, char ** argv)
{
  string name = "";
  stl_plane plane = stl_plane(0, 0, 1, 0);
  int op;
  struct options Option = { false, false, true };
  const struct option longopts[] = 
  {
    { "silent",         no_argument,       NULL, 'S' },
    { "help",           no_argument,       NULL, 'H' },
    { "error-recovery", optional_argument, NULL, 'E' },
    { NULL,         no_argument,       NULL,  0  }
  };

  if (argc >= 6 )
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

