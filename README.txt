StlCut uses plane to cut stl file and generates two new stl files.
The hole made with cut is triangulated.
Use with 2--manifold triangular models.
Plane is defined as normal (a,b,c) and distance from origin (d).
If you have problem with cut, try slightly changing a distance.


StlCut requires g++, admesh and poly2tri installed.

To install stlcut you can chose between:

-------------------------------
1.

Compile stlcut with ./build and use it with ./stlcut file.stl
If your distribution does not provide poly2tri installation, download it from github.com/greenm01/poly2tri/tree/master/poly2tri and place the folder with it to the folder with stlcut and name it poly2tri.
then use LOCAL=1 ./build

Now you can use ./stlcut file.stl a b c d
For aditional info use ./stlcut --help

2.

Disclaimer: This part of installation is based on stlsplit installation (https://github.com/admesh/stlsplit)

This installation additionaly needs Premake 4 (http://premake.github.io/)

On Linux, you would use:

    premake4 gmake
    make

If you encounter error while loading shared libraries use:
export LD_LIBRARY_PATH=.

To install stlcut, put the compiled binary and the shared library into your `$PATH`:

    sudo cp build/stlcut /usr/bin
    sudo cp build/libstlcut.so.1 /usr/lib # or lib64

If you intend to build something with stlcut as a library, you'll also need the header file and .so symlink:

    sudo cp stlcut.h /usr/include
    sudo ln -s libstlcut.so.1 /usr/lib/libstlcut.so # or lib64
    sudo ldconfig

Running the command line tool:

    stlcut file.stl a b c d 
    For aditional info use
    stlcut --help

-------------------------------

You can test if everything is working properly with:
./cuttests ./stl_files/*.stl
For aditional info use
./cuttests --help

------------------------------

How to use stlcut library in your program:

include stlcut.h
Create new instance of Mesh class.
Use openStl(char * name) or setStl(stl_file f) to set your mesh.
(Optional) setOptions(bool silent, bool errorRecovery) First argument, if true, will disable text output. Second argument, if false, will disable recovery from errors during cut.
Use cut(stl_plane plane) To cut through your object. Takes stl_plane which defines cutting plane as parameter.
If cut() returns true, use save(string name) to save new meshes or getFinalStls() to get pointer to their stl_files.

Example:

  Mesh mesh;
  mesh.openStl("bunny.stl");
  if (mesh.cut(stl_plane(0,0,1,0)))
      mesh.save("half_of_bunny");   


If you want to use ADMeshGUI you have to use installation number 2.
