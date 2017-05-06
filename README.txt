StlCut uses plane to cut stl file and generates two new stl files.
The hole made with cut is triangulated.
Use with 2--manifold models.
Plane is defined as normal (a,b,c) and distance from origin (d).
If you have problem with cut, try slightly changing a distance.


StlCut requires g++, admesh and poly2tri installed.

To install stlcut you can chose between:

<<<<<<< HEAD
Additional info:
You should encounter problems with kolo.stl and trubka.stl as these are non-manifold meshes.
=======
-------------------------------
1.

Compile stlcut with ./build and use it with ./stlcut file.stl
If your distribution does not provide poly2tri installation, download from github.com/greenm01/poly2tri/tree/master/poly2tri and place the folder with it to the folder with stlcut and name it poly2tri.
then use LOCAL=1 ./build

Now you can use ./stlcut file.stl a b c d
For aditional info use ./stlcut help

2.

Disclaimer: This part of installation is based on stlsplit installation (https://github.com/admesh/stlsplit)

This installation additionaly needs Premake 4 (http://premake.github.io/)

On Linux, you would use:

    premake4 gmake
    make
If you encounter problem erro while loading shared libraries use
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
    stlcut help

-------------------------------
As a result, you'll get `Cut_Mesh_1.stl`, `Cut_Mesh_2.stl` if the cut is succesful.

To use test script you have to have Python and pytest.
Run it with pytest -v test_stlcut.py
Test script uses stl files from stl_files folder.
You should encounter problems with kolo.stl and trubka.stl as they are non-manifold meshes.

------------------------

If you want to use ADMeshGUI you have to use installation number 2.
If you cant for some reason, download stlcut and change file meshobject.cpp
Change line #include <stlcut> to #include <stlcut.h> 
>>>>>>> development
