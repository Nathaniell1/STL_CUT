compile with ./build.sh
If you dont have poly2tri, download it and place the folder with it to the STL_CUT folder and name it poly2tri
then use LOCAL=1 ./build

You need to have admesh installed.

tests is a script to compile and test stlcut on a file with random planes (numbers from 0-1), cuts 10 times
how to use tests: ./tests file_to_cut , or just ./tests *.stl to cut every file in current folder


Additional info:
You should encounter problems with kolo.stl and trubka.stl as they are non-manifold meshes.
