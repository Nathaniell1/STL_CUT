compile with: g++ -std=c++11 -o cut stlcut.cpp ./poly2tri/sweep/*.o ./poly2tri/common/*.o -ladmesh

You need to have admesh installed.

tests is a script to compile and test stlcut on a file with random planes (numbers from 0-1), cuts 10 times
how to use tests: ./tests file_to_cut , or just ./tests *.stl to cut every file in current folder


Additional info:
You should encounter problems with kolo.stl and trubka.stl as these are non-manifold meshes.
