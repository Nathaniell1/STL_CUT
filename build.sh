#!/bin/bash -e
if [[ "$LOCAL" == "1" ]]; then
  flags='./poly2tri/sweep/*.o ./poly2tri/common/*.o'
else
  flags='-lpoly2tri'
fi
g++ -g -std=c++11 -o stlcut prgstlcut.cpp stlcut.cpp $flags -ladmesh
g++ -g -std=c++11 -o cuttests tests.cpp stlcut.cpp $flags -ladmesh
