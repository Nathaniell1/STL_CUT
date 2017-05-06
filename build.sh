#!/bin/bash -e
if [[ "$LOCAL" == "1" ]]; then
  flags='./poly2tri/sweep/*.o ./poly2tri/common/*.o'
else
  flags='-lpoly2tri'
fi
g++ -std=c++11 -o stlcut prgstlcut.cpp stlcut.cpp $flags -lpoly2tri -ladmesh