// c++ -o unit06 `root-config --glibs --cflags` -lm h2Factory.cc hChain.cc h2Chain.cc unit06.cpp

#include "h2Factory.h"
#include "hFunctions.h"

int main ()
{
  h2Factory factory ("test2.root", true) ;
  factory.add_h2 ("test", "jet tags eta prod", 6, -3., 3., 4, -2., 2., 2) ;
  factory.Fill ("test", 0, 1., 1.) ;
  factory.Fill ("test", 1, 2., 1.) ;
  int colors[2] = {9,98} ;
  factory.applyToAll (setH2Colors (std::vector<int> (colors, colors + 2))) ;
}