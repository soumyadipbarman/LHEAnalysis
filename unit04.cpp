// c++ -o unit04 `root-config --glibs --cflags` -lm hFactory.cc hChain.cc unit04.cpp

#include "hFactory.h"
#include "hFunctions.h"

int main ()
{
  hFactory factory ("test.root", true) ;
  factory.add_h1 ("test", "jet tags eta prod", 200, -10, 10, 3) ;
  factory.Fill ("test", 0, 1.) ;
  int colors[3] = {92,98,0} ;
  factory.applyToAll (setH1Colors (std::vector<int> (colors, colors + 3))) ;
}