#ifndef h2Factory_h
#define h2Factory_h

#include <map>
#include <string>
#include <functional>
#include "hChain.h"
#include "h2Chain.h"
#include "TString.h"

class h2Factory
{
  public :
    h2Factory (std::string fileName = "no", bool print = true) ;
    ~h2Factory () ;

    void add_h2 (TString baseName, TString baseTitle, 
                 int nbinsx, double minx, double maxx,
                 int nbinsy, double miny, double maxy, int NUM) ;

    template <class UnaryFunction>
    void applyToAll (UnaryFunction function) 
      {
        for (std::map <TString, h2Chain *>::iterator mapIt = m_H2content.begin () ;
             mapIt != m_H2content.end () ;
             ++mapIt)
          {
            function (mapIt->second) ;
          }
      }

    void Fill (const TString & name, int i, double valx, double valy) ;
    h2Chain * operator[] (const TString& name) ;
    void Print (int isLog = 0, int rebin = 1) ;

  private :
  
    std::string m_fileName ;
    bool m_print ;
    std::map <TString, h2Chain *> m_H2content ;

} ;

#endif
