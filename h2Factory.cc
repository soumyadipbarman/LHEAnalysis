#include "h2Factory.h"
#include "TString.h"
#include <iostream>


h2Factory::h2Factory (std::string fileName, 
                    bool print) :
 m_fileName (fileName) ,
 m_print (print)                     
{}


//PG --------------------------------------------------------   


h2Factory::~h2Factory () 
{
//  if (m_H2content.size () && m_print)
//    for (std::map <TString, h2Chain *>::iterator mapIt = m_H2content.begin () ;
//         mapIt != m_H2content.end () ;
//         ++mapIt)
//      {
//        mapIt->second->Print () ; 
//        //PG FIXME generalizzare le opzioni di stampa    
//      }
  if (m_H2content.size () && m_fileName != "no")
    {
      TFile histosFile (m_fileName.c_str (),"recreate") ;
      histosFile.cd () ;
      for (std::map <TString, h2Chain *>::iterator mapIt = m_H2content.begin () ;
           mapIt != m_H2content.end () ;
           ++mapIt)
        {
          mapIt->second->Write (histosFile) ; 
        }
      histosFile.Close () ;
    }  
}


//PG --------------------------------------------------------   


void
h2Factory::add_h2 (TString baseName, TString baseTitle, 
                   int nbinsx, double minx, double maxx,
                   int nbinsy, double miny, double maxy, int NUM)
{
  if (m_H2content.find (baseName) != m_H2content.end ())
    {
      std::cerr << "ERROR : histograms series " << baseName
                << " already existing, NOT replaced" << std::endl ;
      return ;                
    }
  h2Chain * dummy = new h2Chain (baseName, baseTitle, 
                                 nbinsx, minx, maxx,
                                 nbinsy, miny, maxy, NUM) ;
  m_H2content[baseName] = dummy ;     
  return ;                          
}
                 

//PG --------------------------------------------------------   


void 
h2Factory::Fill (const TString & name, int i, double valx, double valy) 
  {
    if (m_H2content.find (name) != m_H2content.end ())
      m_H2content[name]->Fill (i,valx, valy) ;
    return ;
  }


//PG --------------------------------------------------------   


h2Chain * 
h2Factory::operator[] (const TString& name)
  {
    if (m_H2content.find (name) == m_H2content.end ())
      {
        std::cerr << "ERROR : histograms series " << name
                  << " not existing, returning the first one" << std::endl ;
        return m_H2content.begin ()->second ;                
      }
    return m_H2content[name] ;
  
  }


//PG --------------------------------------------------------   


void 
h2Factory::Print (int isLog, int rebin)
{
  if (!m_H2content.size ()) return ;
  for (std::map <TString, h2Chain *>::iterator mapIt = m_H2content.begin () ;
       mapIt != m_H2content.end () ;
       ++mapIt)
    {
      mapIt->second->Print () ; 
    }
}

