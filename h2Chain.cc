#include "h2Chain.h"

  h2Chain::h2Chain (TString baseName, TString baseTitle, 
                    int nbinsx, double minx, double maxx,
                    int nbinsy, double miny, double maxy, int NUM) :
    m_baseName (baseName)
  {
    TString name ;
    TString title ;
    for (int i=0 ; i<NUM ; ++i)
      {
        name = TString ("h2_") ; 
        name += i ;
        name += TString ("_") + m_baseName ;
        title = baseTitle + TString (" ") ;
        title += i ;
        TH2F * dummy = new TH2F (name,title,nbinsx,minx,maxx,nbinsy,miny,maxy) ;
        dummy->SetStats (0) ;
        m_histos.push_back (dummy) ;
      }
  }


//PG --------------------------------------------------------   


h2Chain::~h2Chain () 
{
  for (int i=0 ; i<m_histos.size () ; ++i)
    delete m_histos.at (i) ;
}  


//PG --------------------------------------------------------   


void 
h2Chain::SetColors (std::vector<int> colors)
  {
    //PG this is weak, assumes correct number of elements
    for (int i=0 ; i<m_histos.size () ; ++i)
      {
        m_histos.at (i)->SetMarkerColor (colors.at (i)) ;
        m_histos.at (i)->SetLineColor (colors.at (i)) ;
        m_histos.at (i)->SetFillColor (colors.at (i)) ;
        m_histos.at (i)->SetFillStyle (3001) ;
      }
    //PG this is a hack
//    m_histos.back ()->SetLineColor (kBlack) ;
//    m_histos.back ()->SetLineWidth (2) ;
  } 


//PG --------------------------------------------------------   


void 
h2Chain::Fill (int i, double valx, double valy) 
  {
    m_histos.at (i)->Fill (valx, valy) ;
    return ;
  }


//PG --------------------------------------------------------   


void 
h2Chain::Print () 
  {  
    TCanvas c1 ;
    c1.cd () ;
    m_histos.at (0)->Draw ("box") ;
    for (int i=0 ; i<m_histos.size () - 1 ; ++i)
      {
        m_histos.at (i)->SetStats (0) ;
        m_histos.at (i)->Draw ("boxsame") ;
      }  
    //PG this is a hack
//    m_histos.back ()->Draw ("same") ;
    TString name = m_baseName + TString ("_h2.gif") ;
    c1.Print (name,"gif") ;

    m_histos.at (0)->DrawNormalized ("box") ;
    for (int i=0 ; i<m_histos.size () - 1 ; ++i)
      m_histos.at (i)->DrawNormalized ("boxsame") ;
    //PG this is a hack
//    m_histos.back ()->DrawNormalized ("same") ;
    name = TString ("N") + m_baseName + TString ("_h2.gif") ;
    c1.Print (name,"gif") ;
  }


//PG --------------------------------------------------------   


void 
h2Chain::PrintEach () 
  {  
    TCanvas c1 ;
    c1.cd () ;
    for (int i=0 ; i<m_histos.size () ; ++i)
      {
        m_histos.at (i)->SetStats (0) ;
        m_histos.at (i)->Draw () ;
        TString name = m_baseName ;
                name += TString ("_") ;
                name += i ;
                name += TString ("_h2.gif") ;
        c1.Print (name,"gif") ;
      }        
  }


//PG --------------------------------------------------------   


void 
h2Chain::Write (TFile & outputFile)
  {
    outputFile.cd () ;
    for (int i=0 ; i<m_histos.size () ; ++i)
      m_histos.at (i)->Write () ;
    return ;
  }


//PG --------------------------------------------------------   


void 
h2Chain::Scale (int index, double factor)
{
  m_histos.at (index)->Scale (factor) ;
}

