#ifndef SusySkimMaker_LeptonVetoWeight_h
#define SusySkimMaker_LeptonVetoWeight_h

// ROOT include(s)
#include "TFile.h"
#include "TString.h"
#include "TH2.h"
#include "TH1.h"

// C++ includes
#include <fstream>
#include <iostream>
#include <stdarg.h>
#include <map>

/*
 
  Written by: Max Baugh, based on template from Matthew Gignac

*/


class LeptonVetoWeight
{

 public:
   LeptonVetoWeight();

   ~LeptonVetoWeight(){}

   bool Initialize();
   double GetWeight(int f, float pt, float eta, float phi, float mu);
   double GetUncertainty(int f, float pt, float eta, float phi, float mu);

 private:

   TH2* m_muonRates;
   TH1* m_muonPtRates;
   TH1* m_muonEtaRates;

   TH2* m_electronRates;
   TH1* m_electronPtRates;
   TH1* m_electronEtaRates;

};


#endif
