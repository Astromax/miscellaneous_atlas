#ifndef SusySkimMaker_SmearLeptonPt_h
#define SusySkimMaker_SmearLeptonPt_h

// ROOT include(s)
#include "TFile.h"
#include "TString.h"
#include "TF1.h"
#include "TLorentzVector.h"

// C++ includes
#include <fstream>
#include <iostream>
#include <stdarg.h>
#include <map>
#include <vector>
/*
 
  Written by: Max Baugh

*/

class SmearLeptonPt
{

 public:
   SmearLeptonPt();

   ~SmearLeptonPt(){}
   
   bool Initialize();
   std::vector<float> SmearPt(int f, float pt);
   float SmearZmass(float probePt_smeared, float probeEta, float probePhi, float probeE, float tagPt, float tagEta, float tagPhi, float tagE);

   struct ParseResult{
     std::vector<float> pt_low;
     std::vector<float> pt_high;
     std::vector<TF1*> smears;
     std::vector<TF1*> smears_MeanUp;
     std::vector<TF1*> smears_MeanDown;
     std::vector<TF1*> smears_SigmaUp;
     std::vector<TF1*> smears_SigmaDown;
     std::vector<TF1*> smears_AlphaUp;
     std::vector<TF1*> smears_AlphaDown;
   };

 private:

   void initSmearing();
   ParseResult parse(const std::string& filename);
   ParseResult m_elparseResult;
   ParseResult m_muparseResult;
   double crystallBall(double* x, double *par);

};


#endif
