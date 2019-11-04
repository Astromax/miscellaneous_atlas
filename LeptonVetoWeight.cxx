#include "SusySkimMaker/LeptonVetoWeight.h"
#include "PathResolver/PathResolver.h"
#include "SusySkimMaker/MsgLog.h"

// ROOT includes
#include "TSystem.h"

LeptonVetoWeight::LeptonVetoWeight()
{


}
// ------------------------------------------------------------------------ //
bool LeptonVetoWeight::Initialize()
{

  TFile* file = TFile::Open("/gpfs/slac/atlas/fs1/u/mbaugh/RooFit/RateMaps/ZTAP_muon_efficiencies_data_BW+L_09-25.root", "READ");
  m_muonRates = dynamic_cast<TH2*>( file->Get("Post_Processed_Efficiency_TransferMap_vs_probeEta_and_probePt") );
  m_muonRates->SetDirectory(0);
  file->Close();

  /*
  TFile* fileB = TFile::Open("/gpfs/slac/atlas/fs1/u/mbaugh/RooFit/RateMaps/ZTAP_muon_efficiencies_data_BW+L_09-29.root", "READ");
  m_muonCVRates = dynamic_cast<TH2*>( fileB->Get("Post_Processed_Efficiency_NoMSTrack_TransferMap_vs_probeEta_and_probePt") );
  m_muonCVRates->SetDirectory(0);
  fileB->Close();
  */

  TFile* file2 =TFile::Open("/gpfs/slac/atlas/fs1/u/mbaugh/RooFit/RateMaps/ZTAP_electron_efficiencies_data_BW+L_09-25.root", "READ");
  m_electronRates = dynamic_cast<TH2*>( file2->Get("Post_Processed_Efficiency_TransferMap_vs_probeEta_and_probePt") );
  m_electronRates->SetDirectory(0);
  file2->Close();

  /*
  TFile* file2B = TFile::Open("/gpfs/slac/atlas/fs1/u/mbaugh/RooFit/RateMaps/ZTAP_electron_efficiencies_data_BW+L_09-29.root", "READ");
  m_electronCVRates= dynamic_cast<TH2*>( fileB->Get("Post_Processed_Efficiency_CaloVeto_TransferMap_vs_probeEta_and_probePt") );
  m_electronCVRates->SetDirectory(0);
  file2B->Close();
  */

  return true;
}
// ------------------------------------------------------------------------ // 
double LeptonVetoWeight::GetWeight(int f, float pt, float eta, float phi, float mu)
{

  double w1(0.0), w2(0.0);

  //Info("LeptonVetoWeight::GetWeight","Getting appropriate weight...");

  //If Pt out of range, use value in the highest reliable Pt bin
  if(pt>100) pt = 90;

  //If Pt below 20, set to zero
  if(pt<20) return 0.0;

  // Check for muons first, then electrons, if neither set to 0 (nothing else is faking it)
  if (f==1){
    w1 = m_muonRates->GetBinContent(m_muonRates->FindBin(eta, pt) );
    //w2 = m_muonCVRates->GetBinContent(m_muonCVRates->FindBin(eta, phi) );
  } else if(f==2){
    w1 = m_electronRates->GetBinContent(m_electronRates->FindBin(eta, pt) );
    //w2 = m_electronCVRates->GetBinContent(m_electronCVRates->FindBin(eta, pt) );
  } 
  //Info("LeptonVetoWeight::GetWeight","Lepton Veto Weight for Pt %f, eta %f,  and flavor %i is %.10f",pt,eta,f,w);

  //double w = w1 * w2;

  return w1;
}
// ------------------------------------------------------------------------ // 
double LeptonVetoWeight::GetUncertainty(int f, float pt, float eta, float phi, float mu)
{

  double e1(0.0), e2(0.0), w1(0.0), w2(0.0);

  if(pt>100) pt = 90;

  //If Pt below 20, set to zero                                                                                                                                                                                                           
  if(pt<20) return 0.0;

  if(f==1){
    e1 = m_muonRates->GetBinError(m_muonRates->FindBin(eta, pt) );
    //e2 = m_muonCVRates->GetBinError(m_muonCVRates->FindBin(eta, phi) );
    w1 = m_muonRates->GetBinContent(m_muonRates->FindBin(eta, pt) );
    //w2 = m_muonCVRates->GetBinContent(m_muonCVRates->FindBin(eta, phi) );
      
  }else if(f==2){
    e1 = m_electronRates->GetBinError(m_electronRates->FindBin(fabs(eta), pt) );
    //e2 = m_electronCVRates->GetBinError(m_electronCVRates->FindBin(eta, pt) );
    w1 = m_electronRates->GetBinContent(m_electronRates->FindBin(eta, pt) );
    //w2 = m_electronCVRates->GetBinContent(m_electronCVRates->FindBin(eta, pt) );
  }

  //double e = sqrt((e1*w2) ** 2 + (e2*w1) ** 2);

  return e1;
}
