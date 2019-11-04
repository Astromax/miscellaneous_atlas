#include "SusySkimMaker/SmearLeptonPt.h"
#include "PathResolver/PathResolver.h"
#include "SusySkimMaker/MsgLog.h"

// ROOT includes
#include "TSystem.h"

// C++ includes
#include <iostream>
#include <fstream>
#include <array>
#include <vector>

SmearLeptonPt::SmearLeptonPt()
{


}
// ------------------------------------------------------------------------ //
bool SmearLeptonPt::Initialize()
{
  initSmearing();

  return true;
}
// ------------------------------------------------------------------------------------------ //
void SmearLeptonPt::initSmearing()
{

  //Fancy methods
  m_elparseResult = parse("/gpfs/slac/atlas/fs1/u/mbaugh/SusySkimAna/dummybranch/source/SusySkimLongLived/dummy_electron_v0.txt");
  m_muparseResult = parse("/gpfs/slac/atlas/fs1/u/mbaugh/SusySkimAna/dummybranch/source/SusySkimLongLived/dummy_muon_v0.txt");

}
// ------------------------------------------------------------------------------------------ //
/**/
SmearLeptonPt::ParseResult SmearLeptonPt::parse(const std::string& filename)
{
  std::vector<float> pt_low;
  std::vector<float> pt_high;
  std::vector<float> Mean;
  std::vector<float> mean_low;
  std::vector<float> mean_high;
  std::vector<float> Sigma;
  std::vector<float> sigma_low;
  std::vector<float> sigma_high;
  std::vector<float> Alpha;
  std::vector<float> alpha_low;
  std::vector<float> alpha_high;
  std::vector<float> Theta_MS;
  std::vector<float> Theta_MA;
  std::vector<float> Theta_SA;

  TF1* qoverpt;
  TF1* qoverpt_MeanUp;
  TF1* qoverpt_MeanDown;
  TF1* qoverpt_SigmaUp;
  TF1* qoverpt_SigmaDown;
  TF1* qoverpt_AlphaUp;
  TF1* qoverpt_AlphaDown;
  std::vector<TF1*> smears;
  std::vector<TF1*> smears_MeanUp;
  std::vector<TF1*> smears_MeanDown;
  std::vector<TF1*> smears_SigmaUp;
  std::vector<TF1*> smears_SigmaDown;
  std::vector<TF1*> smears_AlphaUp;
  std::vector<TF1*> smears_AlphaDown;

  std::ifstream file(filename);
  if(file.is_open()){
    std::string line;
    std::getline(file, line); //Skip the first line 
    while(std::getline(file, line)){
      std::stringstream ss(line);
      float ptlow,pthigh,mean,meanlow,meanhigh,sig,siglow,sighigh,alpha,alphalow,alphahigh,thetaMS,thetaMA,thetaSA;
      if(ss >> ptlow >> pthigh >> mean >> meanlow >> meanhigh >> sig >> siglow >> sighigh >> alpha >> alphalow >> alphahigh >> thetaMS >> thetaMA >> thetaSA){
        pt_low.push_back(ptlow);
        pt_high.push_back(pthigh);
        Mean.push_back(mean);
        mean_low.push_back(meanlow);
        mean_high.push_back(meanhigh);
        Sigma.push_back(sig);
        sigma_low.push_back(siglow);
        sigma_high.push_back(sighigh);
        Alpha.push_back(alpha);
        alpha_low.push_back(alphalow);
        alpha_high.push_back(alphahigh);
	Theta_MS.push_back(thetaMS);
	Theta_MA.push_back(thetaMA);
	Theta_SA.push_back(thetaSA);
        Info("SmearLeptonPt::initSmearing", "This line in the file has pt_low %0.1f and pt_high %0.1f",ptlow,pthigh);
      }
    }
  }
    
  for(unsigned int i=0; i<pt_low.size(); i++){
    float delta_meanUp = mean_high[i] - Mean[i];
    float delta_meanDown = mean_low[i] - Mean[i];
    float delta_SigmaUp = sigma_high[i] - Sigma[i];
    float delta_SigmaDown = sigma_low[i] - Sigma[i];
    float delta_AlphaUp = alpha_high[i] - Alpha[i];
    float delta_AlphaDown = alpha_low[i] - Alpha[i];

    qoverpt = new TF1(TString::Format("FittedSmearingFunction_%i",i), [&](double*x, double*p){ return crystallBall(x,p); }, -1000, 1000, 4);
    qoverpt->SetParameter(0, 1.0);
    qoverpt->SetParameter(1, Mean[i]);
    qoverpt->SetParameter(2, Sigma[i]);
    qoverpt->SetParameter(3, Alpha[i]);
    smears.push_back(qoverpt);

    qoverpt_MeanUp = new TF1(TString::Format("FittedSmearingFunction_MeanUp_%i",i), [&](double*x, double*p){ return crystallBall(x,p); }, -1000, 1000, 4);                                                                                             
    qoverpt_MeanUp->SetParameter(0, 1.0);
    qoverpt_MeanUp->SetParameter(1, mean_high[i]);
    qoverpt_MeanUp->SetParameter(2, Sigma[i] + Theta_MS[i] * delta_meanUp);
    qoverpt_MeanUp->SetParameter(3, Alpha[i] + Theta_MA[i] * delta_meanUp);
    smears_MeanUp.push_back(qoverpt_MeanUp);

    qoverpt_MeanDown = new TF1(TString::Format("FittedSmearingFunction_MeanDown_%i",i), [&](double*x, double*p){ return crystallBall(x,p); }, -1000, 1000, 4);
    qoverpt_MeanDown->SetParameter(0, 1.0);
    qoverpt_MeanDown->SetParameter(1, mean_low[i]);
    qoverpt_MeanDown->SetParameter(2, Sigma[i] + Theta_MS[i] * delta_meanDown);
    qoverpt_MeanDown->SetParameter(3, Alpha[i] + Theta_MA[i] * delta_meanDown);
    smears_MeanDown.push_back(qoverpt_MeanDown);

    qoverpt_SigmaUp = new TF1(TString::Format("FittedSmearingFunction_SigmaUp_%i",i), [&](double*x, double*p){ return crystallBall(x,p); }, -1000, 1000, 4);
    qoverpt_SigmaUp->SetParameter(0, 1.0);
    qoverpt_SigmaUp->SetParameter(1, Mean[i] + Theta_MS[i] * delta_SigmaUp);
    qoverpt_SigmaUp->SetParameter(2, sigma_high[i]);
    qoverpt_SigmaUp->SetParameter(3, Alpha[i] + Theta_SA[i] * delta_SigmaUp);
    smears_SigmaUp.push_back(qoverpt_SigmaUp);

    qoverpt_SigmaDown = new TF1(TString::Format("FittedSmearingFunction_SigmaDown_%i",i), [&](double*x, double*p){ return crystallBall(x,p); }, -1000, 1000, 4);
    qoverpt_SigmaDown->SetParameter(0, 1.0);
    qoverpt_SigmaDown->SetParameter(1, Mean[i] + Theta_MS[i] * delta_SigmaDown);
    qoverpt_SigmaDown->SetParameter(2, sigma_low[i]);
    qoverpt_SigmaDown->SetParameter(3, Alpha[i] + Theta_SA[i] * delta_SigmaDown);
    smears_SigmaDown.push_back(qoverpt_SigmaDown);

    qoverpt_AlphaUp = new TF1(TString::Format("FittedSmearingFunction_AlphaUp_%i",i), [&](double*x, double*p){ return crystallBall(x,p); }, -1000, 1000, 4);
    qoverpt_AlphaUp->SetParameter(0, 1.0);
    qoverpt_AlphaUp->SetParameter(1, Mean[i] + Theta_MA[i] * delta_AlphaUp);
    qoverpt_AlphaUp->SetParameter(2, Sigma[i] + Theta_SA[i] * delta_AlphaUp);
    qoverpt_AlphaUp->SetParameter(3, alpha_high[i]);
    smears_AlphaUp.push_back(qoverpt_AlphaUp);

    qoverpt_AlphaDown = new TF1(TString::Format("FittedSmearingFunction_AlphaDown_%i",i), [&](double*x, double*p){ return crystallBall(x,p); }, -1000, 1000, 4);
    qoverpt_AlphaDown->SetParameter(0, 1.0);
    qoverpt_AlphaDown->SetParameter(1, Mean[i] + Theta_MA[i] * delta_AlphaDown);
    qoverpt_AlphaDown->SetParameter(2, Sigma[i] + Theta_SA[i] * delta_AlphaDown);
    qoverpt_AlphaDown->SetParameter(3, alpha_low[i]);
    smears_AlphaDown.push_back(qoverpt_AlphaDown);
  }

  return ParseResult {pt_low, pt_high, smears, smears_MeanUp, smears_MeanDown, smears_SigmaUp, smears_SigmaDown, smears_AlphaUp, smears_AlphaDown};
}
// ------------------------------------------------------------------------------------------ //
double SmearLeptonPt::crystallBall(double* x, double *par)
{
  double constant = par[0];
  double mean = par[1];
  double sigma = par[2];
  double alpha = par[3];//*sigma;

  if (sigma < 0.)     return 0.;
  if (alpha < 0.)     return 0.;
  double z = (x[0] - mean)/sigma;
  alpha = std::abs(alpha);
  double norm1 = sigma*sqrt(2*M_PI)*erf(alpha/sqrt(2));
  double norm2 = sigma*exp(-alpha*alpha/2)/alpha;
  double norm3 = norm2;
  constant /= (norm1 + norm2 + norm3);
  if (z  < - alpha){
    return constant * std::exp( alpha * (z + 0.5 * alpha));
  }else if (z  > + alpha){
    return constant * std::exp(-alpha * (z - 0.5 * alpha));
  }else{
    return constant * std::exp(- 0.5 * z * z);
  }

}
// ------------------------------------------------------------------------ // 
std::vector<float> SmearLeptonPt::SmearPt(int f, float pt)
{

  float smearedPt(0.0), smearedPt_MeanUp(0.0), smearedPt_MeanDown(0.0), smearedPt_SigmaUp(0.0), smearedPt_SigmaDown(0.0), smearedPt_AlphaUp(0.0), smearedPt_AlphaDown(0.0);
  std::vector<float> smearedPts;

  //Info("SmearLeptonPt::SmearPt","Smearing the lepton Pt...");
  auto qoverpt = 1.0 / pt; // GeV^-1

  if(f==2 && m_elparseResult.smears.size()==0){
    std::cout<<"Electron smearing functions are missing! \n";
    smearedPts = {pt, pt, pt, pt, pt, pt, pt};
    return smearedPts;
  }

  if(f==2 && m_muparseResult.smears.size()==0){
    std::cout<<"Muon smearing functions are missing! \n";
    smearedPts = {pt, pt, pt, pt, pt, pt, pt};
    return smearedPts;
  }

  int target_index(-1);
  auto smearedQoverpt = qoverpt;
  auto smearedQoverpt_MeanUp = qoverpt;
  auto smearedQoverpt_MeanDown = qoverpt;
  auto smearedQoverpt_SigmaUp = qoverpt;
  auto smearedQoverpt_SigmaDown = qoverpt;
  auto smearedQoverpt_AlphaUp = qoverpt;
  auto smearedQoverpt_AlphaDown = qoverpt;

  if(f==1){
    for(unsigned int i=0; i<m_muparseResult.pt_low.size(); i++){
      if(i==0 && pt < m_muparseResult.pt_low[i]){
	target_index = i;
	break;
      }
      if(pt >= m_muparseResult.pt_low[i] && pt < m_muparseResult.pt_high[i]){
	target_index = i;
	break;
      }
    }
    if(target_index == -1) target_index = m_muparseResult.pt_low.size() - 1;
    //Info("SmearLeptonPt::SmearPt", "The target index for muon pt %0.2f is %i",pt,target_index);
    smearedQoverpt += m_muparseResult.smears[target_index]->GetRandom()*(1e-3);//TeV^-1 -> GeV^-1
    //float delta = m_muparseResult.smears[target_index]->GetRandom()*(1e-3);
    //smearedQoverpt += delta;
    //Info("SmearLeptonPt::SmearPt", "Original pt: %0.2f, delta qopt: %0.2f",pt, delta);
    smearedQoverpt_MeanUp += m_muparseResult.smears_MeanUp[target_index]->GetRandom()*(1e-3);
    smearedQoverpt_MeanDown += m_muparseResult.smears_MeanDown[target_index]->GetRandom()*(1e-3);
    smearedQoverpt_SigmaUp += m_muparseResult.smears_SigmaUp[target_index]->GetRandom()*(1e-3);
    smearedQoverpt_SigmaDown += m_muparseResult.smears_SigmaDown[target_index]->GetRandom()*(1e-3);
    smearedQoverpt_AlphaUp += m_muparseResult.smears_AlphaUp[target_index]->GetRandom()*(1e-3);
    smearedQoverpt_AlphaDown += m_muparseResult.smears_AlphaDown[target_index]->GetRandom()*(1e-3);
  }else if(f==2){
    for(unsigned int i=0; i<m_elparseResult.pt_low.size(); i++){
      if(i==0 && pt < m_elparseResult.pt_low[i]){
	target_index = i;
	break;
      }
      if(pt >= m_elparseResult.pt_low[i] && pt < m_elparseResult.pt_high[i]){
	target_index = i;
	break;
      }
    }
    if(target_index == -1) target_index = m_elparseResult.pt_low.size() - 1;
    //Info("SmearLeptonPt::SmearPt", "The target index for electron pt %0.2f is %i",pt,target_index);
    smearedQoverpt += m_elparseResult.smears[target_index]->GetRandom()*(1e-3);//TeV^-1 -> GeV^-1
    smearedQoverpt_MeanUp += m_elparseResult.smears_MeanUp[target_index]->GetRandom()*(1e-3);
    smearedQoverpt_MeanDown += m_elparseResult.smears_MeanDown[target_index]->GetRandom()*(1e-3);
    smearedQoverpt_SigmaUp += m_elparseResult.smears_SigmaUp[target_index]->GetRandom()*(1e-3);
    smearedQoverpt_SigmaDown += m_elparseResult.smears_SigmaDown[target_index]->GetRandom()*(1e-3);
    smearedQoverpt_AlphaUp += m_elparseResult.smears_AlphaUp[target_index]->GetRandom()*(1e-3);
    smearedQoverpt_AlphaDown += m_elparseResult.smears_AlphaDown[target_index]->GetRandom()*(1e-3);
  }
  smearedPt = fabs(1 / smearedQoverpt);
  smearedPt_MeanUp = fabs(1 / smearedQoverpt_MeanUp);
  smearedPt_MeanDown = fabs(1 / smearedQoverpt_MeanDown);
  smearedPt_SigmaUp = fabs(1 / smearedQoverpt_SigmaUp);
  smearedPt_SigmaDown = fabs(1 / smearedQoverpt_SigmaDown);
  smearedPt_AlphaUp = fabs(1 / smearedQoverpt_AlphaUp);
  smearedPt_AlphaDown = fabs(1 / smearedQoverpt_AlphaDown);

  smearedPts.push_back(smearedPt);
  smearedPts.push_back(smearedPt_MeanUp);
  smearedPts.push_back(smearedPt_MeanDown);
  smearedPts.push_back(smearedPt_SigmaUp);
  smearedPts.push_back(smearedPt_SigmaDown);
  smearedPts.push_back(smearedPt_AlphaUp);
  smearedPts.push_back(smearedPt_AlphaDown);

  /*
  Info("SmearLeptonPt::SmearPt","Lepton Pt %f has been smeared to %f",pt,smearedPt);
  Info("SmearLeptonPt::SmearPt","Lepton Pt %f has smeared_MeanUp of %f",pt,smearedPt_MeanUp);
  Info("SmearLeptonPt::SmearPt","Lepton Pt %f has smeared_MeanDown of %f",pt,smearedPt_MeanDown);
  Info("SmearLeptonPt::SmearPt","Lepton Pt %f has smeared_SigmaUp of %f",pt,smearedPt_SigmaUp);
  Info("SmearLeptonPt::SmearPt","Lepton Pt %f has smeared_SigmaDown of %f",pt,smearedPt_SigmaDown);
  Info("SmearLeptonPt::SmearPt","Lepton Pt %f has smeared_AlphaUp of %f",pt,smearedPt_AlphaUp);
  Info("SmearLeptonPt::SmearPt","Lepton Pt %f has smeared_AlphaDown of %f",pt,smearedPt_AlphaDown);
  */

  return smearedPts;
}
//------------------------------------------------------------------------------------//
float SmearLeptonPt::SmearZmass(float probePt_smeared, float probeEta, float probePhi, float probeE, float tagPt, float tagEta, float tagPhi, float tagE)
{

  TLorentzVector probe;
  TLorentzVector tag;

  probe.SetPtEtaPhiE(probePt_smeared, probeEta, probePhi, probeE);
  tag.SetPtEtaPhiE(tagPt, tagEta, tagPhi, tagE);

  float Zmass_smeared = (probe + tag).M();

  return Zmass_smeared;

}
