// Usage:
// > root -b doAll.C

// C++
#include <iostream>
#include <vector>
#include <map>

// ROOT
#include "TBenchmark.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH3D.h"
#include "TROOT.h"
#include "TTreeCache.h"
#include "TLorentzVector.h"
#include <sstream>
#include <iostream>
#include <fstream>

// CMS3
//#include "CMS3_old20150505.cc"
//#include "CMS3_fuckingsync.cc"
//#include "CMS3_Moriond17.cc"
#include "CMS3.cc"

//MT2 variants

using namespace std;
using namespace tas;

class Cut{
public:
  string variable;
  float cut1;
  float cut2;
};

float dRbetweenVectors(LorentzVector vec1,LorentzVector vec2 ){                                                        
  float dphi = std::min(::fabs(vec1.Phi() - vec2.Phi()), 2 * M_PI - fabs(vec1.Phi() - vec2.Phi()));
  float deta = vec1.Eta() - vec2.Eta();
  return sqrt(dphi*dphi + deta*deta);
}

float calculateMt(LorentzVector p4, LorentzVector met){
  float phi1 = p4.Phi();
  float phi2 = met.Phi();
  float Et1  = p4.Et();
  float Et2  = met.Et();

  return sqrt(2*Et1*Et2*(1.0 - cos(phi1-phi2)));
}

float dPhiBetweenVectors(float phi1,float phi2 ){                                                                                                              
  return fabs(std::min(float(fabs(phi1-phi2)),float(2*M_PI-fabs(phi1-phi2))));
}

float fill_mt_met_lep(){
  return  mt_met_lep();
}

float fill_n_of_first_bjet(){
  float fill = 0;
  for (unsigned int i=0; i<ak4pfjets_passMEDbtag().size(); ++i){
    if (ak4pfjets_passMEDbtag()[i]){
      fill = i+1;
      break;
    }
  }
  return fill;
}

float fill_n_bjets(){
  float count = 0;
  for (unsigned int i=0; i<ak4pfjets_p4().size(); ++i){
    if(ak4pfjets_passMEDbtag()[i]){
      count = count + 1;
    }
  }
  return count; 
}

float fill_pfmet(){
  return pfmet();
}

float fill_pt_lep1_p4(){
  return lep1_p4().pt();
}

float fill_pt_ISR(){
  return ak4pfjets_p4()[0].Pt();
  }

float fill_dphi_ISR_any(){
  float phi_ISR = ak4pfjets_p4()[0].Phi();
  float dphi = 4;
  for (unsigned int i = 1; i<ak4pfjets_p4().size(); ++i){
    if(dPhiBetweenVectors(phi_ISR, ak4pfjets_p4()[i].Phi())<dphi){
      dphi = dPhiBetweenVectors(phi_ISR, ak4pfjets_p4()[i].Phi());
    }
  }
  if (dphi > dPhiBetweenVectors(phi_ISR, lep1_p4().Phi())){
    dphi = dPhiBetweenVectors(phi_ISR, lep1_p4().Phi());
  }
  if (dphi_ak4pfjet_met()[0]<dphi){
    dphi = dphi_ak4pfjet_met()[0];
  }
  return dphi;
}

float fill_dphi_ISR_jet(){
  float phi_ISR = ak4pfjets_p4()[0].Phi();
  float dphi = 4;
  for (unsigned int i = 1; i<ak4pfjets_p4().size(); ++i){
    if(dPhiBetweenVectors(phi_ISR, ak4pfjets_p4()[i].Phi())<dphi){
      dphi = dPhiBetweenVectors(phi_ISR, ak4pfjets_p4()[i].Phi());
    }
  }
  if(dphi>3.5){
    dphi = 137;
  }
  return dphi;
}

float fill_dphi_ISR_bjet() {
  float phi_ISR = ak4pfjets_p4()[0].Phi();
  float dphi = 4;
  for (unsigned int i = 1; i<ak4pfjets_p4().size(); ++i){
    if(dPhiBetweenVectors(phi_ISR, ak4pfjets_p4()[i].Phi())<dphi&&ak4pfjets_passMEDbtag()[i]){
      dphi = dPhiBetweenVectors(phi_ISR, ak4pfjets_p4()[i].Phi());
    }
  }
  if (dphi >3.5){
    dphi = 137;
  }
  return dphi;

}

float fill_dphi_ISR_lep(){
  float phi_ISR = ak4pfjets_p4()[0].Phi();
  float dphi = dPhiBetweenVectors(phi_ISR, lep1_p4().Phi());
  return dphi;
}

float fill_dphi_ISR_ptmiss(){
  return dphi_ak4pfjet_met()[0];
}

float fill_dphi_bjet_lep(){
  float dphi = 4;
  for (unsigned int i = 0; i<ak4pfjets_p4().size(); ++i){
    if(dPhiBetweenVectors(ak4pfjets_p4()[i].Phi(), lep1_p4().Phi())<dphi&&ak4pfjets_passMEDbtag()[i]){
      dphi = dPhiBetweenVectors(ak4pfjets_p4()[i].Phi(), lep1_p4().Phi());
    }
  }
  if(dphi > 3.5){
    dphi = 137;
  }
  return dphi;
}

float fill_dphi_bjet_ptmiss(){
  float dphi = 4;
  for (unsigned int i = 0; i<dphi_ak4pfjet_met().size(); ++i){
    if(dphi>dphi_ak4pfjet_met()[i]&&ak4pfjets_passMEDbtag()[i]){
      dphi = dphi_ak4pfjet_met()[i];
    }
  }
  if(dphi>3.5){
    dphi = 137;
  }
  return dphi;
}

float fill_dphi_lep_ptmiss(){
  return dPhiBetweenVectors(lep1_p4().Phi(), pfmet_phi());
}

float fill_mt(){
  float dphi = dPhiBetweenVectors(ak4pfjets_p4()[0].Phi(), pfmet_phi());
  float mt = sqrt(2*pfmet()*ak4pfjets_p4()[0].Pt()*(1-cos(dphi)));
  return mt;
}

float fill_R(){
  return abs(pfmet()/(ak4pfjets_p4()[0].Pt()));
}

float fill_r_pt(){
  float r_pt = 137;
  for (unsigned int i = 0; i<ak4pfjets_p4().size(); ++i){
    if(ak4pfjets_passMEDbtag()[i]){
      r_pt = (ak4pfjets_p4()[i].Pt()-lep1_p4().pt())/(ak4pfjets_p4()[i].Pt()+lep1_p4().pt());
      break;
    }
  }
  return r_pt;
}

float fill_r_pt_prime(){
  float r_pt_prime = (ak4pfjets_p4()[0].Pt()-lep1_p4().pt())/(ak4pfjets_p4()[0].Pt()+lep1_p4().pt());
  return r_pt_prime;
}

float fill_R_bE(){
  float pt_b = 137;
  for (unsigned int i = 0; i<ak4pfjets_p4().size(); ++i){
    if(ak4pfjets_passMEDbtag()[i]){
      pt_b = pt_b + ak4pfjets_p4()[i].Pt();
    }
  }
  return pt_b/pfmet();
}

float fill_R_lE(){
  return lep1_p4().pt()/pfmet();
}

float fill_M_jj(){
  float M_jj = 137;
  float M_jj_diff;
  float M_W = 80;
  for(unsigned int i = 1; i<ak4pfjets_p4().size(); ++i){
    if(ak4pfjets_passMEDbtag()[i]) continue;
    for(unsigned int j = i+1; j<ak4pfjets_p4().size(); ++j){
      if(ak4pfjets_passMEDbtag()[j]) continue;
      else if(abs((ak4pfjets_p4()[i]+ak4pfjets_p4()[j]).M()-M_W)<M_jj_diff||(i==1&&j==2)){
	M_jj_diff = abs((ak4pfjets_p4()[i]+ak4pfjets_p4()[j]).M()-M_W);
	M_jj = (ak4pfjets_p4()[i]+ak4pfjets_p4()[j]).M();
      }
    }    
  }
  return M_jj;
}

float fill_csv_sum(){
  float M_jj = 137;
  float M_W = 80;
  float a;
  float b;
  for(unsigned int i = 1; i<ak4pfjets_p4().size(); ++i){
    if(ak4pfjets_passMEDbtag()[i]) continue;
    for(unsigned int j = i+1; j<ak4pfjets_p4().size(); ++j){
      if(ak4pfjets_passMEDbtag()[j]) continue;
      else if(abs((ak4pfjets_p4()[i]+ak4pfjets_p4()[j]).M()-M_W)<M_jj||(i==1&&j==2)){
	M_jj = abs((ak4pfjets_p4()[i]+ak4pfjets_p4()[j]).M()-M_W);
	a = i;
	b = j;
      }
    }
  }
    return ak4pfjets_CSV()[a] + ak4pfjets_CSV()[b];
}

 float fill_csv_min(){
   float M_jj = 137;
   float M_W = 80;
   float a;
   float b;
   for(unsigned int i = 1; i<ak4pfjets_p4().size(); ++i){
     if(ak4pfjets_passMEDbtag()[i]) continue;
     for(unsigned int j = i+1; j<ak4pfjets_p4().size(); ++j){
       if(ak4pfjets_passMEDbtag()[j]) continue;
       else if(abs((ak4pfjets_p4()[i]+ak4pfjets_p4()[j]).M()-M_W)<M_jj||(i==1&&j==2)){
	 M_jj = abs((ak4pfjets_p4()[i]+ak4pfjets_p4()[j]).M()-M_W);
	 a = i;
	 b = j;
       }
     }
   }
   return min(ak4pfjets_CSV()[a], ak4pfjets_CSV()[b]);
 }

  float fill_csv_max(){
    float M_jj = 137;
    float M_W = 80;
    float a;
    float b;
    for(unsigned int i = 1; i<ak4pfjets_p4().size(); ++i){
      if(ak4pfjets_passMEDbtag()[i]) continue;
      for(unsigned int j = i+1; j<ak4pfjets_p4().size(); ++j){
	if(ak4pfjets_passMEDbtag()[j]) continue;
	else if(abs((ak4pfjets_p4()[i]+ak4pfjets_p4()[j]).M()-M_W)<M_jj||(i==1&&j==2)){
	  M_jj = abs((ak4pfjets_p4()[i]+ak4pfjets_p4()[j]).M()-M_W);
	  a = i;
	  b = j;
	}
      }
    }
    return max(ak4pfjets_CSV()[a], ak4pfjets_CSV()[b]);
  }

float fill_n_loose_bjets(){
  float M_jj = 137;
  float M_W = 80;
  float a;
  float b;
  for(unsigned int i = 1; i<ak4pfjets_p4().size(); ++i){
    if(ak4pfjets_passMEDbtag()[i]) continue;
    for(unsigned int j = i+1; j<ak4pfjets_p4().size(); ++j){
      if(ak4pfjets_passMEDbtag()[j]) continue;
      else if(abs((ak4pfjets_p4()[i]+ak4pfjets_p4()[j]).M()-M_W)<M_jj||(i==1&&j==2)){
	M_jj = abs((ak4pfjets_p4()[i]+ak4pfjets_p4()[j]).M()-M_W);
	a = i;
	b = j;
      }
    }
  }
  int n = 0;
  if(ak4pfjets_CSV()[a]>0.5426) ++n;
  if(ak4pfjets_CSV()[b]>0.5426) ++n;
  return n;
}

float fill_n_loose_minus_good(){
  return nloosebtags()-ngoodbtags();
}

int ScanChain( TChain* chain, bool fast = true, int nEvents = -1, string skimFilePrefix = "test") {

  vector<string> histonames;
  vector<int> hbins;
  vector<float> hlow;
  vector<float> hup;
  vector<Cut> fullCuts;
  vector<vector<Cut> > cuts;

  map<string, TH1F*>histos;
  map<string, set<int> > Nsignalpoints;
  map<string, float> funcMap;

  Cut n_loose_minus_good;
  n_loose_minus_good.variable = "n_loose_minus_good";
  n_loose_minus_good.cut1 = 0;
  n_loose_minus_good.cut2 = 1;
  fullCuts.push_back(n_loose_minus_good);

  /*Cut csv_sum;
  csv_sum.variable = "csv_sum";
  csv_sum.cut1 = 0;
  csv_sum.cut2 = 0.7;
  fullCuts.push_back(csv_sum);*/

  Cut dphi_ISR_ptmiss;
  dphi_ISR_ptmiss.variable = "dphi_ISR_ptmiss";
  dphi_ISR_ptmiss.cut1 = 0.5;//2
  dphi_ISR_ptmiss.cut2 = 4; 
  fullCuts.push_back(dphi_ISR_ptmiss);

  Cut dphi_lep_ptmiss;
  dphi_lep_ptmiss.variable = "dphi_lep_ptmiss";
  dphi_lep_ptmiss.cut1 = 0;
  dphi_lep_ptmiss.cut2 = 2; //1.5
  fullCuts.push_back(dphi_lep_ptmiss);

  Cut n_of_first_bjet;
  n_of_first_bjet.variable = "n_of_first_bjet";
  n_of_first_bjet.cut1 = 0;
  n_of_first_bjet.cut2 = 10;
  fullCuts.push_back(n_of_first_bjet);

  Cut pfmet; 
  pfmet.variable = "pfmet";
  pfmet.cut1 = 500;
  pfmet.cut2 = 1500;
  fullCuts.push_back(pfmet);

  Cut pt_ISR;
  pt_ISR.variable = "pt_ISR";
  pt_ISR.cut1 = 250;
  pt_ISR.cut2 = 1500;
  fullCuts.push_back(pt_ISR);

  Cut pt_lep1_p4;
  pt_lep1_p4.variable = "pt_lep1_p4";
  pt_lep1_p4.cut1 = 0;
  pt_lep1_p4.cut2 = 150;
  fullCuts.push_back(pt_lep1_p4);

  Cut R;
  R.variable = "R";
  R.cut1 = 0.5;
  R.cut2 = 2;
  fullCuts.push_back(R);

  //histonames.push_back("csv_sum");
  //hbins.push_back(20);      hlow.push_back(0);       hup.push_back(2);

  histonames.push_back("n_loose_minus_good");
  hbins.push_back(10);   hlow.push_back(0); hup.push_back(10);  

  histonames.push_back("dphi_ISR_ptmiss");
  hbins.push_back(20);     hlow.push_back(0);       hup.push_back(4);

  histonames.push_back("dphi_lep_ptmiss");
  hbins.push_back(20);     hlow.push_back(0);       hup.push_back(4);

  histonames.push_back("n_of_first_bjet");
  hbins.push_back(10);   hlow.push_back(0); hup.push_back(10);

  histonames.push_back("pfmet");
  hbins.push_back(27);     hlow.push_back(150);     hup.push_back(1500);

  histonames.push_back("pt_ISR");
  hbins.push_back(20);      hlow.push_back(0);       hup.push_back(1500);

  histonames.push_back("pt_lep1_p4");
  hbins.push_back(20);      hlow.push_back(0);       hup.push_back(500);

  histonames.push_back("R");
  hbins.push_back(30);     hlow.push_back(0);       hup.push_back(3);
  
  int numHistos = histonames.size();

  for(int i = 0; i<numHistos; ++i){
    int stop = 0;
    for (unsigned int j = 0; j<fullCuts.size(); ++j){
      for (unsigned int k = j+1; k<fullCuts.size(); ++k){
	for (unsigned int l = k+1; l<fullCuts.size(); ++l){
	  if (j == k|| j==l || l==k || histonames[i]==fullCuts[j].variable || histonames[i]==fullCuts[k].variable || histonames[i] == fullCuts[l].variable){//don't want repeated cut combinations, don't want to graph variables that aren't part of the cut
	    stop = 1;
	  continue;
	  }
	  vector<Cut> cutSet;
	  cutSet.push_back(fullCuts[j]);
	  cutSet.push_back(fullCuts[k]);
	  cutSet.push_back(fullCuts[l]);
	  
	  //declare cut histograms
	  histonames.push_back(histonames[i] + "_cut_" + cutSet[0].variable + "_" + cutSet[1].variable + "_" + cutSet[2].variable);
	  hbins.push_back(hbins[i]);
	  hlow.push_back(hlow[i]);
	  hup.push_back(hup[i]);
	}
	if(stop==1) continue;
      }
      if(stop==1) continue;
    }
  }

 for (unsigned int j = 0; j<fullCuts.size(); ++j){
   for (unsigned int k = j+1; k<fullCuts.size(); ++k){
     for (unsigned int l = k+1; l<fullCuts.size(); ++l){
       if (j == k || j==l || l==k){//don't want repeated cut combinations, don't want to graph variables that aren't part of the cut
	 continue;
       }
       vector<Cut> cutSet;
       cutSet.push_back(fullCuts[j]);
       cutSet.push_back(fullCuts[k]);
       cutSet.push_back(fullCuts[l]);
       cuts.push_back(cutSet);
     }
   }
 }
 //make cut combinations
  for(unsigned int i = 0; i<histonames.size(); ++i){
    for(unsigned int b = 0; b<4; ++b){//why 4?
      
      string samplename = skimFilePrefix;

      if(samplename.find(string("LostLeptonAndTop")) != string::npos){
	
	if(b==0) samplename = "LostLepton";
	else if (b==1) samplename = "TT1l";
	else continue;
	
      }

      else if(samplename.find(string("Signal_T2tt")) != string::npos){

	if(b==0) samplename = "Signal_T2tt_Wcorridor"; 
	else if(b==1) samplename = "Signal_T2tt_topcorridor";
        else if(b==2) samplename = "Signal_T2tt_betweencorridor";
        else if(b==3) samplename = "Signal_T2tt_highDM";
	else continue; //we don't consider the other signals
        set<int> dummy;
        Nsignalpoints[samplename] = dummy;

      }
      else if(samplename.find(string("SignalGen_T2tt")) != string::npos){

        if(b==0) samplename = "SignalGen_T2tt_Wcorridor";      
        else if(b==1) samplename = "SignalGen_T2tt_topcorridor";    
        else if(b==2) samplename = "SignalGen_T2tt_betweencorridor";
        else if(b==3) samplename = "SignalGen_T2tt_highDM";         
        set<int> dummy;
        Nsignalpoints[samplename] = dummy;

      }

      else if(b>=1) continue;//for all other samples, we don't need splitting
      string mapname = histonames[i] + "_" + samplename;
      if(histos.count(mapname) == 0 ){
	if(histonames[i]=="pfmet_cut_n_loose_minus_good_dphi_lep_ptmiss_pt_lep1_p4"){
	  float xbins3[6] = {250, 300, 400, 500, 650, 900};
	  histos[mapname] = new TH1F(mapname.c_str(), "", 5, xbins3);
	  continue; 
	}

	histos[mapname] = new TH1F(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
      }
      histos[mapname]->Sumw2(); 
    }
  }
   
     unsigned int nEventsRunning = 0;
     //Goal: Loop over events to Analyze
     unsigned int nEventsTotal = 0;
     unsigned int nEventsChain = chain->GetEntries();
     if( nEvents >= 0 ) nEventsChain = nEvents;
     TObjArray *listOfFiles = chain->GetListOfFiles();
     TIter fileIter(listOfFiles);
     TFile *currentFile = 0;
     TH3D* counterhistSig;
     TH1D* counterhist;

     // File Loop - loop over all files in the TChain
     while ( (currentFile = (TFile*)fileIter.Next()) ) {
       cout << "Running over file " << currentFile->GetTitle() << endl;
       // Get File Content - load counter histograms and the tree containing all events
       TFile *file = new TFile( currentFile->GetTitle() );
       if(skimFilePrefix.find(string("Signal")) != string::npos){
	 counterhistSig = (TH3D*)file->Get("h_counterSMS");
	 counterhistSig->SetDirectory(0);
       } else { 
	 counterhist = (TH1D*)file->Get("h_counter");
	 counterhist->SetDirectory(0);
       }

       TTree *tree = (TTree*)file->Get("t");
       cms3.Init(tree);
    
       // Loop over Events from the tree of the current file
       if( nEventsTotal >= nEventsChain ) continue;
       unsigned int nEventsTree = tree->GetEntriesFast();
       for( unsigned int event = 0; event < nEventsTree; ++event) {
 
	 // Get Event Content
	 if( nEventsTotal >= nEventsChain ) continue;
	 if(fast) tree->LoadTree(event);
	 cms3.GetEntry(event);//load branches (needed in C++)
	 ++nEventsTotal;

	 //Progress
	 CMS3::progress( nEventsTotal, nEventsChain );

	 float weight = 1;
	 
	 
	 //define what type of sample it is
	 string samplename = skimFilePrefix;
	 //cout<<samplename<<endl;
	 if(samplename.find(string("LostLeptonAndTop")) != string::npos){

	   if(is1lepFromTop()) samplename = "TT1l";
	   else samplename = "LostLepton";

	 }

	 else if(samplename.find(string("Signal_T2tt")) != string::npos){
	   
	   if(mass_lsp()<100.) continue;
	   if(     (mass_stop()-mass_lsp())<  98.) samplename = "Signal_T2tt_Wcorridor";
	   else if((mass_stop()-mass_lsp())> 165.&&(mass_stop()-mass_lsp())< 185.) samplename = "Signal_T2tt_topcorridor";
	   else if((mass_stop()-mass_lsp())>= 99.&&(mass_stop()-mass_lsp())<=150.) samplename = "Signal_T2tt_betweencorridor";
	   else if((mass_stop()-mass_lsp())>=199.&&(mass_stop()-mass_lsp())<=250.) samplename = "Signal_T2tt_highDM";
	   else continue; //not of interest for our study

	 }
    
	 else if(samplename.find(string("SignalGen_T2tt")) != string::npos){
	   if(mass_lsp()<100.) continue;
	   if((mass_stop()-mass_lsp())<  98.) samplename = "SignalGen_T2tt_Wcorridor";
	   else if((mass_stop()-mass_lsp())> 165.&&(mass_stop()-mass_lsp())< 185.) samplename = "SignalGen_T2tt_topcorridor";
	   else if((mass_stop()-mass_lsp())>= 99.&&(mass_stop()-mass_lsp())<=150.) samplename = "SignalGen_T2tt_betweencorridor";
	   else if((mass_stop()-mass_lsp())>=199.&&(mass_stop()-mass_lsp())<=250.) samplename = "SignalGen_T2tt_highDM";
	   else continue; //not of interest for our study
      }

	 //get correct event weight!
	 //expected event number (for real data) N = cross section * Luminosity
	 //if we generate M simulated events, we need to normalize M to the correctly expected event yield: N = weight * M
	 //--> weight = N / M = cross section * luminosity / M
	 //in our ntuples: scale1fb = cross section * 1 fb^-1 / M --> can scale to any expected data set size
	 if(!is_data() && skimFilePrefix.find(string("Signal")) == string::npos) {//this is all background
	   int Nevts = counterhist->GetBinContent(counterhist->FindBin(22));//this is the total number of events in the unskimmed/generated sample
	   double nevts = double(Nevts);
	   weight = scale1fb()*150.;//this is the event weight
	   double lepSFweight = weight_lepSF()*nevts/counterhist->GetBinContent(counterhist->FindBin(28));
	   if(lepSFweight>=0&&!std::isinf(lepSFweight)&&!std::isnan(lepSFweight)) weight *= lepSFweight;
	 } 
	 else if(!is_data()){ //data
	   int Nevts = counterhistSig->GetBinContent(counterhistSig->FindBin(mass_stop(),mass_lsp(),36));
	   Nsignalpoints[samplename].insert(counterhistSig->FindBin(mass_stop(),mass_lsp(),36));//this is a unique identifier for the signal point
	   double nevts = double(Nevts);
	   double lepSFweight = weight_lepSF()*nevts/counterhistSig->GetBinContent(counterhistSig->FindBin(mass_stop(),mass_lsp(),27));
	   double lepFSweight = weight_lepSF_fastSim()*nevts/counterhistSig->GetBinContent(counterhistSig->FindBin(mass_stop(),mass_lsp(),33));//the signal is produced with FastSimulation - needs an additional scale factor between fast simulation and full simulation
        weight = xsec()*150000/nevts;//xsec is given in pb, not fb
        if(lepSFweight>=0&&!std::isinf(lepSFweight)&&!std::isnan(lepSFweight)) weight *= lepSFweight;
        if(lepFSweight>=0&&!std::isinf(lepFSweight)&&!std::isnan(lepFSweight)) weight *= lepFSweight;
        if(!filt_fastsimjets() )     continue;
	}
      TString currentfilename = currentFile->GetTitle();

      //have 2 ttbar-->2l2nu samples that model same thing - combined event weight must be take into account
      if(currentfilename.Contains("ttbar_diLept_madgraph_pythia8_25ns")) weight *= 5.77109e+06/(5.77109e+06 + 2.34556e+07);
      if(currentfilename.Contains("ttbar_diLept_madgraph_pythia8_ext1_25ns")) weight *= 2.34556e+07/(5.77109e+06 + 2.34556e+07);

      //inclusive WJets sample contains also events with nu-pT > 200, but we have an extra sample for that - veto the double counting 
      if(skimFilePrefix.find(string("WJets")) != string::npos){
        if(currentfilename.Contains("nupT200")){ if(nupt()<200) continue; }
        else { if(nupt()>=200) continue; }
      }

      //preselection
      if(nvtxs()<0)               continue;
      if(ngoodleps()!=1)          continue;
      if(nvetoleps()!=1)          continue;
      if(!PassTrackVeto())        continue;
      if(!PassTauVeto())          continue;
      if(ngoodjets()<5)           continue;
      if(ngoodbtags()<1)          continue;
      if(fill_pfmet()<250)             continue;
      if(mt_met_lep()<150)        continue;
      if(mindphi_met_j1_j2()<0.5) continue;
      if(ak4pfjets_passMEDbtag()[0]) continue;

	 //here begins the core of this looper - we have no our baseline events - let's investigate those
      funcMap["csv_sum"] = fill_csv_sum();
      funcMap["dphi_ISR_ptmiss"] = fill_dphi_ISR_ptmiss();
      funcMap["dphi_lep_ptmiss"] = fill_dphi_lep_ptmiss();
      funcMap["n_of_first_bjet"] = fill_n_of_first_bjet();
      funcMap["n_loose_minus_good"] = fill_n_loose_minus_good();
      funcMap["pfmet"] = fill_pfmet();
      funcMap["pt_ISR"] = fill_pt_ISR();
      funcMap["pt_lep1_p4"] = fill_pt_lep1_p4();
      funcMap["R"] = fill_R();
      funcMap["R_lE"] = fill_R_lE();    
      funcMap["r_pt"] = fill_r_pt();
      
	 //now we fill histograms
      string mapname;
      //fill uncut histograms
      for (int i = 0; i<numHistos; ++i){
	mapname = histonames[i] + "_" + samplename;
	if (funcMap[histonames[i] ] == 137.) continue;
	/*else if(funcMap["dphi_lep_ptmiss"]>(-(funcMap["csv_sum"]-1)/0.35)){
	  continue;
	  }*/
	/*else if((funcMap["pt_lep1_p4"]>50&&funcMap["dphi_lep_ptmiss"]>2)||(funcMap["dphi_lep_ptmiss"]<2&&funcMap["pt_lep1_p4"]>-funcMap["dphi_lep_ptmiss"]*100+250)){
	  continue;
	  }*/ 
	else histos[mapname] -> Fill(funcMap[histonames[i] ], weight);
      }

      //fill cut histograms
      for (int i=0; i<numHistos; ++i){
	if (funcMap[histonames[i] ] == 137.) continue;
	for (unsigned int j = 0; j<cuts.size(); ++j){
	  int flag = 0;
	  if (histonames[i]== cuts[j][0].variable||histonames[i]== cuts[j][1].variable||histonames[i] == cuts[j][2].variable) continue;
	  for(unsigned int k = 0; k<cuts[j].size(); ++k){
	    if (funcMap[cuts[j][k].variable ] < cuts[j][k].cut1){
	      flag = 1;
	      break;
	    }
	    else if (funcMap[cuts[j][k].variable ] > cuts[j][k].cut2){
	      flag = 1;
	      break;
	    }
	    else if (cuts[j][k].variable == "n_of_first_bjet" && (funcMap["n_of_first_bjet"] == 1 || (funcMap["n_of_first_bjet"] == 2))){
	      flag = 1;
	      break;
	    }
	    /*else if(funcMap["dphi_lep_ptmiss"]>(-(funcMap["csv_sum"]-1)/0.35)){
	      flag = 1; 
	      break;
	      }*/
	    /*else if((funcMap["pt_lep1_p4"]>50&&funcMap["dphi_lep_ptmiss"]>2)||(funcMap["dphi_lep_ptmiss"]<2&&funcMap["pt_lep1_p4"]>-funcMap["dphi_lep_ptmiss"]*100+250)){
	      flag = 1; 
	      break;
	      }*/
	  }
	  if (flag == 1) continue;
	  mapname = histonames[i] + "_cut_" + cuts[j][0].variable + "_" + cuts[j][1].variable + "_" + cuts[j][2].variable + "_" + samplename;
	  histos[mapname] -> Fill(funcMap[histonames[i] ], weight);
	}
      }
       }//event loop
  
       // Clean Up
       delete tree;
       file->Close();
       delete file;
     }//file loop

     for(map<string,TH1F*>::iterator h=histos.begin(); h!=histos.end();++h){
       //add overflow
       h->second->SetBinContent(h->second->GetNbinsX(), h->second->GetBinContent(h->second->GetNbinsX() )+ h->second->GetBinContent(h->second->GetNbinsX()+1) );
       h->second->SetBinError(h->second->GetNbinsX(), sqrt(pow(h->second->GetBinError(h->second->GetNbinsX() ),2)+pow(h->second->GetBinError(h->second->GetNbinsX()+1),2) ) );
       //add underflow
       h->second->SetBinContent(1, h->second->GetBinContent(1)+ h->second->GetBinContent(0) );
       h->second->SetBinError(1, sqrt(pow(h->second->GetBinError(1),2)+pow(h->second->GetBinError(0),2) ) );
       //for signal only - normalize the signal yield to what we expect from a single event point (i.e. get averaged yield per point)
       if(h->first.find("Signal_T2tt_Wcorridor")      !=string::npos) h->second->Scale(1./Nsignalpoints["Signal_T2tt_Wcorridor"]      .size());
       if(h->first.find("Signal_T2tt_topcorridor")    !=string::npos) h->second->Scale(1./Nsignalpoints["Signal_T2tt_topcorridor"]    .size());
       if(h->first.find("Signal_T2tt_betweencorridor")!=string::npos) h->second->Scale(1./Nsignalpoints["Signal_T2tt_betweencorridor"].size());
       if(h->first.find("Signal_T2tt_highDM")         !=string::npos) h->second->Scale(1./Nsignalpoints["Signal_T2tt_highDM"]         .size());
       if(h->first.find("SignalGen_T2tt_Wcorridor")      !=string::npos) h->second->Scale(1./Nsignalpoints["SignalGen_T2tt_Wcorridor"]      .size());
       if(h->first.find("SignalGen_T2tt_topcorridor")    !=string::npos) h->second->Scale(1./Nsignalpoints["SignalGen_T2tt_topcorridor"]    .size());
       if(h->first.find("SignalGen_T2tt_betweencorridor")!=string::npos) h->second->Scale(1./Nsignalpoints["SignalGen_T2tt_betweencorridor"].size());
       if(h->first.find("SignalGen_T2tt_highDM")         !=string::npos) h->second->Scale(1./Nsignalpoints["SignalGen_T2tt_highDM"]         .size());
  }

     TFile *samples = new TFile("samples.root", "UPDATE");
     samples->cd();
     for(unsigned int i=0; i<histonames.size(); ++i){
       for(unsigned int b = 0; b<4; ++b){
	 string samplename = skimFilePrefix;
	 if(samplename.find(string("Top")) != string::npos){
	   if(b==0) samplename = "LostLepton";
	   else if(b==1) samplename = "TT1l";
	   else continue;
	 }
	 else if(samplename.find(string("Signal_T2tt")) != string::npos){
	  if(b==0) samplename = "Signal_T2tt_Wcorridor";
	  if(b==1) samplename = "Signal_T2tt_topcorridor";
	  if(b==2) samplename = "Signal_T2tt_betweencorridor";
	  if(b==3) samplename = "Signal_T2tt_highDM"; 
	 }
	 else if(samplename.find(string("SignalGen_T2tt")) != string::npos){
	   if(b==0) samplename = "SignalGen_T2tt_Wcorridor";
	   if(b==1) samplename = "SignalGen_T2tt_topcorridor";
	   if(b==2) samplename = "SignalGen_T2tt_betweencorridor";
	   if(b==3) samplename = "SignalGen_T2tt_highDM";
	 }
	 else if(b==1) continue;
	 string mapname = histonames[i] + "_" + samplename;
	 histos[mapname] -> Write(mapname.c_str(), TObject::kOverwrite);
       }
     }
     samples->Close();
     cout << "Saved histos in " << samples->GetName()<<endl;
     
     return 0;



}//function loop
