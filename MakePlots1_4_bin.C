#include "TList.h"
#include "TMath.h"
#include "TString.h"
#include "TObject.h"
#include "TFile.h"
#include "TPad.h"
#include "TDirectory.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TH3F.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TH3D.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TSystem.h"
#include "TMap.h"
#include "TStopwatch.h"
#include "TColor.h"
#include "TLegend.h"
#include "THStack.h"
#include "TLatex.h"
#include "TLine.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <map>
#include <string>
#include <cmath>

using namespace std;

void MakePlots1_4_bin(){

  
  bool  logy   = true; //should the y axis being plotted in logarithmic scale (true) or linear scale (false)
  bool  data   = false;//set this to false - we don't look at data
  float lumi   = 150.; //check what luminosity was put into the weight in ExampleLooper
  bool  addgen = true; //average reco and genMET signal samples



  vector<string> histonames;
  vector<string> histox;
  vector<string> selecttitle;
  vector<string> bgnames;
  vector<string> signames;
  vector<string> bgleg;
  vector<string> sigleg;
  vector<Color_t> bgcol;
  vector<Color_t> sigcol;
  vector<string> cutnames;

  map<string, TH1F*> hist;
  map<string, THStack*> stack;

  // cutnames.push_back("csv_sum");
  cutnames.push_back("n_loose_minus_good");
  cutnames.push_back("dphi_ISR_ptmiss");
  cutnames.push_back("dphi_lep_ptmiss");
  //cutnames.push_back("n_loose_bjets");
  cutnames.push_back("n_of_first_bjet");
  cutnames.push_back("pfmet");
  cutnames.push_back("pt_ISR");
  cutnames.push_back("pt_lep1_p4");
  cutnames.push_back("R");

  bgnames.push_back("Znunu");        
  bgcol.push_back(kMagenta-5);
  bgleg.push_back("Z#rightarrow#nu#bar{#nu}");

  bgnames.push_back("WJets");        
  bgcol.push_back(kOrange-2);
  bgleg.push_back("1#font[12]{l} not from top");

  bgnames.push_back("TT1l");       
  bgcol.push_back(kRed-7);
  bgleg.push_back("1#font[12]{l} from top");

  bgnames.push_back("LostLepton");   
  bgcol.push_back(kCyan-3);
  bgleg.push_back("Lost Lepton");

  histonames.push_back("n_loose_minus_good");  
  // histonames.push_back("csv_sum");
  histonames.push_back("dphi_ISR_ptmiss");
  histonames.push_back("dphi_lep_ptmiss");
  histonames.push_back("n_of_first_bjet");
  histonames.push_back("pfmet");
  histonames.push_back("pt_ISR");
  histonames.push_back("pt_lep1_p4");
  histonames.push_back("R");

  int numHistos = histonames.size();

  //histox.push_back("Max(CSV_{a}, CSV_{b})");
  //histox.push_back("Min(CSV_{a}, CSV_{b})");
  histox.push_back("N Loose BJets - N Good BJets");
  //histox.push_back("CSV_{a} + CSV_{b}");
  histox.push_back("Minimum #phi [ISR, p_{T}^{miss}]");
  histox.push_back("Minimum #phi [lepton, p_{T}^{miss}]");
  //histox.push_back("jet_motherid");
  histox.push_back("Position of 1st BJet");
  histox.push_back("p_{T}^{miss}");
  histox.push_back("p_{T} [ISR]");
  histox.push_back("p_{T} [l]");
  histox.push_back("R");

  vector<float> met_bin_fill;
  met_bin_fill.push_back(3);
  met_bin_fill.push_back(7);
  met_bin_fill.push_back(15);
  met_bin_fill.push_back(25);
  met_bin_fill.push_back(35);
  met_bin_fill.push_back(50);

  for(int i = 0; i<numHistos; ++i){
    int stop = 0;
    for (unsigned int j = 0; j<cutnames.size(); ++j){
      for(unsigned int k = j+1; k<cutnames.size(); ++k){
	for(unsigned int l = k+1; l<cutnames.size(); ++l){
	  for(unsigned int m = l+1; m<cutnames.size(); ++m){
	    //if (histonames[i]==cutnames[j]||histonames[i]==cutnames[k]||histonames[i]==cutnames[l]||histonames[i]==cutnames[m]){
	    //stop = 1;
	    //continue;
	    //}
	    //declare cut histograms
	    histonames.push_back(histonames[i] + "_cut_" + cutnames[j] + "_" + cutnames[k] + "_" + cutnames[l] + "_" + cutnames[m]);
	    histox.push_back(histox[i] + " with cuts " + cutnames[j] + ", " + cutnames[k] + ", " + cutnames[l] + "_" + cutnames[m]);
	  }
	  if(stop==1) continue;
	}
	if(stop==1) continue;
      }
      if(stop==1) continue;
    }
  }
  
  signames.push_back("Signal_T2tt_Wcorridor");       
  sigleg.push_back("#tilde{t}#rightarrowt#tilde{#chi}^{0}_{1} (W corridor)");     sigcol.push_back(kGreen+2);

  signames.push_back("Signal_T2tt_topcorridor");     
  sigleg.push_back("#tilde{t}#rightarrowt#tilde{#chi}^{0}_{1} (t corridor)");     sigcol.push_back(kBlue+1);

  signames.push_back("Signal_T2tt_betweencorridor"); 
  sigleg.push_back("#tilde{t}#rightarrowt#tilde{#chi}^{0}_{1} (between corridor)"); 
  sigcol.push_back(kYellow+1);

  signames.push_back("Signal_T2tt_highDM");          
  sigleg.push_back("#tilde{t}#rightarrowt#tilde{#chi}^{0}_{1} (above corridor)");   
  sigcol.push_back(kMagenta+1);

  gStyle->SetOptFit(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TFile *f = new TFile("samples.root","READ");

  for (int i = 0; i < histonames.size(); ++i){
    for (int b = 0; b < bgnames.size(); ++b){
      string mapname =  histonames[i] + "_" + bgnames[b];
      string bgname = "bgsum_" + histonames[i];
			 
      if (b == 0){

	hist[mapname] = (TH1F*)f -> Get(mapname.c_str());
	hist[bgname] = (TH1F*)hist[mapname] -> Clone(bgname.c_str());
	hist[mapname]->SetLineColor(bgcol[b]);
	hist[mapname]->SetMarkerColor(bgcol[b]);
	hist[mapname]->SetFillColor(bgcol[b]);
	hist[mapname]->GetXaxis()->SetTitle(histox[i].c_str());
	hist[bgname]->GetXaxis()->SetTitle(histox[i].c_str());
      }

      else{
	
	hist[mapname] = (TH1F*)f -> Get(mapname.c_str());
	hist[bgname]-> Add(hist[mapname]);
	hist[mapname]->SetLineColor(bgcol[b]);
	hist[mapname]->SetMarkerColor(bgcol[b]);
	hist[mapname]->SetFillColor(bgcol[b]);
	hist[mapname]->GetXaxis()->SetTitle(histox[i].c_str());

	}
    }

    for(unsigned int b = 0; b<signames.size();++b){
      string mapname = histonames[i] + "_" + signames[b];
      hist[mapname] = (TH1F*)f->Get(mapname.c_str());//get histogram
      if(addgen){
        string signamegen = signames[b];
        signamegen.replace(0,6,"SignalGen");
        string mapname2 =  histonames[i] + "_" + signamegen;
        hist[mapname2] = (TH1F*)f->Get(mapname2.c_str());//get histogram
        hist[mapname]->Add(hist[mapname2],1.);
        hist[mapname]->Scale(0.5);
      }
      hist[mapname]->SetLineWidth(3);
      hist[mapname]->SetLineStyle(7);
      hist[mapname]->SetLineColor(sigcol[b]);
      hist[mapname]->SetMarkerColor(sigcol[b]);
      hist[mapname]->GetXaxis()->SetTitle(histox[i].c_str());
      }
    
  }


for(map<string,TH1F*>::iterator h=    hist.begin(); h!=    hist.end();++h) {
    //here do common styles
    hist[h->first]->GetXaxis()->SetLabelFont(42);
    hist[h->first]->GetXaxis()->SetLabelSize(0.04);
    hist[h->first]->GetXaxis()->SetTitleSize(0.05);
    hist[h->first]->GetXaxis()->SetTitleOffset(0.9);
    hist[h->first]->GetXaxis()->SetTitleFont(42);
    hist[h->first]->GetXaxis()->SetNdivisions(505);
    
    if(hist[h->first]->GetYaxis()->GetBinWidth(1)==1){
      hist[h->first]->GetYaxis()->SetTitle("events");
    }
    hist[h->first]->GetYaxis()->SetLabelFont(42);
    hist[h->first]->GetYaxis()->SetLabelSize(0.04);
    hist[h->first]->GetYaxis()->SetTitleSize(0.05);
    hist[h->first]->GetYaxis()->SetTitleOffset(1.2);
    hist[h->first]->GetYaxis()->SetTitleFont(42);
    hist[h->first]->GetZaxis()->SetLabelFont(42);
    hist[h->first]->GetZaxis()->SetLabelSize(0.035);
    hist[h->first]->GetZaxis()->SetTitleSize(0.035);
    hist[h->first]->GetZaxis()->SetTitleFont(42);
    //if(logy) hist[h->first]->SetMaximum(2.5*hist[h->first]->GetMaximum());
    //else     hist[h->first]->SetMaximum(1.25*hist[h->first]->GetMaximum());
    //cout << h->first << endl;
  }

  for(unsigned int i = 0; i<histonames.size(); ++i) {

    string bgname = "bgsum_" + histonames[i];
 
    float maximum = 0; float minimum = 0;
    string stackname = histonames[i];
    string axisname = histonames[i] + "_axis";
    stack[stackname] = new THStack();
    stack[stackname] -> SetName(stackname.c_str());

    for (int b = 0; b<bgnames.size(); ++b){
      string mapname =  histonames[i] + "_" + bgnames[b];
      stack[stackname] -> Add(hist[mapname]);
      if(hist[mapname] -> GetMinimum()>0) minimum += hist[mapname]->GetMinimum();
      if(hist[mapname] -> GetMaximum()>0) maximum += hist[mapname]->GetMaximum();

    }
 
    if(logy){//for plotting - define minimum maximum of y axis range
      minimum *=0.9;
      maximum *=100.;
      minimum = pow(10.0, floor(log10(minimum)));
      maximum = pow(10.0, ceil(log10(maximum)));
      if(minimum==0) minimum = 0.02;
      if(minimum>1&&minimum<=5) minimum = 0.2;
      if(minimum>5) minimum = 2;
    } else {
      minimum *=0.;
      maximum *=2.;
    }
    hist[bgname]->SetMaximum(maximum);
    hist[bgname]->SetMinimum(minimum);
    stack[stackname]->SetMaximum(maximum);
    stack[stackname]->SetMinimum(minimum);
    stack[stackname]->SetHistogram(hist[bgname]);
    
    if(histonames[i]=="pt_ISR_cut_dphi_lep_ptmiss_pt_ISR_pt_lep1_p4_R"){
      int numBins = hist[bgname] -> GetNbinsX();
      string signame = histonames[i] + "_Signal_T2tt_topcorridor";
      double err;
      double sg_err;
      double yield = hist[bgname] -> IntegralAndError(6, numBins, err);
      double sg_yield = hist[signame] -> IntegralAndError(6, numBins, sg_err);

      cout<<"---> "<<bgname<<">250 <---"<<endl;
      cout<<"BG Yield: "<<yield<<endl;
      cout<<"SG Yield: "<<sg_yield<<endl;
      cout<<"BG Error: "<<err<<endl;
      cout<<"SG Error: "<<sg_err<<endl;
      
      yield = hist[bgname] -> IntegralAndError(8, numBins, err);
      sg_yield = hist[signame] -> IntegralAndError(8, numBins, sg_err);
      cout<<"---> "<<bgname<<">350 <---"<<endl;
      cout<<"BG Yield: "<<yield<<endl;
      cout<<"SG Yield: "<<sg_yield<<endl;
      cout<<"BG Error: "<<err<<endl;
      cout<<"SG Error: "<<sg_err<<endl;

      yield = hist[bgname] -> IntegralAndError(10, numBins, err);
      sg_yield = hist[signame] -> IntegralAndError(10, numBins, sg_err);
      cout<<"---> "<<bgname<<">450 <---"<<endl;
      cout<<"BG Yield: "<<yield<<endl;
      cout<<"SG Yield: "<<sg_yield<<endl;
      cout<<"BG Error: "<<err<<endl;
      cout<<"SG Error: "<<sg_err<<endl;

      yield = hist[bgname] -> IntegralAndError(12, numBins, err);
      sg_yield = hist[signame] -> IntegralAndError(12, numBins, sg_err);
      cout<<"---> "<<bgname<<">550 <---"<<endl;
      cout<<"BG Yield: "<<yield<<endl;
      cout<<"SG Yield: "<<sg_yield<<endl;
      cout<<"BG Error: "<<err<<endl; 
      cout<<"SG Error: "<<sg_err<<endl;

      yield = hist[bgname] -> IntegralAndError(15, numBins, err);
      sg_yield = hist[signame] -> IntegralAndError(15, numBins, sg_err);
      cout<<"---> "<<bgname<<">700 <---"<<endl;
      cout<<"BG Yield: "<<yield<<endl;
      cout<<"SG Yield: "<<sg_yield<<endl;
      cout<<"BG Error: "<<err<<endl; 
      cout<<"SG Error: "<<sg_err<<endl;
      
    }

    /*if (histonames[i]=="pfmet_cut_n_loose_minus_good_dphi_lep_ptmiss_pt_ISR_pt_lep1_p4"||histonames[i]=="pfmet_cut_dphi_lep_ptmiss_n_of_first_bjet_pt_ISR_pt_lep1_p4"||histonames[i]=="pfmet_cut_dphi_lep_ptmiss_pt_ISR_pt_lep1_p4_R"){
      vector<int> met_bounds;
      met_bounds.push_back(1500);
      int numBins = hist[bgname] -> GetNbinsX();
      for(unsigned int j = 0; j<met_bin_fill.size(); ++j){
	//cout<<j<<endl;
	int flag = 0;
	for(unsigned int k = (met_bounds[met_bounds.size()-1]-150)/50; k>1; --k){
	  double err;
	  double yield = hist[bgname] -> IntegralAndError(k, (met_bounds[met_bounds.size()-1]-150)/50, err);
	  string signame = histonames[i] + "_Signal_T2tt_topcorridor";
	  double sg_err;
	  double sg_yield = hist[signame] -> IntegralAndError(k, (met_bounds[met_bounds.size()-1]-150)/50, sg_err);
	  //cout<<"help "<<hist[bgname] -> IntegralAndError(1, numBins, err);
	  if(yield<met_bin_fill[j]){
	    continue;
	  }
	  else if((met_bounds[met_bounds.size()-1]-150)/50==k){
	    continue;
	  }
	  else{
	    cout<<"--> "<<bgname<<" "<<(k-1)*50+150<<" to "<<met_bounds[met_bounds.size()-1]<<" <--"<<endl;
	    //cout << (met_bounds[met_bounds.size()-1]-150)/50 << " " << k << endl;
	    cout<<"BG Yield: "<<yield<<endl;
	    cout<<"SG Yield: "<<sg_yield<<endl;
	    cout<<"BG Error: "<<err<<endl;
	    cout<<"SG Error: "<<sg_err<<endl;
	    met_bounds.push_back((k-1)*50+150);

	    if(j==met_bin_fill.size()-1){
	      yield = hist[bgname] -> IntegralAndError(1, (met_bounds[met_bounds.size()-1]-150)/50, err);
	      sg_yield = hist[signame] -> IntegralAndError(1, (met_bounds[met_bounds.size()-1]-150)/50, sg_err);
	    cout << (met_bounds[met_bounds.size()-1]-150)/50 << " " << 1 << endl;
	      cout<<"--> "<<bgname<<" "<<250<<" to "<<met_bounds[met_bounds.size()-1]<<" <--"<<endl;
	      cout<<"BG Yield: "<<yield<<endl;
	      cout<<"SG Yield: "<<sg_yield<<endl;
	      cout<<"BG Error: "<<err<<endl;
	      cout<<"SG Error: "<<sg_err<<endl;
	      //for(int l = numBins; l>0; --l){
	      //cout<<"bg " << hist[bgname]->GetXaxis()->GetBinLowEdge(l) << "-" <<  hist[bgname]->GetXaxis()->GetBinLowEdge(numBins+1) << ": " << hist[bgname]->IntegralAndError(l,numBins,err) << " +/- " << err << endl;
	      //}
	    }

	    flag = 1;
	    break;;
	  }
	}
	if(flag==1) continue;
      }
      }*/

  }

  /*    for(unsigned int i = 0; i<histonames.size(); ++i){

      //now we are ready to draw the pretty picture - first define all needed quantities and then Draw
      TCanvas *c1 = new TCanvas("c1", "",334,192,600,600);//plots are done on a canvas
      c1->SetFillColor(0);
      c1->SetBorderMode(0);
      c1->SetBorderSize(2);
      //if(logy) c1->SetLogy();    // Log y
      c1->SetTickx(1);
      c1->SetTicky(1);
      c1->SetLeftMargin(0.18);
      c1->SetRightMargin(0.05);
      c1->SetTopMargin(0.07);
      c1->SetBottomMargin(0.15);
      c1->SetFrameFillStyle(0);
      c1->SetFrameBorderMode(0);
      c1->SetFrameFillStyle(0);
      c1->SetFrameBorderMode(0);
    
      TPad *plotpad = new TPad("plotpad", "Pad containing the overlay plot",0,0.165,1,1);//0,0.18,1,1);
      plotpad->Draw();
      plotpad->cd();
      plotpad->Range(-85.71429,-3.864499,628.5714,6.791402);//(133.1169,-3.101927,782.4675,0.7583922);
      plotpad->SetFillColor(0);
      plotpad->SetBorderMode(0);
      plotpad->SetBorderSize(2);
      if(logy) plotpad->SetLogy();
      plotpad->SetTickx(1);
      plotpad->SetTicky(1);
      plotpad->SetLeftMargin(0.12);
      plotpad->SetRightMargin(0.04);
      plotpad->SetTopMargin(0.05);
      // plotpad->SetBottomMargin(0.15);
      plotpad->SetFrameFillStyle(0);
      plotpad->SetFrameBorderMode(0);
      plotpad->SetFrameFillStyle(0);
      plotpad->SetFrameBorderMode(0);
    
      plotpad->cd();
      //TLatex *tLumi = new TLatex(0.95,0.944,"36.6 fb^{-1} (13 TeV)");
      string ls = Form("%f",lumi);//luminosity of the samples - see in ExampleLooper --> weight
      ls.erase ( ls.find_last_not_of('0') + 1, std::string::npos );
      ls.erase ( ls.find_last_not_of('.') + 1, std::string::npos );
      TLatex *tLumi = new TLatex(0.95,0.954,Form("%s fb^{-1} (13 TeV)",ls.c_str()));
      //TLatex *tLumi = new TLatex(0.95,0.944,"(13 TeV)");
      tLumi->SetNDC();
      tLumi->SetTextAlign(31);
      tLumi->SetTextFont(42);
      tLumi->SetTextSize(0.042);
      tLumi->SetLineWidth(2);
      TLatex *tECM = new TLatex(0.95,0.954,"(13 TeV)");
      //TLatex *tECM = new TLatex(0.95,0.944,"(13 TeV)");
      tECM->SetNDC();
      tECM->SetTextAlign(31);
      tECM->SetTextFont(42);
      tECM->SetTextSize(0.042);
      tECM->SetLineWidth(2);
      //tLumi->Draw();
      TLatex *tCMS = new TLatex(0.12,0.954,"CMS");
      tCMS->SetNDC();
      tCMS->SetTextAlign(11);
      tCMS->SetTextFont(61);
      tCMS->SetTextSize(0.0525);
      tCMS->SetLineWidth(2);
      //tCMS->Draw();
      TLatex *tSim = new TLatex(0.225,0.954,"Supplementary");
      tSim->SetNDC();
      tSim->SetTextAlign(11);
      tSim->SetTextFont(52);
      tSim->SetTextSize(0.042);
      tSim->SetLineWidth(2);
      TLatex *tPrel = new TLatex(0.225,0.954,"Preliminary");
      tPrel->SetNDC();
      tPrel->SetTextAlign(11);
      tPrel->SetTextFont(52);
      tPrel->SetTextSize(0.042);
      tPrel->SetLineWidth(2);
      TLegend *leg1 = new TLegend(0.2,0.735,0.5,0.925,NULL,"brNDC");//legend is always important so that reader knows what color/samples belong together
      leg1->SetBorderSize(0);
      leg1->SetTextSize(0.035);
      leg1->SetLineColor(1);
      leg1->SetLineStyle(1);
      leg1->SetLineWidth(2);
      leg1->SetFillColor(0);
      leg1->SetFillStyle(1001);

      for(unsigned int i = 0; i<bgnames.size(); ++i){
	leg1->AddEntry(hist[histonames[0]+"_"+bgnames[i] ], bgnames[i].c_str(),"f");
      }
      string stackname = histonames[i];
      stack[stackname] -> Draw("hist");

      TLegend *leg2 = new TLegend(0.5,0.735,0.85,0.925,NULL,"brNDC");//have 2 legends - one for background, another for signal+data(if available)
      leg2->SetBorderSize(0);
      leg2->SetTextSize(0.035);
      leg2->SetLineColor(1);
      leg2->SetLineStyle(1);
      leg2->SetLineWidth(2);
      leg2->SetFillColor(0);
      leg2->SetFillStyle(1001);

      for(unsigned int j = 0; j<signames.size(); ++j){

	hist[histonames[i] + "_" + signames[j] ] -> Draw("histsame");
	leg2 -> AddEntry(hist[histonames[i] + "_" + signames[j] ], sigleg[j].c_str(),"f");
      }
    
      //hist[stackname]->Draw("sameaxis");
      leg1->Draw();
      leg2->Draw();
      tCMS->Draw();
      tPrel->Draw();
      tLumi->Draw();
      
      c1->cd();
      TPad *ratiopad = new TPad("ratiopad", "Pad containing the ratio",0,0,1,0.16); //0,0,1,0.26);
      ratiopad->Draw();
      ratiopad->cd();
      //ratiopad->Range(-85.71429,-0.4,628.5714,2.266667);  //(133.1169,0.06923079,782.4675,1.607692);
      ratiopad->SetFillColor(0);
      ratiopad->SetBorderMode(0);
      ratiopad->SetBorderSize(2);
      ratiopad->SetTickx(1);
      ratiopad->SetTicky(1);
      ratiopad->SetLeftMargin(0.12);
      ratiopad->SetRightMargin(0.04);
      //  ratiopad->SetTopMargin(0.04);
      ratiopad->SetBottomMargin(0.2);
      ratiopad->SetFrameFillStyle(0);
      ratiopad->SetFrameBorderMode(0);
      ratiopad->SetFrameFillStyle(0);
      ratiopad->SetFrameBorderMode(0);
      float ratiomax = 0;
      for(unsigned int j = 0; j<signames.size(); ++j){
	hist[histonames[i] + "_" + signames[j]+"_Ratio" ]  = (TH1F*)hist[histonames[i] + "_" + signames[j] ]->Clone((histonames[i] + "_" + signames[j]+"_Ratio").c_str());
	hist[histonames[i] + "_" + signames[j]+"_Ratio" ]->Divide(hist["bgsum_" + histonames[i] ]);
	hist[histonames[i] + "_" + signames[j]+"_Ratio" ]->SetMinimum(0);
	if(hist[histonames[i] + "_" + signames[j]+"_Ratio" ]->GetMaximum()>ratiomax) ratiomax = hist[histonames[i] + "_" + signames[j]+"_Ratio" ]->GetMaximum();
	hist[histonames[i] + "_" + signames[j]+"_Ratio" ]->GetYaxis()->SetLabelSize(0.14);
	hist[histonames[i] + "_" + signames[j]+"_Ratio" ]->GetYaxis()->SetNdivisions(504);
	hist[histonames[i] + "_" + signames[j]+"_Ratio" ]->GetYaxis()->SetTitle("sig/bg");
	hist[histonames[i] + "_" + signames[j]+"_Ratio" ]->GetYaxis()->SetTitleSize(0.14);
	hist[histonames[i] + "_" + signames[j]+"_Ratio" ]->GetYaxis()->SetTitleOffset(0.28);
	hist[histonames[i] + "_" + signames[j]+"_Ratio" ]->GetXaxis()->SetLabelSize(0.0);
	hist[histonames[i] + "_" + signames[j]+"_Ratio" ]->GetXaxis()->SetTitleSize(0.0);

      }
      for(unsigned int j = 0; j<signames.size(); ++j){
	hist[histonames[i] + "_" + signames[j]+"_Ratio" ]->SetMaximum(TMath::Min(ratiomax*1.2,1.));
	if(j==0) hist[histonames[i] + "_" + signames[j]+"_Ratio" ]->Draw("hist");
	else     hist[histonames[i] + "_" + signames[j]+"_Ratio" ]->Draw("histsame");
      }
      string outname = "plots/" + stackname + ".pdf";
      c1->cd();
      c1->SaveAs(outname.c_str());//save the pretty picture
      c1->Clear();
        
     
  }*/

} //end function
