//#ifndef FUNC_H
//#define FUNC_H

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TCanvas.h"
#include "TH1D.h"
#include "TFile.h"
#include "TH2D.h"
#include "TLatex.h"
#include "Math/Vector4D.h"
#include "TStyle.h"
#include <map>

#include "TDavixFile.h"


#include "TMath.h"
#include "TH1F.h"
#include "TROOT.h"
#include "THStack.h"
#include "TLegend.h"
#include "TF1.h"
#include "TSystem.h"

#include <math.h> // for “fabs”

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <vector>

//RooFit include
#include "RooDataHist.h"
#include "RooFit.h"
#include "RooRealVar.h"
#include "RooBreitWigner.h"
#include "RooCBShape.h"
#include "RooFFTConvPdf.h"
#include "RooAddPdf.h"
#include "RooChebychev.h"
#include "RooPlot.h"
#include "RooFitResult.h"

using namespace ROOT::VecOps;
using RNode = ROOT::RDF::RNode;

using namespace ROOT;                                                                                                                                           
using Vec_t = const ROOT::RVec<float>&;

using namespace RooFit;


float ComputeInvariantMass(Vec_t px, Vec_t py, Vec_t pz, Vec_t e) {                                                                                              
  ROOT::Math::PxPyPzEVector p1(px[0], py[0], pz[0], e[0]);                                                                                                     
  ROOT::Math::PxPyPzEVector p2(px[1], py[1], pz[1], e[1]);                                                                                                     
  return (p1 + p2).mass();                                                                                                                                     
}


RVecF ComputeEnergy(Vec_t px, Vec_t py, Vec_t pz, double m_e) {                                                                                                  
  RVecF energy;                                                                                                                                                
  energy[0]=ROOT::Math::sqrt(px[0]*px[0]+py[0]*py[0]+pz[0]*pz[0]+m_e*m_e);                                                                                     
  energy[1]=ROOT::Math::sqrt(px[1]*px[1]+py[1]*py[1]+pz[1]*pz[1]+m_e*m_e);                                                                                     
  energy[0]=gRandom->Gaus(energy[0],2);                                                                                                                        
  energy[1]=gRandom->Gaus(energy[1],2);                                                                                                                        
  return energy;                                                                                                                                               
}             
        


std::vector<double> myGetFitParameters(TH1D map_histog, double mean_bw, double input_width, double input_sigma, TString path, int m_sf){
                    
  std::vector<double> parameters;
  parameters.clear();


  double fit_lowcut = 80.;
  double fit_highcut = 100.;
  
  double r_fit_lowcut = 80.;
  double r_fit_highcut = 100.;
  

  //change the range of the plots to get the mean!!!
  RooRealVar x( "x", "x", r_fit_lowcut, r_fit_highcut);//84,98//80-100
  x.setBins(10000,"cache") ;
  x.setMin("cache",64.) ;
  x.setMax("cache",118.) ;

  //std::cout<<"mean_bw = "<<mean_bw<<std::endl;
  // Breit-Wigner                        
  RooRealVar m0( "m0", "m0", mean_bw, fit_lowcut, fit_highcut);//80-100
  RooRealVar width( "width", "width", input_width, 1., 4.);//2.49,1.,4.
  RooBreitWigner bw( "bw", "bw", x, m0, width);

  // Crystal-Ball                                                                                                                         
  RooRealVar mean( "mean", "mean", 0. );
  RooRealVar sigma( "sigma", "sigma", input_sigma, 1., 5.);//2.6,1.,5.
  RooRealVar alpha( "alpha", "alpha", 1.3 );
  RooRealVar n( "n", "n", 5.1 );
  RooCBShape cb( "cb", "cb", x, mean, sigma, alpha, n );
  // convolution                                                                                                                         
  RooFFTConvPdf pdf_sig( "pdf_sig", "pdf_sig", x, bw, cb );//pdf_sig
  //add polynomial bkg
  // Build Chebychev polynomial p.d.f.
  RooRealVar coef0("c0","coefficient #0",1.0,-.01,0.01);
  RooRealVar coef1("c1","coefficient #1",-0.1,-.01,0.01);
  RooRealVar coef2("c2","coefficient #2",-0.1,-.01,0.01);
  RooChebychev bkg1("bkg1","bkg1",x,RooArgList(coef0,coef1,coef2));
  RooRealVar fsig("fsig","signal fraction",0.9,0.,1.);
  RooAddPdf pdf("pdf","pdf",RooArgList(pdf_sig,bkg1),RooArgList(fsig));
  RooDataHist histo("histo","histo",x,Import(map_histog));

  x.setRange("signal",fit_lowcut,fit_highcut) ;
  
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
  ROOT::Math::MinimizerOptions::SetDefaultTolerance(0.0000001);
  ROOT::Math::MinimizerOptions::SetDefaultPrecision(0.0000001);
  
  pdf.fitTo(histo,SumW2Error(kTRUE), Range("signal")) ;//pdf
  
  TCanvas canv( "canv", "canv", 800., 600. );
  RooPlot* frame1 = x.frame(Bins(100),Title("Convolution of a Breit-Wigner and a Crystal-Ball, Chebychev pol. bkg")) ;
  histo.plotOn(frame1,Name("Data")) ;
  pdf.plotOn(frame1,Name("pdf"),LineColor(kRed)) ;
  pdf.paramOn(frame1,Layout(0.60));
  pdf.plotOn(frame1,Components("bkg1"),LineStyle(kDotted),LineColor(kBlue));

  TCanvas* canvas = new TCanvas("canvas","canvas",800,600) ;
  canvas->cd() ; 
  TPad *pad1 =  new TPad("pad1","pad1name",0.01,0.31,0.99,0.99);
  TPad *pad2 =  new TPad("pad2","pad2name",0.01,0.01,0.99,0.41);
  pad1->Draw();pad2->Draw();
  pad1->cd();
  pad1->SetBottomMargin(0.16);
  pad2->SetBottomMargin(0.24);
  frame1->GetYaxis()->SetTitleOffset(1.4) ;
  frame1->GetXaxis()->SetTitle((Form("m_{ee} [GeV], sf_%i ", m_sf)));
  frame1-> Draw();
  pad1->Modified();
  pad1->RedrawAxis();
  pad1->Update();
  pad2->cd();

  TF1* tf1_model = pdf.asTF(x);

  TH1* clone_data = (TH1*)histo.createHistogram("clone_data",x,Binning(50));
  RooDataHist *pdfHisto_data = pdf.generateBinned(x,1000000);
  TH1* clone_fit_data = (TH1*)pdfHisto_data->createHistogram("clone_fit_data",x,Binning(50));
  clone_fit_data->Scale(clone_data->Integral()/clone_fit_data->Integral());

  RooDataHist *pdfHisto = pdf_sig.generateBinned(x,1000000);
  TH1* clone_fit = (TH1*)pdfHisto->createHistogram("clone_fit",x,Binning(50));
  clone_fit->Scale(clone_data->Integral()/clone_fit->Integral());

  clone_data->Divide(clone_fit_data);

  int x1 = fit_lowcut;
  int x2 = fit_highcut;
  int bin1 = clone_data->FindBin(x1);
  int bin2 = clone_data->FindBin(x2);

  for(int i = 0; i < clone_data->GetNbinsX()+1; i++){
    if(i<bin1)
      clone_data->SetBinContent(i,0.);
    if(i>bin2)
      clone_data->SetBinContent(i,0.);
  }

  clone_data->GetXaxis()->SetTitle((Form("m_{ee} [GeV], sf_%i", m_sf)));
  clone_data->GetYaxis()->SetTitle("DATA / FIT");
  clone_data->GetXaxis()->SetRangeUser(fit_lowcut, fit_highcut);//80-100
  clone_data->GetYaxis()->SetRangeUser(0.7,1.3);
  clone_data->GetXaxis()->SetLabelSize(0.1);
  clone_data->GetYaxis()->SetLabelSize(0.08);
  clone_data->GetXaxis()->SetTitleSize(0.08);
  clone_data->GetYaxis()->SetTitleSize(0.09);
  clone_data->GetYaxis()->SetTitleOffset(0.6);
  clone_data->GetXaxis()->SetTitleOffset(1.2);
  clone_data->Draw("E1");
  pad2->Modified();
  pad2->SetGridy();
  pad2->RedrawAxis();
  pad2->Update();

  TString output_folder = path+"FitPlots";
  struct stat sb;
  if (stat(output_folder, &sb) == 0 && S_ISDIR(sb.st_mode)) {
    gSystem->cd(output_folder);    
  } 
  else {
    gSystem->mkdir(path+"FitPlots", kTRUE);
  }

  canvas->SaveAs(path+Form("FitPlots/sf_%i.pdf",m_sf));
  canvas->SaveAs(path+Form("FitPlots/sf_%i.root",m_sf));

  double modelMean = tf1_model->GetMaximumX();

  parameters.push_back(modelMean);//0
  parameters.push_back(m0.getError());//1
  parameters.push_back(m0.getVal());//2
  parameters.push_back(sigma.getVal());//3
  parameters.push_back(sigma.getError());//4
  parameters.push_back(mean.getVal());//5
  parameters.push_back(mean.getError());//6
  parameters.push_back(alpha.getVal());//7
  parameters.push_back(alpha.getError());//8
  parameters.push_back(n.getVal());//9
  parameters.push_back(n.getError());//10
  parameters.push_back(width.getVal());//11
  parameters.push_back(width.getError());//12

  cout << "width = "<<width.getVal() << endl;

  delete frame1;
  delete clone_data;
  delete clone_fit;
  delete pad1;
  delete pad2;
  delete canvas;

  return parameters;

}

