//This in particular compile using  g++ Analyzer.cxx Spotsize.cxx -o spotsize.exe `root-config --libs --cflags`
//Then use as ./spotsize.exe path_to_rootfile

#include <iostream>
#include <string>
#include <vector>
#include "Analyzer.h"
#include <TTree.h>
#include <TFile.h>
#include <TArrayF.h>
#include <TH1.h>
#include <TMarker.h>
#include <TH2.h>
#include <TF1.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TGraphErrors.h>
#include <TSystem.h>
#include <TLine.h>
#include <TEllipse.h>
#include <TMath.h>

using namespace std;



void ScIndicesElem(int nSc,int* ScNelements, vector<int>* B, vector<int> *E){
  B->clear();
  E->clear();
  
  int parcount=0;
 
  for(int i=0;i<nSc;i++){
    B->push_back(parcount);
    E->push_back(parcount+ScNelements[i]);

    parcount+=ScNelements[i];
  }
}



int main(int argc, char** argv){

  ////////////////////////////////////Get File //////////////////////////////////////
  TFile* f = TFile::Open(Form("%s",argv[1]));
  TTree* tree = (TTree*)f->Get("Events");
  
  ///////////////////////////////////Set Branches and Define Variables////////////////////////////////////
  int nmax=250000;
  int nscmax=50;
  int npixel=2304;
  int npixelsmall=250;
  float slimnesslimit=0.6;
  
  unsigned int nSc;
  int run;
  int event;
             //Pixels
  
  int scID[nmax];
  int scIDall[nmax];
  int ScNpixels[nscmax];
  float XPix[nmax];
  float YPix[nmax];
  float ZPix[nmax];
  int ScNpixelsall[nscmax];
  float XPixall[nmax];
  float YPixall[nmax];
  float ZPixall[nmax];
  float width[nscmax];
  float length[nscmax];
  
  tree->SetBranchAddress("run",&run);
  tree->SetBranchAddress("event",&event);           
  tree->SetBranchAddress("nSc",&nSc);
  tree->SetBranchAddress("sc_ID",scID);
  tree->SetBranchAddress("sc_nintpixels",ScNpixels);
  tree->SetBranchAddress("sc_ypixelcoord",XPix);
  tree->SetBranchAddress("sc_xpixelcoord",YPix);
  tree->SetBranchAddress("sc_zpixel",ZPix);
  tree->SetBranchAddress("sc_IDall",scIDall);
  tree->SetBranchAddress("sc_nallintpixels",ScNpixelsall);
  tree->SetBranchAddress("sc_yallpixelcoord",XPixall);
  tree->SetBranchAddress("sc_xallpixelcoord",YPixall);
  tree->SetBranchAddress("sc_zallpixel",ZPixall);
  tree->SetBranchAddress("sc_width",width);
  tree->SetBranchAddress("sc_length",length);


  /////////////////////////////////Analysis Variables ////////////////////////////////////////////////
  vector<int> BeginScPix;
  vector<int> EndScPix;
  vector<int> BeginScallPix;
  vector<int> EndScallPix;
  
  TGraphErrors *integral_vs_radius=new TGraphErrors();
  TGraphErrors *height_vs_radius=new TGraphErrors();
  TGraphErrors *res_vs_radius=new TGraphErrors();
  TGraphErrors *integral_vs_radius_all=new TGraphErrors();
  TGraphErrors *integraloverradius_vs_radius_all=new TGraphErrors();
  TGraphErrors *height_vs_radius_all=new TGraphErrors();
  TGraphErrors *res_vs_radius_all=new TGraphErrors();
  TH1F *Original_radii=new TH1F("Original_radii","Original_radii",50,30,150);
  
  int points=20;
  double radii[points]={1.,0.9,0.8,0.7,0.6,0.5,0.45,0.4,0.38,0.35,0.3,0.25,0.22,0.2,0.19,0.18,0.17,0.15,0.1,0.05};
  TH1F* hintegral[points];
  for(int i=0;i<points;i++) hintegral[i]=new TH1F(Form("Histo_integral_rad%f",radii[i]),Form("Histo_integral_rad%f",radii[i]),100,0,50000);
  TH1F* hheight[points];
  for(int i=0;i<points;i++) hheight[i]=new TH1F(Form("Histo_height_rad%f",radii[i]),Form("Histo_height_rad%f",radii[i]),200,0,500);
  TH1F* hintegralall[points];
  for(int i=0;i<points;i++) hintegralall[i]=new TH1F(Form("Histo_integralall_rad%f",radii[i]),Form("Histo_integralall_rad%f",radii[i]),100,0,50000);
  TH1F* hheightall[points];
  for(int i=0;i<points;i++) hheightall[i]=new TH1F(Form("Histo_heightall_rad%f",radii[i]),Form("Histo_heightall_rad%f",radii[i]),200,0,500);
  
  /////////////////////////////////Analysis //////////////////////////////////////////////////////////
  int counter=0;
  int counterall=0;
  TH2F* Track=new TH2F("Track","Track",npixel,0,npixel,npixel,0,npixel);					
  TH2F* Trackall=new TH2F("Trackall","Trackall",npixel,0,npixel,npixel,0,npixel);
  TFile* fout = new TFile(Form("Tracks.root"),"recreate");

  for(int k=0;k<tree->GetEntries();k++)
  {
    tree->GetEntry(k);
    cout << "Nev: "<< k << endl;
   
    ScIndicesElem(nSc,ScNpixels,&BeginScPix,&EndScPix);
    ScIndicesElem(nSc,ScNpixelsall,&BeginScallPix,&EndScallPix);
	
	//Start the cycle on the supercluster of the event
    for(int i=0;i<nSc;i++)
    {   
	   if(scID[BeginScPix[i]]>=0)
	   {
		   if(width[i]/length[i]>0.6 && width[i]/length[i]<1)
		   {
			for(int j=BeginScPix[i];j<EndScPix[i];j++)
			{
				Track->SetBinContent(XPix[j],YPix[j],ZPix[j]);
			}  //chiudo for j (fill histos)
			
			//Track->SetName(Form("Track%i_run%i_evt%i",counter,run,event));
			//Track->SetTitle(Form("Track%i_run%i_evt%i",counter,run,event));
			//Track->Write();
			
			Analyzer Traccia("Track",npixelsmall,Track,npixel);
			Original_radii->Fill(Traccia.GetRadius());
			//Traccia.SavetoFile(Form("Track_small%i_run%i_evt%i",counter,run,event));
			
			for(int j=0;j<points;j++)
			{
				if(Traccia.GetRadius()<90)
				{
					Traccia.IntegralAnalysis(radii[j]);
					hintegral[j]->Fill(Traccia.GetIntegral());
					hheight[j]->Fill(Traccia.GetHeight());
				}
				
			}
			
			//Drawing the event
			if(event%40==0 || Traccia.GetRadius()>90)
			{
				TCanvas *c1=new TCanvas(Form("Canvas_run%i_evt%i",run,event),Form("Canvas_run%i_evt%i",run,event),1000,1000);
				TH2F *test=Traccia.GetHisto();
				double x,y,r;
				Traccia.GetCentre(x,y);
				r=Traccia.GetRadius();
				TEllipse la1(x,y,r*0.2,r*0.2);
				la1.SetLineWidth(2);
				la1.SetLineColor(kRed);
				la1.SetFillStyle(0);
				TLine *ltest=new TLine(x,y,x+r/sqrt(2),y+r/sqrt(2));
				//cout<<x<<"   "<<y<<"     "<<Traccia.GetRadius()<<endl;
				ltest->SetLineWidth(2);
				ltest->SetLineColor(kRed);
				//Traccia.Barycenter(x,y);
				double trash=0;
				double phi=Traccia.AngleLineMaxRMS(trash);
				TMarker m(x,y,20);
				m.SetMarkerColor(kBlack);
				m.SetMarkerSize(2);
				TF1 *mainline=new TF1("mainline","tan([0])*(x-[1])+[2]",0,250);
				mainline->SetParameters(phi,x,y);
				mainline->SetLineWidth(2);
				mainline->SetLineColor(kBlue);
				test->Draw("colz");
				ltest->Draw("same");
				m.Draw("same");
				mainline->Draw("same");
				la1.Draw("");
				c1->Write();
				delete test;
				delete mainline;
				delete ltest;
				delete c1;
			}
			////
			
			
			counter++;
			Track->Reset();
			//cout<<counter<<endl;
		   }	
		}		//end if on overthreshold sc
	   
	   
	   /////////NOW THE TRACK WITH ALL THE PIXELS
	   
	   if(scIDall[BeginScallPix[i]]>=0)
	   {
		   if(width[i]/length[i]>0.6 && width[i]/length[i]<1)
		   {
			for(int j=BeginScallPix[i];j<EndScallPix[i];j++)
			{
				Trackall->SetBinContent(XPixall[j],YPixall[j],ZPixall[j]);
			}//chiudo for j (fill histos)
			
			//Trackall->SetName(Form("Trackall%i_run%i_evt%i",counter,run,event));
			//Trackall->SetTitle(Form("Trackall%i_run%i_evt%i",counter,run,event));
			//Trackall->Write();
			
			Analyzer Traccia("Trackall",npixelsmall,Trackall,npixel);
			//Traccia.SavetoFile(Form("Trackall_small%i_run%i_evt%i",counterall,run,event));
			
			//Drawing the event
			if(event%40==0)
			{
				TCanvas *c1=new TCanvas(Form("Canvas_all_run%i_evt%i",run,event),Form("Canvas_all_run%i_evt%i",run,event),1000,1000);
				TH2F *test=Traccia.GetHisto();
				double x,y,r;
				Traccia.GetCentre(x,y);
				r=Traccia.GetRadius();
				TEllipse la1(x,y,r*0.2,r*0.2);
				la1.SetLineWidth(2);
				la1.SetLineColor(kRed);
				la1.SetFillStyle(0);
				TLine *ltest=new TLine(x,y,x+r/sqrt(2),y+r/sqrt(2));
				//cout<<x<<"   "<<y<<"     "<<Traccia.GetRadius()<<endl;
				ltest->SetLineWidth(2);
				ltest->SetLineColor(kRed);
				//Traccia.Barycenter(x,y);
				double trash=0;
				double phi=Traccia.AngleLineMaxRMS(trash);
				TMarker m(x,y,20);
				m.SetMarkerColor(kBlack);
				m.SetMarkerSize(2);
				TF1 *mainline=new TF1("mainline","tan([0])*(x-[1])+[2]",0,250);
				mainline->SetParameters(phi,x,y);
				mainline->SetLineWidth(2);
				mainline->SetLineColor(kBlue);
				test->Draw("colz");
				ltest->Draw("same");
				m.Draw("same");
				mainline->Draw("same");
				la1.Draw("");
				c1->Write();
				delete test;
				delete mainline;
				delete ltest;
				delete c1;
			}
			////
			for(int j=0;j<points;j++)
			{
				if(Traccia.GetRadius()<90)
				{
					Traccia.IntegralAnalysis(radii[j]);
					hintegralall[j]->Fill(Traccia.GetIntegral());
					hheightall[j]->Fill(Traccia.GetHeight());
				}
			}
			
			counterall++;
			Trackall->Reset();
			//cout<<counter<<endl;
		   }	
		}		//end if on allhit sc
	
	}	//end if on supercluster number
	
  }//chiudo for k (loop on the events)
  
  //Now I have the histograms full of info from the tracks and I do the analysis
  Original_radii->Fit("gaus","QR","",60,100);
  double avg_rad=Original_radii->GetFunction("gaus")->GetParameter(1);
  
  for(int i=0;i<points;i++)
  {
	  hintegral[i]->Fit("gaus","Q");
	  double a,b,c,d;
	  hintegral[i]->Write();
	  a=hintegral[i]->GetFunction("gaus")->GetParameter(1);
	  b=hintegral[i]->GetFunction("gaus")->GetParameter(2);
	  c=hintegral[i]->GetFunction("gaus")->GetParError(1);
	  d=hintegral[i]->GetFunction("gaus")->GetParError(2);
	  integral_vs_radius->SetPoint(i,radii[i],a);
	  integral_vs_radius->SetPointError(i,0,c);
	  
	  res_vs_radius->SetPoint(i,radii[i],b/a);
	  res_vs_radius->SetPointError(i,0,b/a*sqrt((c/a)*(c/a) +(d/b)*(d/b) ));
	  
	  hheight[i]->Fit("gaus","Q");
	  hheight[i]->Write();
	  height_vs_radius->SetPoint(i,radii[i],hheight[i]->GetFunction("gaus")->GetParameter(1));
	  height_vs_radius->SetPointError(i,0,hheight[i]->GetFunction("gaus")->GetParError(1));
	  
	  hintegralall[i]->Fit("gaus","Q");
	  hintegralall[i]->Write();
	  a=hintegralall[i]->GetFunction("gaus")->GetParameter(1);
	  b=hintegralall[i]->GetFunction("gaus")->GetParameter(2);
	  c=hintegralall[i]->GetFunction("gaus")->GetParError(1);
	  d=hintegralall[i]->GetFunction("gaus")->GetParError(2);
	  integral_vs_radius_all->SetPoint(i,radii[i],a);
	  integral_vs_radius_all->SetPointError(i,0,c);
	  integraloverradius_vs_radius_all->SetPoint(i,radii[i],a/(2*radii[i]*avg_rad));
	  integraloverradius_vs_radius_all->SetPointError(i,0,c/(2*radii[i]*avg_rad));
	  
	  res_vs_radius_all->SetPoint(i,radii[i],b/a);
	  res_vs_radius_all->SetPointError(i,0,b/a*sqrt((c/a)*(c/a) +(d/b)*(d/b) ));
	  
	  hheightall[i]->Fit("gaus","Q");
	  hheightall[i]->Write();
	  height_vs_radius_all->SetPoint(i,radii[i],hheightall[i]->GetFunction("gaus")->GetParameter(1));
	  height_vs_radius_all->SetPointError(i,0,hheightall[i]->GetFunction("gaus")->GetParError(1));
  }
  
  integral_vs_radius->SetMarkerStyle(20);
  integral_vs_radius->SetName(Form("integral_vs_radius"));
  integral_vs_radius->GetXaxis()->SetTitle(Form("Radius perc"));
  integral_vs_radius->GetYaxis()->SetTitle(Form("integral"));
  
  res_vs_radius->SetMarkerStyle(20);
  res_vs_radius->SetName(Form("Res_vs_radius"));
  res_vs_radius->GetXaxis()->SetTitle("Radius perc");
  res_vs_radius->GetYaxis()->SetTitle("Res");
  
  height_vs_radius->SetMarkerStyle(20);
  height_vs_radius->SetName(Form("Height_vs_radius"));
  height_vs_radius->GetXaxis()->SetTitle(Form("Radius perc"));  
  height_vs_radius->GetYaxis()->SetTitle(Form("height"));
  
  integral_vs_radius_all->SetMarkerStyle(20);
  integral_vs_radius_all->SetName(Form("integral_vs_radius_all"));
  integral_vs_radius_all->GetXaxis()->SetTitle(Form("Radius perc"));
  integral_vs_radius_all->GetYaxis()->SetTitle(Form("integral"));
  
  integraloverradius_vs_radius_all->SetMarkerStyle(20);
  integraloverradius_vs_radius_all->SetName(Form("integraloverradius_vs_radius_all"));
  integraloverradius_vs_radius_all->GetXaxis()->SetTitle(Form("Radius perc"));
  integraloverradius_vs_radius_all->GetYaxis()->SetTitle(Form("integral/radius (# / pixels)"));
  
  res_vs_radius_all->SetMarkerStyle(20);
  res_vs_radius_all->SetName(Form("Res_vs_radius_all"));
  res_vs_radius_all->GetXaxis()->SetTitle(Form("Radius perc"));
  res_vs_radius_all->GetYaxis()->SetTitle(Form("Res"));
  
  height_vs_radius_all->SetMarkerStyle(20);
  height_vs_radius_all->SetName(Form("Height_vs_radius_all"));
  height_vs_radius_all->GetXaxis()->SetTitle("Radius perc");  
  height_vs_radius_all->GetYaxis()->SetTitle("height");
  
  
  integral_vs_radius->Write();
  res_vs_radius->Write();
  height_vs_radius->Write();
  
  integral_vs_radius_all->Write();
  integraloverradius_vs_radius_all->Write();
  res_vs_radius_all->Write();
  height_vs_radius_all->Write();
  Original_radii->Write();
  
  TCanvas *resolutions=new TCanvas("resol","resol",1000,700);
  res_vs_radius->GetYaxis()->SetRangeUser(0.1,0.2);
  res_vs_radius->Draw("apl");
  res_vs_radius_all->SetMarkerColor(kRed);
  res_vs_radius_all->SetLineColor(kRed);
  res_vs_radius_all->Draw("samepl");
  TLine *old=new TLine(0.18,0.1,0.18,0.2);
  old->SetLineWidth(2);
  old->SetLineColor(kBlue);
  TLine *our=new TLine(0.38,0.1,0.38,0.2);
  our->SetLineWidth(2);
  our->SetLineColor(kGreen+2);
  old->Draw("same");
  our->Draw("same");
  resolutions->Write();
  
  TCanvas *ints=new TCanvas("ints","ints",1000,700);
  integral_vs_radius->Draw("apl");
  integral_vs_radius_all->SetMarkerColor(kRed);
  integral_vs_radius_all->SetLineColor(kRed);
  integral_vs_radius_all->Draw("samepl");
  TLine *old2=new TLine(0.18,0,0.18,33000);
  old2->SetLineWidth(2);
  old2->SetLineColor(kBlue);
  TLine *our2=new TLine(0.38,0,0.38,33000);
  our2->SetLineWidth(2);
  our2->SetLineColor(kGreen+2);
  old2->Draw("same");
  our2->Draw("same");
  ints->Write();
  
  fout->Save();
  fout->Close();
  delete Track;
  
  return 0;
}

