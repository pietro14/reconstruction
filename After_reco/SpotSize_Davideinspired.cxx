//This in particular compile using  g++ Analyzer.cxx Spotsize.cxx -o davedesize.exe `root-config --libs --cflags`
//Then use as ./davedesize.exe path_to_rootfile

#include <iostream>
#include <fstream>
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
#include <TF2.h>
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
  int nmax=7250000;
  int nscmax=450;
  int npixel=2304;
  int npixelsmall=250;
  
  float slimnesslimit=0.63;
  int xlimitmin=900;
  int xlimitmax=1400;
  int ylimitmin=850;
  int ylimitmax=1550;
  
  
  
  
  unsigned int nSc;
  int run;
  int event;
             //Pixels
  
  int *scID=new int[nmax];
  int *scIDall=new int[nmax];
  int ScNpixels[nscmax];
  float *XPix=new float[nmax];
  float *YPix=new float[nmax];
  float *ZPix=new float[nmax];
  int ScNpixelsall[nscmax];
  float *XPixall=new float[nmax];
  float *YPixall=new float[nmax];
  float *ZPixall=new float[nmax];
  float width[nscmax];
  float length[nscmax];
  float xpos[nscmax];
  float ypos[nscmax];
  float integral[nscmax];
  
  for(int con=0;con<nscmax;con++)
  {
	  ScNpixels[con]=0;
	  ScNpixelsall[con]=0;
	  width[con]=0;
	  length[con]=0;
	  xpos[con]=0;
	  ypos[con]=0;
  }  
  for(int con=0;con<nmax;con++)
  {
	  scID[con]=0;
	  scIDall[con]=0;
	  XPix[con]=0;
	  YPix[con]=0;
	  ZPix[con]=0;
	  XPixall[con]=0;
	  YPixall[con]=0;
	  ZPixall[con]=0;
  }  
  
  
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
  tree->SetBranchAddress("sc_ymean",xpos);
  tree->SetBranchAddress("sc_xmean",ypos);
  tree->SetBranchAddress("sc_integral",integral);

  /////////////////////////////////Analysis Variables ////////////////////////////////////////////////
  vector<int> BeginScPix;
  vector<int> EndScPix;
  vector<int> BeginScallPix;
  vector<int> EndScallPix;
  
  
  TH2F *Spotssumall=new TH2F("Spotssumall","Spotssumall",200,0,200,200,0,200);
  TH2F *Tracksmall=NULL;
  double xcentre;
  double ycentre;
  int rebinning=1;
  int nclustertot=0;
  /////////////////////////////////Analysis //////////////////////////////////////////////////////////
  
  ////Necessary for the class Analyzer to work
  TH2F* Track=new TH2F("Track","Track",npixel,0,npixel,npixel,0,npixel);					
  TH2F* Trackall=new TH2F("Trackall","Trackall",npixel,0,npixel,npixel,0,npixel);
  
  tree->GetEntry(0);
  ofstream foutsig("Sigma_values.txt",ios_base::app);
  foutsig<<"Run block "<<run<<"-";
  TFile* fout = new TFile(Form("Tracks%d.root",run),"recreate");
  ////
  
  
  for(int k=0;k<tree->GetEntries();k++)
  {
	  //Necessary structure to read the superclusters from reco_run tree
    
	tree->GetEntry(k);
	
	if(ScNpixelsall[0]>2304*1500 || nSc>40) 
	{
		cout << "Nev: "<< k << " skipped" << endl;
		continue;
	}
	cout << "Nev: "<< k << endl;
	
    
    ScIndicesElem(nSc,ScNpixels,&BeginScPix,&EndScPix);
    ScIndicesElem(nSc,ScNpixelsall,&BeginScallPix,&EndScallPix);
	
	//Start the cycle on the supercluster of the event
    for(int i=0;i<nSc;i++)
    {   
	   
	   /////////NOW THE TRACK WITH ALL THE PIXELS
	   
	   if(scIDall[BeginScallPix[i]]>=0)
	   {
		   if(width[i]/length[i]>slimnesslimit && width[i]/length[i]<1 && xpos[i]>xlimitmin && xpos[i]<xlimitmax && ypos[i]>ylimitmin && ypos[i]<ylimitmax /*&& integral[i]>5000 && integral[i]<9500*/)
		   {
			for(int j=BeginScallPix[i];j<EndScallPix[i];j++)
			{
				Trackall->SetBinContent(XPixall[j],YPixall[j],ZPixall[j]);
			}//chiudo for j (fill histos)
			
			
			Analyzer Traccia("Trackall",npixelsmall,Trackall,npixel);
			
			Tracksmall=(TH2F*)Traccia.GetHisto();
			//////choose the centre
			//Traccia.GetCentre(xcentre,ycentre);
			
			TH1F *hx=(TH1F*)Tracksmall->ProjectionX("hx",1,200/rebinning);
			TH1F *hy=(TH1F*)Tracksmall->ProjectionY("hy",1,200/rebinning);
			hx->Fit("gaus","q");
			xcentre=hx->GetFunction("gaus")->GetParameter(1)/rebinning;
			hy->Fit("gaus","q");
			ycentre=hy->GetFunction("gaus")->GetParameter(1)/rebinning;
				
			/*TF2 *f2 = new TF2("f2","[0]*TMath::Gaus(x,[1],[2])*TMath::Gaus(y,[3],[4])",0,100,0,100);
			f2->SetParameters(Tracksmall->GetBinContent(Tracksmall->GetMaximumBin()),xcentre,hx->GetFunction("gaus")->GetParameter(2),ycentre,hy->GetFunction("gaus")->GetParameter(2));
			Tracksmall->Fit("f2","q0");
			xcentre=Tracksmall->GetFunction("f2")->GetParameter(1);
			ycentre=Tracksmall->GetFunction("f2")->GetParameter(3);*/
			
			delete hx;
			delete hy;
			
			//////to rebin earlier
			int edgexmin=100/rebinning;
			int edgeymin=100/rebinning;
			int intxcentre=(int)xcentre;
			int intycentre=(int)ycentre;
			if(intxcentre-edgexmin<1) edgexmin=intxcentre;
			if(intycentre-edgeymin<1) edgeymin=intycentre;
			
			if(k%400==0)
			{
				double xcen,ycen;
				TCanvas *c1=new TCanvas(Form("Track%d_%d",k,i),Form("Track%d_%d",k,i),1000,700);
				
				Traccia.GetCentre(xcen,ycen);
				TMarker baricentre(xcen,ycen,22);
				baricentre.SetMarkerColor(kMagenta);
				baricentre.SetMarkerSize(1);
				
				TH1F *hx=(TH1F*)Tracksmall->ProjectionX("hx",1,200/rebinning);
				TH1F *hy=(TH1F*)Tracksmall->ProjectionY("hy",1,200/rebinning);
				hx->Fit("gaus","");
				xcen=hx->GetFunction("gaus")->GetParameter(1)/rebinning;
				hy->Fit("gaus","");
				ycen=hy->GetFunction("gaus")->GetParameter(1)/rebinning;
				TMarker singlegaus(xcen,ycen,24);
				singlegaus.SetMarkerColor(kRed);
				singlegaus.SetMarkerSize(1);
				
				TF2 *f2 = new TF2("f2","[0]*TMath::Gaus(x,[1],[2])*TMath::Gaus(y,[3],[4])",0,100,0,100);
				f2->SetParameters(Tracksmall->GetBinContent(Tracksmall->GetMaximumBin()),xcentre,hx->GetFunction("gaus")->GetParameter(2),ycentre,hy->GetFunction("gaus")->GetParameter(2));
				Tracksmall->Fit("f2","0");
				xcen=Tracksmall->GetFunction("f2")->GetParameter(1);
				ycen=Tracksmall->GetFunction("f2")->GetParameter(3);
				TMarker doublegaus(xcen,ycen,23);
				doublegaus.SetMarkerColor(kBlue);
				doublegaus.SetMarkerSize(1);
				delete f2;
				delete hx;
				delete hy;
				
				Tracksmall->Draw("COLZ");
				baricentre.Draw("same");
				singlegaus.Draw("same");
				doublegaus.Draw("same");
				
				c1->Write();
			}
			
			//cout<<k<<"   "<<i<<endl;
			
			//cout<<1+100-edgexmin<<"   "<<intxcentre-100+1<<"    "<<Tracksmall->GetNbinsX()<<endl;
			
			for(int lax=1+100/rebinning-edgexmin;lax<=200/rebinning && intxcentre-100/rebinning+lax<=Tracksmall->GetNbinsX();lax++)
			{
				for(int lay=1+100/rebinning-edgeymin;lay<=200/rebinning && intycentre-100/rebinning+lay<=Tracksmall->GetNbinsY();lay++)
				{
					Spotssumall->SetBinContent(lax,lay, Spotssumall->GetBinContent(lax,lay) + Tracksmall->GetBinContent( intxcentre-100/rebinning+lax , intycentre-100/rebinning+lay ) );
					//cout<< lax<<"   "<<lay<<"   "<<Spotssumall->GetBinContent(lax,lay) <<"     "<<Tracksmall->GetBinContent( intxcentre-100+lax , intycentre-100+lay )<<endl;
				}
			}
			
			nclustertot++;
			delete Tracksmall;
			Trackall->Reset();
		   }	
		}		//end if on allhit sc
	
	}	//end if on supercluster number
	
  }//chiudo for k (loop on the events)
  
  //Spotssumall->Scale(1./nclustertot);		//used to compare different algorithm
  Spotssumall->Write();
  
  TH1F *Spotssumall_x=(TH1F*)Spotssumall->ProjectionX("Spotssumall_x",1,200/rebinning);
  TH1F *Spotssumall_y=(TH1F*)Spotssumall->ProjectionY("Spotssumall_y",1,200/rebinning);
 
 
 
  /////////////versione doppia gaussiana unbinned
  double base, jump, mean, sx, sy, sx2, sy2, sigma, sigmaerr, sigma2, sigmaerr2, ratio, dratio, rx, ry;
  TF1 *fitfunc=new TF1("fitfunc","[0]*exp(-0.5*((x-[1])/[2])^2)+[3]+[5]*exp(-0.5*((x-[1])/[4])^2)",0,200);
  int centre;
  
  base=Spotssumall_x->GetBinContent(1);
  centre=Spotssumall_x->GetMaximumBin();
  jump= Spotssumall_x->GetBinContent(centre)-base;
  fitfunc->SetParameters( jump, centre , 10, base, 30,jump/10);
  fitfunc->SetParLimits(0,2000,1e7);
  fitfunc->SetParLimits(1,85/rebinning,115/rebinning);
  fitfunc->SetParLimits(2,0,25);
  fitfunc->SetParLimits(4,5,40);
  fitfunc->SetParLimits(5,0,1e6);
  Spotssumall_x->Fit("fitfunc","l");
  sx=fitfunc->GetParameter(2);
  sx2=fitfunc->GetParameter(4);
  if(sx>sx2)
  {
	  fitfunc->SetParameters( fitfunc->GetParameter(5), centre , sx2, base, sx,fitfunc->GetParameter(0));
	  Spotssumall_x->Fit("fitfunc","l");
	  sx=fitfunc->GetParameter(2);
	  sx2=fitfunc->GetParameter(4);
  }
  rx=fitfunc->GetParameter(0)/fitfunc->GetParameter(5);
  Spotssumall_x->Write();
  
  base=Spotssumall_y->GetBinContent(1);
  centre=Spotssumall_y->GetMaximumBin();
  jump= Spotssumall_y->GetBinContent(centre)-base;
  fitfunc->SetParameters( jump, centre, 10, base,30,jump/10);
  Spotssumall_y->Fit("fitfunc","l");
  sy=fitfunc->GetParameter(2);
  sy2=fitfunc->GetParameter(4);
  if(sy>sy2)
  {
	  fitfunc->SetParameters( fitfunc->GetParameter(5), centre , sy2, base, sy,fitfunc->GetParameter(0));
	  Spotssumall_y->Fit("fitfunc","l");
	  sy=fitfunc->GetParameter(2);
	  sy2=fitfunc->GetParameter(4);
  }
  ry=fitfunc->GetParameter(0)/fitfunc->GetParameter(5);
  Spotssumall_y->Write();
  
  sigma=(sx+sy)/2*rebinning;
  sigmaerr=abs(sx-sy)/2*rebinning;
  sigma2=(sx2+sy2)/2*rebinning;
  sigmaerr2=abs(sx2-sy2)/2*rebinning;
  ratio=(rx+ry)/2;
  dratio=abs(rx-ry)/2;
  
  foutsig<<run<<"\t\t\t"<<sigma<<" +/- "<<sigmaerr<<"\t\t\t\t\t"<<sigma2<<" +/- "<<sigmaerr2<<"\t\t\t\t\t"<<ratio<<" +/- "<<dratio<<endl;
  cout<<"Dimension pixels: "<<sigma<<" +/- "<<sigmaerr<<endl;
  cout<<"Dimension2 pixels: "<<sigma2<<" +/- "<<sigmaerr2<<endl;
  cout<<"Second gaussian (tales) is: ("<<100./(ratio+1)<<" +/- "<<100.*dratio/pow(ratio+1,2)<<") % of the total"<<endl;
  
  foutsig.flush();
  foutsig.close();
  
  fout->Save();
  fout->Close();
  delete Track;
  delete[] scID;
  delete[] XPix;
  delete[] YPix;
  delete[] ZPix;
  delete[] scIDall;
  delete[] XPixall;
  delete[] YPixall;
  delete[] ZPixall;
  
  return 0;
}

