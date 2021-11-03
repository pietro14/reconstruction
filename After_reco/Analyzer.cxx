#include <cstdio>
#include <cmath>
#include <Riostream.h>
#include <TFile.h>
#include "Analyzer.h"
#include "TAxis.h"
#include <TH1F.h>
#include <TF1.h>
#include <TH2F.h>

//ClassImp(Analyzer)
double RMSOnLine(double XBar, double YBar, double Phi);

//costruttore di default
Analyzer::Analyzer():
fintegral(0.),
fradius(0.),
fheight(0.),
fnpixelx(0),
fnpixely(0),
fxcentr(0),
fycentr(0),
fTrack(NULL)
{
}

//costruttore standard
Analyzer::Analyzer(const char* nometh2, int npixel, TH2F *Tracklarge, int npixelorig) 
{
	fintegral=0.;
	fradius=0.;
	fheight=0.;
	fycentr=Tracklarge->GetMaximumBin()/(npixelorig+2);
	fxcentr=Tracklarge->GetMaximumBin()%(npixelorig+2);
	int minx=npixelorig+100;
	int miny=npixelorig+100;
	int maxx=-1;
	int maxy=-1;
	for(int i=fxcentr-npixel/2;i<fxcentr+npixel/2;i++)
	{
		for(int j=fycentr-npixel/2;j<fycentr+npixel/2;j++)
		{
			double x=Tracklarge->GetBinContent(i,j);
			if(i<minx && x>0)  minx=i;
			if(i>maxx && x>0)  maxx=i;
			if(j<miny && x>0)  miny=j;
			if(j>maxy && x>0)  maxy=j;
		}
	}

	fnpixelx=maxx-minx+1+10;
	fnpixely=maxy-miny+1+10;
	fTrack=new TH2F(nometh2,nometh2,fnpixelx,0,fnpixelx,fnpixely,0,fnpixely);
	for(int j=1;j<=fnpixelx;j++)
	{
		for(int k=1;k<=fnpixely;k++)
		{
			fTrack->SetBinContent(j,k,Tracklarge->GetBinContent(minx+j-1-5,miny+k-1-5));
		}
	}
	
	Barycenter(fxcentr,fycentr);
	
	for(int j=1;j<=fnpixelx;j++)
	{
		for(int k=1;k<=fnpixely;k++)
		{
			if(fTrack->GetBinContent(j,k)>0)
			{
				double x=sqrt((fycentr-k)*(fycentr-k) + (fxcentr-j)*(fxcentr-j) );
				if(x>fradius)  fradius=x;
			}
		}
	}
	
}

//copyconstructor
Analyzer::Analyzer(const Analyzer& source):
fintegral(source.fintegral),
fradius(source.fradius),
fheight(source.fheight),
fnpixelx(source.fnpixelx),
fnpixely(source.fnpixely),
fxcentr(source.fxcentr),
fycentr(source.fycentr)
{
  if(source.fTrack==NULL)
   {
	   fTrack=NULL;
   }
  else
  {
   TAxis *xb=source.fTrack->GetXaxis();
   double xmin=xb->GetXmin();
   double xmax=xb->GetXmax();
   int nbinsx=xb->GetNbins();
   TAxis *yb=source.fTrack->GetYaxis();
   double ymin=yb->GetXmin();
   double ymax=yb->GetXmax();
   int nbinsy=yb->GetNbins();
   fTrack = new TH2F("TheTrack'","Trackcopy'",nbinsx,xmin,xmax,nbinsy,ymin,ymax);
   for(Int_t i=1;i<=nbinsx;i++) 
   {
	   for(int j=1;j<=nbinsy;j++)	   fTrack->SetBinContent(i,source.fTrack->GetBinContent(i));
   }
  }

}

//distruttore
Analyzer::~Analyzer()
{
	if(fTrack!=NULL)
	{
	 delete fTrack;
	}
}



//restituisce the integral
double Analyzer::GetIntegral()
{
	return fintegral;
}


//restituisce il raggio
double Analyzer::GetRadius()
{
	return fradius;
}

//restituisce the peak height
double Analyzer::GetHeight()
{
	return fheight;
}

//returns centre coordinates
void Analyzer::GetCentre(double &x, double &y)
{
	x=fxcentr;
	y=fycentr;
}

//Saves histogram to root file
void Analyzer::SavetoFile(const char* nometh2)
{
	fTrack->SetName(nometh2);
	fTrack->Write();
	return;
}

//returns histogram pointer
TH2F* Analyzer::GetHisto()
{
	TH2F *copy=(TH2F *)fTrack->Clone();
	return copy;
}

//Resets the variables
void Analyzer::Reset()
{
	fintegral=0;
	fheight=0;
	fradius=0;
	fxcentr=0;
	fycentr=0;
	delete fTrack;
	fTrack=NULL;
}

//Calculates the barycentre of the track
void Analyzer::Barycenter(double &XBar, double &YBar)
{
  double Xb=0;
  double Yb=0;
  double Z=0;
  double ChargeTot=0;
  double min=0.;
  
  for(int i=1;i<fnpixelx;i++)
  {
	for(int j=1;j<fnpixely;j++)
	{
		if(fTrack->GetBinContent(i,j)<min) min=fTrack->GetBinContent(i,j);
	}
  }
  
  if(min>0) min=0;  
  
  for(int i=1;i<fnpixelx;i++)
  {
	for(int j=1;j<fnpixely;j++)
	{  
	  Z=fTrack->GetBinContent(i,j)+min;
	  if(Z!=0)
	  {
		  Xb+=(Z*i);
		  Yb+=(Z*j);
		  ChargeTot+=Z;
	  }
	}
  }
  
  Xb/=ChargeTot;
  Yb/=ChargeTot;

  XBar=Xb;
  YBar=Yb;
  
  return;
}


//Find main line (maximize RMS on this line)
double Analyzer::AngleLineMaxRMS(double &RMSOnLineVal)			//RMSOnLineVal is currently useless in thin code but seems useful for Samuele's code
{
  double XBar,YBar;
  double Sum1=0;
  double Sum2=0;
  double Z=0.;
  double Phi;
  double RmsAng;
  double RmsAngPerp;
  
  Barycenter(XBar,YBar);
  
  for(int i=1;i<fnpixelx;i++)
  {
	for(int j=1;j<fnpixely;j++)
	{  
	  Z=fTrack->GetBinContent(i,j);
	  if(Z!=0)
	  {
		  Sum1+= Z*(i-XBar)*(j-YBar);
		  Sum2+= Z*( (j-YBar)*(j-YBar) - (i-XBar)*(i-XBar)  );
	  }
	}
  }

  Phi=-0.5*TMath::ATan(2*Sum1/Sum2);
  
  RmsAng=RMSOnLine(XBar,YBar,Phi);
  RmsAngPerp=RMSOnLine(XBar,YBar,Phi+TMath::Pi()/2);
  
  if( RmsAng > RmsAngPerp )
  {  
      RMSOnLineVal=RmsAng;
      return Phi;
  } 
  else 
  {
      RMSOnLineVal=RmsAngPerp;
      if(Phi+TMath::Pi()/2>TMath::Pi()/2)     return Phi+TMath::Pi()/2-TMath::Pi();
      else     return Phi+TMath::Pi()/2;
  }
  
}
//Called by AngleLineMaxRMS
double Analyzer::RMSOnLine(double XBar, double YBar, double Phi)
{
  double RMS=0;
  double ChargeTot=0;
  double Z=0.;


  for(int i=1;i<fnpixelx;i++)
  {
	for(int j=1;j<fnpixely;j++)
	{  
	  Z=fTrack->GetBinContent(i,j);
	  if(Z!=0)
	  {
		  RMS+= Z*( (i-XBar)*cos(Phi) + (j-YBar)*sin(Phi) )*( (i-XBar)*cos(Phi) + (j-YBar)*sin(Phi) );
		  ChargeTot+=Z;
	  }
	}
  }
  
  return RMS/=ChargeTot;
}



//Calculates integral and height given the percentage of the radius to be used
void Analyzer::IntegralAnalysis(double percradius)
{
	fheight=0;
	fintegral=0;
	double Z=0.;
	TH1D *tax=fTrack->ProjectionX();
	tax->Fit("gaus","QR","",fxcentr-percradius*fradius,fxcentr+percradius*fradius);
	double localxcentr=tax->GetFunction("gaus")->GetParameter(1);
	TH1D *tay=fTrack->ProjectionY();
	tay->Fit("gaus","QR","",fycentr-percradius*fradius,fycentr+percradius*fradius);
	double localycentr=tay->GetFunction("gaus")->GetParameter(1);
	delete tax;
	delete tay;
	
	//calculates height
	for(int i=localxcentr-2;i<=localxcentr+2;i++)
	{
		for(int j=localycentr-2;j<=localycentr+2;j++)
		{
			fheight+=fTrack->GetBinContent(i,j);
		}
	}
	fheight/=25.;
	
	for(int i=1;i<fnpixelx;i++)
	{
		for(int j=1;j<fnpixely;j++)
		{  
			Z=fTrack->GetBinContent(i,j);
			if(Z!=0)
			{
				double dist=sqrt((fycentr-j)*(fycentr-j) + (fxcentr-i)*(fxcentr-i));
				if(dist<percradius*fradius)
				{
					fintegral+=fTrack->GetBinContent(i,j);  
				}
			}
		}
	}
	
	return;
}

//Function that calculates the RMS on the direction orthogonal to the main direction
double Analyzer::MaxRMS()
{
	double RMS=0.;
	
	return RMS;
	
}
