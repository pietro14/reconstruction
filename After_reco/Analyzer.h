#ifndef ANALYZER_H
#define ANALYZER_H


#include <vector>


class TH2F;


class Analyzer {

	public:
	
	//constructors and destructor
	Analyzer();
	Analyzer(const Analyzer& source);
	Analyzer(const char* nometh2, int npixel ,TH2F* Tracklarge, int npixelorig );
	//Analyzer(const char* nometh2, int npixel ,vector<int> B, vector<int> E, int index);
	virtual ~Analyzer();

	//Simple returners
	double GetIntegral();
	double GetHeight();
	double GetRadius();
	TH2F* GetHisto();
	void GetCentre(double &x, double &y);

	//More complex operations
	void Reset();
	void SavetoFile(const char* nometh2);
	void Barycenter(double &XBar, double &YBar);
	double AngleLineMaxRMS(double &RMSOnLineVal);
	double RMSOnLine(double XBar, double YBar, double Phi);
	void IntegralAnalysis(double percrad);
	double MaxRMS();

	private:
	double fintegral;
	double fradius;
	double fheight;
	double fxcentr;
	double fycentr;	
	TH2F *fTrack;
	int fnpixelx;
	int fnpixely;

	//ClassDef(Analyzer,0)

};

#endif

