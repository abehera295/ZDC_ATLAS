#ifndef Tool_H_
#define Tool_H_

//#include "MulAna.h"
//#include "Rule.h"
#include "TGraphErrors.h"
#include "TFile.h"

class Tool
{
	private:
		TGraphErrors* gTrkEffScale[20][5];
		TGraph* gTrkEff[10][5];
		
	public:
		Tool();
		~Tool();
		double detTrkEff(double cent, double eta, double pt); // determine tracking efficiency
};

#endif



