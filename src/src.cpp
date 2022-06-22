
#include <iostream>
#include "JSL.h"

#include "ChunkArray.h"
#include "simpleProbability.h"


void SyntheticTest(int BinWidth, double ShatterProbability, int GenomeLength, bool RelativeCoords)
{
	ChunkArray ca(BinWidth);
	ca.SyntheticInitalise(GenomeLength,ShatterProbability);


	auto pp = ca.Optimise();
	auto vp = ca.BinCounts(RelativeCoords);
	double eVal = vp.x[vp.x.size()-1]*1.1;
	double mVal = vp.x[0]*0.9;


	double synthMax = eVal;
	if (RelativeCoords)
	{
		synthMax *= GenomeLength;
	}
	synthMax = std::min((int)synthMax +1, GenomeLength);
	auto predX = JSL::Vector::intspace(1,synthMax,std::max((int)(synthMax/10000),BinWidth));
	
	auto predY = simplemodel::AnalyticalProbability(predX,GenomeLength,pp.Value,BinWidth);
	predX = (predX-1) + (double)BinWidth/2;
	
	if (RelativeCoords)
	{
		predX = (predX / GenomeLength); 
	}
	
	auto cumsum = predY;
	for (int i = 1; i < predY.Size(); ++i)
	{
		cumsum[i] += cumsum[i-1];	
	}
	cumsum /= cumsum[cumsum.Size()-1];

	JSL::gnuplot gp;
	gp.SetSuperTitle("Synthetic Testing",20);
	// gp.SetFontSize(JSL::Fonts::SuperTitle,20);
	namespace lp = JSL::LineProperties;
	
	gp.SetMultiplot(2,1);
	gp.SetAxis(0);
	gp.Scatter(vp.x,vp.y,lp::Legend("Simulation"));
	gp.Plot(predX,predY,lp::Legend("Theory p=" + pp.Print()),lp::PenSize(2));
	gp.SetXLog(true);
	gp.SetYLog(true);
	// gp.SetYRange(1e-5,10);
	gp.SetXRange(mVal,eVal);
	std::string correctLabel = "Fragment Length";
	if (RelativeCoords)
	{
		correctLabel = "Fractional " + correctLabel;
	}
	gp.SetXLabel(correctLabel);
	gp.SetYLabel("Counts");
	gp.SetTitle("Absolute Occurence");

	gp.SetLegend(true);
	auto vp2 = vp.CumulativeSum();
	gp.SetAxis(1);
	gp.Scatter(vp2.x,vp2.y,lp::Legend("Simulation"));
	gp.Plot(predX,cumsum,lp::Legend("Theory p=" + pp.Print()),lp::PenSize(2));

	gp.SetXLog(true);
	gp.SetYLog(false);
	gp.SetXLabel(correctLabel);
	gp.SetYLabel("Cumulative Occurence");
	gp.SetXRange(mVal,eVal);
	gp.SetLegend(true);
	gp.Show();
	


}

int main(int argc, char ** argv)
{
	JSL::Argument<double> GenomeLength(100,"L",argc,argv);
	JSL::Argument<double> ShatterProbability(0.1,"p",argc,argv);
	JSL::Argument<int> seed(time(0),"s",argc,argv);
	JSL::Argument<std::string> SourceFile("__null__","f",argc,argv);
	JSL::Argument<int> BinWidth(1,"b",argc,argv);
	JSL::Argument<bool> RelativeCoords(true,"r",argc,argv);
	srand(seed);
	
	SyntheticTest(BinWidth,ShatterProbability,GenomeLength,RelativeCoords);

	return 0;
}
