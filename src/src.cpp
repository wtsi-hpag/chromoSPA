#include <unistd.h>
#include <iostream>

// #define GNUPLOT_NO_TIDY
#include "JSL.h"

#include "ChromosomeArray.h"
#include "simpleProbability.h"


void LengthAnalysis(ChromosomeArray & ca, int binWidth, bool coords, JSL::gnuplot & gp, std::string titleBase)
{
	namespace lp = JSL::LineProperties;
	auto pp = ca.Optimise();
	std::cout << pp << std::endl;
	gp.SetAxis(0,0);
	auto vp = ca.BinCounts(coords);
	double eVal = vp.x[vp.x.size()-1]*1.1;
	double mVal = vp.x[0]*0.9;
	double synthMax = eVal;
	int GenomeLength = ca.UnbrokenLength;
	if (coords)
	{
		synthMax *= GenomeLength;
	}
	synthMax = std::min((int)synthMax +1, GenomeLength);

	
	auto predX = JSL::Vector::logintspace(1,synthMax,10000);
	auto predY = simplemodel::AnalyticalProbability(predX,GenomeLength,pp.Value,binWidth);
	auto cumsum = simplemodel::CumulativeProbability(predX,GenomeLength,pp.Value,binWidth);
	predX = (predX-1) + (double)binWidth/2;

	if (coords)
	{
		predX = (predX / GenomeLength); 
	}
	
	


	std::string correctLabel = "Fragment Length";
	if (coords)
	{
		correctLabel = "Fractional " + correctLabel;
	}
	gp.SetXLabel(correctLabel);
	gp.SetYLabel("Counts");

	auto plotY = predY;
	for (int i = 0; i < predY.Size(); ++i)
	{
		if (predY[i] < 1e-2)
		{
			plotY[i] = 0;
		}
	}

	gp.Scatter(vp.x,vp.y,lp::Legend(titleBase + " observed"));
	gp.Plot(predX,plotY,lp::Legend(titleBase + " bestFit"),lp::PenType(JSL::Dash),lp::Colour("hold"));
	gp.SetXLog(true);
	gp.SetYLog(true);
	gp.SetLegend(true);
	


	gp.SetAxis(1,0);
	auto vp2 = vp.CumulativeSum();
	gp.Plot(vp2.x,vp2.y,lp::Legend(titleBase + " observed"));
	gp.Plot(predX,cumsum,lp::Legend(titleBase + " bestFit"),lp::PenType(JSL::Dash),lp::Colour("hold"));
	// gp.Scatter(predX,cumsum,lp::Legend(titleBase + " bestFit"),lp::Colour("hold"));
	gp.SetXLog(true);
	gp.SetXLabel(correctLabel);
	gp.SetYLabel("Cumulative Occurence");
	// gp.SetYLog(true);

	// gp.SetLegend(true);
}

void PositionAnalysis(ChromosomeArray & ca,JSL::gnuplot & gp,std::string titleBase)
{

}

void SyntheticTest(int BinWidth, double ShatterProbability, int GenomeLength, bool RelativeCoords, int posBins)
{
	ChromosomeArray ca(BinWidth,posBins);
	ca.SyntheticInitalise(GenomeLength,ShatterProbability);
	JSL::gnuplot gp;
	gp.SetMultiplot(2,1);
	LengthAnalysis(ca,BinWidth,RelativeCoords,gp,"Synthetic");
	gp.Show();
}

void RealTest(int binWidth, std::string file,bool coords, int quality, int posBins)
{
	JSL::gnuplot gp;
	gp.SetMultiplot(2,1);
	

	std::vector<int> studies = {1,6};

	for (auto chromosome : studies)
	{
		ChromosomeArray ca(binWidth,posBins);
		
		ca.FileInitialise(file,chromosome,quality);
		std::string base = "Chrom" + std::to_string(chromosome);
		LengthAnalysis(ca,binWidth,coords,gp,base);
	}
	gp.Show();
}

int main(int argc, char ** argv)
{

	JSL::Argument<double> GenomeLength(100,"L",argc,argv);
	JSL::Argument<double> ShatterProbability(0.001,"p",argc,argv);
	JSL::Argument<int> seed(time(0),"s",argc,argv);
	JSL::Argument<std::string> SourceFile("__null__","f",argc,argv);
	JSL::Argument<int> BinWidth(1,"b",argc,argv);
	JSL::Argument<bool> RelativeCoords(true,"r",argc,argv);
	JSL::Argument<int> Quality(3,"q",argc,argv);	
	JSL::Argument<int> PositionBins(100,"x",argc,argv);
	srand(seed);
	
	
	// SyntheticTest(synthBins,ShatterProbability,GenomeLength,RelativeCoords,PositionBins);
	RealTest(BinWidth,SourceFile,RelativeCoords,Quality,PositionBins);
	
	return 0;
}
