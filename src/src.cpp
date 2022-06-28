#include <unistd.h>
#include <iostream>

// #define GNUPLOT_NO_TIDY
#include "JSL.h"

#include "ChunkArray.h"
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
	auto predX = JSL::Vector::logintspace(1,synthMax,100000);
	std::cout << "Size print " << synthMax << "  " << predX.Size() << std::endl;
	auto predY = simplemodel::AnalyticalProbability(predX,GenomeLength,pp.Value,binWidth);
	predX = (predX-1) + (double)binWidth/2;
	
	if (coords)
	{
		predX = (predX / GenomeLength); 
	}
	
	auto cumsum = predY;
	for (int i = 1; i < predY.Size(); ++i)
	{
		cumsum[i] += cumsum[i-1];	
	}
	cumsum /= cumsum[cumsum.Size()-1];

	std::string correctLabel = "Fragment Length";
	if (coords)
	{
		correctLabel = "Fractional " + correctLabel;
	}
	gp.SetXLabel(correctLabel);
	gp.SetYLabel("Counts");


	gp.Scatter(vp.x,vp.y,lp::Legend(titleBase + " observed"));
	gp.Plot(predX,predY,lp::Legend(titleBase + " bestFit"),lp::PenType(JSL::Dash));
	gp.SetXLog(true);
	// gp.SetYLog(true);
	gp.SetLegend(true);
	double b = 0.5;
	gp.SetYRange(b,std::max(predY[0],b+3));
	gp.SetAxis(1,0);
	auto vp2 = vp.CumulativeSum();
	gp.Plot(vp2.x,vp2.y,lp::Legend(titleBase + " observed"));
	gp.Plot(predX,cumsum,lp::Legend(titleBase + " bestFit"),lp::PenType(JSL::Dash));
	gp.SetXLog(true);
	gp.SetXLabel(correctLabel);
	gp.SetYLabel("Cumulative Occurence");
	// gp.SetYLog(false);

	// gp.SetLegend(true);
}

void SyntheticTest(int BinWidth, double ShatterProbability, int GenomeLength, bool RelativeCoords)
{
	ChromosomeArray ca(BinWidth);
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
	auto predX = JSL::Vector::intspace(1,synthMax,std::max((int)(synthMax/1000),BinWidth));
	std::cout << predX << std::endl;
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
	// gp.Show();

}

void RealTest(int binWidth, std::string file,bool coords, int quality)
{
	JSL::gnuplot gp;
	gp.SetMultiplot(2,1);
	
	int chromosome = 1;

	for (int chromosome = 1; chromosome < 2; ++chromosome)
	{
		ChromosomeArray ca(binWidth);
		
		ca.FileInitialise(file,chromosome,quality);
		std::string base = "Chrom" + std::to_string(chromosome);
		LengthAnalysis(ca,binWidth,coords,gp,base);
	}
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
	JSL::Argument<int> Quality(3,"q",argc,argv);	
	srand(seed);
	
	RealTest(BinWidth,SourceFile,RelativeCoords,Quality);
	
	return 0;
}
