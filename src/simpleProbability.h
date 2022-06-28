#pragma once
#include <math.h>
#include <vector>
#include "JSL.h"
#include "compositeData.h"
namespace simplemodel
{
	// double BreakProbability(double factor,int position, int length, int chromosomeLength);

	double AnalyticalProbability(int n,int N, double p);
	double AnalyticalProbability(int n, int N, double p, int width);
	JSL::Vector AnalyticalProbability(JSL::Vector ns, int N, double p, int width);
	double LogLikelihood(const std::vector<int> & nk, double p, int binWidth);

	double LogLikelihood(const std::vector<BinCounter> & nk, double p);

	UncertainValue OptimalProbability(const std::vector<int> & nk, int binWidth);
	UncertainValue OptimalProbability(const std::vector<BinCounter> & nk);
};