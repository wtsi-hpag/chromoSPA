#pragma once
#include <string>
#include <iostream>
#include <vector>
#include "compositeData.h"
#include "JSL.h"
#include "simpleProbability.h"

class ChromosomeArray
{
	public:
		ChromosomeArray(int binWidth);

		void SyntheticInitalise(int genomeLength, double shatterProb);
		void FileInitialise(std::string targetFile, int chromosome, int qualityThreshold);
		
		void Print();
		void OrderedPrint();
		double Likelihood(double p);
		UncertainValue Optimise();

		VectorPair BinCounts(bool relativeX);
		int UnbrokenLength;
	private:
		std::vector<BinCounter> Data;
		std::vector<int> DataIndex;
		int Width;
		int TotalCatches;
		int MaximumBindex;
		
		void IncrementCounter(int length);
};