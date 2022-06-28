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
		ChromosomeArray(int binWidth, int positionCounts);

		void SyntheticInitalise(int genomeLength, double shatterProb);
		void FileInitialise(std::string targetFile, int chromosome, int qualityThreshold);
		


		UncertainValue Optimise();

		VectorPair BinCounts(bool relativeX);
		int UnbrokenLength;
		std::vector<BinCounter> SizeData;
	private:
		
		std::vector<int> DataIndex;
		std::vector<int> PositionData;

		int PositionBinCount;
		int SizeBinWidth;
		int TotalCatches;
		int MaximumBindex;

		void LogSize(int length);
		void LogPosition(int pos);
};