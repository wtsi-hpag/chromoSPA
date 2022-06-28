#include "ChromosomeArray.h"

ChromosomeArray::ChromosomeArray(int sizeBinWidth, int positionBinCount)
{
	SizeBinWidth = sizeBinWidth;
	PositionBinCount = positionBinCount;
	PositionData.resize(PositionBinCount,0);

}

void ChromosomeArray::SyntheticInitalise(int length, double prob)
{
	UnbrokenLength = length;
	MaximumBindex = (length)/SizeBinWidth-1;
	// namespace model = simplemodel;
	TotalCatches = 0;
	int currentLength = 0;

	JSL::ProgressBar pb(length);
	pb.SetName("Initialising Synthetic Fragments");
	// int probCut = RAND_MAX * prob;
	for (int i = 1; i <= length; ++i)
	{
		++currentLength;
		double r = (double)rand()/RAND_MAX;
		// double p = model::BreakProbability(prob,i,currentLength,length);
		if (r < prob)
		{
			LogSize(currentLength);
			currentLength = 0;
		}
		pb.Update(i);
	}
	if (currentLength > 0)
	{
		LogSize(currentLength);
	}
	sort(SizeData.begin(), SizeData.end(), [](const BinCounter  & lhs, const  BinCounter & rhs) {
     return lhs.Start < rhs.Start;});
}

void ChromosomeArray::FileInitialise(std::string file, int chromosome, int qualityThreshold)
{
	//recover the total lengths of all chromosomes from file (inefficient to do every time, but only runs a few times)
	std::vector<long int> lengths;
	std::string lengthFiles = "Data/chromLengths.dat";
	forLineVectorIn(lengthFiles,' ',
		long int length = stoi(FILE_LINE_VECTOR[2]);
		lengths.push_back(length);
	)
	UnbrokenLength = lengths[chromosome - 1];

	MaximumBindex = (UnbrokenLength)/SizeBinWidth-1;
	//extract all breaks associated with this chromosome
	std::vector<BreakData> breaks;
	forLineVectorIn(file,' ',
		int source = stoi(FILE_LINE_VECTOR[1]);
		int quality = stoi(FILE_LINE_VECTOR[5]);
		if (source == chromosome && quality >= qualityThreshold)
		{
			breaks.push_back(BreakData(FILE_LINE_VECTOR));
		}
	);

	//sort the breaks based on their position in the unbroken chromosome - necessary because the datafile saves them in a slightly different order
	sort(breaks.begin(), breaks.end(), [](const BreakData  & lhs, const  BreakData & rhs) { return lhs.SourcePosition< rhs.SourcePosition;});


	int prevIdx = 0;
	for (int i = 0; i < breaks.size(); ++i)
	{
		int newIdx = breaks[i].SourcePosition;
		int length = newIdx - prevIdx;
		if (length > 0)
		{
			LogSize(length);
		}
		prevIdx = newIdx;
	}
	int finalLength = lengths[chromosome - 1] - prevIdx;
	if (finalLength > 0)
	{
		LogSize(finalLength);
	}


	sort(SizeData.begin(), SizeData.end(), [](const BinCounter  & lhs, const  BinCounter & rhs) {
     return lhs.Start < rhs.Start;});


}


void ChromosomeArray::LogSize(int length)
{
	int bindex = (length-1)/SizeBinWidth;
	if (bindex >= DataIndex.size())
	{
		DataIndex.resize(bindex + 1,-1);
	}

	int relocate = DataIndex[bindex];
	if (relocate != -1)
	{
		++SizeData[relocate].Count;
	}
	else
	{
		int binStart = bindex * SizeBinWidth + 1;
		SizeData.push_back(BinCounter(binStart,binStart +SizeBinWidth-1,1));
		DataIndex[bindex] = SizeData.size() - 1;

		if (bindex == MaximumBindex)
		{
			std::cout << length << " is in bindex " << bindex << ", max is " << MaximumBindex << std::endl;
			SizeData[SizeData.size()-1].isFinalBin = true;
		}
	}
	++TotalCatches;
}

void ChromosomeArray::LogPosition(int position)
{
	double pos = (double)position/UnbrokenLength;

	int bin = pos * PositionBinCount;
	++PositionData[bin];
}

UncertainValue ChromosomeArray::Optimise()
{
	return simplemodel::OptimalProbability(SizeData);
}

VectorPair ChromosomeArray::BinCounts(bool relative_x)
{
	
	VectorPair p = DecompileBinArray(SizeData);

	if (relative_x)
	{
		for (int i = 0; i < p.x.size(); ++i)
		{
			p.x[i]/=UnbrokenLength;
		}
	}
	return p;
}

