#include "ChunkArray.h"

ChunkArray::ChunkArray(int width)
{
	Width = width;
}

void ChunkArray::SyntheticInitalise(int length, double prob)
{
	UnbrokenLength = length;
	MaximumBindex = (length)/Width-1;
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
			IncrementCounter(currentLength);
			currentLength = 0;
		}
		pb.Update(i);
	}
	if (currentLength > 0)
	{
		IncrementCounter(currentLength);
	}
}


void ChunkArray::IncrementCounter(int length)
{
	int bindex = (length-1)/Width;
	if (bindex >= DataIndex.size())
	{
		DataIndex.resize(bindex + 1,-1);
	}

	int relocate = DataIndex[bindex];
	if (relocate != -1)
	{
		++Data[relocate].Count;
	}
	else
	{
		int binStart = bindex * Width + 1;
		Data.push_back(BinCounter(binStart,binStart +Width-1,1));
		DataIndex[bindex] = Data.size() - 1;

		if (bindex == MaximumBindex)
		{
			Data[Data.size()-1].isFinalBin = true;
		}
	}
	++TotalCatches;
}

double ChunkArray::Likelihood(double p)
{
	return simplemodel::LogLikelihood(Data,p);
}
UncertainValue ChunkArray::Optimise()
{
	return simplemodel::OptimalProbability(Data);
}

VectorPair ChunkArray::BinCounts(bool relative_x)
{
	sort(Data.begin(), Data.end(), [](const BinCounter  & lhs, const  BinCounter & rhs) {
     return lhs.Start < rhs.Start;});
	VectorPair p = DecompileBinArray(Data);

	if (relative_x)
	{
	for (int i = 0; i < p.x.size(); ++i)
	{
		p.x[i]/=UnbrokenLength;
	}
	}
	return p;
}

