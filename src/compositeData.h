#pragma once
#include <ostream>
#include <iomanip>
#include <sstream>
#include <vector>
//Declare as structs to imply that they are 'simple' objects which are merely souped-up containers

struct VectorPair
{
	std::vector<double> x;
	std::vector<double> y;
	
	VectorPair CumulativeSum();
};

struct BreakData
{
	int SourceChromo;
	int BreakChromo;
	int SourcePosition;
	int BreakPosition;
	int Quality;
	BreakData(std::vector<std::string> inputLine);
};

struct UncertainValue
{
	double Value;
	double Error;

	UncertainValue(double v, double e);
	bool IsConsistent(double v);
	double DeviationsFromNorm(double x);
	std::string Print();
};

std::ostream & operator<<(std::ostream & os, const UncertainValue & uv);

struct BinCounter
{
	int Start;
	int End;
	int Count;
	bool isFinalBin;
	BinCounter(int start, int end);
	BinCounter(int start, int end, int pop);
};

VectorPair DecompileBinArray(const std::vector<BinCounter> & input);