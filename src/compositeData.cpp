#include "compositeData.h"

UncertainValue::UncertainValue(double v, double e) : Value(v), Error(e)
{

}
std::string UncertainValue::Print()
{
	std::ostringstream os;
	os << *this;
	return os.str();
}
std::ostream & operator<<(std::ostream & os, const UncertainValue & uv)
{
	os << std::setprecision(2)<< uv.Value << "±" << uv.Error;
	return os;
}
bool UncertainValue::IsConsistent(double v)
{
	return ( DeviationsFromNorm(v) <= 1.0);
}
double UncertainValue::DeviationsFromNorm(double x)
{
	return abs(x - Value)/Error;
}

BinCounter::BinCounter(int start, int end)
{
	Start = start;
	End = end;
	Count = 0;
}
BinCounter::BinCounter(int start, int end, int pop)
{
	Start = start;
	End = end;
	Count = pop;
}

VectorPair DecompileBinArray(const std::vector<BinCounter> & input)
{
	int n = input.size();
	VectorPair vp;
	vp.x.resize(n);
	vp.y.resize(n);
	for (int i = 0; i < n; ++i)
	{
		vp.x[i] = (double)(input[i].Start + input[i].End -1)/2;
		vp.y[i] = input[i].Count;
	}
	return vp;
}

VectorPair VectorPair::CumulativeSum()
{
	double sum = 0;
	for (int i = 0; i < y.size(); ++i)
	{
		sum += y[i];
	}
	VectorPair p;
	p.x = x;
	p.y = y;
	for (int i = 1; i < p.y.size(); ++i)
	{
		p.y[i] += p.y[i-1];
		p.y[i-1] /= sum;
	}
	p.y[p.y.size()-1]/=sum;
	return p;
}