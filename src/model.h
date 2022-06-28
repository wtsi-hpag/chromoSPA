#pragma once
#include <math.h>
#include <vector>
#include "JSL.h"
#include "compositeData.h"
template <int ModelDimension=1>
class Model
{
	public:
		virtual double Probability(int length,const std::vector<double> &  modelParameters)
		{
			JSL::Assert("Model dimensions must be consistent",modelParameters.size() == ModelDimension);

			return 0;
		}

		virtual double BinnedProbability(int binStart, int binWidth, const std::vector<double> & modelParameters)
		{
			JSL::Assert("Model dimensions must be consistent",modelParameters.size() == ModelDimension);
			double s = 0;
			int top = std::min(binWidth,1000-binStart);
			for (int i = 0; i < binWidth; ++i)
			{
				s+=Probability(binStart + i,modelParameters);
			}
			return s;
		}
		virtual JSL::Vector BinnedProbability(const JSL::Vector & lengths, int binWidth, const std::vector<double> &  modelParameters)
		{
			JSL::Assert("Model dimensions must be consistent",modelParameters.size() == ModelDimension);
			int prevBin = -1;
			double prevValue = 0;
			JSL::Vector output(lengths.Size());
			for (int i = 0; i < lengths.Size(); ++i)
			{
				int bindex = (lengths[i]-1)/binWidth;
				int binStart = bindex * binWidth + 1;
				if (binStart != prevBin)
				{
					prevValue = BinnedProbability(binStart,binWidth, modelParameters);
					prevBin = binStart;
				}
				
				output[i] = prevValue;
				
			}
			return output;
		}
		virtual JSL::Vector CumulativeBinnedProbability(const JSL::Vector & lengths, int binWidth,  const std::vector<double> &  modelParameters)
		{
			JSL::Assert("Model dimensions must be consistent",modelParameters.size() == ModelDimension);
			
			JSL::Vector out(lengths.Size());
			int precomputedBins = -1;
			for (int i = 0; i < lengths.Size(); ++i)
			{
				int containingBin = (lengths[i] -1)/binWidth;
				if (i > 0)
				{
					out[i] = out[i-1];
				}
				if (containingBin > precomputedBins)
				{
					int tmp_bin = precomputedBins + 1;
					while (tmp_bin <= containingBin)
					{
						out[i] += BinnedProbability(tmp_bin*binWidth+1,binWidth,modelParameters);
						++tmp_bin;
					}
					precomputedBins = containingBin;
				}
			}
			return out;
		}

		virtual double LogLikelihoodInstance(int binStart, int dataCount, int binWidth, const std::vector<double> &  modelParameters)
		{
			// double p = BinnedProbability(binStart,binWidth,modelParameters);
			// return dataCount * log(p);
			
			int i = (int)round((double)(binStart-1)/binWidth);
			double Q = 1.0 - modelParameters[0];
			double p1 = i * binWidth * log(Q);
			double p2 = log(1.0 - pow(Q,binWidth));
			// std::cout << i << "  " << binStart << std::endl;
			return dataCount * (p1 + p2);
		}
		virtual double LogPrior(const std::vector<double> & modelParameters)
		{
			return 0;
		}
		virtual double LogLikelihood(const std::vector<BinCounter> & data, const std::vector<double> & modelParameters)
		{
			double s = LogPrior(modelParameters);
		
			for (int i = 0; i < data.size(); ++i)
			{
				s+= LogLikelihoodInstance(data[i].Start,data[i].Count,data[i].End-data[i].Start+1,modelParameters);
			}
			return s;
		}
		
};



class UniformBreakModel : public Model<1>
{
	public:
		int MaxLength;
		UniformBreakModel(int totalLength)
		{
			MaxLength = totalLength;
		}
		double Probability(int length,const std::vector<double> &  modelParameters)
		{

			double p = modelParameters[0];
			double Q = pow(1.0 - p,length-1);
			if (length > MaxLength)
			{
				return 0;
			}
			if (length < MaxLength)
			{
				return Q * p;
			}
			return Q;
		}

		double LogLikelihoodClever(const std::vector<BinCounter> & data, const std::vector<double> & modelParameters)
		{
			S1=0;S2=0;
			for (int i = 0; i < data.size(); ++i)
			{
				int W = (data[i].End - data[i].Start + 1);
				int bindex = (int)round( (double)(data[i].Start - 1)/W);
				S1+= data[i].Count * bindex;
				S2 += data[i].Count;
			}
			
			double p = modelParameters[0];
			double W  = data[0].End - data[0].Start + 1;
			double Q = pow(1.0 - p,W);
			
			int s = LogPrior(modelParameters) + S1 * W * log(1.0 - p) + S2 * log(1.0-Q);

			return s;
		}

	private:
		double S1;
		double S2;
};