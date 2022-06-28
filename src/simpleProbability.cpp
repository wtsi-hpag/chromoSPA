#include "simpleProbability.h"
#include "JSL.h"
namespace simplemodel
{

	double AnalyticalProbability(int n,int N,double p)
	{
		JSL::Assert("Probability meaningless for lengths < 1 and > N", n > 0 && n <= N);
		double prob = pow(1.0-p,n-1);
		
		if (n < N)
		{
			prob *= p;
		}
		return prob;
	}
	double AnalyticalProbability(int n,int N,double p, int width)
	{
		JSL::Assert("Probability meaningless for length " + std::to_string(n), n > 0 && n <= N);
		int bindex = n/width;
		int start = bindex * width;
		int endPoint = std::min(start+width,N);
		double Q = 1.0 - p;

		double prob = pow(Q,start -1) - pow(Q,endPoint - 1);
		if (endPoint >= N)
		{
			prob += pow(Q,N-1);
		}
		// std::cout << n << "  " << bindex << "  " << start << "  " << prob << std::endl;
		return prob;
	}
	JSL::Vector AnalyticalProbability(JSL::Vector ns,int N,double p, int width)
	{
		JSL::Vector out = ns;
		for (int i = 0; i < ns.Size(); ++i)
		{
			out[i] = AnalyticalProbability(ns[i],N,p,width);
			// double q= AnalyticalProbability(ns[i],N,p);
			
			
			out[i] *= (double)(N*p);
	
		
		}
		return out;
	}

	JSL::Vector CumulativeProbability(JSL::Vector ns,int N,double p, int width)
	{
		
		JSL::Vector out = JSL::Vector(ns.Size());
		double s = 0;
		double s2 = 0;
		int previousBindex = -1;
		for (int i = 0; i < ns.Size(); ++i)
		{
			
			int bindex = round(ns[i] / width);
			// std::cout << "New loop for length " << ns[i] << "  " << bindex << std::endl;
			if (i > 0)
			{
				out[i] = out[i-1];
			
			}
			if (bindex > previousBindex)
			{
				// std::cout << "am here" << std::endl;
				int tmp_bindex = previousBindex + 1;
				// std::cout << tmp_bindex << "  " << bindex << std::endl;
				while (tmp_bindex <= bindex)
				{
					int start = std::max(1,tmp_bindex * width);
					out[i] += AnalyticalProbability(start,N,p,width);
					
					// std::cout << "\t" << tmp_bindex << "  " << tmp_bindex * width + 1 << "  " << out[i] << std::endl;
					++tmp_bindex;
				}
				previousBindex = tmp_bindex -1;
			}
			// exit(10);
		
		}
		return out;
	}

	double Prior(double p)
	{
		if (p <= 0 || p >=1)
		{
			return -99999999999;
		}
		else
		{
			return 0;
		}
	}
	double likelihoodContribution(int nk, int k, double p, int binWidth)
	{
		return nk * ( (k-1) * log(1.0 - p) + log(binWidth * p) - 0.5*(binWidth - 1));
	}
	double LogLikelihood(const std::vector<int> & nk, double p, int binWidth)
	{
		double L = Prior(p);
		for (int i = 0; i < nk.size(); ++i)
		{
			int binStart = i*binWidth + 1;
			// int binEnd = binStart + binWidth;

			L += likelihoodContribution(nk[i],binStart,p,binWidth);		
		}
		return L;
	}

	double LogLikelihood(const std::vector<BinCounter> & nk, double p)
	{
		double L = Prior(p);
		for (int i = 0; i < nk.size(); ++i)
		{
			L += likelihoodContribution(nk[i].Count,nk[i].Start,p,nk[i].End - nk[i].Start + 1);
		}
		return L;
	}


	UncertainValue OptimalProbability(const std::vector<int> & nk, int binWidth)
	{
		double zeroMoment = 0;
		double oneMoment = 0;
		for (int i = 0; i < nk.size(); ++i)
		{
			int binStart = i*binWidth + 1;
			zeroMoment += nk[i];
			oneMoment += binStart * nk[i];
		}
		double pOpt = (2 * zeroMoment)/ (2 * oneMoment + (binWidth - 1) * zeroMoment);
		double errCompute = 0;
		for (int i = 0; i < nk.size(); ++i)
		{
			int binStart = i*binWidth + 1;
			errCompute += nk[i] * (1.0 + binStart * pOpt * pOpt - 2 * pOpt);
		}
		double err = pOpt * (1.0 - pOpt)/errCompute;
		
		return UncertainValue(pOpt,err);
	}
	UncertainValue OptimalProbability(const std::vector<BinCounter> & nk)
	{
		double zeroMoment = 0;
		double oneMoment = 0;
		double binWidth = nk[0].End - nk[0].Start + 1;
		for (int i = 0; i < nk.size(); ++i)
		{
			int binStart = nk[i].Start;
			int binEnd = nk[i].End;
		
			double idx = (binStart - 1)/binWidth;

			if (!nk[i].isFinalBin)
			{
				zeroMoment += nk[i].Count;
			}

			oneMoment += idx * nk[i].Count;
		
		}
		double qOptPow = oneMoment/(oneMoment + zeroMoment);
		double qOpt = pow(qOptPow,1.0/binWidth);
		double pOpt = 1.0 - qOpt;

		
		double errDenom = oneMoment * pow(1.0 - qOptPow,2) + zeroMoment * (qOptPow + binWidth - 1);
		double err = qOpt * (1.0 - qOptPow)/sqrt(binWidth * errDenom);
		
		return UncertainValue(pOpt,err);
	}
	
}