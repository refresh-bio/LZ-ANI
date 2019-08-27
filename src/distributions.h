#pragma once
#include <cmath>
#define M_SQRT1_2  0.707106781186547524401  // 1/sqrt(2)


// ************************************************************************************
class IDistribution
{
public:
	virtual double pdf(double x) const = 0;
	virtual double cdf(double x) const = 0;
	virtual double invcdf(double x) const = 0;
};


// ************************************************************************************
// Normal distribution wrapper
class NormalDistribution : public IDistribution
{
protected:
	double meanX;
	double devX;
	double meanX2;

public:
	virtual double mean() const { return meanX; }
	virtual double dev() const { return devX; }

	/// Initializes mean and standard deviation
	NormalDistribution(double mean, double dev) : meanX(mean), devX(dev), meanX2(mean*mean) {}

	/// Probability density function
	virtual double pdf(double x) const override { return 0; }

	/// Cumulative distribution function (from Handbook of Mathematical Functions)
	virtual double cdf(double x) const override {
		// get standardized argument
		double sx = (x - meanX) / devX;
		double out = 0.5 * std::erfc(-sx * M_SQRT1_2);
		return out;
	}

	/// Inversed cumulative distribution function (from Handbook of Mathematical Functions)
	virtual double invcdf(double p) const override {

		if (p <= 0.0 || p >= 1.0) {
			throw std::invalid_argument("Invalid invcdf() argument");
		}

		double X = p < 0.5
			? -rationalApproximation(std::sqrt(-2.0 * std::log(p))) // F^-1(p) = - G^-1(p)
			: rationalApproximation(std::sqrt(-2.0 * std::log(1 - p))); 	// F^-1(p) = G^-1(1-p)

																			// de-standardize
		X = X * devX + meanX;
		return X;
	}

protected:
	double rationalApproximation(double t) const
	{
		// Abramowitz and Stegun formula 26.2.23.
		// The absolute value of the error should be less than 4.5 e-4.
		const double c[] = { 2.515517, 0.802853, 0.010328 };
		const double d[] = { 1.432788, 0.189269, 0.001308 };
		return t - ((c[2] * t + c[1])*t + c[0]) /
			(((d[2] * t + d[1])*t + d[0])*t + 1.0);
	}
};


// ************************************************************************************
// Binomial distribution - normal approximation
class BinomialDistributionApproximation : public NormalDistribution {
public:
	BinomialDistributionApproximation(uint32_t n, double p_success) : NormalDistribution(
		n * p_success, std::sqrt(n * p_success * (1.0 - p_success))
	) {}

};