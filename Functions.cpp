////////////////////////////////////////////////////
// Functions.cpp
// Sample functions implementation file for example AmiBroker plug-in
//
// Copyright (C)2001 Tomasz Janeczko, amibroker.com
// All rights reserved.
//
// Last modified: 2001-09-24 TJ
// 
// You may use this code in your own projects provided that:
//
// 1. You are registered user of AmiBroker
// 2. The software you write using it is for personal, noncommercial use only
//
// For commercial use you have to obtain a separate license from Amibroker.com
//
////////////////////////////////////////////////////

#define _SCL_SECURE_NO_WARNINGS 1
#define _USE_MATH_DEFINES

#include <math.h>
#include <functional>
#include <numeric>

#include "StdAfx.h"
#include "Plugin.h"
#include "OptionPricing.h"

//=============================================================================
// Private Data
//

// Default Values (all zero)
static float ZeroValues[] =
{
	0.0,				// dividend date 1 (0 => no dividend)
	0.0,				// dividend amount 1
	0.0,				// dividend date 2 (0 => no dividend)
	0.0,				// dividend amount 2
	1.0,				// Early Exercise? (non-zero => yes)
	0.0 				// Minimum periods (zero => use the
						// number of days to expiry)
};

static float DefaultValuesVol[] =
{
	256.0,				// default to 1 year volatilityLookBack period
};

//=============================================================================
// Helper functions
//

// Copy "empty" values from Src to Dst until we get a "real" value
static int SkipEmptyValues(int nSize, float *Src, float *Dst)
{
	int i;

	for (i = 0; i < nSize && IS_EMPTY(Src[i]); i++)
	{
		Dst[i] = EMPTY_VAL;
	}

	return i;
}

///////////////////////////////////////////
// Each AFL function has the following prototype:
//
// AmiVar FunctionName(int NumArgs, AmiVar *ArgsTable)
//
// You can define as many functions as you want.
// To make them visible you should add them to the function
// table at the bottom of this file
//
///////////////////////////////////////////////


AmiVar TVV(int NumArgs, AmiVar* ArgsTable)
{
	AmiVar				R;				// return value
	float*				Rv;				// Result vector
	float*				Sv;				// Share price vector
	float*				Tv;				// on This date vector
	float				K;				// striKe price
	long				E;				// Expiry date
	bool				CorP;			// true for Call, false for Put
	float				v;				// annualised Volatility
	float				r;				// Risk free return rate
	int					nDividends = 0;	// no. elements in ...
	DiscreteDividend	vDividends[2];	// the discrete dividends for the share
	bool				ee;				// is Early Exercise possible?
	int					minPeriods;		// minimum number of periods to use

	R = gSite.AllocArrayResult();

	Rv = R.array;
	Sv = ArgsTable[0].array;					// array 1
	Tv = ArgsTable[1].array;					// array 2
	K  = ArgsTable[2].val;						// float 1
	E  = abDateToJDN((long) ArgsTable[3].val);	// float 2
	CorP = (ArgsTable[4].val != 0);				// float 3
	v = ArgsTable[5].val;						// float 4
	r = ArgsTable[6].val;						// float 5
	if (ArgsTable[7].val != 0)					// optional floats 1 & 2
	{
		vDividends[nDividends].date = abDateToJDN((long) ArgsTable[7].val);
		vDividends[nDividends].amount = ArgsTable[8].val;
		++nDividends;
	}
	if (ArgsTable[9].val != 0)					// optional floats 3 & 4
	{
		vDividends[nDividends].date = abDateToJDN((long) ArgsTable[9].val);
		vDividends[nDividends].amount = ArgsTable[10].val;
		++nDividends;
	}
	ee = (ArgsTable[11].val != 0);				// optional float 5
	minPeriods = (long) ArgsTable[12].val;		// optional float 6

	int					n = gSite.GetArraySize();
	int					i;
	int					j;

	j = SkipEmptyValues(n, Sv, Rv);

    for (i = j; i < n; ++i)
    {
		Rv[i] = TV1day(Sv[i],
					   abDateToJDN((long) Tv[i]),
					   K,
					   E,
					   CorP,
					   v,
					   r,
					   nDividends,
					   vDividends,
					   ee,
					   minPeriods);
	}

    return R;
}

struct Greeks {
	float delta;
	float gamma;
	float vega;
	float theta;
	float rho;	
};

AmiVar execFormula(int NumArgs, AmiVar* ArgsTable, std::function<float (float, float, long, bool, float, float, Greeks&)> f)
{
	AmiVar				R;				// return value
	float*				Rv;				// Result vector
	float*				Sv;				// Share price vector
	float*				Tv;				// on This date vector
	float				K;				// striKe price
	long				E;				// Expiry date
	bool				CorP;			// true for Call, false for Put
	float				v;				// annualised Volatility
	float				r;				// Risk free return rate

	AmiVar				delta;			// greeks arrays
	AmiVar				gamma;			
	AmiVar				vega;			
	AmiVar				theta;
	AmiVar				rho;

	delta	= gSite.AllocArrayResult();
	gamma	= gSite.AllocArrayResult();
	vega	= gSite.AllocArrayResult();
	theta	= gSite.AllocArrayResult();
	rho		= gSite.AllocArrayResult();

	R		= gSite.AllocArrayResult();

	Rv		= R.array;
	Sv		= ArgsTable[0].array;					// array 1
	Tv		= ArgsTable[1].array;					// array 2
	K		= ArgsTable[2].val;						// float 1
	E		= abDateToJDN((long) ArgsTable[3].val);	// float 2
	CorP	= (ArgsTable[4].val != 0);				// float 3
	v		= ArgsTable[5].val;						// float 4
	r		= ArgsTable[6].val;						// float 5

	int n	= gSite.GetArraySize();
	int	i;
	int	j;

	j = SkipEmptyValues(n, Sv, Rv);


    for (i = j; i < n; ++i)
    {
		auto currentPrice		= Sv[i];
		auto currentDate		= abDateToJDN( (long) Tv[i]);
		auto daysToExpiration	= E - currentDate;
		Greeks greeks;
		Rv[i]					= f (currentPrice, K, daysToExpiration, CorP, v, r, greeks);

		delta.array[i]	= greeks.delta;
		gamma.array[i]	= greeks.gamma;
		vega.array[i]	= greeks.vega;
		theta.array[i]	= greeks.theta;
		rho.array[i]	= greeks.rho;
	}

	gSite.SetVariable("bsDelta",	delta);
	gSite.SetVariable("bsGamma",	gamma);
	gSite.SetVariable("bsVega",		vega);
	gSite.SetVariable("bsTheta",	theta);
	gSite.SetVariable("bsRho",		rho);

    return R;
}

inline float w(float l) {
	static const float a1 = 0.31938153f;
	static const float a2 = -0.356563782f;
	static const float a3 = 1.781477937f;
	static const float a4 = -1.821255978f;
	static const float a5 = 1.330274429f;
	static const float b = 0.2316419f;
	static const float c = 1.0f / (float)std::sqrt(2.0f * M_PI);

    static const float k = 1.0f / (1.0f + b * 1.0f);
    static const float k2 = k * k;
    static const float k4 = k2 * k2;
    return c * std::exp(-0.5f * l * l) * (a1 * k + a2 * k2 + a3 * k * k2 + (a4 * k4 + a5 * k * k4));
}

inline float cnd(float x) {
    if(x < 0.0f)
        return w (-x);
    else
        return 1.0f - w (x); 
}

inline float snd(float x) { return float (std::exp(- (std::pow(x,2.0f)) / 2.0f) / std::sqrt(2.0f * M_PI));}

auto blackScholesEuro(float price, float strike, long days, bool CorP, float v, float r, Greeks& greeks) -> float {
      float s			= price, x = strike; long t = days;
      float sqrtt		= (float) std::sqrt((float)t);
      float d1			= (std::log(s / x) + (r + v * v / 2.0f) * t) / (v * sqrtt);
      float d2			= d1 - v * sqrtt;
      x					= x * std::exp(-r * t);
      float ert			= std::exp(- r * t);
      float snd1		= snd(d1);

      greeks.gamma		= snd1 / (s * v * sqrtt);
      greeks.vega		= s * snd1 * sqrtt;
      
	  float cnd1		= cnd (d1);
      float cnd2		= cnd (d2);
      
	  if(CorP) {
		greeks.delta = cnd1;
		greeks.theta	= - (s * snd1 * v) / (2.0f * sqrtt) - r * x * ert * cnd2;
		greeks.rho		= x * t * ert * cnd2;
        return			s * cnd1 - x * cnd2;
      } else {
        float cndm2		= cnd (-d2);
		greeks.delta	= cnd1 - 1.0f;
		greeks.theta	= - (s * snd1 * v) / (2.0f * sqrtt) + r * x * ert * cndm2;
		greeks.rho		= -x * t * ert * cndm2;
        return			x * cndm2 - s * cnd (-d1);
      };
}

AmiVar bsValueEuropean(int NumArgs, AmiVar* ArgsTable)
{
	return execFormula(NumArgs, ArgsTable, blackScholesEuro);
}


float calcHistVol(float* prices, int currentBar, int lookBack) {
	auto startPrice	= currentBar - lookBack;
	if(startPrice < 0) return EMPTY_VAL;

	auto sum	= std::accumulate(&prices[startPrice], &prices[currentBar], 0.0f);
	auto mean	= sum / lookBack;
	
	auto sq_sum	= std::inner_product(&prices[startPrice], &prices[currentBar], &prices[startPrice], 0.0f);
	auto stdev	= std::sqrt(sq_sum / lookBack - mean * mean);
	return stdev;
}

AmiVar bsHistVol(int NumArgs, AmiVar* ArgsTable)
{
	AmiVar				R;				// return value
	float*				Rv;				// Result vector
	float*				Sv;				// Share price vector
	float				lookBack;		// Lookback period in days

	R = gSite.AllocArrayResult();

	Rv = R.array;
	Sv = ArgsTable[0].array;					// array 1
	lookBack = ArgsTable[1].val;				// float 1

	int					n = gSite.GetArraySize();
	int					i;
	int					j;

	j = SkipEmptyValues(n, Sv, Rv);

    for (i = j; i < n; ++i)
    {
		Rv[i] = calcHistVol(Sv, i, (int) lookBack);
	}

    return R;
}

AmiVar IVV(int NumArgs, AmiVar *ArgsTable)
{
	AmiVar				R;				// return value
	float*				Rv;				// Result vector
	float*				Ov;				// Option price vector
	float*				Sv;				// Share price vector
	float*				Tv;				// on This date vector
	float				K;				// striKe price
	long				E;				// Expiry date
	bool				CorP;			// true for Call, false for Put
	float				v;				// annualised Volatility
	float				r;				// Risk free return rate
	int					nDividends = 0;	// no. elements in ...
	DiscreteDividend	vDividends[2];	// the discrete dividends for the share
	bool				ee;				// is Early Exercise possible?
	int					minPeriods;		// minimum number of periods to use

	R = gSite.AllocArrayResult();

	Rv = R.array;
	Ov = ArgsTable[0].array;					// array 1
	Sv = ArgsTable[1].array;					// array 2
	Tv = ArgsTable[2].array;					// array 3
	K  = ArgsTable[3].val;						// float 1
	E  = abDateToJDN((long) ArgsTable[4].val);	// float 2
	CorP = (ArgsTable[5].val != 0);				// float 3
	v = ArgsTable[6].val;						// float 4
	r = ArgsTable[7].val;						// float 5
	if (ArgsTable[8].val != 0)					// optional floats 1 & 2
	{
		vDividends[nDividends].date = abDateToJDN((long) ArgsTable[8].val);
		vDividends[nDividends].amount = ArgsTable[9].val;
		++nDividends;
	}
	if (ArgsTable[10].val != 0)					// optional floats 3 & 4
	{
		vDividends[nDividends].date = abDateToJDN((long) ArgsTable[10].val);
		vDividends[nDividends].amount = ArgsTable[11].val;
		++nDividends;
	}
	ee = (ArgsTable[12].val != 0);				// optional float 5
	minPeriods = (long) ArgsTable[13].val;		// optional float 6

	int					n = gSite.GetArraySize();
	int					i;
	int					j;

	j = SkipEmptyValues(n, Sv, Rv);

    for (i = j; i < n; ++i)
    {
		Rv[i] = IV1day(Ov[i],
					   Sv[i],
					   abDateToJDN((long) Tv[i]),
					   K,
					   E,
					   CorP,
					   v,
					   r,
					   nDividends,
					   vDividends,
					   ee,
					   minPeriods);
	}

    return R;
}


/////////////////////////////////////////////
// Function table now follows
//
// You have to specify each function that should be
// visible for AmiBroker.
// Each entry of the table must contain:
// "Function name",
//	{ FunctionPtr,
//			<no. of array args>,
//			<no. of string args>,
//			<no. of float args>,
//			<no. of default args>,
//			<pointer to default values table float *>

FunctionTag gFunctionTable[] =
{
	//"mdaTVV",	{ TVV, 2, 0, 5, 6, ZeroValues },
	//"mdaIVV",	{ IVV, 3, 0, 5, 6, ZeroValues },
	"bsVEuro",	{ bsValueEuropean, 2, 0, 5, 0, NULL },
	"bsHistVol",{ bsHistVol, 1, 0, 0, 1, DefaultValuesVol }
};

int gFunctionTableSize = sizeof(gFunctionTable)/sizeof(FunctionTag);
