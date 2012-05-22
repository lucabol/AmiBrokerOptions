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

#define _CRT_SECURE_NO_WARNINGS
#define _SCL_SECURE_NO_WARNINGS 1
#define _USE_MATH_DEFINES

#include <math.h>
#include <functional>
#include <numeric>

#include "boost/date_time/gregorian/gregorian.hpp"

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

// Convert abDate to boost date
inline boost::gregorian::date abDateToDate(long abDate)
{
    int					dd;
    int					mm;
    int					yy;

	// Extract the day, month, and year from the given date
    dd = (abDate % 100);
    mm = (abDate / 100) % 100;
    yy = (abDate / 10000) + 1900;
	return boost::gregorian::date(yy, mm, dd);
}

// Convert from boost date to abDate
inline long dateToAbDate(boost::gregorian::date d) {
	return 10000 * (d.year() - 1900) + 100 * d.month() + d.day();
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


AmiVar execFormula(int NumArgs, AmiVar* ArgsTable, std::function<double (double, double, double, bool, double, double, Greeks&)> f)
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
		auto daysToExpiration	= (E - currentDate);
		Greeks greeks;
		Rv[i]					= (float) f (currentPrice, K, daysToExpiration, CorP, v, r, greeks);

		delta.array[i]	= (float)greeks.delta;
		gamma.array[i]	= (float)greeks.gamma;
		vega.array[i]	= (float)greeks.vega;
		theta.array[i]	= (float)greeks.theta;
		rho.array[i]	= (float)greeks.rho;
	}

	gSite.SetVariable("bsDelta",	delta);
	gSite.SetVariable("bsGamma",	gamma);
	gSite.SetVariable("bsVega",		vega);
	gSite.SetVariable("bsTheta",	theta);
	gSite.SetVariable("bsRho",		rho);

    return R;
}

inline double w(double l) {
	static const auto a1    = 0.31938153;
	static const auto a2    = -0.356563782;
	static const auto a3    = 1.781477937;
	static const auto a4    = -1.821255978;
	static const auto a5    = 1.330274429;
	static const auto b     = 0.2316419;
	static const auto c     = 1.0 / std::sqrt(2.0 * M_PI);

    auto k                  = 1.0 / (1.0 + b * l);
    auto k2                 = k * k;
    auto k4                 = k2 * k2;
    auto value              = c * std::exp(-0.5 * l * l) * (a1 * k + a2 * k2 + a3 * k * k2 + (a4 * k4 + a5 * k * k4));
    return value;
}

inline double cnd(double x) {
    if(x < 0.0)
        return w (-x);
    else
        return 1.0 - w (x); 
}

inline double snd(double x) { return std::exp(- (std::pow(x,2.0)) / 2.0) / std::sqrt(2.0 * M_PI);}

// Black Scholes formula, to calculate volatility surface read doc below formula 12. Just use terms to 3rd power with p.34 table 4 values middle table
// http://finance.business.queensu.ca/psfile/DaglishHullSuoRevised.pdf, F = forward value of S for a contract maturying at T -> F = S*exp(r * (T -t)) , K = strike price, S = asset price
extern "C"
auto blackScholesEuro(double price, double strike, double days, bool CorP, double v, double r, Greeks& greeks) -> double {
      auto s			= price, x = strike; auto t = days / 365.0;
      auto sqrtt		= std::sqrt(t);
      auto d1			= (std::log(s / x) + (r + v * v / 2.0) * t) / (v * sqrtt);
      auto d2			= d1 - v * sqrtt;
      x					= x * std::exp(-r * t);
      auto ert			= std::exp(- r * t);
      auto snd1		    = snd(d1);

      greeks.gamma		= snd1 / (s * v * sqrtt);
      greeks.vega		= s * snd1 * sqrtt;
      greeks.vega       = greeks.vega / 100.0;
      
	  auto cnd1		    = cnd (d1);
      auto cnd2		    = cnd (d2);
      
	  if(CorP) {
		greeks.delta = cnd1;
		greeks.theta	= - (s * snd1 * v) / (2.0 * sqrtt) - r * x * ert * cnd2;
		//greeks.theta	= (s * snd1 * v) / (2.0 * sqrtt) + r * x * ert * cnd2;
		greeks.theta	= greeks.theta / 365.0;
		greeks.rho		= x * t * ert * cnd2;
		greeks.rho		= greeks.rho / 100.0;
        auto value      = s * cnd1 - x * cnd2;
        return value;		
      } else {
        auto cndm2		= cnd (-d2);
		greeks.delta	= cnd1 - 1.0;
		greeks.theta	= - (s * snd1 * v) / (2.0 * sqrtt) + r * x * ert * cndm2;
		greeks.theta	= greeks.theta / 365.0;
		greeks.rho		= -x * t * ert * cndm2;
		greeks.rho		= greeks.rho / 100.0;
        return			x * cndm2 - s * cnd (-d1);
      };
}

AmiVar bsValueEuropean(int NumArgs, AmiVar* ArgsTable)
{
	return execFormula(NumArgs, ArgsTable, blackScholesEuro);
}



void print(boost::gregorian::date d) {
	std::stringstream ss;
	ss << d;
	std::cout << ss.str().c_str() << std::endl;
}

// Returns the expiry date some months out. 0 = this month, 1 = next month, 2 = two months out, etc...
boost::gregorian::date calcExpiry(boost::gregorian::date d, int expiries) {
	using namespace boost::gregorian;

	typedef nth_day_of_the_week_in_month nth_dow;
	months single(expiries);
	date dTarget = d + single;
	nth_dow ndm(nth_dow::third, Friday, dTarget.month());
	auto value = ndm.get_date(dTarget.year());
	print(value);
	return  value;
}

extern "C"
char* testExpiry(int year, int month, int day, int expiries) {
	auto d		= boost::gregorian::date(year, month, day);
	auto r		= calcExpiry(d, expiries);

	std::stringstream ss;
	ss << r;
	static char tmp[10];
	strcpy(tmp, ss.str().c_str());
	return tmp;
}

inline long round(double number) {     return number < 0.0 ? (long) ceil(number - 0.5) : (long) floor(number + 0.5); }

AmiVar bsExpiry(int NumArgs, AmiVar* ArgsTable)
{
	AmiVar			R;				// return value
	float*			Rv;				// Result vector
	float*			Tv;				// on This date vector
	float			expiries;		// 1 = firstExpirty, 2 = secondExpiry, ...

	R				= gSite.AllocArrayResult();
	Rv				= R.array;
	Tv				= gSite.CallFunction("DateNum", 0, NULL).array;
	
	int	n			= gSite.GetArraySize();
	expiries		= ArgsTable[0].val;
	int	i;
	int j;

	j				= SkipEmptyValues(n, Tv, Rv);

    for (i = j; i < n; ++i)
    {
		auto currentDate	= abDateToDate(round(Tv[i]));
		auto d				= calcExpiry(currentDate, round(expiries));
		Rv[i]				= (float) dateToAbDate(d);
	}
	return R;
}

bool condor(
			boost::gregorian::date now,	// date to evaluate 
			double spot,				// spot price underlying
			double v,					// ATM volatility
			double r,					// risk free rate
			double step,				// % of spot price to keep as distance between wings, also distance between high open interest strikes
			int minCallShortDist,		// min distance from the short call strike in steps
			int minPutShortDist,		// min distance from the short put strike in steps
			int minDays,				// min number of days to expiry
			int maxDays,				// max number of days to expiry
			double maxDelta,			// max acceptable delta value for shorts in steps
			double minPremium,			// min accepted premium as % of step
			Condor& ret					// return value
			)
{
	// convert params to dollar signs
	auto stepPr			= round(step * spot);
	auto toUSD			= [stepPr] (double x) { return round(stepPr * x);};
	auto minCpr			= toUSD( minCallShortDist );
	auto minPpr			= toUSD( minPutShortDist );
	auto premiumPr		= toUSD( minPremium );

	// calc strike values for short legs
	auto atm			= round(spot / stepPr) * (long) stepPr;
	auto callShort		= atm + minCpr;
	auto putShort		= atm - minPpr;
	
	auto addDays		= [](boost::gregorian::date d, int dys) -> boost::gregorian::date {
		using namespace boost::gregorian;

		auto toAdd		= days(dys);
		auto dTarget	= d + toAdd;
		print(dTarget);
		return  dTarget;
	};

	// calc min & max allowed expiry dates
	auto minDate		= addDays(now, minDays);
	auto maxDate		= addDays(now, maxDays);
	auto expiry			= calcExpiry(now, 0);

	// find first good expiry
	while(expiry < minDate)
		expiry			= calcExpiry(expiry, +1);

	Greeks g;
	auto scholes		= [spot, v, r, &g] (double strike, int days, bool CorP) {
		return blackScholesEuro(spot, strike, days, CorP, v, r, g);
	};

	// find a condor that works at this expiry
	auto findCondor		= [callShort, putShort, spot, stepPr, v, r, minCpr, minPpr, scholes, maxDelta, premiumPr, &g, &ret] (int days) -> bool {
		auto cStrike	= callShort;
		auto pStrike	= putShort;

		while(true) {
			// it would be faster to 'continue' after each one of these conditions is not meet, but less clear
			auto shCall			= scholes(cStrike, days, true);
			auto cDelta			= g.delta;
			auto shPut			= scholes(pStrike, days, false);
			auto pDelta			= g.delta;
			auto lgCall			= scholes(cStrike + stepPr, days, true);
			auto lgPut			= scholes(pStrike - stepPr, days, false);
			auto net			= shCall + shPut - lgCall - lgPut;

			if(cDelta <= maxDelta && (-pDelta) <= maxDelta && net >= premiumPr) {
				// it's a good condor
				ret.longCallStrike		= lgCall;
				ret.shortCallStrike		= shCall;
				ret.shortPutStrike		= shPut;
				ret.longPutStrike		= lgPut;
				ret.netPremium			= net;
				return true;
			} else {
				// it has too big deltas, but the premium is good enough to continue the search at this expiry
				if(net > premiumPr) {
					cStrike					= cStrike + stepPr;
					pStrike					= pStrike - stepPr;
				} else {
					// premium less than required, no point continuing at this expiry, it just gets worse (true?)
					return false;
				}
			}
		}
	};

	// increases the expiry until it finds a condor or the expiry is too far out
	while (expiry < maxDate) {
		auto days		= (expiry - now).days();
		if(findCondor(days))
			return true;
		expiry			= calcExpiry(expiry, +1);
	}

	return false;
}

extern "C"
bool testCondor(int year, int month, int day,double spot,double v,double r,double step,int minCallShortDist,int minPutShortDist,
				int minDays,int maxDays,double maxDelta,double minPremium,Condor& ret) {
	return condor(boost::gregorian::date(year, month, day), spot, v, r, step, minCallShortDist, minPutShortDist, minDays, maxDays, maxDelta, minPremium, ret);
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
	"bsVEuro",	{ bsValueEuropean, 2, 0, 5, 0, NULL },
	"bsExpiry",	{ bsExpiry, 0, 0, 1, 0, NULL }
};

int gFunctionTableSize = sizeof(gFunctionTable)/sizeof(FunctionTag);
