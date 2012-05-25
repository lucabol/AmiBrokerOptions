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


AmiVar execFormula(int NumArgs, AmiVar* ArgsTable, std::function<double (double, double, double, bool, double, double, Greeks&, bool)> f)
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
	bool				useSmile;		// true -> use it, false -> standard bs

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
	useSmile= (ArgsTable[7].val != 0);				// bool 7

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
		Rv[i]					= (float) f (currentPrice, K, daysToExpiration, CorP, v, r, greeks, useSmile);

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

// Eq 12 from paper below, gives inconsistent results
extern "C"
double volatilitySmile1(double spot, double strike, double time, double r) {
	static const auto	c1	= -0.2547240; 
	static const auto	c2	=  0.1319096; 
	static const auto	c3	=  0.5595131;

	auto F					= spot * std::exp(r * time);
	auto K					= std::log(strike / F) / std::sqrt(time);
	auto ret				= c1 * K + c2 * std::pow(K, 2) + c3 * std::pow(std::abs(K), 3);

	return ret;
}


// Eq 11 from paper below
extern "C"
double volatilitySmile(double spot, double strike, double time, double r, bool CorP) {
	static const auto	b0	=  0.0058480;
	static const auto	b1	= -0.2884075; 
	static const auto	b2	=  0.0322727;
	static const auto	b3	= -0.0075740;
	static const auto	b4	=  0.0015705;
	static const auto	b5	=  0.0414902;

	auto F					= spot * std::exp(r * time);
	auto K					= std::log(strike / F);
	auto ret				= b0 + b1 * K + b2 * std::pow(K, 2) + b3 * time + b4 * std::pow(time, 2) + b5 * K * time;

	auto crazy				= CorP ? 3.5 : 2.0; // derived empirically from volatility smile one day
	ret						= crazy * ret;

	return ret;
}

// Eq 14 from paper below gives inconsistent results
extern "C"
double volatilitySmile2(double spot, double strike, double time, double r) {
	static const auto	d0	=  0.4419942;
	static const auto	d1	= -0.2394926; 
	static const auto	d2	=  0.0912235;
	static const auto	d3	=  0.2581095;

	auto F					= spot * std::exp(r * time);
	auto K					= std::log(strike / F) / std::pow(time, d0);
	auto ret				= d1 * K + d2 * std::pow(K, 2) + d3 * std::pow(K, 3);
	return ret;
}

// Black Scholes formula, to calculate volatility surface read doc below formula 11. 
// http://finance.business.queensu.ca/psfile/DaglishHullSuoRevised.pdf, F = forward value of S for a contract maturying at T -> F = S*exp(r * (T -t)) , K = strike price, S = asset price
extern "C"
auto blackScholesEuro(double price, double strike, double days, bool CorP, double v, double r, Greeks& greeks, bool useSmile)
					  -> double {

      auto s			= price;
	  auto x			= strike;
	  auto t			= days / 365.0;
	  if(useSmile) {
		  v				= v + volatilitySmile(s, x, t, r, CorP);
	  }

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


#ifdef _DEBUG
void print(boost::gregorian::date d) {
	std::stringstream ss;
	ss << d;
	std::cout << ss.str().c_str() << std::endl;
}
#else
void print(boost::gregorian::date d) {}
#endif

// Returns the expiry date some months out. 0 = this month, 1 = next month, 2 = two months out, etc...
boost::gregorian::date calcExpiry(boost::gregorian::date d, int expiries) {
	using namespace boost::gregorian;

	typedef nth_day_of_the_week_in_month nth_dow;
	months single(expiries);
	date dTarget = d + single;
	nth_dow ndm(nth_dow::third, Friday, dTarget.month());
	auto value = ndm.get_date(dTarget.year());
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
		return blackScholesEuro(spot, strike, days, CorP, v, r, g, true);
	};

	// find a condor that works at this expiry
	auto findCondor		= [callShort, putShort, spot, stepPr, v, r, minCpr, minPpr, scholes, maxDelta, premiumPr, &g, &ret] (int days) -> bool {
		ret.shortCallStrike				= callShort;
		ret.shortPutStrike				= putShort;
		auto shCallPremium				= 0.0;
		auto shPutPremium				= 0.0;

		// find short call strike price < maxDelta
		while(true) {
			shCallPremium				= scholes(ret.shortCallStrike, days, true);
			if(g.delta <= maxDelta)
				break;
			else
				ret.shortCallStrike		+= stepPr;
		}

		// find short put strike price < maxDelta
		while(true) {
			shPutPremium				= scholes(ret.shortPutStrike, days, false);
			if( (- g.delta) <= maxDelta)
				break;
			else
				ret.shortPutStrike		-= stepPr;
		}

		// check premium is adeguate
		ret.longCallStrike				= ret.shortCallStrike + stepPr;
		ret.longPutStrike				= ret.shortPutStrike  - stepPr;
		auto lgCall						= scholes(ret.longCallStrike, days, true);
		auto lgPut						= scholes(ret.longPutStrike,  days, false);
		ret.netPremium					= shCallPremium + shPutPremium - lgCall - lgPut;

		return ret.netPremium > premiumPr;
	};

	// increases the expiry until it finds a condor or the expiry is too far out
	while (expiry < maxDate) {
		auto days		= (expiry - now).days();
		if(findCondor(days)) {
			ret.year	= expiry.year();
			ret.month	= expiry.month();
			ret.day		= expiry.day();
			return true;
		}
		expiry			= calcExpiry(expiry, +1);
	}

	return false;
}

AmiVar bsCondor(int NumArgs, AmiVar* ArgsTable)
{
	using namespace std;

	// results
	AmiVar results						= gSite.AllocArrayResult();			
	AmiVar longCallStrike				= gSite.AllocArrayResult();			
	AmiVar shortCallStrike				= gSite.AllocArrayResult();				
	AmiVar shortPutStrike				= gSite.AllocArrayResult();				
	AmiVar longPutStrike				= gSite.AllocArrayResult();				
	AmiVar netPremium					= gSite.AllocArrayResult();				
	AmiVar expiries						= gSite.AllocArrayResult();

	int n								= gSite.GetArraySize();

	auto dates							= vector<boost::gregorian::date>(n);
	float* nows							= ArgsTable[0].array;
	for(int i = 0; i < n; ++i)			dates.push_back(abDateToDate((long)nows[i]));

	float* spots						= ArgsTable[1].array;
	float* vols							= ArgsTable[2].array;
	float* rates						= ArgsTable[3].array;
	float  steps						= ArgsTable[4].val;
	float  minCall						= ArgsTable[5].val;
	float  minPut						= ArgsTable[6].val;
	float  minDays						= ArgsTable[7].val;
	float  maxDays						= ArgsTable[8].val;
	float  maxDelta						= ArgsTable[9].val;
	float  minPremium					= ArgsTable[10].val;
	
	int	i;
	int	j;

	// add empty values at the start of results if there is no spot
	for (j = 0; j < n && IS_EMPTY(spots[j]); j++)
	{
		longCallStrike.array[j]			= EMPTY_VAL;
		shortCallStrike.array[j]		= EMPTY_VAL;
		longPutStrike.array[j]			= EMPTY_VAL;
		shortPutStrike.array[j]			= EMPTY_VAL;
		netPremium.array[j]				= EMPTY_VAL;
		expiries.array[j]				= EMPTY_VAL;
	}

	Condor c;
    for (i = j; i < n; ++i)
    {
		results.array[i]				= condor(dates[i], spots[i], vols[i], rates[i], steps, round(minCall), round(minPut), round(minDays),
										  round(maxDays), maxDelta, minPremium, c);
		longCallStrike.array[i]			= (float) c.longCallStrike;
		shortCallStrike.array[i]		= (float) c.shortCallStrike;
		longPutStrike.array[i]			= (float) c.longPutStrike;
		shortPutStrike.array[i]			= (float) c.shortPutStrike;
		netPremium.array[i]				= (float) c.netPremium;
		expiries.array[i]				= (float) dateToAbDate(boost::gregorian::date(c.year, c.month, c.day));
	}

	gSite.SetVariable("bsLongCallStrike",  longCallStrike);
	gSite.SetVariable("bsShortCallStrike", shortCallStrike);
	gSite.SetVariable("bsLongPutStrike",   longPutStrike);
	gSite.SetVariable("bsShortPutstrike",  shortPutStrike);
	gSite.SetVariable("bsNetPremium",	   netPremium);
	gSite.SetVariable("bsExpiry",		   expiries);

    return results;
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
	"bsVEuro",	{ bsValueEuropean, 2, 0, 6, 0, NULL },
	"bsExpiry",	{ bsExpiry,		   0, 0, 1, 0, NULL },
	"bsCondor",	{ bsCondor,		   4, 0, 7, 0, NULL }
};

int gFunctionTableSize = sizeof(gFunctionTable)/sizeof(FunctionTag);
