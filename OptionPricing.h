//============================================================================
// Functions exported by this DLL plus data structures required to use them
//
 
#include	"Plugin.h"

// This represents a discrete dividend
struct	DiscreteDividend
{
	long				date;
	float				amount;
};

extern "C" long abDateToJDN(long abDate);

extern "C" PLUGINAPI float
TV1day(float				S,			// Share price ...
	   long					T,			// on This date
	   float				K,			// option striKe price ...
	   long					E,			// on this Expiry date
	   bool					CorP,		// true for Call, false for Put
	   float				v,			// annualised Volatility
	   float				r,			// Risk free return rate
	   int					nDividends,	// no. elements in ...
	   DiscreteDividend*	vDividends,	// the discrete dividends for the share
	   bool					ee,			// is Early Exercise possible?
	   int					minPeriods);// minimum number of periods to use

extern "C" PLUGINAPI float
TV1CAN(float				S,			// Share price ...
	   long					T,			// on This date
	   float				K,			// option striKe price ...
	   long					E,			// on this Expiry date
	   float				v,			// annualised Volatility
	   float				r,			// Risk free return rate
	   int					minPeriods);// minimum number of periods to use

extern "C" PLUGINAPI float
TV1CAD(float				S,			// Share price ...
	   long					T,			// on This date
	   float				K,			// option striKe price ...
	   long					E,			// on this Expiry date
	   float				v,			// annualised Volatility
	   float				r,			// Risk free return rate
	   int					nDividends,	// no. elements in ...
	   DiscreteDividend*	vDividends,	// the discrete dividends for the share
	   int					minPeriods);// minimum number of periods to use

extern "C" PLUGINAPI float
TV1PAN(float				S,			// Share price ...
	   long					T,			// on This date
	   float				K,			// option striKe price ...
	   long					E,			// on this Expiry date
	   float				v,			// annualised Volatility
	   float				r,			// Risk free return rate
	   int					minPeriods);// minimum number of periods to use

extern "C" PLUGINAPI float
TV1PAD(float				S,			// Share price ...
	   long					T,			// on This date
	   float				K,			// option striKe price ...
	   long					E,			// on this Expiry date
	   float				v,			// annualised Volatility
	   float				r,			// Risk free return rate
	   int					nDividends,	// no. elements in ...
	   DiscreteDividend*	vDividends,	// the discrete dividends for the share
	   int					minPeriods);// minimum number of periods to use

extern "C" PLUGINAPI float
TV1CEN(float				S,			// Share price ...
	   long					T,			// on This date
	   float				K,			// option striKe price ...
	   long					E,			// on this Expiry date
	   float				v,			// annualised Volatility
	   float				r,			// Risk free return rate
	   int					minPeriods);// minimum number of periods to use

extern "C" PLUGINAPI float
TV1CED(float				S,			// Share price ...
	   long					T,			// on This date
	   float				K,			// option striKe price ...
	   long					E,			// on this Expiry date
	   float				v,			// annualised Volatility
	   float				r,			// Risk free return rate
	   int					nDividends,	// no. elements in ...
	   DiscreteDividend*	vDividends,	// the discrete dividends for the share
	   int					minPeriods);// minimum number of periods to use

extern "C" PLUGINAPI float
TV1PEN(float				S,			// Share price ...
	   long					T,			// on This date
	   float				K,			// option striKe price ...
	   long					E,			// on this Expiry date
	   float				v,			// annualised Volatility
	   float				r,			// Risk free return rate
	   int					minPeriods);// minimum number of periods to use

extern "C" PLUGINAPI float
TV1PED(float				S,			// Share price ...
	   long					T,			// on This date
	   float				K,			// option striKe price ...
	   long					E,			// on this Expiry date
	   float				v,			// annualised Volatility
	   float				r,			// Risk free return rate
	   int					nDividends,	// no. elements in ...
	   DiscreteDividend*	vDividends,	// the discrete dividends for the share
	   int					minPeriods);// minimum number of periods to use

extern "C" PLUGINAPI void
TVVCAN(int					n,			// No. entries in:
	   float*				vS,			// Share price vector ...
	   long*				vT,			// on This date vector
	   float*				vR,			// Results vector
	   float				K,			// option striKe price ...
	   long					E,			// on this Expiry date
	   float				v,			// annualised Volatility
	   float				r,			// Risk free return rate
	   int					minPeriods);// minimum number of periods to use

extern "C" PLUGINAPI void
TVVCAD(int					n,			// No. entries in:
	   float*				vS,			// Share price vector ...
	   long*				vT,			// on This date vector
	   float*				vR,			// Results vector
	   float				K,			// option striKe price ...
	   long					E,			// on this Expiry date
	   float				v,			// annualised Volatility
	   float				r,			// Risk free return rate
	   int					nDividends,	// no. elements in ...
	   DiscreteDividend*	vDividends,	// the discrete dividends for the share
	   int					minPeriods);// minimum number of periods to use

extern "C" PLUGINAPI void
TVVPAN(int					n,			// No. entries in:
	   float*				vS,			// Share price vector ...
	   long*				vT,			// on This date vector
	   float*				vR,			// Results vector
	   float				K,			// option striKe price ...
	   long					E,			// on this Expiry date
	   float				v,			// annualised Volatility
	   float				r,			// Risk free return rate
	   int					minPeriods);// minimum number of periods to use

extern "C" PLUGINAPI void
TVVPAD(int					n,			// No. entries in:
	   float*				vS,			// Share price vector ...
	   long*				vT,			// on This date vector
	   float*				vR,			// Results vector
	   float				K,			// option striKe price ...
	   long					E,			// on this Expiry date
	   float				v,			// annualised Volatility
	   float				r,			// Risk free return rate
	   int					nDividends,	// no. elements in ...
	   DiscreteDividend*	vDividends,	// the discrete dividends for the share
	   int					minPeriods);// minimum number of periods to use

extern "C" PLUGINAPI void
TVVCEN(int					n,			// No. entries in:
	   float*				vS,			// Share price vector ...
	   long*				vT,			// on This date vector
	   float*				vR,			// Results vector
	   float				K,			// option striKe price ...
	   long					E,			// on this Expiry date
	   float				v,			// annualised Volatility
	   float				r,			// Risk free return rate
	   int					minPeriods);// minimum number of periods to use

extern "C" PLUGINAPI void
TVVCED(int					n,			// No. entries in:
	   float*				vS,			// Share price vector ...
	   long*				vT,			// on This date vector
	   float*				vR,			// Results vector
	   float				K,			// option striKe price ...
	   long					E,			// on this Expiry date
	   float				v,			// annualised Volatility
	   float				r,			// Risk free return rate
	   int					nDividends,	// no. elements in ...
	   DiscreteDividend*	vDividends,	// the discrete dividends for the share
	   int					minPeriods);// minimum number of periods to use

extern "C" PLUGINAPI void
TVVPEN(int					n,			// No. entries in:
	   float*				vS,			// Share price vector ...
	   long*				vT,			// on This date vector
	   float*				vR,			// Results vector
	   float				K,			// option striKe price ...
	   long					E,			// on this Expiry date
	   float				v,			// annualised Volatility
	   float				r,			// Risk free return rate
	   int					minPeriods);// minimum number of periods to use

extern "C" PLUGINAPI void
TVVPED(int					n,			// No. entries in:
	   float*				vS,			// Share price vector ...
	   long*				vT,			// on This date vector
	   float*				vR,			// Results vector
	   float				K,			// option striKe price ...
	   long					E,			// on this Expiry date
	   float				v,			// annualised Volatility
	   float				r,			// Risk free return rate
	   int					nDividends,	// no. elements in ...
	   DiscreteDividend*	vDividends,	// the discrete dividends for the share
	   int					minPeriods);// minimum number of periods to use

extern "C" PLUGINAPI float
IV1day(float				O,			// Option Price ...
	   float				S,			// Share price ...
	   long					T,			// on This date
	   float				K,			// option striKe price ...
	   long					E,			// on this Expiry date
	   bool					CorP,		// true for Call, false for Put
	   float				v,			// annualised Volatility
	   float				r,			// Risk free return rate
	   int					nDividends,	// no. elements in ...
	   DiscreteDividend*	vDividends,	// the discrete dividends for the share
	   bool					ee,			// is Early Exercise possible?
	   int					minPeriods);// minimum number of periods to use

extern AmiVar TVV(int NumArgs, AmiVar* ArgsTable);

extern AmiVar IVV(int NumArgs, AmiVar *ArgsTable);

struct Greeks {
	double delta;
	double gamma;
	double vega;
	double theta;
	double rho;	
};

extern "C" PLUGINAPI double
blackScholesEuro(double price, double strike, double days, bool CorP, double v, double r, Greeks& greeks);