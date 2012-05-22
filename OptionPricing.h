//============================================================================
// Functions exported by this DLL plus data structures required to use them
//
 
#include	"Plugin.h"

extern "C"
long abDateToJDN(long abDate);

struct Greeks {
	double delta;
	double gamma;
	double vega;
	double theta;
	double rho;	
};

struct Condor {
	double longCallStrike;
	double shortCallStrike;
	double shortPutStrike;
	double longPutStrike;

	double netPremium;
};

extern "C" PLUGINAPI double
blackScholesEuro(double price, double strike, double days, bool CorP, double v, double r, Greeks& greeks);

extern "C" PLUGINAPI char*
testExpiry(int year, int month, int day, int expiries);

extern "C" PLUGINAPI bool
testCondor(int year, int month, int day,double spot,double v,double r,double step,int minCallShortDist,int minPutShortDist,
				int minDays,int maxDays,double maxDelta,double minPremium,Condor& ret);


