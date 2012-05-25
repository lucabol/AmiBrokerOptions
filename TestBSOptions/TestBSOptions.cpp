#include "stdafx.h"
#include "..\OptionPricing.h"
#include <iostream>

void banner(const char* testName) {
    std::cout << "Testing " << testName << " ... " << std::endl;
}

void test(bool a) {
	if(a) std::cout << "Ok ..." << std::endl; else std::cout << "ERROR !!!" << std::endl;
}

void test(double a, double b) {
    if (std::abs(a - b) < 0.001)
        std::cout << "Ok ..." << std::endl;
    else
        std::cout << "ERROR !!!" << std::endl;
}

void testG(Greeks g, double value, double v, double d, double ga, double ve, double t, double r) {
    test(value, v);
    test(g.delta, d);
    test(g.gamma, ga);
    test(g.vega, ve);
    test(g.theta, t);
    test(g.rho, r);
}

void testD(char* a, char* b) {
	if(strcmp(a, b) == 0) std::cout << "Ok ..." << std::endl; else std::cout << "ERROR !!!" << std::endl; 
}

int _tmain(int argc, _TCHAR* argv[])
{
    Greeks g;

	auto spot			= 1318.86;
	auto vol			= 0.2233;
	auto days			= 59;
	auto i				= 0.09 / 100;

	banner("black Scholes standard");
    double v			= blackScholesEuro(1314.0, 1330.0, 26, true, 0.1769, 0.0045, g, false);
    testG(g, v, 17.878, 0.411, 0.006, 1.36382, -0.47038, 0.37160);

    v					= blackScholesEuro(1314.0, 1330.0, 26, false, 0.1769, 0.0045, g, false);
    testG(g, v, 33.452, -0.589, 0.006, 1.36382, -0.454005, -0.57548);

	banner("expiry dates");
	auto d				= testExpiry(2012, 5, 22, 0);	testD("2012-May-18", d);
	d					= testExpiry(2012, 12, 31, 1);	testD("2013-Jan-18", d);
	d					= testExpiry(2011, 1, 1, 8);	testD("2011-Sep-16", d);

	banner("volatility smile");
	auto voli			= vol + volatilitySmile(1318.86, 800, 59 / 365.0, 0.09 / 100, false);
	test( voli > vol + 0.10);
	voli				= vol +  volatilitySmile(1318.86, 1000, 59 / 365.0, 0.09 / 100, false);
	test( voli > vol + 0.05);
	voli				= vol +  volatilitySmile(1318.86, 1318.86, 59 / 365.0, 0.09 / 100, true);
	test( voli < vol + 0.05 && voli > vol - 0.05);
	voli				= vol +  volatilitySmile(1318.86, 1400, 59 / 365.0, 0.09 / 100, true);
	test( voli < vol - 0.03);
	voli				= vol +  volatilitySmile(1318.86, 1600, 59 / 365.0, 0.09 / 100, true);
	test( voli < vol - 0.10);

	banner("condors");
	Condor c;
	auto found			= testCondor(2012, 05, 25, 1318.86, 0.2233, 0.0045, 25.0/1314.0, 4, 6, 40, 120, 0.10, 3.0 / 25.0, c);
	test(found);
	test(c.month == 8 && c.day == 17 && c.shortCallStrike == 1450 && c.shortPutStrike == 1075 &&
		 c.longCallStrike == 1475 && c.longPutStrike == 1050 && c.netPremium > 4);

    return 0;
}

