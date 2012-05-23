// TestBSOptions.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "..\OptionPricing.h"
#include <iostream>

void test(double a, double b) {
    if (std::abs(a - b) < 0.001)
        std::cout << "Ok ..." << std::endl;
    else
        std::cout << "ERROR !!!" << std::endl;
}

void testG(Greeks g, double value, double v, double d, double ga, double ve, double t, double r) {
    static int x = 0;
    std::cout << "Test: " << x++ << std::endl;
    test(value, v);
    test(g.delta, d);
    test(g.gamma, ga);
    test(g.vega, ve);
    test(g.theta, t);
    test(g.rho, r);
    std::cout << std::endl;
}

void testD(char* a, char* b) {
	if(strcmp(a, b) == 0) std::cout << "Ok ..." << std::endl; else std::cout << "ERROR !!!" << std::endl; 
}

int _tmain(int argc, _TCHAR* argv[])
{
    Greeks g;
    double v			= blackScholesEuro(1314.0, 1330.0, 26, true, 0.1769, 0.0045, g);
    testG(g, v, 17.878, 0.411, 0.006, 1.36382, -0.47038, 0.37160);

    v					= blackScholesEuro(1314.0, 1330.0, 26, false, 0.1769, 0.0045, g);
    testG(g, v, 33.452, -0.589, 0.006, 1.36382, -0.454005, -0.57548);

	auto d				= testExpiry(2012, 5, 22, 0);	testD("2012-May-18", d);
	d					= testExpiry(2012, 12, 31, 1);	testD("2013-Jan-18", d);
	d					= testExpiry(2011, 1, 1, 8);	testD("2011-Sep-16", d);

	auto vol = 0.2248 + volatilitySmile(1316, 825, 60.0 / 365.0, 0.0045);
	vol = 0.2248 +  volatilitySmile(1316, 1000, 60.0 / 365.0, 0.0045);
	vol = 0.2248 +  volatilitySmile(1316, 1400, 60.0 / 365.0, 0.0045);
	vol = 0.2248 +  volatilitySmile(1316, 1500, 60.0 / 365.0, 0.0045);

	Condor c;
	auto found			= testCondor(2012, 05, 25, 1314.0, 0.2248, 0.0045, 25.0/1314.0, 4, 6, 40, 80, 0.10, 3.0 / 25.0, c);

    return 0;
}

