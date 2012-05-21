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

int _tmain(int argc, _TCHAR* argv[])
{
    Greeks g;
    double v = blackScholesEuro(1314.0, 1330.0, 26, true, 0.1769, 0.0045, g);
    testG(g, v, 17.878, 0.411, 0.006, 1.36382, 0.373649, 0.37160);

    v = blackScholesEuro(1314.0, 1330.0, 26, false, 0.1769, 0.0045, g);
    testG(g, v, 33.452, -0.589, 0.006, 1.36382, -0.454005, -0.57548);

    return 0;
}

