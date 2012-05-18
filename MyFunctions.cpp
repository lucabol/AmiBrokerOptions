////////////////////////////////////////////////////
// MyFunctions.cpp
//
// Copyright (C)2002 Matt Atterbury
// All rights reserved.
//
////////////////////////////////////////////////////

#include <math.h>
#include "OptionPricing.h"

//=============================================================================
// Private Constants
//
static const int		MAX_PERIODS = 1024;

//=============================================================================
// Helper functions
//

//  Convert an AmiBroker date to a JDN (No. days from sometime in antiquity)
//  AB dates or numbers of the form YYYMMDD, where:
//      YYY is the year  (the 10K, 100K, and 1M digits)
//      MM  is the month (the 100 and 1K digits)
//      DD  is the day   (the 1 and 10 digits)
//  JDNs are the absolute number of days from a fixed date in the past.
extern "C"
long abDateToJDN(long abDate)
{
    int					dd;
    int					mm;
    int					yy;

	// Extract the day, month, and year from the given date
    dd = (abDate % 100);
    mm = (abDate / 100) % 100;
    yy = (abDate / 10000) + 1900;

    // Now convert to a JDN
    long				a;
    long				y;
    long				m;

    a = (14 - mm) / 12;
    y = yy + 4800 - a;
    m = mm + 12 * a - 3;

    return (dd +
		    (long) ((153 * m + 2) / 5) +
			y * 365 +
			(long) (y / 4) -
			(long) (y / 100) +
			(long) (y / 400) -
			32045);
}

//=============================================================================
// Option Class and Type Specific functions for 1 day
//

//
//	Calculate the Theoretical Value for the given price. This uses the
//	binomial method, catering for up to two dividends. The number of
//	periods used is the number of days until expiry, multiplied by an
//	integer until is greater than or equal to the given mimimum, up to
//	a maximum of 1024 periods.
//
//	The formula is taken from "Option Volatility & Pricing" by Sheldon
//	Natenberg, with adjustment made for discrete dividends as found on
//	"the fount of all knowledge" (aka the World Wide Web).
//
extern "C" float
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
	   int					minPeriods)	// minimum number of periods to use
{
    // bits-n-pieces
	int					i;
	int					j;
	int					k;
    long				n;
    float				u;
    float				d;
    float				RR;
    float				p;
    long				X;
	float				U[MAX_PERIODS+1];
	float				V[MAX_PERIODS+1];

	// Work out the number of periods to use.
	n = E - T;
	if (minPeriods > 0)
	{
		n = ((long) (minPeriods + n - 1) / n) * n;
	}
	if (n > MAX_PERIODS) n = MAX_PERIODS;

    // Process one day
    u = (float) exp(v*sqrt((double)(E-T)/365/n));	// upward movement mult'r
    d = 1 / u;										// downward movement mult'r
    RR = (1 + (r * (E - T) / 365 / n));				// annual risk-free rate
    p = (RR - d) / (u - d);							// P(upward movement)

    // Initialise the array of final day values
    for (i = 0; i < n+1; ++i)
	{
        U[i] = (float) (S * pow(u, i) * pow(d, n - i));
		for (k = 0; k < nDividends; ++k)
		{
			if (T < vDividends[k].date && E > vDividends[k].date)
			{
				U[i] -= vDividends[k].amount;
			}
		}
        V[i] = 0.0;
        if (CorP)
		{
            if (U[i] > K) V[i] = U[i] - K;
		}
		else
		{
            if (K > U[i]) V[i] = K - U[i];
		}
	}

    // Work back through the days to get the specific days results
    for (i = n; i > 0; --i)
	{
        X = T + (E - T) * i / n;
        for (j = 0; j < i; ++j)
		{
            U[j] = S;
			for (k = 0; k < nDividends; ++k)
			{
				if (T < vDividends[k].date && X > vDividends[k].date)
				{
					U[i] -= vDividends[k].amount;
				}
			}
            U[j] = (float) (U[j] * pow(u, j) * pow(d, i - j - 1));
            V[j] = (p * V[j + 1] + (1 - p) * V[j]) / RR;
			if      (!ee)
			{
				;
			}
            else if (CorP)
			{
                if (V[j] < (U[j] - K)) V[j] = U[j] - K;
			}
            else
			{
                if (V[j] < (K - U[j])) V[j] = K - U[j];
			}
		}
	}

    return V[0];
}

static void
TV1prepare(long		T,			// on This date
		   long		E,			// on this Expiry date
		   float	v,			// annualised Volatility
		   float	r,			// Risk free return rate
		   int		minPeriods,	// minimum number of periods to use
		   long*	p_n,
		   float*	p_u,
		   float*	p_d,
		   float*	p_RR,
		   float*	p_p)
{
	// Work out the number of periods to use.
	*p_n = E - T;
	if (minPeriods > 0)
	{
		*p_n = ((long) (minPeriods + *p_n - 1) / *p_n) * *p_n;
	}
	if (*p_n > MAX_PERIODS) *p_n = MAX_PERIODS;

    // Calculate the loop invariants
    *p_u = (float) exp(v*sqrt((double)(E-T)/365/ *p_n));// upward mov. mult'r
    *p_d = 1 / *p_u;									// downward mov. mult'r
    *p_RR = (1 + (r * (E - T) / 365 / *p_n));			// risk-free rate
    *p_p = (*p_RR - *p_d) / (*p_u - *p_d);				// P(upward movement)
}

// TV1CAN
// Theoretical Value, 1 day, Call, American, No dividends
extern "C" PLUGINAPI float
TV1CAN(float				S,			// Share price ...
	   long					T,			// on This date
	   float				K,			// option striKe price ...
	   long					E,			// on this Expiry date
	   float				v,			// annualised Volatility
	   float				r,			// Risk free return rate
	   int					minPeriods)	// minimum number of periods to use
{
    // bits-n-pieces
	int					i;
	int					j;
    long				n;
    float				u;
    float				d;
    float				RR;
    float				p;
	float				U[MAX_PERIODS+1];
	float				V[MAX_PERIODS+1];

	TV1prepare(T, E, v, r, minPeriods, &n, &u, &d, &RR, &p);

    // Initialise the array of final day values
    for (i = 0; i < n+1; ++i)
	{
        U[i] = (float) (S * pow(u, i) * pow(d, n - i));
        V[i] = 0.0;
        if (U[i] > K) V[i] = U[i] - K;
	}

    // Work back through the days to get the specific days results
    for (i = n; i > 0; --i)
	{
        for (j = 0; j < i; ++j)
		{
            U[j] = (float) (S * pow(u, j) * pow(d, i - j - 1));
            V[j] = (p * V[j + 1] + (1 - p) * V[j]) / RR;
            if (V[j] < (U[j] - K)) V[j] = U[j] - K;
		}
	}

    return V[0];
}

// TV1CAD
// Theoretical Value, 1 day, Call, American, Discrete dividends
extern "C" PLUGINAPI float
TV1CAD(float				S,			// Share price ...
	   long					T,			// on This date
	   float				K,			// option striKe price ...
	   long					E,			// on this Expiry date
	   float				v,			// annualised Volatility
	   float				r,			// Risk free return rate
	   int					nDividends,	// no. elements in ...
	   DiscreteDividend*	vDividends,	// the discrete dividends for the share
	   int					minPeriods)	// minimum number of periods to use
{
    // bits-n-pieces
	int					i;
	int					j;
	int					k;
    long				n;
    float				u;
    float				d;
    float				RR;
    float				p;
    long				X;
	float				U[MAX_PERIODS+1];
	float				V[MAX_PERIODS+1];

	TV1prepare(T, E, v, r, minPeriods, &n, &u, &d, &RR, &p);

    // Initialise the array of final day values
    for (i = 0; i < n+1; ++i)
	{
        U[i] = (float) (S * pow(u, i) * pow(d, n - i));
		for (k = 0; k < nDividends; ++k)
		{
			if (T < vDividends[k].date && E > vDividends[k].date)
			{
				U[i] -= vDividends[k].amount;
			}
		}
        V[i] = 0.0;
        if (U[i] > K) V[i] = U[i] - K;
	}

    // Work back through the days to get the specific days results
    for (i = n; i > 0; --i)
	{
        X = T + (E - T) * i / n;
        for (j = 0; j < i; ++j)
		{
            U[j] = S;
			for (k = 0; k < nDividends; ++k)
			{
				if (T < vDividends[k].date && X > vDividends[k].date)
				{
					U[i] -= vDividends[k].amount;
				}
			}
            U[j] = (float) (U[j] * pow(u, j) * pow(d, i - j - 1));
            V[j] = (p * V[j + 1] + (1 - p) * V[j]) / RR;
            if (V[j] < (U[j] - K)) V[j] = U[j] - K;
		}
	}

    return V[0];
}

// TV1PAN
// Theoretical Value, 1 day, Put, American, No dividends
extern "C" PLUGINAPI float
TV1PAN(float				S,			// Share price ...
	   long					T,			// on This date
	   float				K,			// option striKe price ...
	   long					E,			// on this Expiry date
	   float				v,			// annualised Volatility
	   float				r,			// Risk free return rate
	   int					minPeriods)	// minimum number of periods to use
{
    // bits-n-pieces
	int					i;
	int					j;
    long				n;
    float				u;
    float				d;
    float				RR;
    float				p;
	float				U[MAX_PERIODS+1];
	float				V[MAX_PERIODS+1];

	TV1prepare(T, E, v, r, minPeriods, &n, &u, &d, &RR, &p);

    // Initialise the array of final day values
    for (i = 0; i < n+1; ++i)
	{
        U[i] = (float) (S * pow(u, i) * pow(d, n - i));
        V[i] = 0.0;
        if (K > U[i]) V[i] = K - U[i];
	}

    // Work back through the days to get the specific days results
    for (i = n; i > 0; --i)
	{
        for (j = 0; j < i; ++j)
		{
            U[j] = (float) (S * pow(u, j) * pow(d, i - j - 1));
            V[j] = (p * V[j + 1] + (1 - p) * V[j]) / RR;
            if (V[j] < (K - U[j])) V[j] = K - U[j];
		}
	}

    return V[0];
}

// TV1PAD
// Theoretical Value, 1 day, Put, American, Discrete dividends
extern "C" PLUGINAPI float
TV1PAD(float				S,			// Share price ...
	   long					T,			// on This date
	   float				K,			// option striKe price ...
	   long					E,			// on this Expiry date
	   float				v,			// annualised Volatility
	   float				r,			// Risk free return rate
	   int					nDividends,	// no. elements in ...
	   DiscreteDividend*	vDividends,	// the discrete dividends for the share
	   int					minPeriods)	// minimum number of periods to use
{
    // bits-n-pieces
	int					i;
	int					j;
	int					k;
    long				n;
    float				u;
    float				d;
    float				RR;
    float				p;
    long				X;
	float				U[MAX_PERIODS+1];
	float				V[MAX_PERIODS+1];

	TV1prepare(T, E, v, r, minPeriods, &n, &u, &d, &RR, &p);

    // Initialise the array of final day values
    for (i = 0; i < n+1; ++i)
	{
        U[i] = (float) (S * pow(u, i) * pow(d, n - i));
		for (k = 0; k < nDividends; ++k)
		{
			if (T < vDividends[k].date && E > vDividends[k].date)
			{
				U[i] -= vDividends[k].amount;
			}
		}
        V[i] = 0.0;
        if (K > U[i]) V[i] = K - U[i];
	}

    // Work back through the days to get the specific days results
    for (i = n; i > 0; --i)
	{
        X = T + (E - T) * i / n;
        for (j = 0; j < i; ++j)
		{
            U[j] = S;
			for (k = 0; k < nDividends; ++k)
			{
				if (T < vDividends[k].date && X > vDividends[k].date)
				{
					U[i] -= vDividends[k].amount;
				}
			}
            U[j] = (float) (U[j] * pow(u, j) * pow(d, i - j - 1));
            V[j] = (p * V[j + 1] + (1 - p) * V[j]) / RR;
            if (V[j] < (K - U[j])) V[j] = K - U[j];
		}
	}

    return V[0];
}

// TV1CEN
// Theoretical Value, 1 day, Call, European, No dividends
extern "C" PLUGINAPI float
TV1CEN(float				S,			// Share price ...
	   long					T,			// on This date
	   float				K,			// option striKe price ...
	   long					E,			// on this Expiry date
	   float				v,			// annualised Volatility
	   float				r,			// Risk free return rate
	   int					minPeriods)	// minimum number of periods to use
{
    // bits-n-pieces
	int					i;
	int					j;
    long				n;
    float				u;
    float				d;
    float				RR;
    float				p;
	float				U[MAX_PERIODS+1];
	float				V[MAX_PERIODS+1];

	TV1prepare(T, E, v, r, minPeriods, &n, &u, &d, &RR, &p);

    // Initialise the array of final day values
    for (i = 0; i < n+1; ++i)
	{
        U[i] = (float) (S * pow(u, i) * pow(d, n - i));
        V[i] = 0.0;
        if (U[i] > K) V[i] = U[i] - K;
	}

    // Work back through the days to get the specific days results
    for (i = n; i > 0; --i)
	{
        for (j = 0; j < i; ++j)
		{
            U[j] = (float) (S * pow(u, j) * pow(d, i - j - 1));
            V[j] = (p * V[j + 1] + (1 - p) * V[j]) / RR;
		}
	}

    return V[0];
}

// TV1CED
// Theoretical Value, 1 day, Call, European, Discrete dividends
extern "C" PLUGINAPI float
TV1CED(float				S,			// Share price ...
	   long					T,			// on This date
	   float				K,			// option striKe price ...
	   long					E,			// on this Expiry date
	   float				v,			// annualised Volatility
	   float				r,			// Risk free return rate
	   int					nDividends,	// no. elements in ...
	   DiscreteDividend*	vDividends,	// the discrete dividends for the share
	   int					minPeriods)	// minimum number of periods to use
{
    // bits-n-pieces
	int					i;
	int					j;
	int					k;
    long				n;
    float				u;
    float				d;
    float				RR;
    float				p;
    long				X;
	float				U[MAX_PERIODS+1];
	float				V[MAX_PERIODS+1];

	TV1prepare(T, E, v, r, minPeriods, &n, &u, &d, &RR, &p);

    // Initialise the array of final day values
    for (i = 0; i < n+1; ++i)
	{
        U[i] = (float) (S * pow(u, i) * pow(d, n - i));
		for (k = 0; k < nDividends; ++k)
		{
			if (T < vDividends[k].date && E > vDividends[k].date)
			{
				U[i] -= vDividends[k].amount;
			}
		}
        V[i] = 0.0;
        if (U[i] > K) V[i] = U[i] - K;
	}

    // Work back through the days to get the specific days results
    for (i = n; i > 0; --i)
	{
        X = T + (E - T) * i / n;
        for (j = 0; j < i; ++j)
		{
            U[j] = S;
			for (k = 0; k < nDividends; ++k)
			{
				if (T < vDividends[k].date && X > vDividends[k].date)
				{
					U[i] -= vDividends[k].amount;
				}
			}
            U[j] = (float) (U[j] * pow(u, j) * pow(d, i - j - 1));
            V[j] = (p * V[j + 1] + (1 - p) * V[j]) / RR;
		}
	}

    return V[0];
}

// TV1PEN
// Theoretical Value, 1 day, Put, European, No dividends
extern "C" PLUGINAPI float
TV1PEN(float				S,			// Share price ...
	   long					T,			// on This date
	   float				K,			// option striKe price ...
	   long					E,			// on this Expiry date
	   float				v,			// annualised Volatility
	   float				r,			// Risk free return rate
	   int					minPeriods)	// minimum number of periods to use
{
    // bits-n-pieces
	int					i;
	int					j;
    long				n;
    float				u;
    float				d;
    float				RR;
    float				p;
	float				U[MAX_PERIODS+1];
	float				V[MAX_PERIODS+1];

	TV1prepare(T, E, v, r, minPeriods, &n, &u, &d, &RR, &p);

    // Initialise the array of final day values
    for (i = 0; i < n+1; ++i)
	{
        U[i] = (float) (S * pow(u, i) * pow(d, n - i));
        V[i] = 0.0;
        if (K > U[i]) V[i] = K - U[i];
	}

    // Work back through the days to get the specific days results
    for (i = n; i > 0; --i)
	{
        for (j = 0; j < i; ++j)
		{
            U[j] = (float) (S * pow(u, j) * pow(d, i - j - 1));
            V[j] = (p * V[j + 1] + (1 - p) * V[j]) / RR;
		}
	}

    return V[0];
}

// TV1PED
// Theoretical Value, 1 day, Put, European, Discrete dividends
extern "C" PLUGINAPI float
TV1PED(float				S,			// Share price ...
	   long					T,			// on This date
	   float				K,			// option striKe price ...
	   long					E,			// on this Expiry date
	   float				v,			// annualised Volatility
	   float				r,			// Risk free return rate
	   int					nDividends,	// no. elements in ...
	   DiscreteDividend*	vDividends,	// the discrete dividends for the share
	   int					minPeriods)	// minimum number of periods to use
{
    // bits-n-pieces
	int					i;
	int					j;
	int					k;
    long				n;
    float				u;
    float				d;
    float				RR;
    float				p;
    long				X;
	float				U[MAX_PERIODS+1];
	float				V[MAX_PERIODS+1];

	TV1prepare(T, E, v, r, minPeriods, &n, &u, &d, &RR, &p);

    // Initialise the array of final day values
    for (i = 0; i < n+1; ++i)
	{
        U[i] = (float) (S * pow(u, i) * pow(d, n - i));
		for (k = 0; k < nDividends; ++k)
		{
			if (T < vDividends[k].date && E > vDividends[k].date)
			{
				U[i] -= vDividends[k].amount;
			}
		}
        V[i] = 0.0;
        if (K > U[i]) V[i] = K - U[i];
	}

    // Work back through the days to get the specific days results
    for (i = n; i > 0; --i)
	{
        X = T + (E - T) * i / n;
        for (j = 0; j < i; ++j)
		{
            U[j] = S;
			for (k = 0; k < nDividends; ++k)
			{
				if (T < vDividends[k].date && X > vDividends[k].date)
				{
					U[i] -= vDividends[k].amount;
				}
			}
            U[j] = (float) (U[j] * pow(u, j) * pow(d, i - j - 1));
            V[j] = (p * V[j + 1] + (1 - p) * V[j]) / RR;
		}
	}

    return V[0];
}

//=============================================================================
// Option Class and Type Specific functions for a range of days
//

// TVVCAN
// Theoretical Value, Vector of days, Call, American, No dividends
extern "C" PLUGINAPI void
TVVCAN(int					n,			// No. entries in:
	   float*				vS,			// Share price vector ...
	   long*				vT,			// on This date vector
	   float*				vR,			// Results vector
	   float				K,			// option striKe price ...
	   long					E,			// on this Expiry date
	   float				v,			// annualised Volatility
	   float				r,			// Risk free return rate
	   int					minPeriods)	// minimum number of periods to use
{
	int			i;

	for (i = 0; i < n; ++i)
	{
		vR[i] = TV1CAN(vS[i], abDateToJDN(vT[i]), K, abDateToJDN(E),
					   v, r, minPeriods);
	}
}

// TVVCAD
// Theoretical Value, Vector of days, Call, American, Discrete dividends
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
	   int					minPeriods)	// minimum number of periods to use
{
	int			i;

	for (i = 0; i < n; ++i)
	{
		vR[i] = TV1CAD(vS[i], abDateToJDN(vT[i]), K, abDateToJDN(E),
					   v, r, nDividends, vDividends, minPeriods);
	}
}

// TVVPAN
// Theoretical Value, Vector of days, Put, American, No dividends
extern "C" PLUGINAPI void
TVVPAN(int					n,			// No. entries in:
	   float*				vS,			// Share price vector ...
	   long*				vT,			// on This date vector
	   float*				vR,			// Results vector
	   float				K,			// option striKe price ...
	   long					E,			// on this Expiry date
	   float				v,			// annualised Volatility
	   float				r,			// Risk free return rate
	   int					minPeriods)	// minimum number of periods to use
{
	int			i;

	for (i = 0; i < n; ++i)
	{
		vR[i] = TV1PAN(vS[i], abDateToJDN(vT[i]), K, abDateToJDN(E),
					   v, r, minPeriods);
	}
}

// TVVPAD
// Theoretical Value, Vector of days, Put, American, Discrete dividends
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
	   int					minPeriods)	// minimum number of periods to use
{
	int			i;

	for (i = 0; i < n; ++i)
	{
		vR[i] = TV1PAD(vS[i], abDateToJDN(vT[i]), K, abDateToJDN(E),
					   v, r, nDividends, vDividends, minPeriods);
	}
}

// TVVCEN
// Theoretical Value, Vector of days, Call, European, No dividends
extern "C" PLUGINAPI void
TVVCEN(int					n,			// No. entries in:
	   float*				vS,			// Share price vector ...
	   long*				vT,			// on This date vector
	   float*				vR,			// Results vector
	   float				K,			// option striKe price ...
	   long					E,			// on this Expiry date
	   float				v,			// annualised Volatility
	   float				r,			// Risk free return rate
	   int					minPeriods)	// minimum number of periods to use
{
	int			i;

	for (i = 0; i < n; ++i)
	{
		vR[i] = TV1CAN(vS[i], abDateToJDN(vT[i]), K, abDateToJDN(E),
					   v, r, minPeriods);
	}
}

// TVVCED
// Theoretical Value, Vector of days, Call, European, Discrete dividends
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
	   int					minPeriods)	// minimum number of periods to use
{
	int			i;

	for (i = 0; i < n; ++i)
	{
		vR[i] = TV1CAD(vS[i], abDateToJDN(vT[i]), K, abDateToJDN(E),
					   v, r, nDividends, vDividends, minPeriods);
	}
}

// TVVPEN
// Theoretical Value, Vector of days, Put, European, No dividends
extern "C" PLUGINAPI void
TVVPEN(int					n,			// No. entries in:
	   float*				vS,			// Share price vector ...
	   long*				vT,			// on This date vector
	   float*				vR,			// Results vector
	   float				K,			// option striKe price ...
	   long					E,			// on this Expiry date
	   float				v,			// annualised Volatility
	   float				r,			// Risk free return rate
	   int					minPeriods)	// minimum number of periods to use
{
	int			i;

	for (i = 0; i < n; ++i)
	{
		vR[i] = TV1PAN(vS[i], abDateToJDN(vT[i]), K, abDateToJDN(E),
					   v, r, minPeriods);
	}
}

// TVVPED
// Theoretical Value, Vector of days, Put, European, Discrete dividends
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
	   int					minPeriods)	// minimum number of periods to use
{
	int			i;

	for (i = 0; i < n; ++i)
	{
		vR[i] = TV1PAD(vS[i], abDateToJDN(vT[i]), K, abDateToJDN(E),
					   v, r, nDividends, vDividends, minPeriods);
	}
}

//============================================================================
// Implied Volatility Calculations
//

//
//  Description:
//      Calculate the Implied Volatility for the given option. This uses basic
//      Newton-Raphson to find it to within +/- 0.001.
//
//  Arguments:
//      O       Option price (must be the same OHLC as S)
//      S       Share price (usually the closing price)
//      T       date of This price
//      E       Expiry date
//      K       striKe price
//      CorP    Call (or Put)?
//      v       annualised Volatility (used to start the search)
//      r       Riskless rate of return (zero if no carrying costs)
//		v	Number of dividends in ...
//		v	The dividends due on/before the expiry date.
extern "C" float
IV1day(float				O,			// Option price ...
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
	   int					minPeriods)	// minimum number of periods to use
{
    float	civ;
    float	ctv;
    float	liv;
    float	ltv;
    float	tiv;
    int		n;

    liv = 0.0;
    ltv = 0.0;
    civ = v;
    ctv = TV1day(S, T, K, E, CorP, civ,
				 r, nDividends, vDividends, ee, minPeriods);
    n = 1;

    while ((ctv - O) > 0.001 || (O - ctv) > 0.001)
	{
        tiv = (O - ltv) * (civ - liv) / (ctv - ltv) + liv;
        liv = civ;
        ltv = ctv;
        civ = tiv;
        ctv = TV1day(S, T, K, E, CorP, civ,
					 r, nDividends, vDividends, ee, minPeriods);
        n = n + 1;
    }

    return civ;
}
