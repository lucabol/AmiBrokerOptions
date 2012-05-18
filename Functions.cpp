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

#include <math.h>
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
	"mdaTVV",	{ TVV, 2, 0, 5, 6, ZeroValues },
	"mdaIVV",	{ IVV, 3, 0, 5, 6, ZeroValues },
};

int gFunctionTableSize = sizeof(gFunctionTable)/sizeof(FunctionTag);
