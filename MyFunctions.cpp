////////////////////////////////////////////////////
// MyFunctions.cpp
//
// Copyright (C)2002 Matt Atterbury
// All rights reserved.
//
////////////////////////////////////////////////////

#include <math.h>
#include "OptionPricing.h"

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