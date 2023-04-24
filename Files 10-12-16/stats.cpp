#include <cmath>
#include <math.h>
#include "stats.h"

using namespace std;

Stats::Stats()
{	n=0;
	avg=0.0;	
	sumSq=0.0;	
	nDblSrc=0;
}

void Stats::Add(double x, bool DblSrc)
{	double delta;
	n++;
	delta=x-avg;
	avg+=delta/n;
	sumSq+=delta*(x-avg);
	nDblSrc+=DblSrc;
}

double Stats::getAvg()
{	return avg;
}

double Stats::getStDev()
{	double StDev;
	if (n>1)
		StDev=sqrt(sumSq/(n-1));	
	else // St Dev not defined
		StDev=-1.0;

	return StDev;
}

int Stats::getDblSrc()
{	return nDblSrc;
}