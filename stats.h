#ifndef STATS_H
#define STATS_H

class Stats
{	public:
	Stats();										// constructor
	void Add(double x, bool DblSrc, bool nonOpt);	// add runtime x, double sourcing
	double getAvg();								// get average
	double getStDev();								// get standard deviation
	int	getDblSrc();								// get no of double sourcing periods
	int	getNonOpt();								// get no of non-optimal solutions
	
	private:
	int n;			// number of values
	double avg;		// average
	double sumSq;	// sum of squares
	int nDblSrc;	// no of double sourcing periods
	int nNonOpt;	// no of non-optimal solutions
};

#endif