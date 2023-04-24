#ifndef SOLCHAR_H
#define SOLCHAR_H

class sol_char
{	public:

	sol_char();				// constructor

	double obj;					// optimal objective value
	bool double_sourcing;		// there is a double sourcing period in the optimal solution
	bool double_sourcing_far;	// there is a double sourcing period more than 1 stage away
	int reg_network;			// counts the number of regeneration networks
};

#endif