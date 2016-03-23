/*
 * dfa.h
 */

#ifndef L_H_
#define L_H_

#include <string>
#include <vector>
#include "dfa.h"
using namespace std;

class lstar
{
	private:
		dfa* tgdfaL;

		int n_membership_query;


	public:
		lstar(dfa* tgdfa);

		dfa* main_lstar(bool approximate, string samplestestpath);					// Return infered DFA
};



#endif /* L_H_ */
