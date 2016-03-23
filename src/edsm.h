/*
 * edsm.h
 */

/* NOVITA:
 * Eliminate le funzioni: elenco_stati_edsm, elenco_stati_rb, num_edsm_rb e num_edsm_state
 *
 */

#ifndef EDSM_H_
#define EDSM_H_

#include "bluefringe.h"



using namespace std;


/*! \class edsm
    \brief Class for EDSM inference algorithm.

    Class for EDSM inference algorithm, with all members and methods for start an inference process.
 */
namespace gi
{
class edsm : public bluefringe
{
private:


	/**
	 * Return a score for a dfa. It's calculated by heuristic emerged during Abbadingo Competition (1998).
	 *
	 * @param dfa1
	 * @param positive
	 * @param dim_positive
	 * @param negative
	 * @param dim_negative
	 * @return A score for a DFA, greater is better.
	 */
	int merge_heuristic_score(dfaEDSM* dfa1, vector<SYMBOL>* positive, int dim_positive, vector<SYMBOL>* negative, int dim_negative, int* wp , int* wn);

	double merge_heuristic_score(double error_rate_before, double error_rate_after, int dim_strings, double alpha, int earn_states){};


public:

	/**
	 * Instance an object with all the members and methods for EDSM inference process.
	 * @param path It's the path where find positive and negative samples
	 */
	edsm(const char* path);							// Take the path where find samples

	/**
	 * Destroy and EDSM object, freeing the memory.
	 */
	~edsm();


	/**
	 * Start an EDSM inference process.
	 * @param path It's the base path where create all the output files of EDSM.
	 * @return Inferred DFA
	 */
	dfa* run(string path);						// Argument is the base path where create files


};

}


#endif /* EDSM_H_ */
