/*
 * obt.h
 *
 *  Created on: 01 ott 2015
 *      Author: piero
 */

#ifndef LSTAR_H_
#define LSTAR_H_

#include <string>
#include <map>
#include <vector>
#include <iostream>

#include "dfa.h"
#include "utilities.h"

using namespace std;


enum class obtState
{
	OBT_STATE_RED   = 1,
	OBT_STATE_BLU   = 0
};


/*! \class lstar
    \brief Class for LSTAR inference algorithm.

    Class for LSTAR inference algorithm, with all members and methods for start an inference process through Learner and Teacher.
 */
namespace gi
{
class lstar
{
	private:

	// 3 Strutture: - prefissi con associato il tipo, -risultato delle mq, -lista degli esperimenti
	vector< vector<SYMBOL> > 												experiments;							/*!< Columns of Observation Table */
	map<vector<SYMBOL>, bool> 											mq;											/*!< Local values of results for membership query already asked to Teacher */
	map<vector<SYMBOL>, obtState, vector_int_size_less>	pref;											/*!< Blue and red prefixes of Observation Tables */

	unsigned short int    dim_alfabeto;																					/*!< Alphabet size */
	vector<SYMBOL>*   	  alfabeto;																							/*!< Alphabet */

	dfa*    target;																															/*!< Target DFA of inference process. Used from Teacher. */

	public:

		/**
		 * Instance an object with all the members and methods for LSTAR inference process.
		 * @param path It's the target DFA for the inference process, used from Teacher.
		 */
		lstar(const dfa &targetdfa);

		/**
		 * Destroy an LSTAR object, freeing the memory.
		 */
		~lstar();
		
		/**
		 * Start an LSTAR inference process.
		 * @param approximate If it's true, set an LSTAR inference process using an approximation of target DFA.
		 * @param samplestestpath Samples approximating target DFA.
		 * @return The 
		 */
		dfa*  run_lstar(bool approximate, string samplestestpath);												// TODO: To implement approximate version


		/**
		 * Return a whole "row" of Observation Table for a prefix.
		 * @param prefix Unique prefix labeling a row.
		 * @return Row of observation table
		 */
		vector<obtState> get_row(vector<SYMBOL> prefix);
		

		/**
		 * For a non closed table, fix one or more problems (not necessary all) which keep not closed the table.
		 * @return Return 'true' if the table was modified (DOESEN'T mean table is closed!). Otherwise return 'false', that is:
		 * table was not modified, it was already closed.
		 */
		bool close_obt();

		/**
		 * For a non consistent table, fix one or more problems (not necessary all) which keep not consistent the table.
		 * @return Return 'true' if the table was modified (DOESEN'T mean table is consistent!). Otherwise return 'false', that is:
		 * table was not modified, it was already consistent.
		 */
		bool make_obt_consistent();

		/**
		 * Update the table using the witness.
		 * @param witness Witness is the counterexample coming from Teacher.
		 */
		void update_from_counterexample(vector<SYMBOL> witness);	// Update the OB from the counterexample
		

		/**
		 * Build a DFA from Observation Table.
		 * @return Pointer to builded DFA. 
		 */
		dfa* obt_to_dfa();


		/**
		 * Print the Observation table.
		 */
		void print_obt();
};

}

#endif /* LSTAR_H_ */
