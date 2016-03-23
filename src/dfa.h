/*
 * dfa.h
 */

#ifndef DFA_H_
#define DFA_H_

#include <string>
#include <vector>
#include <map>
#include <limits>

#include "utilities.h"

// Type of states
#define DFA_STATE_NON_ACCEPTING 2						/*!< Non accepting state */
#define DFA_STATE_ACCEPTING 1								/*!< Accepting state */
#define DFA_STATE_REJECTING 0									/*!< Rejecting state */
#define DFA_STATE_UNREACHABLE 3							/*!< Unreachable state from other states. Usually is a state to be delete. */

#define ND numeric_limits<int>::max()					/*!< Usually adopted for a non initialized state. */

typedef unsigned short int SYMBOL;

using namespace std;


/*! \class dfa
    \brief Class DFA.

    Class for DFA representation.
*/

//! "gi" is the library namespace
namespace gi
{
class dfa
{
	protected:
		int** ttable;															/*!< Transition table */
		int	  num_states;													/*!< Number of dfa states */
		int	  start_state;													/*!< Index of start state */

		int	  dim_alphabet;												/*!< Size of alphabet */
		char* alphabet;														/*!< Alphabet symbols */
		map<char, SYMBOL> mapped_alphabet;				/*!< Convert symbols into indices */

		/**
		 * Set the transition table (ttable) reference to an extern transition table passed as argument
		 * @param ext_ttab
		 */
		void 	set_ttable(int** ext_ttab);

		/**
		 * Set the number of states
		 * @param n
		 */
		void 	set_num_state(int n);

		/**
		 * Set the alphabet of dfa to alphabet in input of size "d_alf".
		 * Set also "mapped_alphabet" with an index for every symbol.
		 * @param alf
		 * @param d_alf
		 */
		void 	set_alphabet(const char* alf, const int d_alf);

	public:

		/**
		 * Make an instance of null dfa
		 */
		dfa();

		/**
		 * Make an instance of new dfa with default start state to 0
		 * @param n_state	Number of states
		 * @param dim_alf	Size of alphabet
		 * @param alf		Alphabet symbols
		 */
		dfa(const int n_state, const int dim_alf, const char *alf);

		/**
		 * Make an instance of new dfa. Give possibility to set the start state.
		 * @param n_state	Number of states
		 * @param dim_alf	Size of alphabet
		 * @param alf		Alphabet symbols
		 * @param s_state	Start state
		 */
		dfa(const int n_state, const int dim_alf, const char *alf, const int s_state);

		/**
		 * Make an instance of new dfa.
		 * Set tra transistion table making a copy of "tt_copy" passed as argument.
		 * @param n_state	Number of states
		 * @param dim_alf	Size of alphabet
		 * @param alf		Alphabet symbols
		 * @param s_state	Start state
		 * @param tt_copy	Reference to extern transition table to copy inside current dfa
		 */
		dfa(const int n_state, const int dim_alf, const char *alf, const int s_state, const int** tt_copy );

		/**
		 * Constructor for make a copy of a dfa "d1"
		 * @param d1	Dfa to copy
		 */
		dfa(const dfa &d1);

		/**
		 * Destroy a DFA object, freeing the memory.
		 */
		~dfa();

		/**
		 * Return size of alphabet
		 * @return Size of alphabet
		 */
		int  	get_dim_alphabet();

		/**
		 * Return a pointer to alphabet symbols
		 * @return Pointer to alphabet symbols
		 */
		char*	get_alphabet();

		/**
		 * Get number of states
		 * @return Number of states
		 */
		int   	get_num_states();

		/**
		 * Get index of start state
		 * @return Index of start state
		 */
		int   	get_start_state();

		/**
		 * Get index of arrive state for dfa_string argument
		 * @param dfa_string String executed on dfa
		 * @return Index of arrive state for "dfa_string"
		 */
		int		get_arrive_state(vector<SYMBOL> dfa_string);

		/**
		 * Return a reference to ttable()
		 * @return Pointer to ttable
		 */
		int** 	get_ttable();

		/**
		 * Get value of ttable for index "i", "j"
		 * @param i
		 * @param j
		 * @return
		 */
		int 	get_ttable(int i, int j);

		/**
		 * Set a single value "v" for ttable entry for "i","j"
		 * @param i First coordinate
		 * @param j	Second coordinate
		 * @param v	Value to set
		 */
		void 	set_ttable_entry(int i, int j, int v);


		void 	set_acceptor_state(int state);

		void 	set_rejector_state(int state);


		/**
		 * Make a new dfa from the union of current dfa and "dfa_hp".
		 * The first states are from current dfa, last states from "dfa_hp". The total number of states are sum of the number of states for the 2 dfa.
		 * @param dfa_hp Dfa to add to the current dfa
		 * @return Pointer to the union dfa of two dfa
		 */
		dfa*		unionDFA(dfa* dfa_hp);										// Return a union dfa of (this) and "dfa_hp"

		/**
		 * Minimizes the current dfa
		 * @return Pointer to a new instance of dfa, that is the minimized current dfa
		 */
		dfa*  	minimize_TF();												// Minimization using Table-filing algorithmm


		/** Return a measures of complexity for current DFA.
		 *
		 * @return Default value is the number of states.
		 */
		double		get_complexity();

		/**
		 * Print the transition table of current dfa. Before the transition table print the title passse as parameter.
		 * @param title Title printed before the transition table
		 */
		void 	print_dfa_ttable(string title);

		/**
		 * Print the transition table using the alphabet symbols. Before the transition table print the "title".
		 * @param title Title printed before the transition table
		 */
		void 	print_dfa_ttable_mapped_alphabet(string title);

		/**
		 * Print a dot file for the current dfa, with title "title", in the path "file_path".
		 * @param title	Title printed before the transition table
		 * @param file_path Path where make a dot file
		 */
		void    print_dfa_dot(string title, const char *file_path);

		/**
		 * Print a dot file for the current dfa using the alphabet symbols, with title "title", in the path "file_path".
		 * @param title	Title printed before the transition table
		 * @param file_path Path where make a dot file
		 */
		void 	print_dfa_dot_mapped_alphabet(string title, const char *file_path);

// LSTAR
		/**
		 * Read a dfa from a file
		 * @param file_name Path of the dfa
		 * @return Dfa read from file
		 */
		static dfa read_dfa_file(const string file_name);

		/**
		 * Fills "positive" and "negative" vector with examples inside file in path "path_samples".
		 * Set also the number of positive and negative example in "dim_positive" and "dim_negative".
		 * @param path_samples
		 * @param positive
		 * @param dim_positive
		 * @param negative
		 * @param dim_negative
		 */
		//void 	read_example_file(string path_samples, vector<SYMBOL>* &positive, int* dim_positive, vector<SYMBOL>* &negative,  int *dim_negative);


		/**
		 * Make a membership query to dfa with the "str" string. Return "true" if the arrive state for "str" is acceptor, else "false".
		 * @param str A string to make a membership query.
		 * @return "True" o "false" depending on the arrive state: "true" if acceptor, else if not.
		 */
		bool	membership_query(vector<SYMBOL> str);

		/**
		 * Make an equivalence query, that return an empty witness if the current dfa is equivalent to dfa "dfa_hp", argument of the function.
		 * Otherwise return a witness, that is a counterexample that distinguishes the two dfa.
		 * @param dfa_hp Pointer to dfa to compare for equivalence.
		 * @return A conterexample, empty if the two dfa are equivalent.
		 */
		vector<SYMBOL>	   equivalence_query(dfa* dfa_hp);

		/**
		 * Fills a table with the "Table Filling Algorithm", useful for find the equivalent/distinct states, and also for generate a witness.
		 * The Table considered is only the upper triangular matrix, that can be saved in a linear array.
		 * @return A table saved in a linear array.
		 */
		SYMBOL*		  	   table_filling();

		/**
		 * Create a list of states and corrispective equivalent states
		 * @param distincts A table build with the Table Filling Algorithm
		 * @return A pointer to the first vector<int>. Every vector is a list of equivalent state for the state associate to that vector.
		 */
		vector<int>*  	   equivalent_states_list_from_table(SYMBOL* distincts);

		/**
		 * Make a conterexample from a table build with Table Filling Algorithm using the union dfa of the two dfa to compare.
		 * @param distinct Table build with Table Filling Algorithm
		 * @param start_state_dfa_hp Index of the first state of dfa_hp inside the union dfa
		 * @return A witness that distinguishes the two dfa.
		 */
		vector<SYMBOL>	   witness_from_table(SYMBOL* distinct, int start_state_dfa_hp);

		/**
		 * Like the equivalence query function, but it not generate the witness but
		 * use a limitated set of string where select the conterexample.
		 * @param dfa_hp Pointer to dfa to compare for equivalence.
		 * @param samplestestpath Path of file with a list of counterexample.
		 * @return	A conterexample, empty if the two dfa are equivalent.
		 */
	//	vector<SYMBOL> 	   equivalence_query_approximate(dfa* dfa_hp,string samplestestpath);

		/**
		 * Used when dfa target in LSTAR is not disponibile, instead use a samples
		 * @param dfa_hp Pointer to dfa to compare for equivalence.
		 * @param samplestestpath Path of file with a list of counterexample.
		 * @return A witness that distinguishes the two dfa.
		 */
		//vector<SYMBOL>	   witness_approximate(dfa* dfa_hp, string samplestestpath);
};

}


struct vector_int_size_less
{
    bool operator()(vector<SYMBOL> l, vector<SYMBOL> r) const
    {
        if (l.size() < r.size())
            return true;

        if (l.size() > r.size())
            return false;

        typedef	vector<SYMBOL>::const_iterator It_posneg;

        for(It_posneg i=l.begin(), j=r.begin(); i<l.end(); ++i, ++j){
        	if(*i == *j)
        		continue;
        	else
        		return (*i) < (*j);
        }

        return false;
    };
};



#endif /* DFA_H_ */
