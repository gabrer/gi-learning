/*
 * obt.h
 *
 *  Created on: 14/mag/2015
 *      Author: Gabriele
 */

#ifndef OBT_H_
#define OBT_H_

#include <string>
#include <map>
#include <vector>
#include <iostream>

#include "dfa.h"

using namespace std;

// Mi definisco un predicato che uso nella map per cambiare il modo in cui ordina di default
// Infatti in ordine lessicografico 1 viene dopo 01, mentre a me interessa che si tenga conto
//  anche della lunghezza e quindi 1 venga prima di 01
struct strsize_less
{
    bool operator()(std::string const& l, std::string const& r) const
    {
        if (l.size() < r.size())
            return true;

        if (l.size() > r.size())
            return false;

        return l < r;
    };
};


class obt
{
	private:

		typedef	map<string, char, strsize_less>::const_iterator It;			// Iteratore per le map
		typedef vector<string> string_v;

															// 3 Strutture: - prefissi con associato il tipo, -risultato delle mq, -lista degli esperimenti
		string_v 						exp;				// Contiene i titoli delle colonne
		map<string, char> 				mq;					// Risultato della Mem.query per una stringa (utile anche per sapere se è stata già richiesta al Teacher)
		map<string, char, strsize_less> pref;				// Prefissi: 0->BLU, 1->RED

		char    dim_alfabeto;
		//DA_ELIMINARE char*	symbols;
		string* alfabeto;

		dfa*    target;										// DFA target per la tabella. Usato per le mermbership query



	public:

		obt(dfa* targetdfa, /*char* symbol,*/ char dim_alfabeto, string* alf);

		bool closed();										// Return true if it is closed
		bool consistent();									// Return true if it is consistent
		void update_from_counterexample(string witness);	// Update the OB from the counterexample

		dfa* obt_to_dfa();
		void print_obt();
};


#endif /* OBT_H_ */
