/*
 * lstar2.cpp
 *
 *  Created on: 14/mag/2015
 *      Author: Gabriele
 */

//#include <iostream>
//#include <string>
//#include <vector>
//#include <fstream>
//#include <map>
//#include <cmath>
//#include <stdlib.h>

#include "dfa.h"
#include "lstar.h"
#include "obt.h"

lstar::lstar(dfa* dfatg)
{
	tgdfaL = dfatg->copyDFA();

	//ANGLUIN
	/*tg[0][0]=1;
	tg[0][1]=2;
	tg[0][2]=1;
	 tg[1][0]=0;
	 tg[1][1]=3;
	 tg[1][2]=0;
	tg[2][0]=3;
	tg[2][1]=0;
	tg[2][2]=0;
	 tg[3][0]=2;
	 tg[3][1]=1;
	 tg[3][2]=0;*/
}

dfa* lstar::main_lstar(bool approximate, string samplestestpath)
{
	dfa* dfahp = NULL;

	// Minimize tgdfa
	tgdfaL = tgdfaL->minimize_TF();

	// Create an array of symbols of alphabet
	string* alf = new string[tgdfaL->get_dim_alp()];
	for(int i=0; i<tgdfaL->get_dim_alp(); ++i)
		alf[i] = dfa::intTostring(i);

	// Observation table
	obt* obtable = new obt(tgdfaL, /*tgdfaL->get_symbols(),*/ tgdfaL->get_dim_alp(), alf);

	bool closedTable = false;
	bool consistentTable = false;
	bool finded = false;

	int count_generation = 0;
	while(!finded)
	{
		#ifdef DEBUG_1
		cout << endl << endl << "********************************" << endl;
		cout << count_generation << " GENERATION"<< endl;
		cout << "********************************" << endl;
		#endif

		closedTable = false;
		consistentTable = false;

		// Close and consistent
		while(!closedTable || !consistentTable )
		{
			closedTable = obtable->closed();
			consistentTable = obtable->consistent();
		}

		// Print transition table
		#ifdef DEBUG_2
		obtable->print_obt();
		#endif

		// DFA from transition table
		dfa* dfatmp = obtable->obt_to_dfa();
		#ifdef DEBUG_2
		dfatmp->print_dfa_ttable("DFA_HP dalla TABELLA DI OSSERVAZIONE");
		#endif

		count_generation++;


		// Generate witness if necessary
		string witness = "";
		if(!approximate)
			witness = tgdfaL->equivalence_query(dfatmp);
		else
			witness = tgdfaL->equivalence_query_approximate(dfatmp,samplestestpath);

		/*STAMPA DOT DEI FARI PASSAGGI if(!approximate){
			string nome = "W"+dfa::intTostring(count_generation);
			//merged->print_dfa_with_color(nome);
			string percorso = "/Users/Gabriele/Dropbox/Workspace_eclipse/lstar/lstar/src/intermezzo2/"+nome+".dot";
			dfatmp->print_dfa_dot(nome, percorso.c_str());
		}*/

		// Findend automata when witness is 'e', else update the obb. table
		if(witness.compare("e")){
			obtable->update_from_counterexample(witness);
		}
		else
		{
			cout << "LSTAR: Automa trovato!"<<endl;
			dfatmp = dfatmp->minimize_TF();
			//dfatmp->print_dfa_ttable("DFA TARGET INFERITO");
			dfahp = dfatmp;
			//tgdfaL->print_dfa_ttable("DFA TARGET TEACHER");
			finded = true;
		}
	}

	dfahp->set_n_memb_query(tgdfaL->get_n_memb_query());

	delete tgdfaL;

	return dfahp;
}

