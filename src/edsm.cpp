/*
 * edsm.cpp
 *
 */
 

#include <utilities.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <map>
#include <cmath>
#include <limits>


#include <cstdlib>  // For exit() function

#include "omp.h"		//OpenMP

#include "edsm.h"


#define ND numeric_limits<int>::max()
#define MINF numeric_limits<int>::min()


using namespace std;

/* Convezioni:
 *	 colonna[dim_alfabeto]   -> tipo;
 *	 tipo:   0-> neutro, 1-> accettante, 2-> rigettante, 3-> eliminato
 */
#define COLONNA_TIPO = "dim_alfabeto"


//TODO: cambiare in una define (così puoi attivarlo e disattivarlo in compilazione)
bool db = false;													// Attiva DEBUG


gi::edsm::edsm(const char * path):bluefringe(path){};

//TODO: verificare che venga invoca il distruttore della classe "bluefringe"
gi::edsm::~edsm(){};


int gi::edsm::merge_heuristic_score(dfaEDSM* dfa1, vector<SYMBOL>* positive, int dim_positive, vector<SYMBOL>* negative, int dim_negative, int* wp = NULL, int* wn=NULL)
{
	int* tp, *tn = NULL;
	tp = new int[dfa1->get_num_states()];							//Inizializzo un contatore per ogni stato
	tn = new int[dfa1->get_num_states()];

	for(int i=0; i<dfa1->get_num_states(); ++i){
		tp[i] = 0;
		tn[i] = 0;
	}


	// TODO Implementare COUNTER di OpenMP
	// Controllo quante stringhe POSITIVE vengono riconosciute
	for(int i=0; i<dim_positive; ++i){
		int statoFinale = dfa1->get_arrive_state(positive[i]);
		if(statoFinale != ND){
			if(wp == NULL)
				tp[statoFinale]++;
			else
				tp[statoFinale] += wp[i];
		}
	}

	//cout << endl;
	for(int i=0; i<dfa1->get_num_states(); ++i)
		if(db)
			cout << "tp["<<i<<"]: "<<tp[i]<<" ";

	// TODO Implementare COUNTER di OpenMP
	// Controllo quante stringhe NEGATIVE vengono riconosciute
	for(int i=0; i<dim_negative; ++i){
		int statoFinale = dfa1->get_arrive_state(negative[i]);
		if(statoFinale != ND){
			if(wn == NULL)
				tn[statoFinale]++;
			else
				tn[statoFinale] += wn[i];
		}
	}

	//cout << endl;
	for(int i=0; i<dfa1->get_num_states(); ++i)
		if(db)
			cout << "tn["<<i<<"]: "<<tn[i]<<" ";

	// Calcolo lo SCORE
	int sc = 0;
	for(int i=0; i<dfa1->get_num_states(); ++i)
	{
		if(sc != MINF)
		{
			if(tn[i] > 0){
				if(tp[i] > 0){
					if(db)
						cout <<endl<< "Stato che accetta Positiva&Negativa:"<<i<<endl;
					sc = MINF;
					break;								// E' inutile continuare a ciclare se metto sc = ND
				}
				else
					sc = sc + tn[i] -1;					// Nello pseudo De La Higuera mette -1 (in situazioni complesse porta ad automi finali diversi)
			}else{
				if(tp[i] > 0)
					sc = sc + tp[i] -1;
			}
		}
	}

	if(tp!=NULL)
		delete[] tp;
	if(tn!=NULL)
		delete[] tn;

	return sc;
}


gi::dfa* gi::edsm::run(string base_path)
{
	// Samples from txtfile
	int n_symbols=0;	//number of negative examples
	int dim_positive=0; //number of positive examples
	int dim_negative=0; //number of negative examples
	
	int max_count, *curr_count=NULL;
	
	int	n_red=0, n_blue=0, actual_blue_state=0;

	// example strings
	vector<SYMBOL>* positive=NULL;
	vector<SYMBOL>* negative=NULL;
	char* symbols = NULL;
	int* wp, *wn;
	
	//bool *promotion=NULL;
	bool promoted =false;
	
	dfaEDSM *dfa1 =NULL, /* *dfa_best = NULL,*/ **merged = NULL; // dfa...

	// One dfa_best for every Blue state
	vector<dfaEDSM*> dfa_best;
	vector<int> dfa_score;


	//get positive and negative samples
	read_samples(positive, &dim_positive, negative, &dim_negative, wp, wn);

	// Costruisco PTA
	//dfaEDSM* dfa1 = build_pta(positive,dim_positive);
	
	// Build APTA
	dfa1 = build_apta(positive, dim_positive, negative, dim_negative);

	// Print it!
	if(dfa1->get_num_states() < 1000)
	{
		dfa1->print_dfa_dot("APTA", (base_path+"APTA.dot").c_str() );
		dfa1->print_dfa_dot_mapped_alphabet("APTAALF", (base_path+"APTA_ALF.dot").c_str());
	}else{
		clog<<"APTA too big! I can't print it"<<endl;
	}

	n_blue = dfa1->get_num_blue_states();
	n_red = dfa1->get_num_red_states();

	set_fringe_size(n_red,n_blue);

	cout <<" START Edsm inference process..."<<endl;

	while_count=-1;
	// ESDM
	while(n_blue>0)
	{		
		while_count++;
		
		promoted=false;

		// BLUE ciclo
		for(int i=0; i<n_blue && (!promoted); ++i)
		{	
			actual_blue_state = dfa1->get_blue_states()->at(i);
			// Reset variable for the new run

			// array for the heuristic values of the red group
			if(curr_count != NULL)
				delete[] curr_count;
			curr_count= new int [n_red];
			
			// dfa coming from possible merges
			if(merged != NULL)
				delete[] merged;
			merged = new dfaEDSM*[n_red];
			
			// initialize values
			for(int j=0; j<n_red; ++j){
				curr_count[j] = MINF;
				merged[j] = NULL;
			}			
			

			// RED ciclo
			#pragma omp parallel default(shared)
			{
			#pragma omp for
			for(int j=0; j<n_red; ++j){
				merged[j] = new dfaEDSM(*dfa1);

				merge(merged[j], dfa1->get_red_states()->at(j), actual_blue_state );

				// TODO: Questa riga si può probabilmente eliminare, da fare debug estensivo
				merged[j]->remove_blue_state(actual_blue_state);

				curr_count[j] = merge_heuristic_score(merged[j], positive, dim_positive, negative, dim_negative);
			}
			}
			// end for RED

			// For Statistical purpose
			num_heuristic_merge_valued +=  n_red;
			
			// check if there some merge, else start process for promote
			promoted = true;
			max_count=MINF;
			int j_max=ND;
			for(int j=0; j<n_red; ++j){
				if(curr_count[j]>max_count){
					max_count = curr_count[j];
					j_max = j;
					promoted=false;
				}			
			}

			//cout << "Max_count:"<<max_count<<endl;


			// PROMOTION
			if(promoted){

				promote(dfa1, actual_blue_state);

				//cout << "PROMOZIONE"<<endl;

				#ifdef ALL_DFA_EDSM
				string name = "P"+dfa::intTostring(while_count)+dfa::intTostring(blue[i]);
				dfa1->print_dfa_dot(name, (base_path+name+".dot").c_str());
				#endif


				//Free memory
				typedef	vector<dfaEDSM*>::const_iterator It;
				for(It p1=dfa_best.begin(); p1!=dfa_best.end(); ++p1){
					if(dfa1 == (*p1))
						continue;
					delete (*p1);
				}

				dfa_best.clear();
				dfa_score.clear();
			}
			else	// - Merge accettato come candidato per il merge finale. Lo aggiungo alla lista dei migliori.
			{
				//merged[j_max]->remove_blue_state(actual_blue_state);

				dfa_best.push_back(merged[j_max]);
				dfa_score.push_back(max_count);
			}


			// Free the array with dfa merged for calculate score, leave only the dfa selected as best
			if(merged != NULL){
				for(int j=0; j<n_red; ++j){
					if (j == j_max && (!promoted))			// Leave the best
						continue;
					if(merged[j] != NULL)
						delete merged[j];
				}
				delete[] merged;
				merged = NULL;
			}
		}// end for BLUE	
		
		
		// MERGE
		if(!promoted){ // Do definitive merge, no promotion done. Select best merge between all candidates in "dfa_best"

			// Select the best merge between all the blue states
			int best_score = -1;
			int index_best_dfa = 0;
			for(int t=0; t<dfa_score.size(); ++t)
				if(dfa_score.at(t) > best_score){
					best_score = dfa_score.at(t);
					index_best_dfa = t;
				}

			// Take the blue states before delete the old dfa
			//int colonna_tipo = dfa1->get_dim_alphabet();
			for(int t=0; t<dfa1->get_num_blue_states(); ++t)
				dfa_best.at(index_best_dfa)->add_blue_state( dfa1->get_blue_states()->at(t) );
			if(dfa1 != NULL) delete dfa1;		// Delete old dfa

			// set dfa1 to the new merged dfa
			dfa1 = dfa_best.at(index_best_dfa);
			nuoviBlu(dfa1);
			eliminaStati(dfa1);


			// Print information
			//cout << "MERGE:"<<max_count<<endl;
			//cout <<" ----------------------------------- "<<endl;
			#ifdef ALL_DFA_EDSM
			string name = "M"+dfa::intTostring(while_count);
			dfa1->print_dfa_dot(name, (base_path+name+".dot").c_str());
			#endif

			++num_actual_merge;

			// Free memory
			typedef	vector<dfaEDSM*>::const_iterator It;
			for(It p1=dfa_best.begin(); p1!=dfa_best.end(); ++p1){
				if(dfa1 == (*p1))
					continue;

				if((*p1) != NULL)
					delete (*p1);
			}

			dfa_best.clear();
			dfa_score.clear();
		}
		
		// update values for the dfa
		n_blue = dfa1->get_num_blue_states();
		n_red = dfa1->get_num_red_states();

		set_fringe_size(n_red,n_blue);

	}
	

	if(curr_count != NULL) delete[] curr_count;
	

	// Setto gli stati Accettanti
	int colonna_tipo = dfa1->get_dim_alphabet();
	for(int i=0; i<dim_positive; ++i){
		int accettante = dfa1->get_arrive_state(positive[i]);

//			if(accettante != ND)
//				dfa1->get_ttable()[accettante][colonna_tipo] = DFA_STATE_ACCEPTING;
	}

	// Setto gli stati Rigettanti
	for(int i=0; i<dim_negative; ++i){
		int rigettante = dfa1->get_arrive_state(negative[i]);
		if(rigettante != ND){
			//cout << "Statp di arrivoN: "<<rigettante<<endl;
			dfa1->get_ttable()[rigettante][colonna_tipo] = DFA_STATE_REJECTING;
		}
	}

	// Setto gli stati Eliminati
	eliminaStati(dfa1);
	//dfa1->print_dfa_with_color("AUTOMA FINALE");


	///////////////////////////////////////////////////////////////
	// Delete the unreachable states and insert, if needed, the sink state
	dfaEDSM* finalDFA = dfa1->to_canonical_dfaEDSM_from_red_states();


	//////////////////////////////////////////////////////////////
	// Minimize returna a new dfa, then delete the older
	dfa* finalDFAmin = finalDFA->minimize_TF();

	if(finalDFA) delete finalDFA;

	if(positive) delete[] positive;
	if(negative) delete[] negative;
	if(symbols) delete[] symbols;


	return finalDFAmin;

}

