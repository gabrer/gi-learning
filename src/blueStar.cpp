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

#include "utilities.h"
#include "blueStar.h"

#define STAT_PROMOTION	false

#define COLONNA_TIPO = "dim_alfabeto"

//#define ND 	 numeric_limits<int>::max()
//#define MINF numeric_limits<int>::min()
//#define MAX	 numeric_limits<double>::max()


using namespace std;



gi::blueStar::blueStar(const char * path, double alpha_value, double delta_value):bluefringe(path){
	alpha = alpha_value;
	delta = delta_value;
};

//TODO: verificare che venga invoca il distruttore della classe "bluefringe"
gi::blueStar::~blueStar(){};



double gi::blueStar::error_rate(dfaEDSM* dfa1, vector<SYMBOL>* positive, int dim_positive, vector<SYMBOL>* negative, int dim_negative, int* &wp, int* &wn, const int tot_wp_w, const int tot_wn_w)
{

	int wrong_recognized = 0;
	double e_rate = 0;

	// POSITIVE strings check
	for(int i=0; i<dim_positive; ++i)
		if(!dfa1->membership_query(positive[i])){
				wrong_recognized += wp[i];
			cerr << "EEEEEEEEEEEEEEEE Positive string rejected EEEEEEEEEEEEEEE" << endl;
		}

	// NEGATIVE strings check
	for(int i=0; i<dim_negative; ++i)
		if(dfa1->membership_query(negative[i])){
				wrong_recognized += wn[i];
			//cout << "Negative string accepted" << endl;
		}



	// Total error rate
	e_rate = (double) wrong_recognized/ (double) (tot_wp_w+tot_wn_w);


	#ifdef DEBUG1
	cout << "Stringhe errate: "<<wrong_recognized << " su "<< (double) (tot_wp_w+tot_wn_w) <<endl;
	#endif


	return e_rate;
}



double gi::blueStar::merge_heuristic_score(double error_rate_before, double error_rate_after, int dim_strings, double alpha, int earn_states)
{
	double z_alpha = DMAX;
	double T = error_rate_after - error_rate_before;
	double dev_std_h0 = DMAX;

	double ratio_error_to_n_states = 0;


	// It's "zeta_beta" because is the zeta involved to calculate beta. For our purpose is not necessary to calculate beta
	double zeta_beta = 0;


	// Compute Z_alpha under the H0 hypothesis
	z_alpha = z_alpha_from_u_alpha_two_proportions_test(error_rate_before, error_rate_after, dim_strings, alpha, &dev_std_h0);


	// z_alpha is MINF when the error is too low for normal distribution approximation
	if(z_alpha == DMINF)
		return DMINF * ((1+earn_states)*0.001);					// +1 is for case with earn_states=0; it's needed only the relative value


	// Statistic test
	if( T > z_alpha )
	{
		#ifdef DEBUG1
		cout << "MERGE RIGETTATO"<<endl;
		#endif

		return  DMAX;															// Rejected merge
	}
	else
	{

		zeta_beta = z_alpha - ( (double) delta / (double) dev_std_h0 );


		if(earn_states != 0)
			ratio_error_to_n_states = (-1)*zeta_beta / ((double) earn_states); //(log10(earn_states)/log10(2000));
		else
			ratio_error_to_n_states = (-1)*zeta_beta;


		#ifdef DEBUG1
		cout << "T: " << T << ", z_alpha: "<<z_alpha<<", delta: "<<delta<< "dev_std: "<<dev_std_h0<< endl;
		cout << "zeta-beta: "<< zeta_beta << endl;
		cout << "Stati guadagnati: "<<earn_states<<endl;
		#endif


		return ratio_error_to_n_states;
	}

}



gi::dfa* gi::blueStar::run(string base_path)
{
	// Samples from txtfile
	int n_symbols	 	 = 0;			//number of negative examples
	int dim_positive	 = 0; 		//number of positive examples
	int dim_negative = 0; 		//number of negative examples
	int* wp, *wn;
	int 	tot_wp_w=0, tot_wn_w=0;
	
	double  min_count;
	double *curr_count = NULL;
	double *curr_zalpha = NULL;
	
	int	n_red=0, n_blue=0, actual_blue_state=0;

	double error_rate_before = 0;

	// example strings
	vector<SYMBOL>* positive=NULL;
	vector<SYMBOL>* negative=NULL;
	char* symbols = NULL;

	bool promoted =false;
	bool at_least_one_promotion = false;
	
	dfaEDSM *dfa1 =NULL,  **merged = NULL; // dfa...

	// One dfa_best for every Blue state
	vector<dfaEDSM*> dfa_best;
	vector<double> dfa_score;

	// Blue states candidates to promotion
	vector<int> blue_state_index;
	vector<double> blue_promotion_score;


	//get positive and negative samples
	read_samples(positive, &dim_positive, negative, &dim_negative, wp, wn);

	
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



	for(int i=0; i<dim_positive; ++i)
		tot_wp_w += wp[i];
	for(int i=0; i<dim_negative; ++i)
		tot_wn_w += wn[i];

	cout << "Totale tot_wn_w: "<<tot_wp_w << " e "<< tot_wn_w << endl;


	cout <<" START BlueStar inference process..."<<endl;

	while_count=-1;
	// BLUESTAR
	while(n_blue>0)
	{		
		while_count++;
		
		promoted=false;

		// BLUE ciclo
		for(int i=0; i<n_blue && (!promoted || STAT_PROMOTION); ++i)
		{	

			/*cout << "I: "<<i<< ", n_blue: "<<n_blue << endl;
			cout << "Num blue states: "<<dfa1->get_num_blue_states() << endl;
			cout << "N effettivi: "<<dfa1->get_blue_states()->size() << endl;
			cout << "Old Actual blue state: "<<actual_blue_state << endl;*/
			actual_blue_state = dfa1->get_blue_states()->at(i);


			// Reset variable for the new run

			// array for the heuristic values of the red group, used for select best merge
			if(curr_count != NULL)
				delete[] curr_count;
			curr_count= new double [n_red];
			
			// array for the z_alpha values of the red group, used for select best promotion
			if(STAT_PROMOTION){
				if(curr_zalpha != NULL)
					delete[] curr_zalpha;
				curr_zalpha= new double [n_red];
			}

			// dfa coming from possible merges
			if(merged != NULL)
				delete[] merged;
			merged = new dfaEDSM*[n_red];
			
			// initialize values
			for(int j=0; j<n_red; ++j){
				curr_count[j] = (double) DMAX;
				if(STAT_PROMOTION)
					curr_zalpha[j] = (double) DMAX;
				merged[j] = NULL;
			}			


			// error rate for current dfa
			error_rate_before = error_rate(dfa1, positive, dim_positive, negative, dim_negative, wp, wn, tot_wp_w, tot_wn_w);

			#ifdef DEBUG1
			cout << "BLUE N:"<<n_blue<<endl;
			cout << "Error-rate automa PRE-merging: "<<error_rate_before<<endl;
			#endif

			short 	int error = 0;		// Used for detection of error (OpenMP doesn't allow to handle exception inside a parallel region)

			// RED ciclo
			#pragma omp parallel default(shared)
			{
			#pragma omp for
			for(int j=0; j<n_red; ++j){
				merged[j] = new dfaEDSM(*dfa1);


				merge(merged[j], dfa1->get_red_states()->at(j), actual_blue_state );


				// TODO: Questa riga si può probabilmente eliminare, da fare debug estensivo
				merged[j]->remove_blue_state(actual_blue_state);


				double error_rate_after = error_rate(merged[j], positive, dim_positive, negative, dim_negative, wp, wn, tot_wp_w, tot_wn_w);

				try
				{
					if(STAT_PROMOTION)
						curr_zalpha[j]	=	z_alpha_from_u_alpha_two_proportions_test(error_rate_before, error_rate_after, tot_wp_w+tot_wn_w, alpha);


					curr_count[j]	=	merge_heuristic_score(error_rate_before, error_rate_after, tot_wp_w+tot_wn_w, alpha, dfa1->get_actual_num_states() - merged[j]->get_actual_num_states());
				}
				catch( const char* msg )
				{
					cout << msg << endl;
					++error;
				}

				#ifdef DEBUG1
				cout << "- Blue con rosso -" << endl;
				cout << "Error rate DOPO merge: "<< error_rate_after <<endl;
				cout << "Numero di stati effettivi prima: "<<dfa1->get_actual_num_states()<<", dopo: "<< merged[j]->get_actual_num_states() << endl;
				cout << "Score: "<<curr_count[j]<<endl;
				#endif
			}
			}
			// end for RED

			if(error != 0)
				throw "This inference process is stopped with no DFA";


			// For Statistical purpose
			num_heuristic_merge_valued +=  n_red;
			

			// check if there some merge, else start process for promote
			promoted = true;
			min_count= DMAX;
			int j_min=ND;
			for(int j=0; j<n_red; ++j){
				if(curr_count[j] < min_count){
					min_count = curr_count[j];
					j_min = j;
					promoted=false;
				}			
			}

			//cout << "Max_count:"<<max_count<<endl;


			// EXECUTE Promotion OR Merge
			if(promoted)
			{
				// Immediatly promotion - not adoperated in Blue*
				if(!STAT_PROMOTION){
					promote(dfa1, actual_blue_state);

				}
				else
				{
					// Scelgo il valore minimo tra i curr_zalpha e lo prendo come rappresentante per questa promozione.
					// Quando finirò tutti i controlli, scegliero lo zalpha maggiore e promuoverò il relativo stato
					at_least_one_promotion = true;

					double min_zalpha = DMAX;
					for(int j=0; j<n_red; ++j)
						if(curr_zalpha[j] < min_zalpha)
							min_zalpha = curr_zalpha[j];

					blue_state_index.push_back(actual_blue_state);
					blue_promotion_score.push_back(min_zalpha);

					#ifdef DEBUG1
					cout << "Potenziale promozione, dimensione attuale: "<<blue_promotion_score.size() <<endl;
					for(int t=0; t<blue_promotion_score.size(); ++t)
						cout << "VALORE dentro: " << blue_promotion_score.at(t) << endl;
					#endif
				}



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
			else if(!at_least_one_promotion)	// - Merge accettato come candidato per il merge finale. Lo aggiungo alla lista dei migliori.
			{
				dfa_best.push_back(merged[j_min]);
				dfa_score.push_back(min_count);
			}


			// Free the array with dfa merged for calculate score, leave only the dfa selected as best
			if(merged != NULL){
				for(int j=0; j<n_red; ++j){
					if (j == j_min && (!at_least_one_promotion))			// Leave the best. If there is at least one promotion can delete all
						continue;
					if(merged[j] != NULL)
						delete merged[j];
				}
				delete[] merged;
				merged = NULL;
			}

		}// end for BLUE	


		// MERGE
		if(!at_least_one_promotion && !promoted) 	// Take the best merge as next DFA, no promotion done. Select best merge between all candidates in "dfa_best"
		{

			// Select the best merge between all the blue states (Looking for the lower value)
			double best_score = DMAX;
			int index_best_dfa = ND;
			for(int t=0; t<dfa_score.size(); ++t)
				if(dfa_score.at(t) < best_score){
					best_score = dfa_score.at(t);
					index_best_dfa = t;
				}


			// Take the blue states before delete the old dfa
			for(int t=0; t<dfa1->get_num_blue_states(); ++t)
				dfa_best.at(index_best_dfa)->add_blue_state( dfa1->get_blue_states()->at(t) );


			if(dfa1 != NULL) delete dfa1;		// Delete old dfa


			// set dfa1 to the new merged dfa
			dfa1 = dfa_best.at(index_best_dfa);
			nuoviBlu(dfa1);
			eliminaStati(dfa1);


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

			// Print information
			#ifdef DEBUG1
			cout << endl << "MERGE:"<< min_count<<endl;
			cout << "Indice: "<<index_best_dfa << endl;
			#endif
			//dfa1->print_dfa_with_color("DFA");
			//cout <<" ----------------------------------- "<<endl;
			#ifdef ALL_DFA_EDSM
				string name = "Merged"+intTostring(while_count);
				dfa1->print_dfa_dot(name, (base_path+name+".dot").c_str());
			#endif


		} else if(STAT_PROMOTION){ // Make best promotion, if Statistical Promotion is active

			at_least_one_promotion = false;

			// Select the best promotion between all the blue states (Looking for the upper value of zalpha, minimize the alpha area that is I type erro)
			double best_score = MINF;
			int index_best_blue_state = ND;
			for(int t=0; t<blue_promotion_score.size(); ++t)
				if(blue_promotion_score.at(t) > best_score){
					best_score = blue_promotion_score.at(t);
					index_best_blue_state = blue_state_index.at(t);
				}


			promote(dfa1, index_best_blue_state);


			#ifdef DEBUG1
			cout << "PROMOZIONE - Indice: "<<index_best_blue_state<<", score: "<<best_score<<endl;
			for(int t=0; t<blue_promotion_score.size(); ++t)
				cout << "VALORE: " << blue_promotion_score.at(t) << endl;
			#endif


			blue_promotion_score.clear();
			blue_state_index.clear();

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
		

		// update values for the dfa
		n_blue = dfa1->get_num_blue_states();
		n_red = dfa1->get_num_red_states();

		set_fringe_size(n_red,n_blue);

	}
	

	if(curr_count != NULL) delete[] curr_count;


	// Setto gli stati Rigettanti
	for(int i=0; i<dim_negative; ++i){
		int rigettante = dfa1->get_arrive_state(negative[i]);
		//if(rigettante != ND)
		//	dfa1->set_rejector_state(rigettante);
	}

	// Setto gli stati Accettanti
	for(int i=0; i<dim_positive; ++i){
		int accettante = dfa1->get_arrive_state(positive[i]);
		//if( accettante != ND)
			//dfa1->set_acceptor_state(accettante);
	}

	// Setto gli stati Eliminati
	eliminaStati(dfa1);
	//dfa1->print_dfa_with_color("AUTOMA FINALE");


	///////////////////////////////////////////////////////////////
	// Delete the unreachable states and insert, in needed, the "stato pozzo"
	dfaEDSM* finalDFA = dfa1->to_canonical_dfaEDSM_from_red_states();


	// Final error rate
	error_rate_final_dfa = error_rate(finalDFA, positive, dim_positive, negative, dim_negative,wp, wn, tot_wp_w, tot_wn_w);


	//finalDFA->print_dfa_with_color("EDSM PRE-MINIMIZZAZIONE AUTOMA FINALE");


	//////////////////////////////////////////////////////////////
	// Minimize returna a new dfa, then delete the older
	dfa* finalDFAmin = finalDFA->minimize_TF();

	if(finalDFA) delete finalDFA;

	if(positive) delete[] positive;
	if(negative) delete[] negative;
	if(symbols) delete[] symbols;


	return finalDFAmin;

}


double gi::blueStar::get_error_rate_final_dfa()
{
	return error_rate_final_dfa;
}

