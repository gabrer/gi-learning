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

#include "omp.h"	//OpenMP

#include "edsm.h"


using namespace std;

/* Convezioni:
 *	 colonna[dim_alfabeto]   -> tipo;
 *	 tipo:   0-> neutro, 1-> accettante, 2-> rigettante, 3-> eliminato
 */
#define COLONNA_TIPO = "dim_alfabeto"


//TODO: cambiare in una define (così puoi attivarlo e disattivarlo in compilazione)
//bool db = false;													// Attiva DEBUG

gi::bluefringe::bluefringe(const char * path){
	path_samples = new char[strlen(path)+1];

	dim_alphabet = 0;
	num_actual_merge =0;
	num_heuristic_merge_valued=0;

	fringe_size[0].reserve(1000); 		// Size of blue fringe
	fringe_size[1].reserve(1000);			// Size of red fringe

	strcpy(path_samples, path);

	inverse_mapped_alphabet = NULL;
	while_count = -1;
}


gi::bluefringe::~bluefringe()
{
	if(inverse_mapped_alphabet != NULL)
		delete[] inverse_mapped_alphabet;

	mapped_alphabet.clear();

	delete[] path_samples;

}


void gi::bluefringe::read_samples(vector<SYMBOL>* &positive, int* dim_positive, vector<SYMBOL>* &negative,  int *dim_negative,int* &wp, int* &wn)
{
	cout << "Reading strings from txt file: "<<endl;
	int cp = 0;														// Numero di stringhe positive del linguaggio
	int cn = 0;														//   -    -     -     negative  -      -
	char ch;

	char null_symbol;

	cout << path_samples << endl;

	fstream fin(path_samples, fstream::in);

	if(!fin){
		cerr<<"An error occurred in opening file :("<<endl;
		exit(EXIT_FAILURE);
	}

	while (fin >> noskipws >> ch) {									// Faccio in conteggio previo del numero di stringhe positive e negative presenti nel txt
		if(ch == '\n')
			continue;
		else if(ch == '+')
			cp++;
		else if(ch == '-')
			cn++;
	}
	(*dim_positive) = cp;
	(*dim_negative) = cn;

	positive = new vector<SYMBOL>[cp];
	negative = new vector<SYMBOL>[cn];

	wp	= new int[cp];
	wn	= new int[cn];
	for(int i=0; i<cp; ++i)
		wp[i] = 0;
	for(int i=0; i<cn; ++i)
		wn[i] = 0;

	cout << intTostring(cp) + " positivi"<<endl;
	cout << intTostring(cn) + " negativi"<<endl;

	int flagcp = 0;
	int flagcn = 0;
	bool casopositive = false;
	bool casonegative = false;
	bool primap = true;
	bool priman = true;

	ifstream infile(path_samples);

	bool first = true;
	bool second = false;
	string line;

	while (getline(infile, line))
	{
	    istringstream iss(line);
	    int a;
	    string n;


	    // Read first line for dim alphabet
	    if(first)
	    {
	    	if (!(iss >> a)) { break; } // error
	    	dim_alphabet = a;
	    	//cout << "dimensione alfabeto " << a << endl;
	    	first = false;
	    	second = true;

	    	continue;
	    }


	    // Read second line for alphabet symbol
	    if(second)
	    {
	    	inverse_mapped_alphabet = new char[dim_alphabet];

	    	int counter=-1;
	    	while (iss >> n){
	    		if(counter==-1){
	    			null_symbol = (char) n[0];

	    			++counter;
	    			continue;
	    		}else if(counter>=dim_alphabet)
	    			break;

	    		mapped_alphabet[(char) n[0]] = counter;
	    		inverse_mapped_alphabet[counter++]=n[0];
	    	}

	    	// Alphabet
	    	if(counter!= dim_alphabet){
	    		cerr<<"Error in reading example: number of red alphabet symbols mismatches with the declared one!"<<endl;
	    		cerr<<"Expected symbols: "<<dim_alphabet<<endl;
	    		cerr<<"Red symbols: "<<counter<<endl;


	    		exit(EXIT_FAILURE);
	    	}

	    	// alphabet ok ;)
	    	second= false;
	    }


	    bool weight = true;

	    // Read remaining lines
		while (iss >> n)
		{
			if( !n.compare("+") ){

				weight = true;
				casopositive = true;
				casonegative = false;
				if(primap){												// Se è il primo caso evito l'incremento
					primap =false;
					continue;
				}
				flagcp++;
				continue;

			}else if( !n.compare("-") ){
				weight = true;
				casonegative = true;
				casopositive = false;
				if(priman){												// Se è il primo caso evito l'incremento
					priman =false;
					continue;
				}
				flagcn++;
				continue;
			}

			// se la stringa è vuota, non è necessario aggiungere nulla
			if(((char) n[0]) == null_symbol)
				continue;


			if(weight){
				weight = false;

				if(casopositive)
					wp[flagcp] = stringToint(n);
				else if(casonegative)
					wn[flagcn] = stringToint(n);

				// ELIMINARE DEBUG
				cout << "!T Peso: "<< stringToint(n) << endl;

			} else {

				int tmp = mapped_alphabet[(char) n[0]];

				if(casopositive)
							positive[flagcp].push_back(tmp);
				else if(casonegative)
							negative[flagcn].push_back(tmp);
			}

		}
	}

}


// build apta from sample set
gi::dfaEDSM* gi::bluefringe::build_apta(const vector<SYMBOL>* positive, const int dim_positive, const vector<SYMBOL>* negative, const int dim_negative)
{
	cout << "Costruzione PTA..."<<endl;

	// PREFISSI
	// Calcolo i prefissi e li salvo insieme ad un indice indicatore dello stato
	map<vector<SYMBOL>,int, vector_int_size_less> prefissi;
	typedef	map<vector<SYMBOL>,int,vector_int_size_less>::const_iterator It;

	typedef	vector<SYMBOL>::const_iterator It_posneg;

	// add empty string
	vector<SYMBOL> emptyVec(0);
	prefissi[emptyVec] = 0;

	for(int i=0; i<dim_positive; ++i){
		It_posneg w = positive[i].begin();
		for(It_posneg j=positive[i].begin(); j!=positive[i].end(); ++j){

			vector<SYMBOL>   sub(w,j+1);

			if(prefissi[sub] == 0)
				prefissi[sub] = ND;			// Qui e dopo sono 0 perché non li setto adesso!
		}
	}

	for(int i=0; i<dim_negative; ++i){
		It_posneg w = negative[i].begin();
		for(It_posneg j=negative[i].begin(); j!=negative[i].end(); ++j){
			vector<SYMBOL>   sub(w,j+1);
			if(prefissi[sub] == 0)
				prefissi[sub] = 0;
		}
	}

	// Assegno gli indici di stato corretti per i relativi prefissi
	int indice_stato=0;
	for(It p1=prefissi.begin(); p1!=prefissi.end(); ++p1){
		prefissi[(*p1).first] = indice_stato;
		++indice_stato;
	}

	/* Codice verifica prefissi*/
//	if(false)
//		for(It p1=prefissi.begin(); p1!=prefissi.end(); ++p1){
//			for(It_posneg i=p1->first.begin(); i<p1->first.end(); ++i)
//				cout << (*i);
//			cout<< "; " << p1->second << endl;
//		}

	// *** PTA ***
	dfaEDSM* aptaDFA = new dfaEDSM(prefissi.size(), dim_alphabet, inverse_mapped_alphabet,0);
	int** pta = aptaDFA->get_ttable();

	bool primo_stato=true;
	for(It p=prefissi.begin(); p!=prefissi.end(); ++p)
	{
		vector<SYMBOL> stato = (*p).first;
		//cout << "Stato considerato:"<<stato<<endl;


		// Solo per il Primo Stato -gestisco a parte l'essere red e intanto non accettante-
		// Controllo per ogni elemento dell'alfabeto, per quali è definita una transizione dallo stato iniziale
		if(primo_stato)
		{
			for(int i=0; i<dim_alphabet; ++i){							// Setto il valore delle Transizioni
				vector<SYMBOL> symbol;
				symbol.push_back(i);
				if(prefissi.find(symbol) != prefissi.end())				// Controllo che sia definito uno stato per quando entra dopo lambda un elemento dell'alfabeto
					pta[(*p).second][i] = prefissi[ symbol ];
			}

			pta[(*p).second][dim_alphabet] = 0;							// Primo stato: non accettant
			aptaDFA->add_red_state((*p).second);						// Aggiungo lo stato rosso anche alla gestione tramite VECTOR

			primo_stato=false;
			continue;
		}


		// Stati Accettanti
		for(int i=0; i<dim_positive; ++i){								// Pongo gli stati opportuni ad accettanti: quelli la cui stringa è anche tra le stringhe positive
			if(stato == positive[i]){
				pta[(*p).second][dim_alphabet] = 1;						// Con dim_alfabato arrivo alla prima colonna oltre le lettere dell'alfabeto
				break;
			}
			else
				pta[(*p).second][dim_alphabet] = 0;
		}


		// Transizioni
		for(int i=0; i<dim_alphabet; ++i){						// Verifico per ogni lettera dell'alfabeto che esista la transizione; se esiste memorizzo
			vector<SYMBOL> temp_vect(stato);
			temp_vect.push_back(i);
			if(prefissi.find(temp_vect) != prefissi.end())		// nella riga dello stato considerato ((*p).second) la relativa transizione
				pta[(*p).second][i] = prefissi[temp_vect];
		}

	}

	//Coloro gli stati blu (il primo livello dopo lambda in ESDM):
	for(int i=0; i<dim_alphabet; ++i)
	{
		int statoblu = pta[0][i];
		if(pta[0][i] != ND)
			aptaDFA->add_blue_state(statoblu);							// Lo aggiungo anche al Vector di stati blu
	}

	// Stampo APTA
	/*aptaDFA->print_dfa_with_color("*** APTA DENTRO ***");
	aptaDFA->print_dfa_dot("apta",PATH_DOT_APTA);*/

	cout << "PTA costruito"<<endl;

	return aptaDFA;
}


void gi::bluefringe::merge(dfaEDSM* dfa1, int redstate, int blustate)
{
	int predecessore=0;
	int lettera;
	int num_predecessori=0;

	// Cerco il predecessore dello stato blue
	for(int i=0; i<dfa1->get_num_states(); ++i)
	{
		for(int j=0; j<dfa1->get_dim_alphabet(); ++j)
			if(dfa1->get_ttable()[i][j] == blustate){
				predecessore = i;
				lettera = j;
				num_predecessori++;
			}
	}

	if(num_predecessori != 1)
		cerr << "Num predec:"<< num_predecessori<<" PROBLEMA: "<< num_predecessori << "per lo stato "<<blustate<<", controlla algoritmo!"<<endl;

	// Setto la transizione del precedessore (di bluestate) a redstate,
	// adesso q2 è IRRAGIUNGIBILE nell'albero originale
	dfa1->get_ttable()[predecessore][lettera] = redstate;

	fold(dfa1, redstate, blustate);
}


void gi::bluefringe::fold(dfaEDSM* originale, int redstate, int blustate)
{
	int colonna_tipo = originale->get_dim_alphabet();
	int** current_ttable = originale->get_ttable();

	// Se q2 è accettante, setto ad altrettanto q1
	if(current_ttable[blustate][colonna_tipo] == DFA_STATE_ACCEPTING)
		current_ttable[redstate][colonna_tipo] = DFA_STATE_ACCEPTING;

	// Per ogni lettera, effettuo il controllo se nell'albero originale è definita
	// una transizione presente anche nel subtree:
	// se presente passo agli stati successivi, diversamente effettuo il merge inserendo la transizione.
	for(int i=0; i<originale->get_dim_alphabet(); ++i)
	{
		int statefrom_bluestate = current_ttable[blustate][i];
		int statefrom_redstate = current_ttable[redstate][i];

		// Ormai posso eliminare lo stato blue, lo faccio settando a ND la sua transizione per la lettera corrente
		current_ttable[blustate][i] = ND;

		// TODO: Potrei già eliminare qui lo stato dai blue? penso di si, ma fai debug prima di attivarlo
		//current_ttable[blustate][colonna_tipo+1] = DFA_STATE_WHITE;
		//originale->remove_blue_state(blustate);

		if(statefrom_bluestate != ND){
			if(statefrom_redstate != ND)
				fold(originale, statefrom_redstate, statefrom_bluestate);
			else{
				// Aggiungo la transizione esistente nel subtree ma non nello stato in cui faccio il merge
				current_ttable[redstate][i] = statefrom_bluestate;
				//current_ttable[blustate][i] = ND;
			}
		}
	}
	// Non c'è bisogno che elimini la riga dello stato, sarà irragiungibile da tutti, puoi tenerne conto alla fine.
}


//int gi::bluefringe::edsm_count(dfaEDSM* dfa1, vector<SYMBOL>* positive, int dim_positive, vector<SYMBOL>* negative, int dim_negative)
//{
//	int* tp, *tn = NULL;
//	tp = new int[dfa1->get_num_state()];							//Inizializzo un contatore per ogni stato
//	tn = new int[dfa1->get_num_state()];
//
//	for(int i=0; i<dfa1->get_num_state(); ++i){
//		tp[i] = 0;
//		tn[i] = 0;
//	}
//
//
//	// TODO Implementare COUNTER di OpenMP
//	// Controllo quante stringhe POSITIVE vengono riconosciute
//	for(int i=0; i<dim_positive; ++i){
//		int statoFinale = dfa1->get_arrive_state(positive[i]);
//		if(statoFinale != ND)
//			tp[statoFinale]++;
//	}
//
//	//cout << endl;
//	for(int i=0; i<dfa1->get_num_state(); ++i)
//		if(db)
//			cout << "tp["<<i<<"]: "<<tp[i]<<" ";
//
//	// TODO Implementare COUNTER di OpenMP
//	// Controllo quante stringhe NEGATIVE vengono riconosciute
//	for(int i=0; i<dim_negative; ++i){
//		int statoFinale = dfa1->get_arrive_state(negative[i]);
//		if(statoFinale != ND)
//			tn[statoFinale]++;
//	}
//
//	//cout << endl;
//	for(int i=0; i<dfa1->get_num_state(); ++i)
//		if(db)
//			cout << "tn["<<i<<"]: "<<tn[i]<<" ";
//
//	// Calcolo lo SCORE
//	int sc = 0;
//	for(int i=0; i<dfa1->get_num_state(); ++i)
//	{
//		if(sc != MINF)
//		{
//			if(tn[i] > 0){
//				if(tp[i] > 0){
//					if(db)
//						cout <<endl<< "Stato che accetta Positiva&Negativa:"<<i<<endl;
//					sc = MINF;
//					break;								// E' inutile continuare a ciclare se metto sc = ND
//				}
//				else
//					sc = sc + tn[i] -1;					// Nello pseudo De La Higuera mette -1 (in situazioni complesse porta ad automi finali diversi)
//			}else{
//				if(tp[i] > 0)
//					sc = sc + tp[i] -1;
//			}
//		}
//	}
//
//	if(tp!=NULL)
//		delete[] tp;
//	if(tn!=NULL)
//		delete[] tn;
//
//	return sc;
//}

void gi::bluefringe::promote(dfaEDSM* dfa1, int q)
{
	if(!dfa1->is_inside_blue_states(q))
		cerr << "ATTENZIONE! Lo stato non è BLU!"<< endl;

	// Promuovo a RED lo stato q
	dfa1->add_red_state(q);

	dfa1->remove_blue_state(q);

	// Promuovo a BLU tutti gli stati raggiungibili direttamente da q
	for(int i=0; i<dfa1->get_dim_alphabet(); ++i)
	{
		int transizione = dfa1->get_ttable()[q][i];
		if(transizione != ND)
			dfa1->add_blue_state(transizione);
	}
}



// Promuovo a BLU tutti gli stati raggiungibili in modo diretto da un RED
void gi::bluefringe::nuoviBlu(dfaEDSM* dfa1)
{
	int numred = dfa1->get_num_red_states();

	// Promuovo a BLUE tutti gli stati raggiungibili da un RED in modo diretto
	for(int i=0; i<dfa1->get_dim_alphabet(); ++i)
	{
		for(int j=0; j<numred; ++j)
		{
			int red = dfa1->get_red_states()->at(j);
			int transizione = dfa1->get_ttable()[red][i];
			if(transizione != ND){
				//Pongo lo stato a blu, a meno che non sia rosso
				if(!dfa1->is_inside_red_states(transizione)){
					dfa1->add_blue_state(transizione);

				}
			}
		}
	}
}

void gi::bluefringe::eliminaStati(dfaEDSM* dfa1)
{
	int colonna_tipo = dfa1->get_dim_alphabet();

	// Setto gli stati Eliminati
	for(int k=1; k<dfa1->get_num_states(); ++k)				// Verifico se lo stato K è visitato da qualcuno
	{
		bool visitato = false;
		for(int i=0; i<dfa1->get_num_states(); ++i)
		{
			for(int j=0; j<dfa1->get_dim_alphabet(); ++j)
				if(dfa1->get_ttable()[i][j] == k && dfa1->get_ttable()[i][colonna_tipo] != DFA_STATE_UNREACHABLE)
					visitato = true;
		}

		if(!visitato){
			dfa1->get_ttable()[k][colonna_tipo] = DFA_STATE_UNREACHABLE;		//Setto a Eliminato il tipo dello stato non visitato da alcun altro stato

			// Li elimino anche dal Vector dei Red o dei Blue, se presenti
			dfa1->remove_red_state(k);
			dfa1->remove_blue_state(k);
		}
	}

}


int gi::bluefringe::get_actual_merge(){
	return num_actual_merge;
}

int gi::bluefringe::get_heuristic_merge(){
	return num_heuristic_merge_valued;
}

void gi::bluefringe::set_fringe_size(int r, int b){
	fringe_size[0].push_back(b);		// 0 for Blue, 1 for Red
	fringe_size[1].push_back(r);

}

void gi::bluefringe::print_fringe_size(){

	vector<int>::const_iterator i;
	vector<int>::const_iterator l;

	i=fringe_size[1].begin();			// 0 for Blue, 1 for Red
	l=fringe_size[0].begin();

	while(i!=fringe_size[1].end()){

		cout<<"R:"<<(*i)<<",";
		cout<<"B:"<<(*l) << endl;

		i++;
		l++;
	}

	cout <<"-------"<<endl;
}

int gi::bluefringe::get_while_count(){
	return while_count;
}



//gi::dfa gi::bluefringe::run_edsm(string base_path)
//{
//	// Samples from txtfile
//	int n_symbols=0;	//number of negative examples
//	int dim_positive=0; //number of positive examples
//	int dim_negative=0; //number of negative examples
//
//	int max_count, *curr_count=NULL;
//
//	int	n_red=0, n_blue=0, actual_blue_state=0;
//
//	// example strings
//	vector<SYMBOL>* positive=NULL;
//	vector<SYMBOL>* negative=NULL;
//	char* symbols = NULL;
//
//	//bool *promotion=NULL;
//	bool promoted =false;
//
//	dfaEDSM *dfa1 =NULL, /* *dfa_best = NULL,*/ **merged = NULL; // dfa...
//
//	// One dfa_best for every Blue state
//	vector<dfaEDSM*> dfa_best;
//	vector<int> dfa_score;
//
//
//	//get positive and negative samples
//	read_samples(positive, &dim_positive, negative, &dim_negative);
//
//	// Costruisco PTA
//	//dfaEDSM* dfa1 = build_pta(positive,dim_positive);
//
//	// Build APTA
//	dfa1 = build_apta(positive, dim_positive, negative, dim_negative);
//
//	// Print it!
//	if(dfa1->get_num_state() < 1000)
//	{
//		dfa1->print_dfa_dot("APTA", (base_path+"APTA.dot").c_str() );
//		dfa1->print_dfa_dot_mapped_alphabet("APTAALF", (base_path+"APTA_ALF.dot").c_str());
//	}else{
//		clog<<"APTA too big! I can't print it"<<endl;
//	}
//
//	n_blue = dfa1->get_num_blue_states();
//	n_red = dfa1->get_num_red_states();
//
//	set_fringe_size(n_red,n_blue);
//
//	cout <<" START Edsm inference process..."<<endl;
//
//	while_count=-1;
//	// ESDM
//	while(n_blue>0)
//	{
//		while_count++;
//
//		promoted=false;
//
//		// BLUE ciclo
//		for(int i=0; i<n_blue && (!promoted); ++i)
//		{
//			actual_blue_state = dfa1->get_blue_states()->at(i);
//			// Reset variable for the new run
//
//			// array for the heuristic values of the red group
//			if(curr_count != NULL)
//				delete[] curr_count;
//			curr_count= new int [n_red];
//
//			// dfa coming from possible merges
//			if(merged != NULL)
//				delete[] merged;
//			merged = new dfaEDSM*[n_red];
//
//			// initialize values
//			for(int j=0; j<n_red; ++j){
//				curr_count[j] = MINF;
//				merged[j] = NULL;
//			}
//
//
//			// RED ciclo
//			#pragma omp parallel default(shared)
//			{
//			#pragma omp for
//			for(int j=0; j<n_red; ++j){
//				merged[j] = new dfaEDSM(*dfa1);
//
//				merge(merged[j], dfa1->get_red_states()->at(j), actual_blue_state );
//
//				// TODO: Questa riga si può probabilmente eliminare, da fare debug estensivo
//				merged[j]->remove_blue_state(actual_blue_state);
//
//				curr_count[j] = edsm_count(merged[j], positive, dim_positive, negative, dim_negative);
//			}
//			}
//			// end for RED
//
//			// For Statistical purpose
//			num_heuristic_merge_valued +=  n_red;
//
//			// check if there some merge, else start process for promote
//			promoted = true;
//			max_count=MINF;
//			int j_max=ND;
//			for(int j=0; j<n_red; ++j){
//				if(curr_count[j]>max_count){
//					max_count = curr_count[j];
//					j_max = j;
//					promoted=false;
//				}
//			}
//
//			//cout << "Max_count:"<<max_count<<endl;
//
//
//			// PROMOTION
//			if(promoted){
//
//				edsm_promote(dfa1, actual_blue_state);
//
//				//cout << "PROMOZIONE"<<endl;
//
//				#ifdef ALL_DFA_EDSM
//				string name = "P"+dfa::intTostring(while_count)+dfa::intTostring(blue[i]);
//				dfa1->print_dfa_dot(name, (base_path+name+".dot").c_str());
//				#endif
//
//
//				//Free memory
//				typedef	vector<dfaEDSM*>::const_iterator It;
//				for(It p1=dfa_best.begin(); p1!=dfa_best.end(); ++p1){
//					if(dfa1 == (*p1))
//						continue;
//					delete (*p1);
//				}
//
//				dfa_best.clear();
//				dfa_score.clear();
//			}
//			else	// - Merge accettato come candidato per il merge finale. Lo aggiungo alla lista dei migliori.
//			{
//				//merged[j_max]->remove_blue_state(actual_blue_state);
//
//				dfa_best.push_back(merged[j_max]);
//				dfa_score.push_back(max_count);
//			}
//
//
//			// Free the array with dfa merged for calculate score, leave only the dfa selected as best
//			if(merged != NULL){
//				for(int j=0; j<n_red; ++j){
//					if (j == j_max && (!promoted))			// Leave the best
//						continue;
//					if(merged[j] != NULL)
//						delete merged[j];
//				}
//				delete[] merged;
//				merged = NULL;
//			}
//		}// end for BLUE
//
//
//		// MERGE
//		if(!promoted){ // Do definitive merge, no promotion done. Select best merge between all candidates in "dfa_best"
//
//			// Select the best merge between all the blue states
//			int best_score = -1;
//			int index_best_dfa = 0;
//			for(int t=0; t<dfa_score.size(); ++t)
//				if(dfa_score.at(t) > best_score){
//					best_score = dfa_score.at(t);
//					index_best_dfa = t;
//				}
//
//			// Take the blue states before delete the old dfa
//			//int colonna_tipo = dfa1->get_dim_alphabet();
//			for(int t=0; t<dfa1->get_num_blue_states(); ++t)
//				dfa_best.at(index_best_dfa)->add_blue_state( dfa1->get_blue_states()->at(t) );
//			if(dfa1 != NULL) delete dfa1;		// Delete old dfa
//
//			// set dfa1 to the new merged dfa
//			dfa1 = dfa_best.at(index_best_dfa);
//			nuoviBlu(dfa1);
//			eliminaStati(dfa1);
//
//
//			// Print information
//			//cout << "MERGE:"<<max_count<<endl;
//			//cout <<" ----------------------------------- "<<endl;
//			#ifdef ALL_DFA_EDSM
//			string name = "M"+dfa::intTostring(while_count);
//			dfa1->print_dfa_dot(name, (base_path+name+".dot").c_str());
//			#endif
//
//			++num_actual_merge;
//
//			// Free memory
//			typedef	vector<dfaEDSM*>::const_iterator It;
//			for(It p1=dfa_best.begin(); p1!=dfa_best.end(); ++p1){
//				if(dfa1 == (*p1))
//					continue;
//
//				if((*p1) != NULL)
//					delete (*p1);
//			}
//
//			dfa_best.clear();
//			dfa_score.clear();
//		}
//
//		// update values for the dfa
//		n_blue = dfa1->get_num_blue_states();
//		n_red = dfa1->get_num_red_states();
//
//		set_fringe_size(n_red,n_blue);
//
//	}
//
//
//	if(curr_count != NULL) delete[] curr_count;
//
//
//	// Setto gli stati Accettanti
//	int colonna_tipo = dfa1->get_dim_alphabet();
//	for(int i=0; i<dim_positive; ++i){
//		int accettante = dfa1->get_arrive_state(positive[i]);
//
////			if(accettante != ND)
////				dfa1->get_ttable()[accettante][colonna_tipo] = DFA_STATE_ACCEPTING;
//	}
//
//	// Setto gli stati Rigettanti
//	for(int i=0; i<dim_negative; ++i){
//		int rigettante = dfa1->get_arrive_state(negative[i]);
//		if(rigettante != ND){
//			//cout << "Statp di arrivoN: "<<rigettante<<endl;
//			dfa1->get_ttable()[rigettante][colonna_tipo] = DFA_STATE_REJECTING;
//		}
//	}
//
//	// Setto gli stati Eliminati
//	eliminaStati(dfa1);
//	//dfa1->print_dfa_with_color("AUTOMA FINALE");
//
//
//	//////////////////////////////////////////////////////////////
//	// - Pulisco l'automo dagli stati irragiungibili, aggiorno le transizioni -
//
//	// Conto il numero effettivo di stati finali
//	int colour_column = dfa1->get_dim_alphabet()+1;
//	int n_final_states = 0;
//	for(int i=0; i<dfa1->get_num_state(); ++i)
//		if(dfa1->is_inside_red_states(i))
//			n_final_states++;
//
//
//	// Creo un nuovo automa senza stati irragiungibili
//	int count = 0;
//	dfaEDSM* finalDFA = new dfaEDSM(n_final_states, dfa1->get_dim_alphabet(), dfa1->get_alphabet());
//	map<int,int> updated_transition;
//
//	for(int i=0; i<dfa1->get_num_state(); ++i){
//		if(dfa1->is_inside_red_states(i)){
//
//			for(int j=0; j<dfa1->get_dim_alphabet()+1; ++j){
//
//				// Aggiungo lo stato al nuovo automa
//				finalDFA->get_ttable()[count][j] = dfa1->get_ttable()[i][j];
//
//				updated_transition[i] = count;
//			}
//			++count;
//		}
//	}
//	updated_transition[ND] = ND;
//
//	if(dfa1!=NULL)
//		delete dfa1;
//	dfa1 =NULL;
//
//	bool stato_pozzo = false;
//	// Aggiorno le transizioni
//	for(int i=0; i<finalDFA->get_num_state(); ++i)
//		for(int j=0; j<finalDFA->get_dim_alphabet(); ++j){
//			if(finalDFA->get_ttable()[i][j] == ND)							// Rilevo che c'è una transizione mancante, quindi serve uno stato pozzo
//				stato_pozzo = true;
//			if(updated_transition.find(finalDFA->get_ttable()[i][j]) != updated_transition.end())
//				finalDFA->set_ttable_entry(i, j, updated_transition[ finalDFA->get_ttable()[i][j] ]);
//			else{
//				cerr << "Errore nell'aggiornamento delle stringhe"<<endl;
//				exit(EXIT_FAILURE);
//			}
//		}
//
//
//	// Stampo l'automa prima di applicare il pozzo e la minimizzazione
//	//finalDFA->print_dfa_with_color("AUTOMA FINALE PREPOZZO");
//
//	//finalDFA->print_dfa_dot("FINALEPREPOZZO", percorso.c_str());
//
//	finalDFA->print_dfa_dot_mapped_alphabet("FINALE_PREPOZZO", (base_path+"pulito_pre_pozzo.dot").c_str());
//
//
//	//////////////////////////////////////////////////////////////
//	// Controllo stato pozzo
//	// - Se ci sono transizioni non definite le imposto tutte verso lo stato pozzo
//	if(stato_pozzo)
//	{
//		dfaEDSM* finalDFAPozzo = new dfaEDSM(finalDFA->get_num_state()+1, finalDFA->get_dim_alphabet(), finalDFA->get_alphabet(), 0);
//		int** table = finalDFAPozzo->get_ttable();
//
//		for(int i=0; i<finalDFA->get_num_state(); ++i)
//			for(int j=0; j<finalDFA->get_dim_alphabet()+1; ++j){
//				if(finalDFA->get_ttable()[i][j] == ND)
//					table[i][j] = finalDFA->get_num_state();
//				else
//					table[i][j] = finalDFA->get_ttable()[i][j];
//			}
//
//		for(int j=0; j<finalDFA->get_dim_alphabet(); ++j)
//			table[finalDFA->get_num_state()][j] = finalDFA->get_num_state();
//
//		delete finalDFA;
//		finalDFA = finalDFAPozzo;
//	}
//
//	//finalDFA->print_dfa_with_color("EDSM PRE-MINIMIZZAZIONE AUTOMA FINALE");
//
//
//	//////////////////////////////////////////////////////////////
//	// Minimize returna a new dfa, then delete the older
//	dfa* finalDFA2 = finalDFA->minimize_TF();
//
//	if(finalDFA) delete finalDFA;
//
//	if(positive) delete[] positive;
//	if(negative) delete[] negative;
//	if(symbols) delete[] symbols;
//
//	return *finalDFA2;
//
//}
//
