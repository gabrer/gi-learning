/*
 * DFA IMPLEMENTATION
 */

#include <utilities.h>
#include "dfa.h"
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>	 //isstringstream
#include <map>
#include <cmath>
#include <algorithm>
#include <limits>



using namespace std;

#define BUFFER_SIZE 256

#define DFA_TF_STATE_N numeric_limits<SYMBOL>::max()		// For Table Filling Algorithm
#define DFA_TF_STATE_X numeric_limits<SYMBOL>::max() -1


gi::dfa::dfa(){
	dim_alphabet = 0;
	num_states 	 = 0;
	start_state  = 0;

	ttable = NULL;
	alphabet=NULL;
}


gi::dfa::dfa(const int n_state, const int dim_alf, const char *alf, const int s_state){
	dim_alphabet = dim_alf;
	num_states 	 = n_state;
	start_state  = s_state;

	ttable = new int*[num_states];
	for(int i=0; i<num_states; ++i)
		ttable[i] = new int[dim_alphabet+1];				// "+1" for Type

	for(int i=0; i<num_states; ++i)
		for(int j=0; j<dim_alphabet+1; ++j){
			if(j >= dim_alphabet){
				ttable[i][j] = 0;
				continue;
			}
			ttable[i][j]=ND;
		}


	// Alphabet symbols
	alphabet=NULL;
	set_alphabet(alf, dim_alf);
}

gi::dfa::dfa(const int n_state, const int dim_alf, const char *alf)
:dfa(n_state, dim_alf, alf, 0){}

gi::dfa::dfa(const int n_state, const int dim_alf, const char *alf, const int s_state, const int** tt_copy )
:dfa(n_state, dim_alf, alf, 0){
	for(int i=0; i<n_state; ++i)
		for(int j=0; j<dim_alf+1; ++j)
			ttable[i][j] = tt_copy[i][j];//d1.get(i,j);
}


// Constructor for copy
gi::dfa::dfa(const dfa &d1)
:dfa(d1.num_states, d1.dim_alphabet, d1.alphabet, d1.start_state, (const int**) d1.ttable){}



gi::dfa::~dfa(){
	if(ttable != NULL){
		for(auto i=0; i<num_states; ++i)
			if(ttable[i] != NULL)
				delete[] ttable[i];

		delete [] ttable;
	}

	if(alphabet != NULL){
		delete[] alphabet;
	}

}



gi::dfa* gi::dfa::unionDFA(dfa* dfa_hp)
{
	// DFA UNIONE di dfa_hp e dfa_T
	int count_state = dfa_hp->num_states + num_states;
	//V cout << "Numero di stati totali: "<<dfa_hp->num_state<<" + "<<num_state<<endl;

	dfa* dfa_unione = new dfa(count_state, dim_alphabet, alphabet);

	// Riempo il DFA UNIONE
	// I primi stati sono del target, i restanti dell'hp
	for(int j=0; j<num_states; ++j)												// Automa target
		for(int k=0; k<dim_alphabet+1; ++k)
			dfa_unione->ttable[j][k] = ttable[j][k];

	for(int j=0; j<dfa_hp->num_states; ++j){										// Automa ipotesi
		for(int k=0; k<dim_alphabet+1; ++k){									// Start state dell DFA_hp nel DFA unione pari a num_state del target
			if(k != dim_alphabet)
				dfa_unione->ttable[num_states+j][k] = dfa_hp->ttable[j][k] + num_states;
			else
				dfa_unione->ttable[num_states+j][k] = dfa_hp->ttable[j][k];
		}
	}

	// Stampo DFA unione
	// dfa_unione->print_dfa_ttable("*** DFA UNIONE ***");

	return dfa_unione;
}






// Ritorna l'elenco degli stati equivalenti, se necessario inserisce in "witness" la stringa controesempio
SYMBOL* gi::dfa::table_filling(){
	// The Table considered is only the upper triangular matrix, that can be saved in a linear array of size n(n-1)/2
	// Conversion of index form matrix to array are:
	//
	// From linear index k, to (i,j) for tha matrix (i row, j column)
	// i = n - 2 - floor(sqrt(-8*k + 4*n*(n-1)-7)/2.0 - 0.5)
	// j = k + i + 1 - n*(n-1)/2 + (n-i)*((n-i)-1)/2
	//
	// From (i,j) to k
	// Order such that i<j, then:
	// k = (n*(n-1)/2) - (n-i)*((n-i)-1)/2 + j - i - 1
	//
	// check (http://stackoverflow.com/questions/27086195/linear-index-upper-triangular-matrix)

	// *** TABLE-FILLING ALGORITHM with witness ***
	int i,j,k;
	int n = num_states;
	int tf_l = (num_states*(num_states-1))/2;


	// Tabella delle coppie di stati distinti
	SYMBOL* distinti = new SYMBOL[tf_l];

	// Inizializzo
	//for(k=0; k<tf_l; ++k)
		//distinti[k]=  tableFillingState::TF_STATE_N;

	// Distinguo tra accettanti e non accettanti
	for(i=0; i<(num_states-1); ++i)
		for(j=i+1; j<num_states; ++j){
			// Se uno è accettante mentre l'altro no
			k= (n*(n-1)/2) - (n-i)*((n-i)-1)/2 + j - i - 1;
			if(ttable[i][dim_alphabet] != ttable[j][dim_alphabet]){
				distinti[k] = DFA_TF_STATE_X;
			}else
				distinti[k]=  DFA_TF_STATE_N;
		}


	// Loop  minimizzante
	// Tiene traccia se la tabella è stata modificata durante l'ultimo ciclo
	bool modificata = true;
	while(modificata)
	{
		modificata = false;
		int stato_arrivo_1;
		int stato_arrivo_2;

		for(i=0; i<(num_states-1); ++i){
			for(j=i+1; j<num_states; ++j){
				k= (n*(n-1)/2) - (n-i)*((n-i)-1)/2 + j - i - 1;
				//V cout << endl << "SP1: "<<i << ", SP2: "<<j;
				if(distinti[k] == DFA_TF_STATE_N){

					for(SYMBOL w=0; w<dim_alphabet; ++w)
					{
						stato_arrivo_1 = ttable[i][w];
						stato_arrivo_2 = ttable[j][w];

						// Lo stato di arrivo letto nel dfa potrebbe essere una coppia del tipo (2,2)
						if(stato_arrivo_1 == stato_arrivo_2)
							continue;

						// Facendo la tabella delle coppie di stati visitati noto che sempre j>i
						// quindi le coppie devono avere sempre i<j
						if(stato_arrivo_2 < stato_arrivo_1){
							int tmp = stato_arrivo_1;
							stato_arrivo_1 = stato_arrivo_2;
							stato_arrivo_2 = tmp;
						}

						// Se finisco su me stesso vado avanti
						if(stato_arrivo_1 == i && stato_arrivo_2 == j)
							continue;

						//V cout << endl <<"SA1: "<<stato_arrivo_1 << ", SA2: "<<stato_arrivo_2 << " --> " << (int)distinti[stato_arrivo_1][stato_arrivo_2] << endl;
						// Se la coppia di arrivo è distinta, lo è anche quella originaria
						int i1 = stato_arrivo_1, j1 = stato_arrivo_2;
						int k1 = (n*(n-1)/2) - (n-i1)*((n-i1)-1)/2 + j1 - i1 - 1;

						if(distinti[k1] != DFA_TF_STATE_N){
							//TODO: è corretto?
							//distinti[k1] = w + '0'; //scrive l'indice (in char) della lettera dell'alfabeto per cui i due stati differiscono
							//distinti[k] = tableFillingState::TF_STATE_O;
							distinti[k] = w;
							modificata = true;
							break;											// Necessario!
						}
					}
				}
			}
		}
	}

	// Stampo la tabella degli stati distinti
	#ifdef DEBUG_2
	cout << "--------------------------" << endl;
	cout << "Tabella dei distinti " << endl;
	cout << "--------------------------" << endl;
	for(i=0; i<(num_states-1); ++i){
		for(j=i+1; j<num_states; ++j){
			k=(n*(n-1)/2) - (n-i)*((n-i)-1)/2 + j - i - 1;
			char toprint= (distinti[k]==DFA_TF_STATE_X)?'X':(distinti[k]==DFA_TF_STATE_N)?'@': (char)(distinti[k] + 48);
			cout << "("<< i << "," << j << "):" << toprint << "  ";
		}
		cout << endl;
	}
	cout << "--------------------------" << endl;
	#endif

	return distinti;
}


gi::dfa* gi::dfa::minimize_TF()
{
	const int num_row = num_states;

	// Matrix distinct: table for record distinct states
	bool** distinct = new bool*[num_states-1];								// Tabella delle coppie di stati distinti
	for(int i=0; i<num_row-1; ++i)											// EX:per 6 stati:5 righe e 6 colonne (di queste 6 solo 5 usate effettivamente, ma serve)
		distinct[i] = new bool[num_states];

	// Initialization of distinct table
	for(int i=0; i<num_states-1; ++i)
		for(int j=0; j<num_states; ++j)
			distinct[i][j]=false;

	// Distinguo tra accettanti e non accettanti							// EX: 0 <= i <= 4, j=i+1 (1 <= j <= 5 per la prima iterazione)
	for(int i=0; i<(num_states-1); ++i)
		for(int j=i+1; j<num_states; ++j)
			if(ttable[i][dim_alphabet] != ttable[j][dim_alphabet])				// Se uno è accettante mentre l'altro no
				distinct[i][j] = true;

	// Loop  minimizzante
	bool modificata = true;
	while(modificata)
	{
		modificata = false;

		for(int i=0; i<(num_states-1); ++i){
			for(int j=i+1; j<num_states; ++j){
				if(!distinct[i][j]){

					for(int k=0; k<dim_alphabet; ++k)
					{
						int stato_arrivo_1 = ttable[i][k];
						int stato_arrivo_2 = ttable[j][k];

						if(stato_arrivo_1 == stato_arrivo_2)
							continue;

						if(stato_arrivo_2 < stato_arrivo_1){				// Facendo la tabella delle coppie di stati visitati noto che sempre j>i
							int tmp = stato_arrivo_1;
							stato_arrivo_1 = stato_arrivo_2;
							stato_arrivo_2 = tmp;
						}

						if(distinct[stato_arrivo_1][stato_arrivo_2]){		// Se la coppia di arrivo è distinta, lo è anche quella originaria
							distinct[i][j] = true;
							modificata = true;
						}
					}
				}
			}
		}
	}


	// ** Mi preparo a costruire l'automa minimizzato **
	// Creo un elenco per avere una lista delle coppie di stati equivalenti più direttamente accessibile,
	// 	se è ND allora non ha un equivalente.
	int stato_equivalente[num_states];
	for(int i=0; i<num_states; ++i)
		stato_equivalente[i] = ND;

	// Conto il numero di stati finali, ed assegno ad ogni stato il suo eventuale stato equivalente
	// (devo anche controllare nel caso ho (0,4) (0,5) che quando incontro (4,5) non lo conto come un ulteriore stato in meno)
	int stati_finali = num_states;
	for(int i=0; i<(num_states-1); ++i)
		for(int j=i+1; j<num_states; ++j)
			if(!distinct[i][j] && stato_equivalente[i] == ND && stato_equivalente[j] == ND){
				stato_equivalente[j] = i;
				stati_finali--;
			}

	// Stampo le informazioni riguardo agli stati equivalenti trovati o meno
	/*//V cout << "N di stati finali: " << stati_finali << endl;
	cout << "Equivalenze stati: " << endl;
	for(int i=0; i<num_state; ++i)
		if(stato_equivalente[i] != ND)
			cout << "S:"<<i<<" --> "<<stato_equivalente[i]<<endl;
		else
			cout << "S:"<<i<<endl
	*/

	// STAMPO il vecchio DFA
	#ifdef DEBUG_2
	this->print_dfa_ttable("DFA BEFORE MINIMIZATION");
	#endif

	// Istanzio il nuovo DFA minimizzato
	dfa* dfa_min = new dfa(stati_finali,  dim_alphabet, alphabet, 0);

	int** ttable_min = dfa_min->get_ttable();

	int count = 0;
	for(int i=0; i<num_states; ++i){
		if(stato_equivalente[i] == ND){
			for(int j=0; j<dim_alphabet+1; ++j)
				ttable_min[count][j]=ttable[i][j];
			count++;
		}
	}

	// Aggiorno le transizioni verso stati ormai scomparsi perchè sostituiti con gli equivalenti
	int equivalenze_finora= 0;
	for(int i=0; i<num_states; ++i){
		if(stato_equivalente[i] != ND)
		{
			for(int k=0; k<stati_finali; ++k)
				for(int t=0; t<dim_alphabet; ++t)
					// Sostituisco la transizione allo stato "i" con lo stato equivalente "stato_equivalente[i]"
					if(ttable_min[k][t] == i)
						ttable_min[k][t] = stato_equivalente[i];
		}
	}

	// Qui aggiorno la label, perché se ad esempio ho collassato 2 stati prima dello stato 6,
	//  adesso lo stato 6 si trova nella riga 4, però le transizioni sono rimaste
	//  verso 6 e devono essere sistemate
	for(int i=0; i<num_states; ++i)
	{
		if(stato_equivalente[i] != ND)
			equivalenze_finora++;

		if(equivalenze_finora != 0)
		{
			int nuova_label = i-equivalenze_finora;
			for(int k=0; k<stati_finali; ++k)
				for(int t=0; t<dim_alphabet; ++t)
					if(ttable_min[k][t] == i)
						ttable_min[k][t] = nuova_label;
		}
	}

	// STAMPO il dfa minimizzato
	#ifdef DEBUG_2
	dfa_min->print_dfa_ttable("MINIMIZED DFA");
	#endif

	//libero la memoria allocata
	if(distinct != NULL){
		for(int i=0; i<num_row -1; ++i){
			if(distinct[i] != NULL){ 
				delete[] distinct[i];
			}
		}
		delete [] distinct;
	}

	return dfa_min;
}


void gi::dfa::print_dfa_ttable(string title)
{
	// STAMPA IL DFA
	// Nel numero di colonne è inclusa la colonna finale di stato accettante o meno (accettante -> 1)

	cout << endl<< "--------------------------" << endl;
	cout << title << endl;
	string header = "  ";
	for(int i=0; i<dim_alphabet; ++i)
		header = header + " | "+ intTostring(i);
	header = header + " - A";
	cout << header << endl;

	for(int i=0; i<num_states; ++i){
		cout << "S"<<i<<"  ";

		for(int j=0; j<dim_alphabet+1; ++j)
		{

			if(j < dim_alphabet && ttable[i][j] == ND)			// Valori delle transizioni, o ND o il valore
				cout << " N ";
			else if(j < dim_alphabet)
				cout << " "<< ttable[i][j] <<"  ";

			else if(j == dim_alphabet)								// Tipo dello stato: accettante o meno
			{
				if(ttable[i][j] == DFA_STATE_NON_ACCEPTING )
					cout << "  / ";
				else if(ttable[i][j] == DFA_STATE_ACCEPTING)
					cout << " Ac ";
				else if(ttable[i][j] == DFA_STATE_REJECTING)
					cout << " Ri ";
				else
					cout << "  X ";
			}
		}
		cout << endl;
	}

	cout << "--------------------------" << endl;

}


void gi::dfa::print_dfa_ttable_mapped_alphabet(string title)
{
	// STAMPA IL DFA
	// Nel numero di colonne è inclusa la colonna finale di stato accettante o meno (accettante -> 1)
	// Uso l'alfabeto originale


	// Using Mapped Alphabet

	/*typedef	map<char, char>::const_iterator It;
	cout << "Mappatura inversa correta? "<<endl;
	for(It p1=inverse_mapped_alphabet.begin(); p1!=inverse_mapped_alphabet.end(); ++p1)
		cout << (*p1).first << "; " << (*p1).second << endl;
	cout << "FINE"<<endl;*/

	cout << endl<< "--------------------------" << endl;
	cout << title << endl;
	string header = "  ";
	for(int i=0; i<dim_alphabet; ++i)
		header = header + " | "+ alphabet[i];
	header = header + " - A";
	cout << header << endl;

	for(int i=0; i<num_states; ++i){
		cout << "S"<<i<<"  ";

		for(int j=0; j<dim_alphabet+1; ++j)
		{

			if(j < dim_alphabet && ttable[i][j] == ND)			// Valori delle transizioni, o ND o il valore
				cout << " N ";
			else if(j < dim_alphabet)
				cout << " "<< alphabet[ttable[i][j]] <<"  ";

			else if(j == dim_alphabet)								// Tipo dello stato: accettante o meno
			{
				if(ttable[i][j] == DFA_STATE_NON_ACCEPTING )
					cout << "  / ";
				else if(ttable[i][j] == DFA_STATE_ACCEPTING)
					cout << " Ac ";
				else if(ttable[i][j] == DFA_STATE_REJECTING)
					cout << " Ri ";
				else
					cout << "  X ";
			}
		}
		cout << endl;
	}

	cout << "--------------------------" << endl;

}




//void gi::dfa::print_dfa_with_color(string title)
//{
//	// STAMPA IL DFA
//	// Nel numero di colonne è inclusa la colonna finale di stato accettante o meno (accettante -> 1)
//
//	cout << endl<< "--------------------------" << endl;
//	cout << title << endl;
//	string header = "    ";
//		for(int i=0; i<dim_alphabet; ++i)
//			header = header + " | "+ intTostring(i);
//		header = header + " - A  - C";
//	cout << header << endl;
//
//	for(int i=0; i<num_state; ++i){
//		cout << "S"<<i<<"  ";
//
//		for(int j=0; j<dim_alphabet+1; ++j)
//		{
//
//			if(j < dim_alphabet && ttable[i][j] == ND)				// Valori delle transizioni, o ND o il valore
//				cout << " N ";
//
//			else if(j < dim_alphabet)
//				cout << " "<< ttable[i][j] <<" ";
//
//			else if(j == dim_alphabet)								// Tipo dello stato: accettante o meno
//			{
//				if(ttable[i][j] == DFA_STATE_NON_ACCEPTING )
//					cout << "  / ";
//				else if(ttable[i][j] == DFA_STATE_ACCEPTING)
//					cout << " Ac ";
//				else if(ttable[i][j] == DFA_STATE_REJECTING)
//					cout << " Ri ";
//				else
//					cout << "  X ";
//			}
//		}
//
////		if(!this->is_inside_blue_states(i) && !this->is_inside_red_states(i))
////			cout << " W ";
////		else if(this->is_inside_blue_states(i))
////			cout << " B";
////		else
////			cout << " R";
//
//		cout << endl;
//	}
//
//	cout << "--------------------------" << endl;
//}
//
//void gi::dfa::print_dfa_with_color_mapped_alphabet(string title)
//{
//	// STAMPA IL DFA
//	// Nel numero di colonne è inclusa la colonna finale di stato accettante o meno (accettante -> 1)
//
//	cout << endl<< "--------------------------" << endl;
//	cout << title << endl;
//	string header = "    ";
//		for(int i=0; i<dim_alphabet; ++i)
//			header = header + " | "+ alphabet[i];
//		header = header + " - A  - C";
//	cout << header << endl;
//
//	for(int i=0; i<num_state; ++i){
//		cout << "S"<<i<<"  ";
//
//		for(int j=0; j<dim_alphabet+1; ++j)
//		{
//
//			if(j < dim_alphabet && ttable[i][j] == ND)				// Valori delle transizioni, o ND o il valore
//				cout << " N ";
//			else if(j < dim_alphabet)
//				cout << " "<< alphabet[ttable[i][j]] <<" ";
//
//			else if(j == dim_alphabet)								// Tipo dello stato: accettante o meno
//			{
//				if(ttable[i][j] == DFA_STATE_NON_ACCEPTING )
//					cout << "  / ";
//				else if(ttable[i][j] == DFA_STATE_ACCEPTING)
//					cout << " Ac ";
//				else if(ttable[i][j] == DFA_STATE_REJECTING)
//					cout << " Ri ";
//				else
//					cout << "  X ";
//			}
//		}
//
////		if(!this->is_inside_blue_states(i) && !this->is_inside_red_states(i))
////			cout << " W ";
////		else if(this->is_inside_blue_states(i))
////			cout << " B";
////		else
////			cout << " R";
//
//		cout << endl;
//	}
//
//	cout << "--------------------------" << endl;
//
//}


void gi::dfa::set_ttable(int** ext_ttable){
	ttable = ext_ttable;
}

int** gi::dfa::get_ttable(){
	return ttable;
}

int gi::dfa::get_ttable(int i, int j){

	if(i<num_states && j<dim_alphabet+1)
		return ttable[i][j];
	else{
		cerr<<"dfa::get_ttable: out of bound"<<endl;
		exit(EXIT_FAILURE);
	}

}

void gi::dfa::set_ttable_entry(int i, int j, int v){
	if(i<num_states && j<dim_alphabet+1)
		ttable[i][j]=v;
	else{
			cerr<<"dfa::set_ttable: out of bound"<<endl;
			exit(EXIT_FAILURE);
		}
}


void	gi::dfa::set_acceptor_state(int state)
{
	ttable[state][dim_alphabet] = DFA_STATE_ACCEPTING;
}


void	gi::dfa::set_rejector_state(int state)
{
	ttable[state][dim_alphabet] = DFA_STATE_REJECTING;
}


int gi::dfa::get_dim_alphabet(){
	return dim_alphabet;
}

int gi::dfa::get_num_states(){
	return num_states;
}

int gi::dfa::get_start_state(){
	return start_state;
}


void gi::dfa::print_dfa_dot(string title, const char *file_path)
{
	ofstream myfile;
	myfile.open(file_path);

	string header = "digraph "+title+" {\n";
	string start_state = "__start0 [label=\"\" shape=\"none\"];\n\n";

	start_state = start_state + "rankdir=LR;\nsize=\"8,5\";\n\n";

	//States
	string states = "";
	string shape = "";
	string style="";
	string color="";
	for(int i=0; i<num_states; ++i)
	{
		if(ttable[i][dim_alphabet] == DFA_STATE_UNREACHABLE)
			continue;

		if(ttable[i][dim_alphabet] == DFA_STATE_ACCEPTING){
			shape = "doublecircle";
			style = "rounded,filled";
		}
		else if(ttable[i][dim_alphabet] == DFA_STATE_REJECTING){
			shape = "circle";
			style = "filled";
		} else {
			shape = "circle";
			style = "filled";
		}

		color="white";

		states = states + "s"+intTostring(i)+" [style=\""+style+"\", color=\"black\", fillcolor=\""+color+"\" shape=\""+shape+"\", label=\""+intTostring(i)+"\"];\n";
	}

	// Transizioni
	string transitions = "";
	for(int i=0; i<num_states; ++i){
		for(int j=0; j<dim_alphabet; ++j){
			int arrive_state = ttable[i][j];
			if(arrive_state == ND)
				continue;

			transitions = transitions + "s"+intTostring(i)+" -> s"+intTostring(arrive_state)+" [label=\""+	intTostring(j)+"\"];\n";
		}
	}

	string end = "__start0 -> 0;";
	string footer ="\n}";

	myfile << header << start_state << states << transitions /*<< end*/<<footer;

	myfile.close();
}



void gi::dfa::print_dfa_dot_mapped_alphabet(string title, const char *file_path)
{
	string state_name_prefix = "q";
	ofstream myfile;
	myfile.open(file_path);

	string header = "digraph "+title+" {\n";

	string start_state = "__start0 [style = invis, shape = none, label = \"\", width = 0, height = 0];\n\n";

	start_state = start_state + "rankdir=LR;\nsize=\"8,5\";\n\n";

	string start_arrow = "";
	start_arrow = "subgraph cluster_main { \n\tgraph [pad=\".75\", ranksep=\"0.15\", nodesep=\"0.15\"];\n\t style=invis; \n\t__start0 -> s0 [penwidth=2];\n}\n";

	//States
	string states = "";
	string shape = "";
	string style="";
	string color="";
	for(int i=0; i<num_states; ++i)
	{
		if(ttable[i][dim_alphabet] == DFA_STATE_UNREACHABLE)
			continue;

		if(ttable[i][dim_alphabet] == DFA_STATE_ACCEPTING){
			shape = "doublecircle";
			style = "rounded,filled";
		}
		else if(ttable[i][dim_alphabet] == DFA_STATE_REJECTING){
			shape = "circle";
			style = "filled";
		} else {
			shape = "circle";
			style = "filled";
		}

		color="white";

		states = states + "s"+intTostring(i)+" [style=\""+style+"\", color=\"black\", fillcolor=\""+color+"\" shape=\""+shape+"\", label=\""+state_name_prefix+intTostring(i)+"\"];\n";
	}


	// Transizioni
	string transitions = "";
	//map<string, string> label_for_transiction;					// La chiave individua una coppia di stati tra cui potrebbe esserci una transizione
																	// Il valore è la label da stampare, contenente tutti i simboli per cui avviene quella transizione

	vector< vector<string> > label_for_transiction(num_states, vector<string>(num_states));

	for(int i=0; i<num_states; ++i){
		for(int j=0; j<dim_alphabet; ++j){

			int arrive_state = ttable[i][j];

			if(arrive_state == ND)
				continue;

			string transition_symbol = charToString(alphabet[j]);

			if(label_for_transiction[i][arrive_state].length() == 0)
				label_for_transiction[i][arrive_state] = label_for_transiction[i][arrive_state] + transition_symbol;
			else if(label_for_transiction[i][arrive_state].length() % 9 == 0)			// Inserisce ogni 7 simboli un ritorno a capo nella label
				label_for_transiction[i][arrive_state] = label_for_transiction[i][arrive_state] + "\\n" + transition_symbol;
			else
				label_for_transiction[i][arrive_state] = label_for_transiction[i][arrive_state] + "," +transition_symbol;


			// ORIGINALE un carattere - una transizione
			//transitions = transitions + "s"+intTostring(i)+" -> s"+intTostring(arrive_state)+" [label=\""+inverse_mapped_alphabet[alphabet_symbols[j]]+"\"];\n";
		}
	}


	for(int i=0; i<num_states; ++i)
		for(int j=0; j<num_states; ++j){
			if(label_for_transiction[i][j].compare(""))
				transitions = transitions + "s"+intTostring(i)+" -> s"+intTostring(j)+" [label=\""+label_for_transiction[i][j]+"\"];\n";
		}

	string end = "__start0 -> 0;";
	string footer ="\n}";

	myfile << header << start_state <<  states << start_arrow << transitions /*<< end*/<<footer;

	myfile.close();
}


void gi::dfa::set_num_state(int n){
	num_states = n;
}


void gi::dfa::set_alphabet(const char* alf, const int d_alf)
{
	// Erase the existing alphabet
	if(alphabet!=NULL) delete[] alphabet;

	// Set new size of alphabet
	dim_alphabet = d_alf;

	// Instance for new alphabet
	alphabet = new char[dim_alphabet];

	// Copy alphabet in input to alphabet for current dfa
	// inside mapped_alphabet for every symbol there is the associated index
	for(SYMBOL i=0; i<dim_alphabet; ++i){
		alphabet[i] = alf[i];
		mapped_alphabet[(char)alf[i]] = i;
	}

}


char* gi::dfa::get_alphabet(){
	return alphabet;
}


int gi::dfa::get_arrive_state(vector<SYMBOL> dfa_string)
{
	int state = 0;
	int next_state=ND;
	int lettera = 0;

	typedef	vector<SYMBOL>::const_iterator It_posneg;

	for(It_posneg i=dfa_string.begin(); i!=dfa_string.end(); ++i){
		next_state = ttable[state][(*i)];
		if(next_state == ND){
			state = ND;
			break;
		}
		state = next_state;
	}

	return state;		// Ritorna ND se la stringa non è compatibile, diversamente ritorna lo stato dove termina la stringa
}


gi::dfa gi::dfa::read_dfa_file(const string file_name)
{
	char nameDFA[BUFFER_SIZE];
	char line[BUFFER_SIZE];

	char *alphabet_file=NULL;

	int counter = 0;
	int num_total_line = 0;
	int cstato = 0;
	int calfabeto = 0;
	int ctransizione = 0;
	int current_line = 0;
	int start_state = 0;

	string n;

	ifstream read;
	string template_line;

	dfa res;

	// Open connection to file
	read.open(file_name.c_str());

	if(read.is_open())
		cout << "File " << file_name << " is open."<< endl;
	else{
		cerr << "Error opening " << file_name << endl;
		exit(EXIT_FAILURE);
	}

	// initial state
	start_state = 0;

	// Read first line and set "num states", "dim alf" and "dfa name"
	read.getline(line,BUFFER_SIZE);
	counter = sscanf(line, "%d %d %s", &(res.dim_alphabet), &(res.num_states), nameDFA);

	// Check if the first line is complete
	if(counter != 3){
		cout << "Error in first line of file" << endl;

		exit(EXIT_FAILURE);
	}


	// read the alphabet
	alphabet_file = new char[res.dim_alphabet];
	read.getline(line,BUFFER_SIZE);

	istringstream iss(line);

	counter=0;
	while (iss >> n){
		if(counter >= res.dim_alphabet)
			break;

		alphabet_file[counter++] = n[0];
	}


	// check read alphabet
	if(counter != res.dim_alphabet){
		cout << "Error in second line of file" << endl;

		if(alphabet_file)
			delete [] alphabet_file;

		exit(EXIT_FAILURE);
	}

	res.set_alphabet(alphabet_file, res.dim_alphabet);

	if(alphabet_file)
		delete [] alphabet_file;


	///////////// compute utility values//////////////////
	// expected lines of file
	num_total_line = (res.dim_alphabet+1)*(res.num_states);

	// template of the all lines of the file
	template_line = (string)nameDFA+"[%d][%d] = %d;";
	/////////////////////////////////////////////////////



	// allocate memory for ttable
	res.ttable = new int*[res.num_states];

	// initialize ttable rows
	// "+2" is for algorithm like EDSM with states Type and Colour
	for(int i=0; i<res.num_states; ++i){
		res.ttable[i] = new int[res.dim_alphabet+1];
		for(int j=0; j<res.dim_alphabet+1; ++j)
			res.ttable[i][j]=ND;
	}


	// Parsing the file
	while(!read.eof())
	{
		read.getline(line,BUFFER_SIZE);

		string cline = line;

		// Handler for last line
		string trimmedline = trim(cline);
		if(trimmedline == "")										// Happen only in the last line
			continue;

		++current_line;

		// Integrity check
		if(current_line > num_total_line){

			cerr << "ERROR - Line number greater than max" << endl;
			exit(EXIT_FAILURE);
		}

		// Read line and set transition
		counter = sscanf(line, template_line.c_str(), &cstato, &calfabeto, &ctransizione);
		if(counter != 3)
			cerr << "ERROR in current line"<<current_line<<endl;

		res.ttable[cstato][calfabeto] = ctransizione;
	}

	// Close connection
	read.close();

	return res;
}


//void gi::dfa::read_example_file(string path_samples, vector<SYMBOL>* &positive, int* dim_positive, vector<SYMBOL>* &negative,  int *dim_negative, int* &wp = NULL, int* &wn = NULL)
//{
//	cout << "Reading strings from txt file: "<<endl;
//	int cp = 0;														// Numero di stringhe positive del linguaggio
//	int cn = 0;														//   -    -     -     negative  -      -
//	char ch;
//
//	cout << path_samples << endl;
//
//	fstream fin(path_samples.c_str(), fstream::in);
//
//	if(!fin){
//		cerr<<"An error occurred in opening file :("<<endl;
//		exit(EXIT_FAILURE);
//	}
//
//	while (fin >> noskipws >> ch) {									// Faccio in conteggio previo del numero di stringhe positive e negative presenti nel txt
//		if(ch == '\n')
//			continue;
//		else if(ch == '+')
//			cp++;
//		else if(ch == '-')
//			cn++;
//	}
//	(*dim_positive) = cp;
//	(*dim_negative) = cn;
//
//	positive = new vector<SYMBOL>[cp];
//	negative = new vector<SYMBOL>[cn];
//	wp	= new int[cp];
//	wn	= new int[cn];
//
//	for(int i=0; i<cp; ++i)
//		wp[i] = 0;
//	for(int i=0; i<cn; ++j)
//		wn[i] = 0;
//
//
//	cout << intTostring(cp) + " positivi"<<endl;
//	cout << intTostring(cn) + " negativi"<<endl;
//
//	int flagcp = 0;
//	int flagcn = 0;
//	bool casopositive = false;
//	bool casonegative = false;
//	bool primap = true;
//	bool priman = true;
//
//	ifstream infile(path_samples.c_str());
//
//	bool first = true;
//	bool second = false;
//	string line;
//
//	int local_dim_alphabet;
//
//	while (getline(infile, line))
//	{
//	    istringstream iss(line);
//	    int a;
//	    string n;
//
//
//	    // Read first line for dim alphabet
//	    if(first){
//	    	if (!(iss >> a)) { break; } // error
//	    	local_dim_alphabet = a;
//	    	//cout << "dimensione alfabeto " << a << endl;
//	    	first = false;
//	    	second = true;
//
//	    	continue;
//	    }
//
//	    if(second){
//
//	    	int counter=0;
//	    	while (iss >> n){
//	    		if(counter>=local_dim_alphabet)
//	    			break;
//
//	    		if(mapped_alphabet.find(n[0])!=mapped_alphabet.end()){
//	    			cerr<<"Error in reading example: alphabet different from that of dfa!"<<endl;
//
//	    			exit(EXIT_FAILURE);
//	    		}
//	    	}
//
//	    	// Alphabet
//	    	if(counter!= local_dim_alphabet){
//	    		cerr<<"Error in reading example: number of red alphabet symbols mismatches with the declared one!"<<endl;
//	    		cerr<<"Expected symbols: "<<dim_alphabet<<endl;
//	    		cerr<<"Red symbols: "<<counter<<endl;
//
//
//	    		exit(EXIT_FAILURE);
//	    	}else if(local_dim_alphabet != dim_alphabet){
//	    		cerr<<"Error in reading example: alphabet different from that of dfa!"<<endl;
//
//	    		exit(EXIT_FAILURE);
//	    	}
//
//	    	// alphabet ok ;)
//	    	second= false;
//	    }
//
//	    bool weight = true;
//
//	    // Read lines
//		while (iss >> n)
//		{
//			if( !n.compare("+") ){
//
//				weight = true;
//				casopositive = true;
//				casonegative = false;
//				if(primap){												// Se è il primo caso evito l'incremento
//					primap =false;
//					continue;
//				}
//				flagcp++;
//				continue;
//
//			}else if( !n.compare("-") ){
//				weight = true;
//				casonegative = true;
//				casopositive = false;
//				if(priman){												// Se è il primo caso evito l'incremento
//					priman =false;
//					continue;
//				}
//				flagcn++;
//				continue;
//			}
//
//
//			int tmp = mapped_alphabet[(char) n[0]];
//
//			if(weight){
//				weight = false;
//
//				if(casopositive)
//					wp[flagcp] = stringToint(n);
//				else if(casonegative)
//					wn[flagcn] = stringToint(n);
//
//			} else {
//
//				if(casopositive)
//							positive[flagcp].push_back(tmp);
//				else if(casonegative)
//							negative[flagcn].push_back(tmp);
//			}
//		}
//	}
//}


vector<SYMBOL> gi::dfa::equivalence_query(dfa* dfa_hp){
	vector<SYMBOL> witness;

	#ifdef DEBUG_2
	cout << endl << "--------------------------" << endl;
	cout << "EQUIVALENCE QUERY" << endl;
	cout << "--------------------------" << endl;
	#endif


	// Build union DFA of target dfa (thisone) and dfa_hp
	dfa* dfa_union = this->unionDFA(dfa_hp);
	#ifdef DEBUG_2
	dfa_union->print_dfa_ttable("DFA UNIONE");
	#endif


	// Table-filling algorithm on union dfa
	SYMBOL* distincts_table = dfa_union->table_filling();

	// Extract list of equivalent states from table of distinct states,
	// every vector contain a list of equivalent states for the state that correspond to the vector.
	vector<int>* equivalent_states_list = dfa_union->equivalent_states_list_from_table(distincts_table);

	// Verify if start states of dfas are equivalent:
	// Controllo se fra gli stati equivalenti allo stato 0 (stato iniziale dell'automa corrente)
	// c'è lo stato iniziale dell'automa ipotesi identificato dall'indice "num_state".
	// Nel caso fosse presente, non entro nell'if, e ritorno un vector vuoto come controesempio.
	// Se così non fosse (è il caso dell'"end()") genero un controesempio.
	if(equivalent_states_list[0].end() == find(equivalent_states_list[0].begin(), equivalent_states_list[0].end(), num_states) )
		witness = dfa_union->witness_from_table(distincts_table, num_states);


	// Free allocated memory
	if(distincts_table != NULL)
		delete[] distincts_table;

	delete [] equivalent_states_list;

	delete dfa_union;


	return witness;
}


//vector<SYMBOL>  gi::dfa::equivalence_query_approximate(dfa* dfa_hp, string samplestestpath)
//{
//	#ifdef DEBUG_2
//	cout << endl << "--------------------------" << endl;
//	cout << "EQUIVALENCE QUERY" << endl;
//	cout << "--------------------------" << endl;
//	#endif
//
//
//	// Build union DFA of target dfa (thisone) and dfa_hp
//	dfa* dfa_union = this->unionDFA(dfa_hp);
//	#ifdef DEBUG_2
//	dfa_union->print_dfa_ttable("DFA UNIONE");
//	#endif
//
//
//	// Table-filling algorithm on union dfa
//	SYMBOL* distincts_table = dfa_union->table_filling();
//
//	// Extract list of equivalent states from table of distinct states
//	vector<int>* equivalent_states_list = dfa_union->equivalent_states_list_from_table(distincts_table);
//
//	// Verify if start states of dfas are equivalent
//	bool eq_dfa = false;
//	for(int i=0; i<equivalent_states_list[0].size(); ++i)
//		if(equivalent_states_list[0][i] == this->get_num_states())
//			eq_dfa = true;
//
//
//	// If start states aren't equivalent then dfas aren't equivalent. Generate a witness (counterexample)
//	//  (index of first state of dfa_hp in states list is dfa_hp->num_state)
//	vector<SYMBOL> witness;
//	if( !eq_dfa )
//		witness = dfa_hp->witness_approximate(dfa_hp,samplestestpath);
//
//	// Free allocated memory
//	if(distincts_table != NULL)
//		delete[] distincts_table;
//
//	return witness;
//}


//vector<SYMBOL> gi::dfa::witness_approximate(dfa* dfa_hp, string path_samples)
//{
//	// Test samples
//	int dim_positive;
//	int dim_negative;
//	vector<SYMBOL> *positive=NULL;
//	vector<SYMBOL> *negative=NULL;
//
//	vector<SYMBOL> witness;
//
//	bool negative_bool=true;
//
//	read_example_file(path_samples, positive, &dim_positive, negative,  &dim_negative);
//
//
//	for(int i=0; i<dim_positive; ++i)
//		if(dfa_hp->membership_query(positive[i]) == 0){
//			witness = positive[i];
//			negative_bool=false;
//			break;
//		}
//
//	if(negative_bool)
//		for(int i=0; i<dim_negative; ++i)
//			if(dfa_hp->membership_query(negative[i]) == 1){
//				witness = negative[i];
//				break;
//			}
//
//	if(positive) delete[] positive;
//	if(negative) delete[] negative;
//
//	return witness;
//}


bool gi::dfa::membership_query(vector<SYMBOL> str){

	// Check if arrive_state is ND (for DFA without sink state)
	int arrive_state = get_arrive_state(str);
	if(arrive_state == ND)
		return false;

	if(ttable[arrive_state][dim_alphabet] == 1)
		return true;
	else
		return false;
}


// Call it in a DFA union
vector<SYMBOL> gi::dfa::witness_from_table(SYMBOL* distinct, int start_state_dfa_hp)
{
	// Se automi NON equivalenti, creo la witness
	vector<SYMBOL> wit;

	int icoppia = 0;
	int jcoppia = start_state_dfa_hp;

	SYMBOL input;

	#ifdef DEBUG_2
	cout << "--- Creo il CONTROESEMPIO --- " << endl;
	#endif

	while(1)
	{
		//V cout << "Coppia: "<<icoppia << ","<<jcoppia<<"; vale: "<<distinct[icoppia][jcoppia]<<endl;
		int n = num_states;
		int k = (n*(n-1)/2) - (n-icoppia)*((n-icoppia)-1)/2 + jcoppia - icoppia - 1;

		if(distinct[k] == DFA_TF_STATE_N)
			cout << "PROBLEMA! Richiesta di controesempio con automi equivalenti";

		input = distinct[k];
		// Caso in cui lo start state dei 2, uno è accettante e l'altro no
		if(input == DFA_TF_STATE_X)
			break;

		//V cout << "stinga parziale: "<<wit<<" + " << input << endl;
		// TODO: il push piuttosto che concatenazione è corretto? old: wit = wit + input;
		wit.push_back(input);
		//V cout << "Vale: "<< input << ", ascii: "<< (int)input << ", diff_i:" << (int)input - 48 << endl;

		icoppia = ttable[icoppia][input];
		jcoppia = ttable[jcoppia][input];

		if(distinct[k] == DFA_TF_STATE_X){
			//V cout << "fine stringa" << endl;
			break;
		}
	}

	#ifdef DEBUG_2
	cout << "Controesempio è: "<< wit << endl;
	#endif

	return wit;
}


vector<int>* gi::dfa::equivalent_states_list_from_table(SYMBOL* distincts)
{
	#ifdef DEBUG_2
	cout << endl << "--------------------------" << endl;
	cout << "Lista stati equivalenti:" << endl;
	cout << "--------------------------" << endl;
	#endif

	vector<int>* stati_equivalenti = new vector<int>[num_states];
	int n= num_states;

	for(int i=0; i<(num_states-1); ++i)
		for(int j=i+1; j<num_states; ++j){
			int k= (n*(n-1)/2) - (n-i)*((n-i)-1)/2 + j - i - 1;
			if(distincts[k] == DFA_TF_STATE_N)
				stati_equivalenti[i].push_back(j);
		}

	#ifdef DEBUG_2
	cout << "N di stati finali: " << num_states << endl;
	for(int i=0; i<num_states; ++i)
		cout << "S["<<i<<"] --> "<< stati_equivalenti[i] << endl;

	cout << "--------------------------" << endl;
	#endif

	return stati_equivalenti;
}


double gi::dfa::get_complexity()
{
	return num_states;
}

// Generazione random di samples
//bool gi::dfa::generate_pos_neg_samples(int n_pos_samples,int n_neg_samples, const char * file_path)
//{
//	srand (time(NULL));
//	int rand_trans =0;
//	map<vector<SYMBOL>, int> samples;										// Stringa e indice dello stato
//	typedef	map<string,int>::const_iterator It;
//
//	int freq_to_stop = 1+2*(dim_alphabet);							// The inverse is probability to stop
//
//	string input_simbol = "";
//	bool finded_one_minimun_accepted_str = false;
//	string tmp_sample = "";
//	string sample_pos = "";
//
//	bool go_next = true;
//	int current_state = 0;
//
//	ofstream myfile;
//	myfile.open(file_path);
//
//	int num_attempts = 0;
//
//	// Write positve samples
//	while(samples.size() < n_pos_samples)
//	{
//		rand_trans = rand() % dim_alphabet;									// Input for first transition
//
//		input_simbol = intTostring(rand_trans);								// Int to String
//
//		finded_one_minimun_accepted_str = false;
//		tmp_sample = "";
//		sample_pos = "";
//
//		// Make positive sample
//		go_next = true;
//		current_state = 0;
//		int num_iteration = 0;
//		while(go_next)
//		{
//			//this->print_dfa_ttable("l");
//			current_state = ttable[current_state][rand_trans];
//			input_simbol = intTostring(rand_trans);
////			cout << "state:"<<current_state<< ", input:"<<input_simbol<<endl;
//
//			// Increment the final random string
//			tmp_sample = tmp_sample + " " +input_simbol;
//
//
//			// Visit DFA with input symbol
//			if(ttable[current_state][dim_alphabet] == 1){
//				finded_one_minimun_accepted_str = true;
//				sample_pos = tmp_sample;
//			}
//
//			// TODO Gestire stati pozzo
//
//			// Keep attention if stop because the upper bound about iteration... maybe some problem!
//			//if(num_iteration > 2000)
//				//cout << "! Limitated generation !"<<endl;
//
//			// Stop condition
//			if( ((rand() % freq_to_stop) == 1 && finded_one_minimun_accepted_str) || num_iteration > 5000)	// If I stop i should have some good string
//				go_next = false;
//
//			// Generate new input symbol
//			rand_trans = rand() % dim_alphabet;
//
//			++num_iteration;
//		} // END while(go_next)
//
//		samples[sample_pos] = 1;
//
//		++num_attempts;
//
//		// Num of attempts to do before stop.. could be a problem!
//		if(num_attempts > 60000)
//			break;
//
//
//		// If you want check that every string is inside language
//		/*if( membership_query(sample_pos) == 1)
//			cout << "Correct"<<endl;
//		else{
//			cout << "Wrong"<< endl;
//			exit(1);
//		}*/
//	} // END while(sample.size...)
//
//	if(num_attempts > 60000)
//		cout << "! Poor set of samples !" << endl;
//
//	// Generation of NEGATIVE SAMPLES
//	while(samples.size() - n_pos_samples < n_neg_samples)
//	{
//		//cout << "DimTot: "<<samples.size() <<", Positivi: "<<n_pos_samples <<", Negativi: "<<n_neg_samples<<endl;
//		int mod_type = rand() % 3;											// 0: substituting, 1:inserting, 2:deleting
//		int n_of_editing = getPoisson(3); //Poisson distribution centred on 3
//		#ifdef DEBUG_2
//		cout << "Poissoin ha ritornato "<< n_of_editing << " modifiche"<<endl;
//		#endif
//
//		int pos_to_edit = 0;
//		int n_changes = 0;
//
//		tmp_sample = "";
//
//		int* clean_sample = NULL;
//		for(It p1=samples.begin(); p1!=samples.end(); ++p1)
//		{
//			if((*p1).second == 1)
//				tmp_sample = (*p1).first;
//
//			// Extract sample without white space
//			int dim_clean_sample = tmp_sample.length()/2;
//			if(clean_sample != NULL)
//				delete[] clean_sample;
//			clean_sample = new int[dim_clean_sample];
//
//
//			int c = 0;
//			for(int i=1; i<tmp_sample.length(); i=i+2){
//				clean_sample[c] = tmp_sample[i] - '0';
//				c++;
//			}
//
//			if(dim_clean_sample == 0)
//				continue;
//
//			// Do changes
//			n_changes = 0;
//			while(n_changes < n_of_editing)									// Make "n_changes" for the single string up to "n_of_editing"
//			{
//				pos_to_edit = rand() % dim_clean_sample;
//				//cout << endl << endl <<"Posizione: "<<pos_to_edit<<endl;
//
//				// Round shift of "pos_to_edit" to right
//				int shift = dim_clean_sample-pos_to_edit;
//
//				/*cout << "Originale"<<endl;
//				for(int i=0; i<dim_clean_sample; ++i)
//					cout << clean_sample[i];
//				cout << endl;*/
//
//				rotate(clean_sample, clean_sample + pos_to_edit, clean_sample + dim_clean_sample);
//
//				/*cout << "Rotata"<<endl;
//				for(int i=0; i<dim_clean_sample; ++i)
//					cout << clean_sample[i];
//				cout << endl;*/
//
//				mod_type = rand() % 3;
//
//				++n_changes;
//
//				if(mod_type == 0){											// Substitution
//					//cout << endl<<"SOSTITUZIONE"<<endl;
//
//					int new_char = rand() % dim_alphabet;
//					while(clean_sample[0] == new_char)
//						new_char = rand() % dim_alphabet;
//
//					clean_sample[0] = new_char;
//
//					/*STAMPAcout << "Dopo sostituzione:"<<endl;
//					for(int i=0; i<dim_clean_sample; ++i)
//						cout << clean_sample[i];
//					cout << endl;*/
//				}
//				else if(mod_type == 1){
//					//cout <<endl<< "INSERIMENTO"<<endl;
//					int new_char = rand() % dim_alphabet;
//
//					int* ScleansampleTMP = new int[dim_clean_sample+1];
//					ScleansampleTMP[0] = new_char;
//					for(int i=1; i<dim_clean_sample+1; ++i)
//						ScleansampleTMP[i] = clean_sample[i-1];
//
//					delete[] clean_sample;
//					clean_sample = ScleansampleTMP;
//
//					dim_clean_sample += 1;
//
//					/*STAMPAcout << "Dopo inserimento:"<<endl;
//					for(int i=0; i<dim_clean_sample; ++i)
//						cout << clean_sample[i];
//					cout << endl;*/
//
//				}
//				else if(mod_type == 2){
//					//cout << endl << "CANCELLO"<<endl;
//					if(dim_clean_sample <= 1)
//						continue;
//
//					int* ScleansampleTMP = new int[dim_clean_sample-1];
//					for(int i=0; i<dim_clean_sample-1; ++i){
//						ScleansampleTMP[i] = clean_sample[i+1];
//					}
//
//					delete[] clean_sample;
//					clean_sample = ScleansampleTMP;
//
//					dim_clean_sample -= 1;
//
//					/*STAMPAcout << "Dopo eliminazione:"<<endl;
//					for(int i=0; i<dim_clean_sample; ++i)
//						cout << clean_sample[i];
//					cout << endl;*/
//				}
//
//				/*cout << "Rotata dopo cambiamento"<<endl;
//				for(int i=0; i<dim_clean_sample; ++i)
//					cout << clean_sample[i];
//				cout << endl;*/
//
//				//cout << "Ritorno originale dopo cambiamento:";
//				rotate(clean_sample, clean_sample + (dim_clean_sample-pos_to_edit), clean_sample + dim_clean_sample);
//				/*for(int i=0; i<dim_clean_sample; ++i)
//					cout << clean_sample[i];
//				cout << endl;*/
//			}
//
//			tmp_sample = "";
//			for(int i=0; i<dim_clean_sample; ++i)
//				tmp_sample = tmp_sample+" "+intTostring(clean_sample[i]);
//
//			if(tmp_sample.length() == 0)
//				continue;
//
//			if(this->membership_query(tmp_sample) == 0){
//				//cout << "stringa Non accettata correttamente:"<<tmp_sample<<";"<<endl;
//				samples[tmp_sample] = 0;
//			}
//
//			if(samples.size() - n_pos_samples >=  n_neg_samples)
//				break;
//		} //EDN FOR
//		if(clean_sample != NULL)
//			delete[] clean_sample;
//	}
//
//	myfile << dim_alphabet << "\n";
//
//	for(It p1=samples.begin(); p1!=samples.end(); ++p1)
//		if((*p1).second  == 1)
//			myfile << "+ "+(*p1).first+"\n";
//
//	for(It p1=samples.begin(); p1!=samples.end(); ++p1)
//		if((*p1).second  == 0)
//			myfile << "- "+(*p1).first+"\n";
//
//	myfile.close();
//
//
//	// Verify if the set of samples is structuraly complete
//	bool** scomplete = new bool*[num_state];
//	for(int i=0; i<num_state; ++i)
//		scomplete[i] = new bool[dim_alphabet+2]; 				// Column "dim_alphabet" true if state visited,
//																//   -    "dim_alphabet+1" true if accepting state used
//	for(int i=0; i<num_state; ++i)
//		for(int j=0; j< dim_alphabet+2; ++j)
//			scomplete[i][j] = false;
//	scomplete[0][dim_alphabet] = true;			// In current set often there aren't some transition to firs state
//
//	 // Init every state that is not ACC to true in the "dim_alphabet+1" column
//	for(int i=0; i<num_state; ++i)
//		if(ttable[i][dim_alphabet] != ACC){
//			scomplete[i][dim_alphabet+1] = true;
//		}
//
//	for(It p1=samples.begin(); p1!=samples.end(); ++p1){
//		if((*p1).second  == 1){
//			current_state=0;
//			for(int i=0; i<(*p1).first.length(); ++i)
//			{
//				if((*p1).first[i] == ' ')
//					continue;
//
//				int input_simbol = (*p1).first[i]-'0';
//
//				// Explorated transition
//				scomplete[current_state][input_simbol] = true;
//
//				current_state = ttable[current_state][input_simbol];
//
//				// Accepting state used?
//				if(ttable[current_state][dim_alphabet] == ACC && (i+1 == (*p1).first.length()))
//					scomplete[current_state][dim_alphabet+1] = true;
//
//				// Explorated state
//				scomplete[current_state][dim_alphabet] = true;
//
//
//			}
//		}
//	}
//
//	for(int i=0; i<num_state; ++i)
//		for(int j=0; j< dim_alphabet+2; ++j)
//			if(!scomplete[i][j]){
//				//cout << "Non  completo per " <<i<<","<<j<<endl;
//				return false;
//			}
//
//
//	// Free memory
//	if(scomplete != NULL){
//			for(int i=0; i<num_state; ++i)
//				if(scomplete[i] != NULL)
//					delete[] scomplete[i];
//
//			delete [] scomplete;
//	}
//
//
//	return true;
//}

