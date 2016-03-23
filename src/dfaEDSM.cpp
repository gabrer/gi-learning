#include "dfaEDSM.h"

#include <iostream>
#include <fstream>
#include <sstream> //isstringstream
#include <set>


gi::dfaEDSM::dfaEDSM()
:dfa(){
	blue_states=NULL;
	red_states=NULL;
}


gi::dfaEDSM::dfaEDSM(const int n_state, const int dim_alf, const char *alf, const int s_state)		// Constructor 1
:dfa(n_state, dim_alf, alf, s_state){

	blue_states = new vector<int>();
	red_states = new vector<int>();

}


gi::dfaEDSM::dfaEDSM(const int n_state, const int dim_alf, const char *alf)							// Default start state to 0
:dfaEDSM(n_state, dim_alf, alf, 0){}


gi::dfaEDSM::dfaEDSM(const dfa &d1)
:dfa(d1){

	blue_states = new vector<int>();
	red_states = new vector<int>();
}


gi::dfaEDSM::dfaEDSM(const dfaEDSM &d1)
:dfa(d1.num_states, d1.dim_alphabet, d1.alphabet, d1.start_state, (const int**) d1.ttable){

	blue_states = new vector<int>;
	red_states = new vector<int>;

	copy_blue_states(d1.get_blue_states());

	copy_red_states(d1.get_red_states());

}


gi::dfaEDSM::~dfaEDSM(){

	if(blue_states!=NULL){
		blue_states->clear();
		delete blue_states;
	}

	if(red_states!=NULL){
		red_states->clear();
		delete red_states;
	}
}


gi::dfa* gi::dfaEDSM::to_dfa(){
	return new dfa(this->num_states,this->dim_alphabet, this->alphabet, this->start_state, (const int**) this->ttable);
}


gi::dfaEDSM* gi::dfaEDSM::to_canonical_dfaEDSM_from_red_states(){

	//////////////////////////////////////////////////////////////
	// - Pulisco l'automo dagli stati irragiungibili, aggiorno le transizioni

	// Conto il numero effettivo di stati finali
	int n_final_states = 0;
	for (int i = 0; i < get_num_states(); ++i)
		if (is_inside_red_states(i))
			n_final_states++;

	//TODO: Aggiungi la copia degli stati red e blue nei vectors

	// Creo un nuovo automa senza stati irragiungibili
	int count = 0;
	dfaEDSM* finalDFA = new dfaEDSM(n_final_states, get_dim_alphabet(), get_alphabet());
	map<int, int> updated_transition;

	for (int i = 0; i < num_states; ++i) {
		if (is_inside_red_states(i)) {

			for (int j = 0; j <get_dim_alphabet() + 1; ++j) {

				// Aggiungo lo stato al nuovo automa
				finalDFA->get_ttable()[count][j] = ttable[i][j];

				updated_transition[i] = count;
			}
			++count;
		}
	}
	updated_transition[ND] = ND;


//	if(updated_transition.size() <= 2){
//		cout << "There is only one or zero red state. Returned a copy of originale dfa."<<endl;
//		delete finalDFA;
//		return new dfaEDSM(*this);
//	}


	bool stato_pozzo = false;
	// Aggiorno le transizioni
	for (int i = 0; i < finalDFA->get_num_states(); ++i)
		for (int j = 0; j < finalDFA->get_dim_alphabet(); ++j) {

			if (finalDFA->get_ttable()[i][j] == ND)										// Rilevo che c'è una transizione mancante, quindi serve uno stato pozzo
				stato_pozzo = true;

			if(updated_transition.find(finalDFA->get_ttable()[i][j]) != updated_transition.end())
				finalDFA->set_ttable_entry(i, j, updated_transition[ finalDFA->get_ttable()[i][j] ]);
			else {
				cerr << "Errore nell'aggiornamento delle stringhe"<<endl;
				exit(EXIT_FAILURE);
			}

		}

	// Stampo l'automa prima di applicare il pozzo e la minimizzazione
	//finalDFA->print_dfa_with_color("AUTOMA FINALE PREPOZZO");
	//finalDFA->print_dfa_dot("FINALEPREPOZZO", percorso.c_str());

	//finalDFA->print_dfa_dot_mapped_alphabet("FINALE_PREPOZZO", (base_path + "pulito_pre_pozzo.dot").c_str());



	//////////////////////////////////////////////////////////////
	// Controllo stato pozzo
	// - Se ci sono transizioni non definite le imposto tutte verso lo stato pozzo
	if (stato_pozzo) {
		dfaEDSM* finalDFAPozzo = new dfaEDSM(finalDFA->get_num_states() + 1, finalDFA->get_dim_alphabet(), finalDFA->get_alphabet(), 0);

		int** table = finalDFAPozzo->get_ttable();

		for (int i = 0; i < finalDFA->get_num_states(); ++i)
			for (int j = 0; j < finalDFA->get_dim_alphabet() + 1; ++j) {
				if (finalDFA->get_ttable()[i][j] == ND)
					table[i][j] = finalDFA->get_num_states();
				else
					table[i][j] = finalDFA->get_ttable()[i][j];
			}

		for (int j = 0; j < finalDFA->get_dim_alphabet(); ++j)
			table[finalDFA->get_num_states()][j] = finalDFA->get_num_states();

		delete finalDFA;
		finalDFA = finalDFAPozzo;
	}

	return finalDFA;
}


int gi::dfaEDSM::get_actual_num_states()
{
	int n_reacheable_states = 0;

	set<int> reacheable_states;

	for(int i=0; i<num_states; ++i){
		for(int j=0; j<dim_alphabet; ++j)
			if(ttable[i][j] != i && ttable[i][j] != ND)
				reacheable_states.insert(ttable[i][j]);
	}


	// Insert also the start state, colud be not reacheable from everyone by definition
	reacheable_states.insert(start_state);

	return reacheable_states.size();

}


//TODO: gi::dfaEDSM* gi::dfaEDSM::remove_unreachable_states(){}

//TODO: void gi::dfaEDSM::add_sink_state(){}


void gi::dfaEDSM::add_blue_state(int blue_state_index)
{
	// If is not already inside, insert the new blue state
	if(!this->is_inside_blue_states(blue_state_index))
		blue_states->push_back(blue_state_index);
}

void gi::dfaEDSM::add_red_state(int red_state_index)
{
	// If is not already inside, insert the new blue state
	if(!this->is_inside_red_states(red_state_index))
		red_states->push_back(red_state_index);

	// Reorder the red states list. Actually there is no reason to maintain the input order.
	sort(red_states->begin(), red_states->end());
}

void gi::dfaEDSM::remove_blue_state(int blue_state_index)
{
	// If blue state is inside vector, remove it
	if(this->is_inside_blue_states(blue_state_index))
		blue_states->erase(remove(blue_states->begin(),blue_states->end(), blue_state_index));
}

void gi::dfaEDSM::remove_red_state(int red_state_index)
{
	// If red state is inside vector, remove it
	if(this->is_inside_red_states(red_state_index))
		red_states->erase(remove(red_states->begin(),red_states->end(), red_state_index));
}

vector<int>* gi::dfaEDSM::get_blue_states() const{
	return blue_states;
}

vector<int>* gi::dfaEDSM::get_red_states() const{
	return red_states;
}

int gi::dfaEDSM::get_num_blue_states(){
	return blue_states->size();
}

int gi::dfaEDSM::get_num_red_states(){
	return red_states->size();
}

bool gi::dfaEDSM::is_inside_blue_states(int blue_state_index)
{
	if ( find(blue_states->begin(), blue_states->end(), blue_state_index) != blue_states->end() )
	   return true;
	else
	   return false;
}


bool gi::dfaEDSM::is_inside_red_states(int red_state_index)
{
	if ( find(red_states->begin(), red_states->end(), red_state_index) != red_states->end() )
	   return true;
	else
	   return false;
}


void gi::dfaEDSM::copy_blue_states(vector<int>* new_blue_vector)
{
	if(new_blue_vector != NULL && blue_states != NULL){
		blue_states->clear();

		(*blue_states) = (*new_blue_vector);
	}
}


void gi::dfaEDSM::copy_red_states(vector<int>* new_red_vector)
{
	if(new_red_vector){
		red_states->clear();

		(*red_states) = (*new_red_vector);
	}
}


void gi::dfaEDSM::print_dfa_with_color(string title)
{
	// STAMPA IL DFA
	// Nel numero di colonne è inclusa la colonna finale di stato accettante o meno (accettante -> 1)

	cout << endl<< "--------------------------" << endl;
	cout << title << endl;
	string header = "    ";
		for(int i=0; i<dim_alphabet; ++i)
			header = header + " | "+ intTostring(i);
		header = header + " - A  - C";
	cout << header << endl;

	for(int i=0; i<num_states; ++i){
		cout << "S"<<i<<"  ";

		for(int j=0; j<dim_alphabet+1; ++j)
		{

			if(j < dim_alphabet && ttable[i][j] == ND)				// Valori delle transizioni, o ND o il valore
				cout << " N ";

			else if(j < dim_alphabet)
				cout << " "<< ttable[i][j] <<" ";

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

		if(!this->is_inside_blue_states(i) && !this->is_inside_red_states(i))
			cout << " W ";
		else if(this->is_inside_blue_states(i))
			cout << " B";
		else
			cout << " R";

		cout << endl;
	}

	cout << "--------------------------" << endl;
}


void gi::dfaEDSM::print_dfa_with_color_mapped_alphabet(string title)
{
	// STAMPA IL DFA
	// Nel numero di colonne è inclusa la colonna finale di stato accettante o meno (accettante -> 1)

	cout << endl<< "--------------------------" << endl;
	cout << title << endl;
	string header = "    ";
		for(int i=0; i<dim_alphabet; ++i)
			header = header + " | "+ alphabet[i];
		header = header + " - A  - C";
	cout << header << endl;

	for(int i=0; i<num_states; ++i){
		cout << "S"<<i<<"  ";

		for(int j=0; j<dim_alphabet+1; ++j)
		{

			if(j < dim_alphabet && ttable[i][j] == ND)				// Valori delle transizioni, o ND o il valore
				cout << " N ";
			else if(j < dim_alphabet)
				cout << " "<< alphabet[ttable[i][j]] <<" ";

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

		if(!this->is_inside_blue_states(i) && !this->is_inside_red_states(i))
			cout << " W ";
		else if(this->is_inside_blue_states(i))
			cout << " B";
		else
			cout << " R";

		cout << endl;
	}

	cout << "--------------------------" << endl;

}


void gi::dfaEDSM::print_dfa_dot(string title, const char *file_path)
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

		if(is_inside_red_states(i))
			color="#ff817b";
		else if(is_inside_blue_states(i))
			color="powderblue";
		else
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


void gi::dfaEDSM::print_dfa_dot_mapped_alphabet(string title, const char *file_path)
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

		if(is_inside_red_states(i))
			color="#ff817b";
		else if(is_inside_blue_states(i))
			color="powderblue";
		else
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
