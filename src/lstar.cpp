/*
 * obt.cpp
 *
 *  Created on: 01 ott 2015
 *      Author: piero
 */

/*
 * obt.cpp
 *
 *  Created on: 14/mag/2015
 *      Author: Gabriele
 */

#include <lstar.h>
#include <utilities.h>
#include <vector>



gi::lstar::lstar(const dfa &targetdfa)
{
	target =  new dfa(targetdfa);

	// Inizializzo
	vector<SYMBOL> empty_prefix;
	pref[empty_prefix] = obtState::OBT_STATE_RED;										// Inizializzo lo stato lambda a RED

	// set alphabet size
	dim_alfabeto = target->get_dim_alphabet();

	// allocate vector from the alphabet
	alfabeto = new vector<SYMBOL> [dim_alfabeto];

	// Inizializzo ogni membro dell'alfabeto a BLU
	for(SYMBOL i=0; i<dim_alfabeto; ++i){
		alfabeto[i].push_back(i);
		pref[ alfabeto[i] ] = obtState::OBT_STATE_BLU;
	}




	experiments.push_back(empty_prefix);									// Inizializzo a "lambda" exp

	mq[empty_prefix] = target->membership_query(empty_prefix);
	for(SYMBOL i=0; i<target->get_dim_alphabet(); ++i)
		mq[alfabeto[i]] = target->membership_query(alfabeto[i]);

}

gi::lstar::~lstar(){
	if(alfabeto != NULL)
		delete [] alfabeto;

	if(target != NULL)
		delete target;
}


//bool obt::closed()
bool gi::lstar::close_obt(){
	// Se la tabella non è chiusa, viene chiusa
	// Se la tabella non viene modificata (era chiusa), torna false, altrimenti true

	// *** CHIUSURA ***
	// Per ogni BLU deve esistere una riga uguale tra i RED
	#ifdef DEBUG_2
	cout << "--------------------------" << endl;
	cout << "Verifico la CHIUSURA" << endl;
	cout << "--------------------------" << endl;
	#endif

	// La tabella è chiusa se la riga relativa ad ogni prefisso BLU è identica ad almeno una riga di un prefisso ROSSO

	// Ciclo solo sui BLU
	for(auto p1=pref.begin(); p1!=pref.end(); ++p1){
		if(p1->second == obtState::OBT_STATE_BLU){
			bool redFound = false;


			vector<obtState> rowB = get_row(p1->first);
			#ifdef DEBUG_3
			cout << "BLU da verificare: " << prefB << endl;
			#endif

			for(auto p2=pref.begin(); p2!=pref.end(); ++p2){				// Ciclo tra i ROSSI
				if(p2->second == obtState::OBT_STATE_RED){
					vector<obtState> rowR = get_row(p2->first);

					if( rowB == rowR ){
						redFound = true;
						break;
					}
				}
			}


			// Per lo stato blu non ho trovato nessun corrispettivo stato rosso
			if(!redFound){
				#ifdef DEBUG_2
				cout << "Tabella NON chiusa" << endl;
				cout << "Promuovo a RED "<< p1->first << endl;
				#endif

				vector<SYMBOL> prefB = p1->first;

				// Promuovo a RED il BLU in questione
				pref[prefB] = obtState::OBT_STATE_RED;

				// Aggiungo ai BLU il nuovo RED concatenato ad ogni simbolo dell'alfabeto
				for(int i=0; i<dim_alfabeto; ++i){

					vector<SYMBOL> tmpB = append_vectors(&prefB, alfabeto+i);

					//tmpB.insert(tmpB.end(), alfabeto[i].begin(), alfabeto[i].end());

					//creo un nuovo prefisso
					pref[tmpB] = obtState::OBT_STATE_BLU;

					// Riempio la riga dell'obt per il nuovo prefisso
					for(auto j=experiments.begin(); j!=experiments.end(); ++j){
						vector<SYMBOL> tmpB_row = append_vectors(&tmpB, &(*j));

						mq[tmpB_row] = target->membership_query( tmpB_row );
					}
				}

				// la funzione aggiunge un solo esperimento e poi esce, segnalando che la tabella è stata modificata.
				// N.B.: la tabella NON è necessariamente chiusa, ha solo un problema in meno
				return true;
			}
		}
	}

	#ifdef DEBUG_2
	cout << "--------------------------" << endl;
	cout << "Fine CHIUSURA" << endl;
	cout << "--------------------------" << endl;
	#endif

	return false;
}

bool gi::lstar::make_obt_consistent()
{
	// *** CONSISTENZA ***
	// 2 stati uguali a parità di input devono raggiungere lo stesso stato

	#ifdef DEBUG_2
	cout << "--------------------------" << endl;
	cout << "Verifico la CONSISTENZA"<< endl;
	cout << "--------------------------" << endl;
	#endif

	bool consistent = true;
	for(auto p1=pref.begin(); p1!=pref.end(); ++p1){

		// Ciclo sui ROSSI: s1
		if(p1->second == obtState::OBT_STATE_RED){
			vector<SYMBOL> prefR1 = p1->first;
			vector<obtState> rowR1 = get_row(p1->first);

			// Ciclo sui ROSSI: s2
			for(auto p2=pref.begin(); p2!=pref.end(); ++p2){
				if(p2->second == obtState::OBT_STATE_RED){
					consistent = true;

					// Verifico di lavorare con 2 red distinti (NB: "compare" qui restituisce T se sono diversi)
					if(p1->first == p2->first)
						continue;

					vector<SYMBOL> prefR2 = p2->first;
					vector<obtState> rowR2 = get_row(p2->first);

					// Se sono row uguali, verifico che non esistano inconsistenze
					if(rowR1 == rowR2){
						//V cout << "Hanno righe uguali, quindi potenziali generatori di inconsistenza"<< endl;
						for(int i=0; i<dim_alfabeto; ++i){
							vector<SYMBOL> prefR1A = append_vectors(&prefR1, alfabeto + i);
							vector<SYMBOL> prefR2A = append_vectors(&prefR2, alfabeto + i);

							// Se NON è consistente, aggiorno EXP e la tabella di conseguenza

							vector<obtState> prefR1A_row = get_row(prefR1A);
							vector<obtState> prefR2A_row = get_row(prefR2A);

							if(prefR1A_row != prefR2A_row){
								unsigned int q=0;
								//determino il suffisso per cui le due righe differiscono
								for(auto w=experiments.begin(); w != experiments.end(); ++w, ++q)
									if(prefR1A_row.at(q) != prefR2A_row.at(q)){
										vector<SYMBOL> expc = append_vectors(alfabeto + i, &(*w));

										//aggiungo il nuovo suffisso (esperimento)
										experiments.push_back(expc);

										//riempio la colonna dell'obt per il nuovo suffisso
										for(auto p3=pref.begin(); p3!=pref.end(); ++p3){
											vector<SYMBOL> tmpP = p3->first;
											tmpP = append_vectors(&tmpP, &expc);

											mq[tmpP] = target->membership_query( tmpP );
										}
										#ifdef DEBUG_2
										cout << "--------------------------" << endl;
										cout << "Fine CONSISTENZA"<< endl;
										cout << "--------------------------" << endl;
										#endif

										// la funzione aggiunge un solo suffisso e poi esce, segnalando che la tabella è stata modificata.
										// N.B.: la tabella NON è necessariamente consistente, ha solo un problema in meno
										return true;
									}
							}
						}
					}

				}
			}
		}
	}

	#ifdef DEBUG_2
	cout << "--------------------------" << endl;
	cout << "Fine CONSISTENZA"<< endl;
	cout << "--------------------------" << endl;
	#endif

	return false;
}


gi::dfa* gi::lstar::obt_to_dfa()
{
	// TABELLA DI OSSERVAZIONE --> DFA

	// Creo un elenco degli stati dell'automa finale con la map "stati"
	// ELENCO DEGLI STATI: Label dello stato e 2 interi:
	// il primo per indice numerico dello stato, secondo per accettante o meno
	map<vector<SYMBOL>, int[2], vector_int_size_less> stati;

	int count_state = 0;

	// Qui faccio il controllo delle transizioni guardando le righe della tabella di oss.
	// Ciclo sui ROSSI: s1
	unsigned int a = pref.size();
	unsigned int w=0;
	for(auto i=pref.begin(); i!=pref.end(); ++i)
		cout<<"Pref "<<(w++)<<": "<<i->first<<endl;

	for(auto p1=pref.begin(); p1!=pref.end(); ++p1){
		bool stato_inDFAfinale = true;

		if(p1->second == obtState::OBT_STATE_RED){
			//string prefR1 = p1->first;

			// Ipotizzo che sia uno stato da aggiungere (stato_inDFAfinale=T).
			// Ciclando internamente "uguale" deve diventare sempre F alla fine, se rimane T
			// c'è uno stato già analizzato uguale a quello preso in considerazione quindi
			// non deve appartenere e stato_inDFAfinale=F
			//bool stato_inDFAfinale=true;
			// Ciclo sui ROSSI: s2
			for(auto p2=pref.begin(); p2!=p1; ++p2){
				if(p2->second == obtState::OBT_STATE_RED && get_row(p1->first) == get_row(p2->first)){
						stato_inDFAfinale = false;
						break;
				}
			}

			if(stato_inDFAfinale){
				// Indice dello stato nella tabella di TRANSIZIONE e tipo di stato (accenttante o meno)
				stati[p1->first][0]= count_state;
				stati[p1->first][1]= mq[p1->first];
				++count_state;
			}
		}
	}


	dfa* dfaOBT = new dfa(count_state, target->get_dim_alphabet(), target->get_alphabet(), 0);
	int** dfaOBTtable = dfaOBT->get_ttable();


	// STAMPO elenco stati rossi
	#ifdef DEBUG_2
	cout << "--------------------------" << endl;
	cout << "ELENCO STATI ROSSI (automa finale): " << endl;
	cout << "--------------------------" << endl;

	for(auto p1=stati.begin(); p1!=stati.end(); ++p1)						// Ciclo sugli stati dell'automa finale
		cout << p1->first <<";" << endl;

	cout << "FINE" << endl;
	#endif


	// Calcolo la matrice delle transizioni
	// Ciclo sugli stati dell'automa finale
	for(auto p1=stati.begin(); p1!=stati.end(); ++p1){
		// Stringa label dello stato di partenza
		vector<SYMBOL> stato_partenza = p1->first;
		int indice_stato_partenza = stati[stato_partenza][0];

		for(int i=0; i<dim_alfabeto; ++i){
			// Stringa label dello stato di arrivo
			vector<SYMBOL> stato_arrivo = append_vectors(&stato_partenza, alfabeto + i);

			// Poichè lo stato di arrivo potrebbe essere una row() tra i blue,
			// devo cercare quale sia il corrispondente stato tra i rossi;
			// quindi fissata la label "stato di arrivo" trovo la label
			// tra i rossi che ha la stessa row, e la uso per estrarre l'indice corrispondente
			// nella matrice del dfa finale
			int indice_stato_arrivo = 0;
			// Trovo la riga dello stato di arrivo e cerco tra i RED il primo che ha la stessa riga
			map<vector<SYMBOL>, int[2]>::const_iterator p2=stati.begin();
			for(p2=stati.begin(); p2!=stati.end(); ++p2){

				vector<SYMBOL> p2_pref = p2->first;
				if(get_row(stato_arrivo) == get_row(p2_pref)){
					indice_stato_arrivo = stati[p2_pref][0];
					break;
				}
			}

			// piero: cambiare in assert
			// DA CAMBIARE
			if(p2 == stati.end()){
				cerr<<"Errore, fermate il mondo: indice_stato _arrivo non assegnato!" <<endl;
				return(NULL);
			}

			dfaOBTtable[indice_stato_partenza][i] = indice_stato_arrivo;
		}
		// Setto il tipo di stato: accettante o meno
		dfaOBTtable[indice_stato_partenza][dim_alfabeto] = stati[stato_partenza][1];
	}

	return dfaOBT;
}


void gi::lstar::update_from_counterexample(vector<SYMBOL> witness)
{
	// *** Aggiungo il Controesempio alla tabella di osservazione ***
	#ifdef DEBUG_2
	cout << "Aggiorno la tabella di osservazione attraverso il controesempio"<< endl;
	#endif

	for(auto i=witness.begin() + 1; i!=witness.end(); ++i){
		vector<SYMBOL> pref_witness = prefix_vector(&witness, i);

		// Aggiungo tutti i prefissi del contr.ex. come RED
		pref[pref_witness] = obtState::OBT_STATE_RED;

		// Chiudo i buchi con MQ
		for(auto k=experiments.begin(); k!=experiments.end(); ++k)
			if(mq.count(pref_witness) == 0)
				mq[pref_witness] = target->membership_query( append_vectors(&pref_witness, &(*k)));

		for(int j=0; j<dim_alfabeto; ++j){
			vector<SYMBOL> pref_witness_alf = append_vectors(&pref_witness,  alfabeto +j);

			// Se prefisso+lettera_alfabeto non è presente, la metto tra i blu
			if(pref.count( pref_witness_alf ) == 0){
				// Faccio la chiusura nei prefissi/suffissi (con ogni lettera dell'alfabeto)
				pref[pref_witness_alf] = obtState::OBT_STATE_BLU;

				// Riempo i buchi con MQ
				for(auto k=experiments.begin(); k!=experiments.end(); ++k){
					vector<SYMBOL> pref_witness_alf_exp = append_vectors(&pref_witness_alf,  &(*k));
					if(mq.count( pref_witness_alf_exp ) == 0)
						mq[pref_witness_alf_exp] = target->membership_query( pref_witness_alf_exp );
				}
			}
		}
	}
}

void gi::lstar::print_obt()
{
	cout << "---------------------------------------------";
	cout << endl<< "TABELLA DI OSSERVAZIONE - Elenco righe totali" << endl;
	cout << "---------------------------------------------" << endl;

	cout << "Exp: ";
	for(auto i=experiments.begin(); i!=experiments.end(); ++i)
		cout << (*i) <<", ";
	cout << endl;
	// Ciclo sugli stati dell'automa finale
	for(auto p1=pref.begin(); p1!=pref.end(); ++p1){
		if((int)p1->second == 1)
		{
			cout << "R ";
			vector<SYMBOL> tmpVect = p1->first;
			cout << tmpVect<<" | ";

			for(auto i=experiments.begin(); i<experiments.end(); ++i){

				cout << (int)mq[append_vectors(&tmpVect, &(*i))] <<" ";
			}

			cout << endl;
		}
	}
	cout << endl;

	// Ciclo sugli stati dell'automa finale
	for(auto p1=pref.begin(); p1!=pref.end(); ++p1){
		if((int)p1->second == 0)
		{
			cout << "B ";
			vector<SYMBOL> tmpVect = p1->first;
			cout << tmpVect<<" | ";

			for(auto i=experiments.begin(); i<experiments.end(); ++i){
				cout << (int)mq[append_vectors(&tmpVect, &(*i))] <<" ";
			}

			cout << endl;
		}
	}

	cout << endl<< "---------------------------------------------" << endl;
}

vector<obtState> gi::lstar::get_row(vector<SYMBOL> prefix){
	vector<obtState> row;

	//vector<vector<SYMBOL>>::const_iterator last_index;
	for(auto last_index=experiments.begin(); last_index!=experiments.end(); ++last_index){
		vector<SYMBOL> tmpVect = prefix;
		tmpVect.insert(tmpVect.end(), last_index->begin(), last_index->end());

		row.push_back((mq[tmpVect])?obtState::OBT_STATE_RED:obtState::OBT_STATE_BLU);
	}

	return row;
}

gi::dfa* gi::lstar::run_lstar(bool approximate, string samplestestpath){
	dfa * dfahp = NULL;
	int count_generation = 0;


	// Minimize tgdfa
	dfa * t1 = target->minimize_TF();

	if(target != NULL)
		delete target;

	target = t1;

	// Observation table
	lstar* obtable = new lstar(*target);

	dfa* dfatmp = NULL;

	while(1){
		#ifdef DEBUG_1
		cout << endl << endl << "********************************" << endl;
		cout << count_generation << " GENERATION"<< endl;
		cout << "********************************" << endl;
		#endif

		// Versione lunga del codice delle prossime righe:
		// Per la prima volta controllo che la tabella sia chiusa e consistente;
		// Tutte le altre volte rendo la tabella chiusa e consistenete solo se è stata precedentemente
		// cambiata dalle funzioni obtable->close_obt() e obtable->make_obt_consistent().

		//ciclo fino a quando la tabella non è chiusa e consistente
		//while( obtable->close_obt() || obtable->make_obt_consistent());
		bool modified=false, modified2=true, firstTime=true;

		while(modified2){
			// Rendo la tabella chiusa e memorizzo se l'ho cambiata almeno una volta
			modified=false;
			while(obtable->close_obt())
				modified=true;

			// Rendo la tabella consitente e memorizzo se l'ho cambiata almeno una volta
			// lo faccio solo se la tabella è cambiata dall'ultimo ciclo, ovvero se è stata chiusa,
			// oppure se è il primo ciclo
			modified2=false;
			while((modified || firstTime) && obtable->make_obt_consistent())
				modified2=true;

			firstTime = false;
		}



		// Print transition table
		#ifdef DEBUG_2
		obtable->print_obt();
		#endif

		// DFA from transition table
		if(dfatmp!= NULL)
			delete dfatmp;

		dfatmp = obtable->obt_to_dfa();

		#ifdef DEBUG_2
		dfatmp->print_dfa_ttable("DFA_HP dalla TABELLA DI OSSERVAZIONE");
		#endif

		++count_generation;


		// Generate witness if necessary
		vector<SYMBOL> witness;
		//if(!approximate)
			witness = target->equivalence_query(dfatmp);
		//else
		//	witness = target->equivalence_query_approximate(dfatmp,samplestestpath);

		/*STAMPA DOT DEI FARI PASSAGGI if(!approximate){
			string nome = "W"+dfa::intTostring(count_generation);
			//merged->print_dfa_with_color(nome);
			string percorso = "/Users/Gabriele/Dropbox/Workspace_eclipse/lstar/lstar/src/intermezzo2/"+nome+".dot";
			dfatmp->print_dfa_dot(nome, percorso.c_str());
		}*/

		// Find automaton when witness is 'e', else update the obb. table
		//TODO: problema con "e": come lo si sostituisce?
		// In pratica, a quanto ho capito negli esempi non potrtanno esserci mai stringhe vuote, per cui il problema
		// dell'approssimazione non si pone. Non mi pare corretto, che però, non si possano avere come campioni positivi stringhe nulle
		/*if(witness.compare("e")){
			obtable->update_from_counterexample(witness);
		}*/
		if(0 != witness.size()){
					obtable->update_from_counterexample(witness);
		}else{
			cout << "LSTAR: Automa trovato!"<<endl;
			dfahp = dfatmp->minimize_TF();

			break;
		}

	}

	if(dfatmp != NULL)
		delete dfatmp;

	//dfahp->set_n_memb_query(tgdfaL->get_n_memb_query());

	delete obtable;

	//delete tgdfaL;

	//dfa dfahp2(*dfahp);

	return dfahp;
}



std::ostream& operator<<(std::ostream& os, obtState s){
    return os << static_cast<char>(s) + '0';
}

