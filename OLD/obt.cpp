/*
 * obt.cpp
 *
 *  Created on: 14/mag/2015
 *      Author: Gabriele
 */

#include "obt.h"
#include <vector>


obt::obt(dfa* targetdfa, /*char* symbol,*/ char dim_alf, string* alf)
{
	target =  targetdfa;

	// Inizializzo
	pref[""] = 1;										// Inizializzo lo stato lambda a RED

	for(int i=0; i<targetdfa->get_dim_alp(); ++i)
		pref[ alf[i] ] = 0;								// Inizializzo ogni membro dell'alfabeto a BLU


	exp.push_back("");									// Inizializzo a "lambda" exp

	mq[""] = target->membership_query("");
	for(int i=0; i<targetdfa->get_dim_alp(); ++i)
		mq[alf[i]] = target->membership_query(alf[i]);

	dim_alfabeto = dim_alf;
	// DA_ELIMINARE symbols = symbol;
	alfabeto = alf;
}


bool obt::closed()
{

	// *** CHIUSURA ***
	// Per ogni BLU deve esistere una riga uguale tra i RED
	#ifdef DEBUG_2
	cout << "--------------------------" << endl;
	cout << "Verifico la CHIUSURA" << endl;
	cout << "--------------------------" << endl;
	#endif

	bool closed = true;
	for(It p1=pref.begin(); p1!=pref.end(); ++p1){
		if((*p1).second == 0)											// Ciclo solo sui BLU
		{
			closed = false;

			string prefB = (*p1).first;
			#ifdef DEBUG_3
			cout << "BLU da verificare: " << prefB << endl;
			#endif

			for(It p2=pref.begin(); p2!=pref.end(); ++p2){				// Ciclo tra i ROSSI
				if((*p2).second == 1)
				{
					string prefR = (*p2).first;

					//V cout << "RED confrontato: " << prefR << endl;

					bool uguale = true;
					for(int i=0; i<exp.size(); ++i)						// Verifico che le righe siano uguali
						if( mq[prefB+exp[i]] != mq[prefR+exp[i]] )
							uguale = false;								// Potrei agire giˆ qui dentro, per chiarezza faccio fuori

					if(uguale){
						//V cout << "Riga blu trovata tra i red"<< endl;
						closed = true;
						break;
					}
				}
			}

			if(!closed)
			{
				#ifdef DEBUG_2
				cout << "Tabella NON chiusa" << endl;
				cout << "Promuovo a RED "<<prefB << endl;
				#endif

				// Promuovo a RED il BLU in questione
				pref[prefB] = 1;

				// Aggiungo ai BLU il nuovo RED concatenato ad ogni simbolo dell'alfabeto
				for(int i=0; i<dim_alfabeto; ++i)
				{
					//V cout << "PREFB: "<<prefB<<endl;
					//V cout << "Alfabeto["<<i<<"]: "<<alfabeto[i]<<endl;
					//V cout << "Aggiungo "<< prefB+alfabeto[i] << " ai blu" << endl;
					pref[prefB+alfabeto[i]] = 0;
					//V cout << "exp.size: " << exp.size() << endl;
					for(int j=0; j<exp.size(); ++j){
						mq[prefB+alfabeto[i]+exp[j]] = target->membership_query( prefB+alfabeto[i]+exp[j] );
						//V cout << "    T["<<prefB+alfabeto[i]+exp[j] <<"]:"<<mq[prefB+alfabeto[i]+exp[j]]<<endl;
					}
				}
			}
		}
	}

	#ifdef DEBUG_2
	cout << "--------------------------" << endl;
	cout << "Fine CHIUSURA" << endl;
	cout << "--------------------------" << endl;
	#endif

	return closed;
}

bool obt::consistent()
{
	// *** CONSISTENZA ***
	// 2 stati uguali a paritˆ di input devono raggiungere lo stesso stato

	#ifdef DEBUG_2
	cout << "--------------------------" << endl;
	cout << "Verifico la CONSISTENZA"<< endl;
	cout << "--------------------------" << endl;
	#endif

	bool consistent = true;
	for(It p1=pref.begin(); p1!=pref.end(); ++p1)
	{
		if((*p1).second == 1){											// Ciclo sui ROSSI: s1

			string prefR1 = (*p1).first;
			//cout << "R1: "<<prefR1<<endl;

			for(It p2=pref.begin(); p2!=pref.end(); ++p2)				// Ciclo sui ROSSI: s2
			{
				consistent = true;
				if((*p2).second == 1 && (*p2).first.compare(prefR1)){	// Verifico di lavorare con 2 red distinti (NB: "compare" qui restituisce T se sono diversi)

					string prefR2 = (*p2).first;
					//cout << "R2: "<<prefR2<<endl;

					bool uguale = true;									// Verifico che le righe siano uguali
					for(int i=0; i<exp.size(); ++i)
						if( mq[prefR1+exp[i]] != mq[prefR2+exp[i]] ){
							uguale = false;								// (Potrei agire giˆ qui dentro, per chiarezza faccio fuori)
							//V cout << "Righe diverse, no problem for consistence" << endl;
						}

					if(uguale)											// Se sono raw uguali, verifico che non esistono inconsistenze
					{
						//V cout << "Hanno righe uguali, quindi potenziali generatori di inconsistenza"<< endl;
						for(int i=0; i<dim_alfabeto; ++i)
						{
							string prefR1A = prefR1+alfabeto[i];
							string prefR2A = prefR2+alfabeto[i];
							//V cout << "- Controllo "<<prefR1A<<" e "<<prefR2A<<" -"<<endl;

							string expc = "";
							for(int j=0; j<exp.size(); ++j){							// Verifico che le righe siano uguali
								if( mq[ prefR1A+exp[j] ] != mq[ prefR2A+exp[j] ] ){
									#ifdef DEBUG_2
									cout << "INCONSISTENZA in " << prefR1A+exp[j] << " e " << prefR2A+exp[j] << endl;
									#endif
									consistent = false;									// Potrei agire giˆ qui dentro, per chiarezza faccio fuori
									expc = exp[j];										// Salvo da parte l'esperimento per cui nasce la differenza
									break;
								}
							}

							if(!consistent)												// Se NON  consistente, aggiorno EXP e la tabella di conseguenza
							{
								exp.push_back(alfabeto[i]+expc);
								string new_exp = alfabeto[i]+expc;

								#ifdef DEBUG_2
								cout << "Aggiungo " << new_exp << " agli exp"<<endl;
								#endif

								for(It p3=pref.begin(); p3!=pref.end(); ++p3){			// Riempo i buchi generati dall'inserimento della nuova colonna degli exp
									mq[(*p3).first + new_exp] = target->membership_query( (*p3).first + new_exp );
									//V cout << "MQ["<< (*p3).first + new_exp << "]: "<< mq[(*p3).first + new_exp] << endl;
								}

								#ifdef DEBUG_2
								cout << "--------------------------" << endl;
								cout << "Fine CONSISTENZA"<< endl;
								cout << "--------------------------" << endl;
								#endif

								//break;
								return consistent;	// Ritorno senza testare ora le altre righe, secondo De La Higuera (se mettro break testo subito)
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

	return consistent;
}


dfa* obt::obt_to_dfa()
{
	// TABELLA DI OSSERVAZIONE --> DFA
	//V cout << endl<< "*** From OBB TABLE to TRANSITION TABLE ***" << endl;
	// Creo un elenco degli stati dell'automa finale con la map "stati"
	map<string, int*> stati;											// ELENCO DEGLI STATI: Label dello stato & 2 interi: il primo per indice numerico dello stato, secondo per accettante o meno
	typedef	map<string, int*>::const_iterator It2;

	int count_state = 0;

	// Qui faccio il controllo delle transizioni guardando le righe della tabella di oss.
	for(It p1=pref.begin(); p1!=pref.end(); ++p1)						// Ciclo sui ROSSI: s1
	{
		if((*p1).second == 1){
			string prefR1 = (*p1).first;
			//V cout << "INITprefR1: "<<prefR1 << ", " << endl;

			// Ipotizzo che sia uno stato da aggiungere (stato_inDFAfinale=T). Ciclando internamente "uguale" deve diventare sempre F alla fine, se rimane T
			// c' uno stato giˆ analizzato uguale a quello preso in considerazione quindi non deve appartenere e stato_inDFAfinale=F
			bool stato_inDFAfinale=true;
			for(It p2=pref.begin(); p2!=pref.end(); ++p2)				// Ciclo sui ROSSI: s2
			{
				if((*p2).second ==1)
				{
					if( !(*p2).first.compare(prefR1) ){ 				// Quando arriva a se stesso si ferma, devo controllare solo i RED "predecessori"
						break;
					}
					else												// Controllo solo tra predecessori RED
					{
						string prefR2 = (*p2).first;

						//V cout << "Pref1: "<<prefR1<<", Pref2: "<<prefR2<<";"<<endl;
						stato_inDFAfinale=false;
						for(int i=0; i<exp.size(); ++i)
							if( mq[prefR1+exp[i]] != mq[prefR2+exp[i]] ){
								stato_inDFAfinale = true;									// (Potrei agire giˆ qui dentro, per chiarezza faccio fuori)
								break;
							}
					}
				}
			}

			if(stato_inDFAfinale){
				//V cout << "Aggiungo lo stato: "<<prefR1 << endl;
				stati[prefR1] = new int[2];
				stati[prefR1][0] = count_state;								// Indice dello stato nella tabella di TRANSIZIONE
				stati[prefR1][1] = mq[prefR1];								// Stato accettante o meno
				//V cout << " -> Indice: "<< stati[prefR1][0] <<" ; Accettante: "<< stati[prefR1][1] << endl;
				count_state++;
			}
		}
	}


	dfa* dfaOBT = new dfa(count_state, /*symbols,*/ dim_alfabeto, 0);
	int** dfaOBTtable = dfaOBT->get_ttable();


	// STAMPO elenco stati rossi
	#ifdef DEBUG_2
	cout << "--------------------------" << endl;
	cout << "ELENCO STATI ROSSI (automa finale): " << endl;
	cout << "--------------------------" << endl;

	for(It2 p1=stati.begin(); p1!=stati.end(); ++p1)						// Ciclo sugli stati dell'automa finale
		cout << (*p1).first <<";" << endl;

	cout << "FINE" << endl;
	#endif


	// Calcolo la matrice delle transizioni
	for(It2 p1=stati.begin(); p1!=stati.end(); ++p1)						// Ciclo sugli stati dell'automa finale
	{
		for(int i=0; i<dim_alfabeto; ++i){
			string stato_partenza = (*p1).first;							// Stringa label dello stato di partenza
			string stato_arrivo = (*p1).first + alfabeto[i];				// Stirnga label dello stato di arrivo

			//cout << "stato partenza: "<<stato_partenza<<endl;
			int indice_stato_partenza = stati[stato_partenza][0];
			//cout << "indice stato partenza: "<< indice_stato_partenza << endl;

			// PoichŽ lo stato di arrivo potrebbe essere una raw() tra i blue, devo cercare quale sia il corrispondete stato tra i rossi
			// quindi fissata la label "stato di arrivo" trovo la label tra i rossi che ha la stessa raw, e la uso per estrarre l'indice corrispondente
			// nella matrice del dfa finale
			int indice_stato_arrivo = -1;
			for(It2 p2=stati.begin(); p2!=stati.end(); ++p2)				// Trovo la riga dello stato di arrivo e cerco tra i RED il primo che ha la stessa riga
			{
				bool uguale = true;
				for(int i=0; i<exp.size(); ++i)
					if( mq[ stato_arrivo+exp[i] ] != mq[ (*p2).first+exp[i] ] ){
						uguale = false;
						break;
					}

				if(uguale){
					//cout << "uguale per "<< p2->first << endl;
					indice_stato_arrivo = stati[(*p2).first][0];
					break;
				}
			}
			
			// piero: non dovresti controllare se indice_stato_arrivo Ã¨ stato asseganto? 
			// DA CAMBIARE
			if(indice_stato_arrivo == -1){
				cout<<"Errore, fermate il mondo: indice_stato _arrivo non assegnato!" <<endl;
				return(NULL);
			}
				

			//cout << "stato di arrivo: "<< stato_arrivo <<";"<< endl;
			//cout << "indice stato arrivo: "<< indice_stato_arrivo <<endl << endl;

			int accettante = stati[stato_partenza][1];								// 0-> non accettante; 1-> accettante

			dfaOBTtable[indice_stato_partenza][i] = indice_stato_arrivo;
			dfaOBTtable[indice_stato_partenza][dim_alfabeto] = accettante;			//(Potresti metterlo fuori dal ciclo)
		}
	}

	return dfaOBT;
}


void obt::update_from_counterexample(string witness)
{
	// *** Aggiungo il Controesempio alla tabella di osservazione ***
	#ifdef DEBUG_2
	cout << "Aggiorno la tabella di osservazione attraveros il controesempio"<< endl;
	#endif

	for(int i=1; i<=witness.length(); ++i)
	{
		pref[witness.substr(0,i)] = 1;										// Aggiungo tutti i prefissi del contr.ex. come RED
		for(int k=0; k<exp.size(); ++k)										// Chiudo i buchi con MQ
			if(mq.count(witness.substr(0,i)) == 0)
				mq[witness.substr(0,i)] = target->membership_query( witness.substr(0,i) + exp[k] );

		for(int j=0; j<dim_alfabeto; ++j)
		{
			if(pref.count( witness.substr(0,i)+alfabeto[j] ) == 0){			// Se prefisso+lettera_alfabeto non  presente, la metto tra i blu
				pref[witness.substr(0,i) + alfabeto[j]] = 0;				// Faccio la chiusura nei prefissi/suffissi (con ogni lettera dell'alfabeto)

				// Riempo i buchi con MQ
				for(int k=0; k<exp.size(); ++k)
					if(mq.count( witness.substr(0,i) + alfabeto[j] + exp[k] ) == 0)
						mq[witness.substr(0,i) + alfabeto[j] + exp[k]] = target->membership_query( witness.substr(0,i) + alfabeto[j] + exp[k] );
			}
		}
	}
}

void obt::print_obt()
{
	cout << "---------------------------------------------";
	cout << endl<< "TABELLA DI OSSERVAZIONE - Elenco righe totali" << endl;
	cout << "---------------------------------------------" << endl;

	cout << "Exp: ";
	for(int i=0; i<exp.size(); ++i)
		cout << exp[i] <<", ";
	cout << endl;

	for(It p1=pref.begin(); p1!=pref.end(); ++p1){						// Ciclo sugli stati dell'automa finale
		if((int)(*p1).second == 1)
		{
			cout << "R ";

			cout << (*p1).first<<" | ";

			for(int i=0; i<exp.size(); ++i)
				cout << (int)mq[(*p1).first+exp[i]] <<" ";

			cout << endl;
		}
	}
	cout << endl;

	for(It p1=pref.begin(); p1!=pref.end(); ++p1){						// Ciclo sugli stati dell'automa finale
		if((int)(*p1).second == 0)
		{
			cout << "B ";

			cout << (*p1).first<<" | ";

			for(int i=0; i<exp.size(); ++i)
				cout << (int)mq[(*p1).first+exp[i]] <<" ";

			cout << endl;
		}
	}

	cout << endl<< "---------------------------------------------" << endl;
}
