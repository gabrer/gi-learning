/*
 * minimizzazione.cpp
 *
 *  Created on: 27/apr/2015
 *      Author: Gabriele
 */

#include "minimize.h"
#include <string>
#include <iostream>
#include <limits>
using namespace std;

#define ND numeric_limits<int>::max()

// Invece di una STRUCT puoi crea una CLASSE DFA, con i metodi "minimizza" e "stampa",
//  per gli stati accettanti puoi lasciare una collonna in più nella matrice DFA e aggiungere inoltre un vettore con gli stati accettanti


// STAMPA IL DFA
// Nel numero di colonne è inclusa la colonna finale di stato accettante o meno (accettante -> 1)
void stampa_dfa(int** dfa, int num_righe, int num_colonne){

	cout << endl << "STAMPO DFA"<< endl;
	cout << "    0 | 1 - A"<< endl;
	for(int i=0; i<num_righe; ++i){
		cout << "S"<<i<<"  ";

		for(int j=0; j<num_colonne; ++j)
			cout << dfa[i][j] << "   ";
		cout << endl;
	}

	return;
}

my_dfa* minimize(int** dfa, const int dim_alfabeto, int n_stati)
{
	//const int dim_alfabeto = 2;
	int count_state = n_stati;												// EX: 6 stati

	// MATRICE DFA
	int num_righe = count_state;										// !? VERIFICA CHE NON CI SIA BISOGNO DI -1
	int num_colonne = dim_alfabeto + 1;									//1 finale è la colonna relativa al tipo
	/*int** dfa = new int*[num_righe];
	for(int i=0; i<num_righe; ++i)
		dfa[i] = new int[num_colonne];*/

	// DFA Target (dalla Angluin)
	/*
	dfa[0][0]=1;
	dfa[0][1]=2;
	dfa[0][2]=1;
	dfa[1][0]=0;
	dfa[1][1]=3;
	dfa[1][2]=0;
	dfa[2][0]=3;
	dfa[2][1]=0;
	dfa[2][2]=0;
	dfa[3][0]=2;
	dfa[3][1]=1;
	dfa[3][2]=0;*/

	// Equivalenza di 2 DFA (da Hopcraft)
	/*dfa[0][0]=0;
	dfa[0][1]=1;
	dfa[0][2]=1;
	dfa[1][0]=0;
	dfa[1][1]=1;
	dfa[1][2]=0;

	dfa[2][0]=3;
	dfa[2][1]=4;
	dfa[2][2]=1;
	dfa[3][0]=3;
	dfa[3][1]=4;
	dfa[3][2]=1;
	dfa[4][0]=3;
	dfa[4][1]=4;
	dfa[4][2]=0;*/

	//DFA 3 stati (http://www.informatik.uni-bremen.de/agbs/lehre/ss05/pi2/hintergrund/minimize_dfa.pdf)
	// 1
	/*dfa[0][0]=1;
	dfa[0][1]=0;
	dfa[0][2]=0;
	dfa[1][0]=0;
	dfa[1][1]=1;
	dfa[1][2]=1;
	dfa[2][0]=1;
	dfa[2][1]=0;
	dfa[2][2]=0;*/

	//2
	/*dfa[0][0]=1;
	dfa[0][1]=0;
	dfa[0][2]=0;
	dfa[1][0]=2;
	dfa[1][1]=1;
	dfa[1][2]=1;
	dfa[2][0]=1;
	dfa[2][1]=0;
	dfa[2][2]=0;*/

	//DFA 8 stati da internet (http://www.cs.engr.uky.edu/~lewis/essays/compilers/min-fa.html)
	/*dfa[0][0]=1;
	dfa[0][1]=4;
	dfa[0][2]=0;
	dfa[1][0]=5;
	dfa[1][1]=2;
	dfa[1][2]=0;
	dfa[2][0]=3;
	dfa[2][1]=6;
	dfa[2][2]=1;
	dfa[3][0]=3;
	dfa[3][1]=3;
	dfa[3][2]=0;
	dfa[4][0]=1;
	dfa[4][1]=4;
	dfa[4][2]=0;
	dfa[5][0]=1;
	dfa[5][1]=4;
	dfa[5][2]=0;
	dfa[6][0]=3;
	dfa[6][1]=7;
	dfa[6][2]=0;
	dfa[7][0]=3;
	dfa[7][1]=6;
	dfa[7][2]=1;*/

	/*DFA 8 stati da internet già minimizzato!
	dfa[0][0]=1;
	dfa[0][1]=0;
	dfa[0][2]=0;

	dfa[1][0]=0;
	dfa[1][1]=2;
	dfa[1][2]=0;

	dfa[2][0]=3;
	dfa[2][1]=4;
	dfa[2][2]=1;

	dfa[3][0]=3;
	dfa[3][1]=3;
	dfa[3][2]=0;

	dfa[4][0]=3;
	dfa[4][1]=2;
	dfa[4][2]=0;
	*/

	// DFA ad 8 stati (Hopcraft libro)
	/*dfa[0][0]=1;
	dfa[0][1]=5;
	dfa[0][2]=0;
	dfa[1][0]=6;
	dfa[1][1]=2;
	dfa[1][2]=0;
	dfa[2][0]=0;
	dfa[2][1]=2;
	dfa[2][2]=1;
	dfa[3][0]=2;
	dfa[3][1]=6;
	dfa[3][2]=0;
	dfa[4][0]=7;
	dfa[4][1]=5;
	dfa[4][2]=0;
	dfa[5][0]=2;
	dfa[5][1]=6;
	dfa[5][2]=0;
	dfa[6][0]=6;
	dfa[6][1]=4;
	dfa[6][2]=0;
	dfa[7][0]=6;
	dfa[7][1]=2;
	dfa[7][2]=0;*/

	// DFA a 6 stati  (http://www.cs.nott.ac.uk/~txa/g52mal/LectureNotes/g52mal-notes.pdf)
	/* dfa[0][0]=1;
	dfa[0][1]=4;
	dfa[0][2]=0;
	dfa[1][0]=2;
	dfa[1][1]=3;
	dfa[1][2]=0;
	dfa[2][0]=2;
	dfa[2][1]=2;
	dfa[2][2]=1;
	dfa[3][0]=2;
	dfa[3][1]=3;
	dfa[3][2]=1;
	dfa[4][0]=5;
	dfa[4][1]=4;
	dfa[4][2]=0;
	dfa[5][0]=5;
	dfa[5][1]=4;
	dfa[5][2]=0;*/


	// *** MINIMIZZAZIONE DFA - Table-filling algorithm ***
		bool** distinti = new bool*[count_state-1];								// Tabella delle coppie di stati distinti
		for(int i=0; i<num_righe; ++i)											// EX: 5 righe e 6 colonne (di queste 6 solo 5 usate effettivamente, ma serve)
			distinti[i] = new bool[count_state];

		// Inizializzo
		for(int i=0; i<count_state-1; ++i)
			for(int j=0; j<count_state; ++j)
				distinti[i][j]=false;

		// Distinguo tra accettanti e non accettanti							// EX: 0 <= i <= 4, j=i+1 (1 <= j <= 5 per la prima iterazione)
		for(int i=0; i<(count_state-1); ++i)
			for(int j=i+1; j<count_state; ++j)
				if(dfa[i][dim_alfabeto] != dfa[j][dim_alfabeto])				// Se uno è accettante mentre l'altro no
					distinti[i][j] = true;

		// Loop  minimizzante
		bool modificata = true;
		while(modificata)
		{
			modificata = false;

			for(int i=0; i<(count_state-1); ++i){
				for(int j=i+1; j<count_state; ++j){
					if(!distinti[i][j]){

						for(int k=0; k<dim_alfabeto; ++k)
						{
							int stato_arrivo_1 = dfa[i][k];
							int stato_arrivo_2 = dfa[j][k];

							if(stato_arrivo_1 == stato_arrivo_2)
								continue;

							if(stato_arrivo_2 < stato_arrivo_1){				// Facendo la tabella delle coppie di stati visitati noto che sempre j>i
								int tmp = stato_arrivo_1;
								stato_arrivo_1 = stato_arrivo_2;
								stato_arrivo_2 = tmp;
							}

							if(distinti[stato_arrivo_1][stato_arrivo_2]){		// Se la coppia di arrivo è distinta, lo è anche quella originaria
								distinti[i][j] = true;
								modificata = true;
							}
						}
					}
				}
			}
		}

		// Stampo la tabella degli stati distinti
		cout << "TABELLA DEI DISTINTI " << endl;
		for(int i=0; i<(count_state-1); ++i){
			for(int j=i+1; j<count_state; ++j){
				cout << "("<<i<<","<<j<<"):" << distinti[i][j] << "  ";
			}
			cout << endl;
		}


		// ** Mi preparo a costruire l'automa minimizzato **
		// Creo un elenco per avere una lista delle coppie di stati equivalenti più direttamente accessibile,
		// 	se è ND allora non ha un equivalente.
		int stato_equivalente[count_state];
		for(int i=0; i<count_state; ++i)
			stato_equivalente[i] = ND;

		// Conto il numero di stati finali, ed assegno ad ogni stato il suo eventuale stato equivalente
		// (devo anche controllare nel caso ho (0,4) (0,5) che quando incontro (4,5) non lo conto come un ulteriore stato in meno)
		int stati_finali = count_state;
		for(int i=0; i<(count_state-1); ++i)
			for(int j=i+1; j<count_state; ++j)
				if(!distinti[i][j] && stato_equivalente[i] == ND && stato_equivalente[j] == ND){
					stato_equivalente[j] = i;
					stati_finali--;
				}

		// Stampo le informazioni riguardo agli stati equivalenti trovati o meno
		cout << "N di stati finali: " << stati_finali << endl;
		cout << "Equivalenze stati: " << endl;
		for(int i=0; i<count_state; ++i)
			if(stato_equivalente[i] != ND)
				cout << "S:"<<i<<" --> "<<stato_equivalente[i]<<endl;
			else
				cout << "S:"<<i<<endl;

		// STAMPO il vecchio DFA
		//stampa_dfa(dfa, num_righe, num_colonne);

		// Istanzio il nuovo DFA minimizzato
		int** dfa_min = new int*[stati_finali];
		for(int i=0; i<stati_finali; ++i)
			dfa_min[i] = new int[num_colonne];

		int count = 0;
		for(int i=0; i<count_state; ++i){
			if(stato_equivalente[i] == ND){
				for(int j=0; j<dim_alfabeto+1; ++j)
					dfa_min[count][j]=dfa[i][j];
				count++;
			}
		}

		// Aggiorno le transizioni verso stati ormai scomparsi perché sostituiti con gli equivalenti
		int equivalenze_finora= 0;
		for(int i=0; i<count_state; ++i){
			if(stato_equivalente[i] != ND)
			{
				for(int k=0; k<stati_finali; ++k)
					for(int t=0; t<dim_alfabeto; ++t)
						if(dfa_min[k][t] == i)
							dfa_min[k][t] = stato_equivalente[i];		// Sostituisco la transizione allo stato "i" con lo stato equivalente "stato_equivalente[i]"
			}
		}

		// Qui aggiorno la label, perché se ad esempio ho collassato 2 stati prima dello stato 6,
		//  adesso lo stato 6 si trova nella riga 4, però le transizioni sono rimaste
		//  verso 6 e devono essere sistemate
		for(int i=0; i<count_state; ++i)
		{
			if(stato_equivalente[i] != ND)
				equivalenze_finora++;

			if(equivalenze_finora != 0)
			{
				int nuova_label = i-equivalenze_finora;
				for(int k=0; k<stati_finali; ++k)
					for(int t=0; t<dim_alfabeto; ++t)
						if(dfa_min[k][t] == i)
							dfa_min[k][t] = nuova_label;
			}
		}

		// STAMPO il dfa minimizzato
		cout << endl << "DFA MINIMIZZATO FINALE";
		//stampa_dfa(dfa_min, stati_finali, num_colonne);

		my_dfa* mindfa = new my_dfa();
		mindfa->dim_alfabeto = dim_alfabeto;
		mindfa->num_colonne = num_colonne;
		mindfa->num_righe = stati_finali;
		mindfa->matrixdfa = dfa_min;

		return mindfa;
}



int* elenco_stati_equivalenti(int** dfa_hp, int n_stati_hp, int** dfa_T, int n_stati_T,  int dim_alfabeto)
{
	// MATRICE DFA UNIONE
	int count_state = n_stati_hp + n_stati_T;
	cout << "Numero di stati totali: "<<n_stati_hp<<" + "<<n_stati_T<<endl;

	int num_righe = count_state;
	int num_colonne = dim_alfabeto + 1;										//1 finale è la colonna relativa al tipo

	// DFA Unione di dfa_hp e dfa_T
	int** dfa = new int*[num_righe];
	for(int i=0; i<num_righe; ++i)
		dfa[i] = new int[num_colonne];

	// Riempo il DFA UNIONE
	int c = 0;

	for(int j=0; j<n_stati_T; ++j){									// Automa target
		for(int k=0; k<num_colonne; ++k)
			dfa[c][k] = dfa_T[j][k];
		c++;
	}

	int stato_partenza_dfa_hp = c;
	cout << "Stato di partanze del dfa_HP nel DFA unione è: "<<c<<endl;

	for(int j=0; j<n_stati_hp; ++j){										// Automa ipotesi
		for(int k=0; k<num_colonne; ++k){
			if(k != dim_alfabeto)
				dfa[c][k] = dfa_hp[j][k] + stato_partenza_dfa_hp;
			else
				dfa[c][k] = dfa_hp[j][k];
		}
		c++;
	}

	cout<< endl << " *** DFA UNIONE ***";
	stampa_dfa(dfa, num_righe, num_colonne);
	cout << endl <<endl;


	// *** TABLE-FILLING ALGORITHM ***
	bool** distinti = new bool*[count_state-1];								// Tabella delle coppie di stati distinti
	for(int i=0; i<num_righe; ++i)											// EX: 5 righe e 6 colonne (di queste 6 solo 5 usate effettivamente, ma serve)
		distinti[i] = new bool[count_state];

	// Inizializzo
	for(int i=0; i<count_state-1; ++i)
		for(int j=0; j<count_state; ++j)
			distinti[i][j]=false;

	// Distinguo tra accettanti e non accettanti							// EX: 0 <= i <= 4, j=i+1 (1 <= j <= 5 per la prima iterazione)
	for(int i=0; i<(count_state-1); ++i)
		for(int j=i+1; j<count_state; ++j)
			if(dfa[i][dim_alfabeto] != dfa[j][dim_alfabeto])				// Se uno è accettante mentre l'altro no
				distinti[i][j] = true;

	// Loop  minimizzante
	bool modificata = true;
	while(modificata)
	{
		modificata = false;

		for(int i=0; i<(count_state-1); ++i){
			for(int j=i+1; j<count_state; ++j){
				// cout << endl << "SP1: "<<i << ", SP2: "<<j;
				if(!distinti[i][j]){

					for(int k=0; k<dim_alfabeto; ++k)
					{
						int stato_arrivo_1 = dfa[i][k];
						int stato_arrivo_2 = dfa[j][k];

						// !? AGGIUNTA SISTEMA BUG -- DA VERIFICARE
						if(stato_arrivo_1 == stato_arrivo_2)				// Lo stato di arrivo letto nel dfa potrebbe essere una coppia del tipo (2,2)
							continue;

						if(stato_arrivo_2 < stato_arrivo_1){				// Facendo la tabella delle coppie di stati visitati noto che sempre j>i
							int tmp = stato_arrivo_1;
							stato_arrivo_1 = stato_arrivo_2;
							stato_arrivo_2 = tmp;
						}

						//cout << endl <<"SA1: "<<stato_arrivo_1 << ", SA2: "<<stato_arrivo_2 << " --> ";
						if(distinti[stato_arrivo_1][stato_arrivo_2]){		// Se la coppia di arrivo è distinta, lo è anche quella originaria
							cout << " Distinti!";
							distinti[i][j] = true;
							modificata = true;
						}
					}
				}
			}
		}
	}

	// Stampo la tabella degli stati distinti
	cout << endl << "TABELLA DEI DISTINTI " << endl;
	for(int i=0; i<(count_state-1); ++i){
		for(int j=i+1; j<count_state; ++j){
			cout << "("<<i<<","<<j<<"):" << distinti[i][j] << "  ";
		}
		cout << endl;
	}


	// Creo un ELENCO per avere una lista delle coppie di stati equivalenti più direttamente accessibile,
	// 	se è ND allora non ha un equivalente.
	int* stato_equivalente = new int[count_state];
	for(int i=0; i<count_state; ++i){
		stato_equivalente[i] = ND;
	}

	for(int i=0; i<(count_state-1); ++i){
		for(int j=i+1; j<count_state; ++j){
			if(!distinti[i][j])
				stato_equivalente[i] = j;
		}
	}

	cout << "Devono essere equivalenti 0 e "<< stato_partenza_dfa_hp <<endl;
	return stato_equivalente;
}


// Puoi pensare di includere la creazione della witness all'interno della funzione "elenco_stati_equivalenti" e ritornare una struct
// con l'elenco degli stati equivalenti e la witness relativa
string witness(int** dfa_hp, int n_stati_hp, int** dfa_T, int n_stati_T,  int dim_alfabeto)
{

	// !? PENSA A FARE UNA FUNZIONE CHE RESTITUISCE IL DFA UNIONE e mettere il codice fuori di qui
	// MATRICE DFA UNIONE
	int count_state = n_stati_hp + n_stati_T;
	cout << "Numero di stati totali: "<<n_stati_hp<<" + "<<n_stati_T<<endl;

	int num_righe = count_state;
	int num_colonne = dim_alfabeto + 1;										//1 finale è la colonna relativa al tipo

	// DFA Unione di dfa_hp e dfa_T
	int** dfa = new int*[num_righe];
	for(int i=0; i<num_righe; ++i)
		dfa[i] = new int[num_colonne];

	// Riempo il DFA UNIONE
	int c = 0;																// Conterrà l'indice dello stato iniziale del 2 dfa

	for(int j=0; j<n_stati_T; ++j){											// Automa target
		for(int k=0; k<num_colonne; ++k)
			dfa[c][k] = dfa_T[j][k];
		c++;
	}

	int stato_partenza_dfa_hp = c;
	cout << "Stato di partanze del dfa_HP nel DFA unione è: "<<c<<endl;

	for(int j=0; j<n_stati_hp; ++j){										// Automa ipotesi
		for(int k=0; k<num_colonne; ++k){
			if(k != dim_alfabeto)
				dfa[c][k] = dfa_hp[j][k] + stato_partenza_dfa_hp;			// Le transizioni devono essere all'indice aggiornato dello stato
			else
				dfa[c][k] = dfa_hp[j][k];
		}
		c++;
	}

	cout << endl << " *** DFA UNIONE ***";
	stampa_dfa(dfa, num_righe, num_colonne);
	cout << " ---------------------- " << endl <<endl;


	// *** TABLE-FILLING ALGORITHM with witness ***
	char** distinti = new char*[count_state-1];								// Tabella delle coppie di stati distinti
	for(int i=0; i<num_righe; ++i)											// EX: 5 righe e 6 colonne (di queste 6 solo 5 usate effettivamente, ma serve)
		distinti[i] = new char[count_state];

	// Inizializzo
	for(int i=0; i<count_state-1; ++i)
		for(int j=0; j<count_state; ++j)
			distinti[i][j]=' ';												// ' ' equivale alla cella vuota

	// Distinguo tra accettanti e non accettanti							// EX: 0 <= i <= 4, j=i+1 (1 <= j <= 5 per la prima iterazione)
	for(int i=0; i<(count_state-1); ++i)
		for(int j=i+1; j<count_state; ++j)
			if(dfa[i][dim_alfabeto] != dfa[j][dim_alfabeto])				// Se uno è accettante mentre l'altro no
				distinti[i][j] = 'x';

	// Loop  minimizzante
	cout << " *** Riempo la table degli stati ***" << endl;
	bool modificata = true;													// Tiene traccia se la tabella è stata modificata durante l'ultimo ciclo
	while(modificata)
	{
		modificata = false;

		for(int i=0; i<(count_state-1); ++i){
			for(int j=i+1; j<count_state; ++j){
				// cout << endl << "SP1: "<<i << ", SP2: "<<j;
				if(distinti[i][j] == ' '){

					for(int k=0; k<dim_alfabeto; ++k)
					{
						int stato_arrivo_1 = dfa[i][k];
						int stato_arrivo_2 = dfa[j][k];

						// !? AGGIUNTA SISTEMA BUG -- DA VERIFICARE
						if(stato_arrivo_1 == stato_arrivo_2)				// Lo stato di arrivo letto nel dfa potrebbe essere una coppia del tipo (2,2)
							continue;

						if(stato_arrivo_2 < stato_arrivo_1){				// Facendo la tabella delle coppie di stati visitati noto che sempre j>i
							int tmp = stato_arrivo_1;						// quindi le coppie devono avere sempre i<j
							stato_arrivo_1 = stato_arrivo_2;
							stato_arrivo_2 = tmp;
						}

						//cout << endl <<"SA1: "<<stato_arrivo_1 << ", SA2: "<<stato_arrivo_2 << " --> ";
						if(distinti[stato_arrivo_1][stato_arrivo_2] != ' '){		// Se la coppia di arrivo è distinta, lo è anche quella originaria
							cout << " Distinti!";
							if(!k)											// Entro quando è 0
								distinti[i][j] = '0';
							else											// Entro quando è 1
								distinti[i][j] = '1';
							modificata = true;
						}
					}
				}
			}
		}
	}

	// Stampo la tabella degli stati distinti
	cout << endl << "TABELLA DEI DISTINTI " << endl;
	for(int i=0; i<(count_state-1); ++i){
		for(int j=i+1; j<count_state; ++j){
			cout << "("<<i<<","<<j<<"):" << distinti[i][j] << "  ";
		}
		cout << endl;
	}


	// Creo un ELENCO per avere una lista delle coppie di stati equivalenti più direttamente accessibile,
	// 	se è ND allora non ha un equivalente.
	int* stato_equivalente = new int[count_state];
	for(int i=0; i<count_state; ++i){
		stato_equivalente[i] = ND;
	}

	for(int i=0; i<(count_state-1); ++i){
		for(int j=i+1; j<count_state; ++j){
			if(distinti[i][j] == ' ')
				stato_equivalente[i] = j;
		}
	}

	cout << "Devono essere equivalenti 0 e "<< stato_partenza_dfa_hp <<endl;

	// Se automi NON equivalenti, creo la witness
	string wit = "";
	if(distinti[0][stato_partenza_dfa_hp] != ' ')
	{
		int icoppia = 0;
		int jcoppia = stato_partenza_dfa_hp;
		cout << "Creo il CONTROESEMPIO: " << endl;
		while(1){
			char input = distinti[icoppia][jcoppia];
			cout << "stinga parziale: "<<wit<<" + " << input << endl;
			if(input == 'x')												// Caso in cui lo start state dei 2, uno è accettante e l'altro no
				break;
			wit = wit + input;
			cout << "Vale: "<< input << ", ascii: "<< (int)input << ", diff_i:" << (int)input - 48 << endl;
			cout << "icoppia: "<<icoppia<<endl;
			icoppia = dfa[icoppia][(int)input - 48];
			jcoppia = dfa[jcoppia][(int)input - 48];
			if(distinti[icoppia][jcoppia] == 'x'){
				cout << "fine stringa" << endl;
				break;
			}
		}

		cout << "Controesempio è: "<< wit << endl;
	}

	return wit;
}
