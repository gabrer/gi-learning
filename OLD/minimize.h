/*
 * minimize.h
 *
 *  Created on: 28/apr/2015
 *      Author: Gabriele
 */

#ifndef MINIMIZE_H_
#define MINIMIZE_H_
#include <string>
using namespace std;

void stampa_dfa(int**, int, int);

struct my_dfa{
	int** matrixdfa;
	int num_righe;
	int num_colonne;
	int dim_alfabeto;
} typedef my_dfa;


my_dfa* minimize(int** , const int , int );

int* elenco_stati_equivalenti(int** dfa_hp, int n_stati_hp, int** dfa_T, int n_stati_T,  int dim_alfabeto);

string witness(int** dfa_hp, int n_stati_hp, int** dfa_T, int n_stati_T,  int dim_alfabeto);


#endif /* MINIMIZE_H_ */
