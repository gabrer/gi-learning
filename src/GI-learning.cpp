/*
 ============================================================================
 Name        : GI-learning.cpp
 Author      : Pietro Cottone
 Version     :
 Copyright   : You can do whatever you want with this code, if you cite me ;)
 Description : Hello World in C++,
 ============================================================================
 */

#include <ctime>
#include <iostream>
#include <string>
#include "dfa.h"
#include "bluefringe.h"
#include "edsm.h"
#include "lstar.h"
#include "blueStar.h"
#include "messages.h"
#include <boost/filesystem.hpp>


#define SAMPLE_DIR "examples" 								// training set: put your training data here

#define EDSM_FILE "examples.txt" 							// training set: put your training data here
#define EDSM_FILE_BIG "examples_big.txt" 				// training set: put your training data here

#define LSTAR_FILE "lstar.txt"										// file for lstar
#define LSTAR_FILE_BIG "lstar_big.txt" 						// file for lstar

#define DOT_DIR "results"											// dot file with inferred DFA
#define DOT_FILE_BLUESTAR "inferredDFA_bluestar.dot"
#define DOT_FILE_EDSM "inferredDFA_edsm.dot"
#define DOT_FILE_LSTAR "inferredDFA_lstar.dot"
//#define EDSM_RES_PATH "DFA_dot" 							// by-products of inference

#define MAX_BUFFER_SIZE 256
#define BASE_PATH_EXE 0											// the working dir is where the executable is

#define MAX_ARGC 3
#define MIN_ARGC 1

using namespace std;

namespace fs=boost::filesystem;

void parse_input_args(int, char**, string *, string *);
string check_res_dir(const string);

/**
 *
 * @param argc
 * @param argv
 * @return
 */

int main(int argc, char* argv[]){
	//dfa
	//dfa EDSM_dfa = NULL;
	gi::edsm* edsm_exe = NULL;
	gi::blueStar* bluestar_exe = NULL;


	clock_t tStart;

	//file
	string edsm_example_file="";
	string lstar_dfa_file="";

	// working dir
	string base_path;
	string res_path;

	//parse input arguments
	parse_input_args(argc, argv, &edsm_example_file, &lstar_dfa_file);

	// Setting the locale, in order to choose the correct language
	//string curr_os_locale = setlocale(LC_CTYPE, "");
	//cout<<"Current locale: "<<curr_os_locale<<endl;

	if(BASE_PATH_EXE){
		// complete path of the executable
		fs::path selfpath=fs::system_complete(argv[0]);

		// path of the folder where the executable is: it is assumed as the working directory
		base_path = selfpath.parent_path().c_str();
	}else{
		// the working dir is the current dir
		base_path = fs::current_path().c_str();
	}

	cout<<"Working dir: "<< base_path <<endl;



//	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	// --- Blue* algorithm ---
//	cout << endl<< "************************";
//	cout << endl<< "********  BLUESTAR *********" << endl;
//	cout << "************************" << endl;
//
//	//string example_file = base_path + fs::path::preferred_separator + SAMPLE_DIR + fs::path::preferred_separator + edsm_example_file;
//
//
//	// print current example file
//	cout<<"Example file: "<< edsm_example_file <<endl;
//
//	bluestar_exe = new gi::blueStar(edsm_example_file.c_str(), 0.05, 1000.0);
//
//
//	// Infer DFA from positive and negative samples (the algorithm is feeded by data in input file)
//	res_path = check_res_dir(base_path);
//	res_path = res_path + fs::path::preferred_separator;
//
//
//	// start timer to compute execution time
//	tStart = clock();
//
//	gi::dfa* BLUESTAR_dfa = bluestar_exe->run(res_path);
//
//	// Print execution time
//	cout<<"Time taken to EDSM: " <<(double)(clock() - tStart)/CLOCKS_PER_SEC << endl;
//
//
//	// Print transition table of the inferred automaton
//	BLUESTAR_dfa->print_dfa_ttable("- BlueStar dfa -");
//
//	// Create dot figure for the inferred automaton
//	BLUESTAR_dfa->print_dfa_dot_mapped_alphabet("BlueStar", (base_path + fs::path::preferred_separator + DOT_DIR + fs::path::preferred_separator + DOT_FILE_BLUESTAR).c_str());
//
//	//Error rate
//	cout << "Error-rate: "<< bluestar_exe->get_error_rate_final_dfa() << endl;
//
//	// free allocated memory
//	if(bluestar_exe!=NULL)
//		delete bluestar_exe;


	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// --- EDSM algorithm ---
	cout << endl<< "************************";
	cout << endl<< "********  EDSM *********" << endl;
	cout << "************************" << endl;

	//string example_file = base_path + fs::path::preferred_separator + SAMPLE_DIR + fs::path::preferred_separator + edsm_example_file;


	// print current example file
	cout<<"Example file: "<< edsm_example_file <<endl;

	edsm_exe = new gi::edsm(edsm_example_file.c_str());


	// Infer DFA from positive and negative samples (the algorithm is feeded by data in input file)
	//res_path = check_res_dir(base_path);

	res_path = res_path + fs::path::preferred_separator;


	// start timer to compute execution time
	tStart = clock();

	gi::dfa* EDSM_dfa = edsm_exe->run(res_path);

	// Print execution time
	cout<<"Time taken to EDSM: " <<(double)(clock() - tStart)/CLOCKS_PER_SEC << endl;


	// Print transition table of the inferred automaton
	EDSM_dfa->print_dfa_ttable("- Edsm dfa -");

	// Create dot figure for the inferred automaton
	EDSM_dfa->print_dfa_dot_mapped_alphabet("EDSM", (base_path + fs::path::preferred_separator + DOT_DIR + fs::path::preferred_separator + DOT_FILE_EDSM).c_str());

	cout << "DFA di EDSM stampato in: "<<base_path + fs::path::preferred_separator + DOT_DIR + fs::path::preferred_separator + DOT_FILE_EDSM <<endl;

	// free allocated memory
	if(edsm_exe!=NULL)
		delete edsm_exe;



	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// --- LSTAR algorithm ---
	cout << endl<< "*************************";
	cout << endl<< "********  LSTAR  ********" << endl;
	cout <<  "*************************" << endl;

	gi::lstar* lstar_exe = new gi::lstar(gi::dfa::read_dfa_file(lstar_dfa_file));

	// start timer to compute execution time
	tStart = clock();

	gi::dfa* L_dfa = lstar_exe->run_lstar(false, "");

	// Stat lstar
	//n_memb_query_lstar = L_dfa->get_n_memb_query();

	cout<<"Time taken to LSTAR: "<< (double)(clock() - tStart)/CLOCKS_PER_SEC << endl;

	L_dfa->print_dfa_dot_mapped_alphabet("LSTAR", (base_path + fs::path::preferred_separator + DOT_DIR + fs::path::preferred_separator + DOT_FILE_LSTAR).c_str());


	delete lstar_exe;


	if(EDSM_dfa != NULL)
		delete EDSM_dfa;
//	if(BLUESTAR_dfa != NULL)
//		delete BLUESTAR_dfa;
	if(L_dfa != NULL)
		delete L_dfa;
}








void parse_input_args(int argc, char* argv[], string *bs, string *ls){
	if(argc>MAX_ARGC || argc<MIN_ARGC){
		cerr<<MSG_WRONG_ARGC<<endl;

		exit(EXIT_FAILURE);
	}

	if(argc>=2){
		if(!strcmp("little", argv[1])){
			(*bs) = (*bs) + EDSM_FILE;
		}else if (!strcmp("big", argv[1])){
			(*bs) = (*bs) + EDSM_FILE_BIG;
		}else{
			//cerr<<MSG_WRONG_ARGV<< argv[1] <<endl;
			(*bs) = argv[1];
		}

		if(3==argc){
			if(!strcmp("little", argv[1])){
				(*ls) = (*ls) + LSTAR_FILE;
			}else if (!strcmp("big", argv[1])){
				(*ls) = (*ls) + LSTAR_FILE_BIG;
			}else{
				(*ls) = (*ls) + argv[2];
			}
		}
	}else{
		(*bs) = (*bs) + EDSM_FILE;
		(*ls) = (*ls) + LSTAR_FILE;
	}


}


string check_res_dir(const string  base_path){
	string res_path;

	// Infer DFA from positive and negative samples (the algorithm is feeded by data in input file)
	res_path = base_path + fs::path::preferred_separator + DOT_DIR;

	// move res_dir if it exists
	if(fs::exists(res_path)){
		char append_time[MAX_BUFFER_SIZE];
		time_t rawtime = std::time(NULL);

		struct tm * timeinfo;

		time ( &rawtime );
		timeinfo = localtime ( &rawtime );


		strftime(append_time, MAX_BUFFER_SIZE, "%Y_%m_%d_%H_%M_%S", timeinfo);

		fs::rename(res_path, res_path + "_" + append_time);
	}

	// create res_dir
	fs::create_directory( res_path );

	res_path = res_path + fs::path::preferred_separator;

	return res_path;
}
