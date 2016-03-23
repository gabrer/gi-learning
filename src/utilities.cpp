/*
 * utility.cpp
 *
 *  Created on: 29 set 2015
 *      Author: piero
 */

#include <utilities.h>
#include <sstream>
#include <fstream>
#include <math.h>

#include <boost/accumulators/accumulators.hpp>	// For Accumulators
#include <boost/accumulators/statistics.hpp>			// To calculate mean and variance for an accumulator

// Scores of Zalpha for alpha = {0.05, 0.025, 0.01} from the Zscore table
#define Z_FOR_ALPHA_05 	1.645
#define Z_FOR_ALPHA_025 	1.96
#define Z_FOR_ALPHA_01 	2.33

using namespace std;

string intTostring(int a)
{
    ostringstream temp;
    temp<<a;
    return temp.str();
}

int stringToint(string str)
{
	stringstream stream(str);

	int n;
	while(1) {
	   stream >> n;
	   if(!stream)
	      break;
	}
	return n;
}

string charToString(char c)
{
	stringstream ss;
	string s;
	ss << c;
	ss >> s;

	return s;
}

string trimRight(string str)
{
	while(str[str.length()-1] == ' ')
		str = str.substr(0,str.length()-1);

	return str;
}


double u_alpha_score(double alpha)
{
	double u_alpha = 0;


	// Selection of the Zalpha value
	if(alpha == 0.05)
		u_alpha = Z_FOR_ALPHA_05;

	else if(alpha == 0.025)
		u_alpha = Z_FOR_ALPHA_025;

	else if(alpha == 0.01)
		u_alpha = Z_FOR_ALPHA_01;

	else{
		cerr << "ERR: Non valid value for Alpha: select one between {0.05, 0.025, 0.01}" << endl;
		cerr << "The inference process was terminated with error"<<endl;
		exit(EXIT_FAILURE);
	}

	return u_alpha;
}



double z_alpha_from_u_alpha_two_proportions_test(double prop1, double prop2, int sample_size, double alpha, double* dev_std_h0)
{
	double u_alpha	= 0,	z_alpha = 0;
	double p_est		=	0,	q_est	 = 0;
	double std_err = 0;


	// Selection of the Zalpha value (u_alpha is the value for the standard normal deviate N(0,1))
	u_alpha = u_alpha_score(alpha);


	// - Estimated error rate for population by mean of error rates for samples -
	// p_est is the accumulator set
	boost::accumulators::accumulator_set<double, boost::accumulators::features<boost::accumulators::tag::mean, boost::accumulators::tag::variance> > p_acc;

	// Insert the values inside the accumulator
	p_acc(prop1);
	p_acc(prop2);

	// Compute p estimator and q
	p_est = boost::accumulators::mean(p_acc);
	q_est = 1-p_est;


	// Check if the error rate is zero
	// Without this threashold, if the proportions is too low the approximation conditions from binomial to gaussian distribution will be not satisfied.
	// And, it's not a problem to due with very small error like if it's no error.
	if( p_est <= 0.04 ){

		#ifdef DEBUG1
		cout << "zeta-beta: MINF"<<endl;
		#endif

		return DMINF;
	}


	// Check approximation conditions
	if(!approximation_conditions_binomial_to_gaussian_distr(p_est, q_est, sample_size)){
		cerr << "ERR: The process stop!" << endl;

		throw "Constraints approximation to Normal Gaussian distr not satisfied";
	}


	// Variance and dev. std. under the H0 hypothesis
	double variance_h0 = (double) (2*p_est*q_est) / (double) sample_size;
	*dev_std_h0 = sqrt( variance_h0 );


	// Z-alpha: u_alpha is the z_alpha value in the normal gaussian distribution, multiplicated for std_err became usefull for the specific distribution
	z_alpha = u_alpha * (*dev_std_h0);


	return z_alpha;
}


double z_alpha_from_u_alpha_two_proportions_test(double prop1, double prop2, int sample_size, double alpha)
{
	double tmp = 0;
	return z_alpha_from_u_alpha_two_proportions_test(prop1, prop2, sample_size, alpha, &tmp);
}



bool approximation_conditions_binomial_to_gaussian_distr(double p_est, double q_est, int sample_size)
{
	// Necessary approximation conditions to normal distribution for binomial distribution
	if ( !(sample_size * p_est > 5  && sample_size * q_est > 5) ){
		cerr << endl << "ERR: Approximation conditions to Normal Distributions are not satisfied."<<endl;
		cerr << "ERR: Details: p: "<<p_est <<", q: "<<q_est<< ", N: "<<sample_size<< endl;
		cerr << "ERR: N*p: "<<sample_size * p_est << ", N*q: "<<sample_size * q_est << endl;
		return false;

	} else {

		return true;
	}
}


int getPoisson(double lambda)
{
	//srand(time(NULL));

	double L = exp(-lambda);
	double p = 1.0;
	int k = 0;

	do{
		k++;
		double tmp = ((double) rand() / (RAND_MAX));				// Random between 0 and 1
		p *= tmp;
	}while(p > L);

	return k-1;
}


// Normalized Compression Distance (NCD)
double ncd(double comp_x, double comp_y, double comp_xy)
{
	return ( comp_xy - fmin(comp_x, comp_y) ) / fmax(comp_x, comp_y);
}


