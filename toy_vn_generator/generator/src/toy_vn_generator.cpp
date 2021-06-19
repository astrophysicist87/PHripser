#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <complex>
#include <vector>
#include <iomanip>
#include <omp.h>
#include <string>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

#include "../include/toy_vn_generator.h"

using namespace std;

int main(int argc, char *argv[])
{
	// --------------------------
	// Define parameters for flow_generator.
	// --------------------------
	flow_generator::n_pT                   = 1;
	flow_generator::max_pT                 = 300000000.0;
	flow_generator::min_pT                 = 0.0;
	flow_generator::order                  = 2;
	//flow_generator::eta_low                = 2.0;
	//flow_generator::eta_high               = 10.0;
	// --------------------------
	const int nLoops                       = 100000;
	const int n_events_per_loop            = 1;
	flow_generator::N_particles_per_event  = 100;
	flow_generator::N_events_to_generate   = n_events_per_loop;
	flow_generator::N_total_events         = nLoops*n_events_per_loop;
	//flow_generator::N_total_events         = 100000;
	// --------------------------
	flow_generator::use_seed               = false;
	// --------------------------

	// Start the calculation.
	cout << endl;
	cout << "<<<======================================>>>" << endl;
	cout << "<<<         BEGIN FLOW GENERATOR         >>>" << endl;
	cout << "<<<======================================>>>" << endl;
	cout << endl;

	// if passing in filenames, process these
	// adding some command line arguments to simplify things...
	flow_generator::resultsDirectory = argv[1];
	flow_generator::output_name = flow_generator::resultsDirectory + "/data.dat";

	//flow_generator::include_psi2_PTslope_fluctuations
	//		= bool(std::stoi(argv[2]));
	flow_generator::v2_magnitude = std::stof(argv[2]);

	flow_generator::print_randomly_generated_data = true;

	// generate dataset(s) randomly
	//#pragma omp parallel for
	for (int iLoop = 0; iLoop < nLoops; iLoop++)
	{
		cout << "Loop " << iLoop + 1 << "/" << nLoops << endl;
		if ( not flow_generator::use_seed )
			flow_generator::seed_index = iLoop;

		// generate dataset randomly
		int dataset_index = ( flow_generator::print_randomly_generated_data ) ?
							iLoop : -1;
		flow_generator::generate_random_data( dataset_index );				
	}

	cout << endl;
	cout << "<<<======================================>>>" << endl;
	cout << "<<<          END FLOW GENERATOR          >>>" << endl;
	cout << "<<<======================================>>>" << endl;

	return 0;
}

