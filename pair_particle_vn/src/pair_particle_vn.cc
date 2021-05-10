#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <complex>
#include <vector>
#include <iomanip>
#include <string>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

#include "../include/pair_particle_vn.h"
#include "../include/pair_particle_vn_Plumberg.h"
#include "../include/pair_particle_vn_Heinz.h"
#include "../include/pair_particle_vn_w_errors.h"

using namespace std;

int main(int argc, char *argv[])
{
	constexpr bool use_plumberg = false;
	constexpr bool use_heinz = false;
	constexpr bool use_PCA = true;

	if ( use_plumberg )
	{
		cout << endl;
		cout << "<<<======================================>>>" << endl;
		cout << "<<<            BEGIN PLUMBERG            >>>" << endl;
		// Load dataset.
		const int nLoops = 1;
		for (int iLoop = 0; iLoop < nLoops; iLoop++)
		{
			cout << "Loop " << iLoop + 1 << "/" << nLoops << endl;
	
			// generate dataset randomly
			plumberg::generate_random_data();
			
			// update Q vectors
			plumberg::update_Q_n();
		}
	
		// Finally evaluate ensemble-averaged Vn
		plumberg::get_V_n_and_eigendecomposition();

		cout << "<<<             END PLUMBERG             >>>" << endl;
		cout << "<<<======================================>>>" << endl;
		cout << endl;
	}

	if ( use_heinz )
	{
		cout << endl;
		cout << "<<<======================================>>>" << endl;
		cout << "<<<             BEGIN  HEINZ             >>>" << endl;
		// Load dataset.
		const int nLoops = 10;
		heinz::use_seed = false;
		for (int iLoop = 0; iLoop < nLoops; iLoop++)
		{
			cout << "Loop " << iLoop + 1 << "/" << nLoops << endl;
			if ( not heinz::use_seed )
				heinz::seed_index = iLoop;

			// generate dataset randomly
			heinz::generate_random_data();
			
			// update Q vectors
			heinz::update_Q_n();
		}
	
		// Finally evaluate ensemble-averaged Vn
		heinz::get_V_n_and_eigendecomposition();
		cout << "<<<              END  HEINZ              >>>" << endl;
		cout << "<<<======================================>>>" << endl;
		cout << endl;
	}

	if ( use_PCA )
	{
		// --------------------------
		// Define parameters for PCA.
		// --------------------------
		PCA::Vn_mode                = 0;
		PCA::n_pT                   = 6;
		PCA::max_pT                 = 3.0;
		PCA::min_pT					= 0.3;
		PCA::order                  = 2;
		PCA::eta_low                = 2.0;
		PCA::eta_high               = 10.0;
		// --------------------------
		const int nLoops            = 10;
		const int n_events_per_loop = 10000;
		PCA::N_particles_per_event  = 1000;
		PCA::N_events_to_generate   = n_events_per_loop;
		PCA::N_total_events         = nLoops*n_events_per_loop;
		//PCA::N_total_events         = 100000;
		// --------------------------
		PCA::use_seed = true;
		bool RNG_mode = false;
		// --------------------------

		// Start the calculation.
		cout << endl;
		cout << "<<<======================================>>>" << endl;
		cout << "<<<              BEGIN  PCA              >>>" << endl;

		// if passing in filenames, process these
		if ( ( not RNG_mode ) and argc > 1 )
		{
			// Overwrite default set above from command-line
			PCA::N_total_events = stoi( argv[1] );

			// Do this once
			PCA::initialize_vectors();

			// Set output filename from second command-line argument
			PCA::output_name = argv[2];

			// Remaining command-line arguments are files to analyze
			for ( int iArg = 3; iArg < argc; iArg++ )
			{
				// Load dataset from file.
				string filename = string( argv[iArg] );
				PCA::read_in_data( filename );
	
				// Compute Q vectors.
				PCA::update_Q_n();
			}
		}
		else
		{
			// adding some command line arguments to simplify things...
			PCA::resultsDirectory = argv[1];
			PCA::output_name = PCA::resultsDirectory + "/eigensystem.dat";

			PCA::include_psi2_PTslope_fluctuations
					= bool(std::stoi(argv[2]));

			PCA::Vn_mode = std::stoi(argv[3]);

			if ( PCA::Vn_mode == 0 )
				PCA::print_randomly_generated_data = true;

			// Do this once
			PCA::initialize_vectors();

			// generate dataset(s) randomly
			for (int iLoop = 0; iLoop < nLoops; iLoop++)
			{
				cout << "Loop " << iLoop + 1 << "/" << nLoops << endl;
				if ( not PCA::use_seed )
					PCA::seed_index = iLoop;
		
				// generate dataset randomly
				int dataset_index = ( PCA::print_randomly_generated_data ) ?
									iLoop : -1;
				PCA::generate_random_data( dataset_index );
				
				// update Q vectors
				PCA::update_Q_n();
			}
		}
	
		// Finally evaluate ensemble-averaged Vn,
		// jackknife estimates, and eigen-decompose.
		PCA::get_V_n_with_errors_and_eigendecomposition();
		cout << "<<<               END  PCA               >>>" << endl;
		cout << "<<<======================================>>>" << endl;
		cout << endl;
	}


  return 0;
}

