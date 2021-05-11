#ifndef TOY_VN_GENERATOR_H
#define TOY_VN_GENERATOR_H

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <complex>
#include <vector>
#include <iomanip>
#include <stdio.h>
#include <algorithm>
#include <chrono>
#include <random>
#include <omp.h>
#include "sampled_distribution.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

using namespace std;

namespace flow_generator
{
	// ------------------------------------
	// Initialize default parameters here.
	// ------------------------------------

	int n_pT                               = 6;
	double max_pT                          = 3.0;
	double min_pT				           = 0.0;
	int order                              = 2;
	//double eta_low                         = 2.0;
	//double eta_high                        = 10.0;

	int N_events_to_generate               = 1000;
	int N_particles_per_event              = 100;
	long N_total_events                    = 0;
	long N_running_total_events            = 0;

	bool use_seed                          = true;
	int seed_index                         = -1;

	bool print_randomly_generated_data     = false;

	bool include_psi2_PTslope_fluctuations = false;

	string resultsDirectory                = "./results/";
	string output_name                     = resultsDirectory + "/data.dat";

	// ------------------------------------
	// Define structs and vectors.
	// ------------------------------------

	// struct for event information
	typedef struct
	{
		int eventID;
		int Nparticles;
	
		vector<vector<double> > eta_vec;
		vector<vector<double> > pphi_vec;
		vector<vector<double> > pT_vec;
	} event_data;

	// declare needed vectors
	vector<event_data>                         all_events;
	
	vector<vector<double> >                    pT_vec;
	vector<vector<double> >                    eta_vec;
	vector<vector<double> >                    phi_vec;
	vector<vector<int> >                       event_vec;

	//---------------------------------------------------------
	// Generate data randomly to mimic (non-factorizable) flow.
	void generate_random_data(int dataset=-1)
	{
		double delta_pT = (max_pT-min_pT)/(n_pT);

		pT_vec.clear();
		eta_vec.clear();
		phi_vec.clear();
		event_vec.clear();
		pT_vec = vector<vector<double> >(n_pT, vector<double> (0));
		eta_vec = vector<vector<double> >(n_pT, vector<double> (0));
		phi_vec = vector<vector<double> >(n_pT, vector<double> (0));
		event_vec = vector<vector<int> >(n_pT, vector<int> (0));


		// Set up distributions
		unsigned seed = chrono::system_clock::now().time_since_epoch().count();
		default_random_engine generator;
		if ( use_seed )
			generator = default_random_engine(seed);
		else
			generator = default_random_engine(seed_index);

		uniform_real_distribution<double> mean_pT_distribution(0.90, 1.10);
		//const double mean_pT = mean_pT_distribution(generator);
		const double mean_pT = 1.0;
		exponential_distribution<double> pT_distribution( mean_pT );
		normal_distribution<double> eta_distribution(0.0, 2.0);
		uniform_real_distribution<double> v2_distribution(0.05, 0.45);
		uniform_real_distribution<double> psi2_distribution(0.0, 2.0*M_PI);
		uniform_real_distribution<double> pphi_distribution(0.0, 2.0*M_PI);
		normal_distribution<double> psi2_pTslope_distribution(0.0,0.1);

	
		auto phiFunc = [](double x, double v2) { return x + v2*sin(2.0*x); };
	    std::mt19937 gen;
		if ( use_seed )
			gen = std::mt19937(seed);
		else
			gen = std::mt19937(seed_index);

		double fluctuation_switch_factor
				= double(include_psi2_PTslope_fluctuations);

		for (int iEvent = 0; iEvent < N_events_to_generate; iEvent++)
		{
	
			//double v2 = v2_distribution(generator);
			double v2 = 0.25;
			double psi2 = 0.0*psi2_distribution(generator);
			double psi2_pTslope = fluctuation_switch_factor
									* psi2_pTslope_distribution(generator);
			
			for (int iParticle = 0; iParticle < N_particles_per_event; iParticle++)
			{				
				// ---------------------------------------------
				// Try to mimic different types of flow effects
				// ---------------------------------------------
				const double pT_value = pT_distribution(generator);
			    //Sampled_distribution<> dist(phiFunc, 0.0, 2.0*M_PI, v2*pT_value / (mean_pT+pT_value));
			    Sampled_distribution<> dist(phiFunc, 0.0, 2.0*M_PI, v2);
				const double pphi_value = dist(gen) + psi2 + psi2_pTslope*pT_value;
				//const double pphi_value = dist(gen) + psi2 + 0.5*M_PI*pT_value / (mean_pT+pT_value);
								
				// Sort into bins.
				double pT = pT_value;
				double eta = eta_distribution(generator);
				double phi = pphi_value;
				int event = iEvent;
				
				if(pT > max_pT)
				{
					//cout << "Warning: pT = " << pT << " larger than max_pT = " << max_pT << " in event = " << event << endl;
					continue;
				}
				if(pT < min_pT)
				{
					//cout << "Warning: pT = " << pT << " smaller than min_pT = " << min_pT << " in event = " << event << endl;
					continue;
				}
				else
				{
					int bin = int((pT-min_pT)/delta_pT);
					eta_vec[bin].push_back(eta);
					phi_vec[bin].push_back(phi);
					pT_vec[bin].push_back(pT);
					event_vec[bin].push_back(event);
					//cout << "Check: " << pT << "   " << eta << "   " << phi << "   " << event << " " << bin << endl;
				}
			}
		}
						
		//----------------------------------------------
		all_events.clear();
		all_events = vector<event_data>(N_events_to_generate);
		for ( auto & this_event : all_events )
		{
			this_event.eta_vec = vector<vector<double> >( n_pT, vector<double>(0) );
			this_event.pphi_vec = vector<vector<double> >( n_pT, vector<double>(0) );
			this_event.pT_vec = vector<vector<double> >( n_pT, vector<double>(0) );
		}
		for (int bin = 0; bin < n_pT; bin++)
		{
			const int n_particles_in_this_bin = event_vec[bin].size();
			for (int iParticle = 0; iParticle < n_particles_in_this_bin; iParticle++)
			{
				int event_ID_of_this_particle = event_vec[bin][iParticle];
				
				all_events[event_ID_of_this_particle].eta_vec[bin].push_back( eta_vec[bin][iParticle] );
				all_events[event_ID_of_this_particle].pphi_vec[bin].push_back( phi_vec[bin][iParticle] );
				all_events[event_ID_of_this_particle].pT_vec[bin].push_back( pT_vec[bin][iParticle] );
				all_events[event_ID_of_this_particle].eventID = event_ID_of_this_particle;
			}
		}
		//----------------------------------------------


		//----------------------------------------------
		// Tally number of completed events
		N_running_total_events += N_events_to_generate;
		cout << "Generated " << N_particles_per_event * N_events_to_generate
				<< " particles belonging to "
				<< N_events_to_generate << " events."
				<< endl;
		cout << "Currently analyzed " << N_running_total_events
				<< " / " << N_total_events << " events." << endl;


		//----------------------------------------------
		// Dump simulated data to file if desired
		if ( print_randomly_generated_data )
		{
			string dataset_stem = ( dataset >= 0 ) ?
									std::string(10 - dataset.length(), '0')
									+ std::to_string( dataset ) : "";
			string RNG_filename = resultsDirectory + "/dataset_" + dataset_stem + ".dat";
			ofstream RNG_output( RNG_filename.c_str() );
			for ( auto & this_event : all_events )
			for (int bin = 0; bin < n_pT; bin++)
			{
				const int n_particles_in_this_bin = this_event.eta_vec[bin].size();
				const auto & this_pT_bin = this_event.pT_vec[bin];
				const auto & this_pphi_bin = this_event.pphi_vec[bin];
				const auto & this_eta_bin = this_event.eta_vec[bin];
				for (int iParticle = 0; iParticle < n_particles_in_this_bin; iParticle++)
					RNG_output << this_pT_bin[iParticle] << "   "
                               << this_eta_bin[iParticle] << "   "
                               << this_pphi_bin[iParticle] << "   "
                               << this_event.eventID << endl;
			}
			RNG_output.close();
		}
		
		return;
	}

}


#endif
