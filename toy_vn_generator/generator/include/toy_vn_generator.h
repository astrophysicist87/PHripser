#ifndef PAIR_PARTICLE_VN_W_ERRORS_H
#define PAIR_PARTICLE_VN_W_ERRORS_H

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

namespace PCA
{
	// ------------------------------------
	// Initialize default parameters here.
	// ------------------------------------

	// Vn_mode chooses which PCA formula to use
	// 0 - Eq. (9) in 1906.08915 (without event-wise pair normalization)
	// 1 - Eq. (8) in 1906.08915 (with event-wise pair normalization)
	int Vn_mode                 = 0;

	int n_pT                    = 6;
	double max_pT               = 3.0;
	double min_pT				= 0.0;
	int order                   = 2;
	double eta_low              = 2.0;
	double eta_high             = 10.0;

	int N_events_to_generate    = 1000;
	int N_particles_per_event   = 100;
	long N_total_events         = 0;
	long N_running_total_events = 0;

	bool use_seed = true;
	int seed_index = -1;

	bool print_randomly_generated_data = false;

	bool include_psi2_PTslope_fluctuations = true;

	string resultsDirectory = "./results/";
	//string resultsFilename = "eigensystem.dat";
	string output_name = resultsDirectory + "/eigensystem.dat";

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

	vector<vector<long> >                      Npairs;
	vector<vector<double> >                    average_Npairs_per_event;
	vector<vector<complex<double> > >          V_n_vec;

	vector<vector<complex<double> > >          Q_i;
	vector<vector<complex<double> > >          Q_j;
	vector<vector<complex<double> > >          Q_ij;

	// Need to store all events separately for error estimation
	vector<vector<vector<complex<double> > > > Q_i_jackknife_estimates;
	vector<vector<vector<complex<double> > > > Q_j_jackknife_estimates;
	vector<vector<vector<complex<double> > > > Q_ij_jackknife_estimates;

	vector<vector<vector<complex<double> > > > jackknife_V_n_estimates;
	vector<vector<vector<double> > >           jackknife_Npairs;

	//-----------------------
	// Load data from a file.
	void read_in_data(string filename)
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


		// load data here
		ifstream input( filename );
		string line = "";

		int N_events_in_this_file = 0, N_particles = 0;
		int initial_event = 0;
		while ( getline( input, line ) )
		{
			// cout << "Entered Loop" << endl;
			
			N_particles++;
			string value;
			vector<double> value_vector(0);
			for ( istringstream is (line); is >> value; )
				value_vector.push_back( stod( value ) );
			
			// Sort into bins.
			double pT = value_vector[0];
			double eta = value_vector[1];
			double phi = value_vector[2];
			int event = value_vector[3];
			if (N_particles==1) initial_event = event;
				
			if ( pT >= max_pT )
				continue;
			if ( pT <= min_pT )
				continue;
			else
			{
				int bin = int((pT-min_pT)/delta_pT);
				eta_vec[bin].push_back(eta);
				phi_vec[bin].push_back(phi);
				pT_vec[bin].push_back(pT);
				event_vec[bin].push_back(event);
			}

			N_events_in_this_file = event - initial_event + 1;
		}
		

		// Update running total number of events studied so far
		N_running_total_events += N_events_in_this_file;
		cout << "Read in " << N_particles
				<< " particles belonging to "
				<< N_events_in_this_file << " events from "
				<< filename << "."
				<< endl;
		cout << "Currently analyzed " << N_running_total_events
				/*<< " / " << N_total_events */<< " events." << endl;
		//N_total_events = N_running_total_events;

		//----------------------------------------------
		// Rearrange to make data more conveniently accessible in evaluating observables
		all_events.clear();
		all_events = vector<event_data>(N_events_in_this_file);
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
				
				all_events[event_ID_of_this_particle - initial_event].eta_vec[bin].push_back( eta_vec[bin][iParticle] );
				all_events[event_ID_of_this_particle - initial_event].pphi_vec[bin].push_back( phi_vec[bin][iParticle] );
				all_events[event_ID_of_this_particle - initial_event].pT_vec[bin].push_back( pT_vec[bin][iParticle] );
				all_events[event_ID_of_this_particle - initial_event].eventID = event_ID_of_this_particle;
			}
		}
		//----------------------------------------------
		
		return;
	}


	//-----------------------
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
			double psi2 = psi2_distribution(generator);
			double psi2_pTslope = fluctuation_switch_factor
									* psi2_pTslope_distribution(generator);
			
			for (int iParticle = 0; iParticle < N_particles_per_event; iParticle++)
			{				
				// ---------------------------------------------
				// Try to mimic different types of flow effects
				// ---------------------------------------------
				const double pT_value = pT_distribution(generator);
			    Sampled_distribution<> dist(phiFunc, 0.0, 2.0*M_PI, v2*pT_value / (mean_pT+pT_value));
				const double pphi_value = dist(gen) + psi2 + psi2_pTslope*pT_value;
			    //Sampled_distribution<> dist(phiFunc, 0.0, 2.0*M_PI, v2*pT_value / (mean_pT+pT_value));
				//const double pphi_value = dist(gen) + psi2 + 0.5*M_PI*pT_value / (mean_pT+pT_value);
								
				// Sort into bins.
				double pT = pT_value;
				double eta = eta_distribution(generator);
				double phi = pphi_value;
				int event = iEvent;
				
				if(pT >= max_pT)
				{
					//cout << "Warning: pT = " << pT << " larger than max_pT = " << max_pT << " in event = " << event << endl;
					continue;
				}
				if(pT <= min_pT)
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

		N_running_total_events += N_events_to_generate;
		cout << "Generated " << N_particles_per_event * N_events_to_generate
				<< " particles belonging to "
				<< N_events_to_generate << " events."
				<< endl;
		cout << "Currently analyzed " << N_running_total_events
				<< " / " << N_total_events << " events." << endl;

		if ( print_randomly_generated_data )
		{
			string dataset_stem = ( dataset >= 0 ) ?
									"_" + std::to_string( dataset ) : "";
			string RNG_filename = resultsDirectory + "/dataset" + dataset_stem + ".dat";
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
	
	
	//-----------------------
	// Perform eigendecomposition on each flow matrix in V_n_vec.
	void get_eigen( 
	  vector<vector<vector<complex<double> > > > & V_n_vec,
	  vector<vector<vector<complex<double> > > > & eigenvectors,
	  vector<vector<double> > & eigenvalues )
	{
		
		int number_of_events = V_n_vec.size();
		int dim = V_n_vec[0][0].size();
		gsl_vector *eval = gsl_vector_alloc(dim);
		gsl_matrix_complex *evec = gsl_matrix_complex_alloc(dim,dim);

		for (int iEvent = 0; iEvent < number_of_events; iEvent++)
		{
			// Solves eigensystem
			gsl_matrix_complex * m = gsl_matrix_complex_alloc (dim, dim);
			for (int i=0; i<dim; i++)
			for (int j=0; j<dim; j++)
			{
				complex<double> value = V_n_vec[iEvent][i][j];
				gsl_complex z_value;
				GSL_SET_COMPLEX(&z_value, real(value), imag(value));
				gsl_matrix_complex_set (m, i, j, z_value);
			}
			gsl_eigen_hermv_workspace *w = gsl_eigen_hermv_alloc(dim);
			gsl_eigen_hermv (m, eval, evec, w);
			gsl_eigen_hermv_free(w);
			gsl_eigen_hermv_sort(eval, evec, GSL_EIGEN_SORT_VAL_DESC);
			
			// Store results
			for (int i=0; i<dim; i++)
			{
				double lambda = gsl_vector_get(eval, i);
				gsl_vector_complex_view evec_i = gsl_matrix_complex_column(evec, i);
		
			    eigenvalues[iEvent][i] = lambda;
			
				for (int j = 0; j < dim; ++j)
			    {
					gsl_complex z = gsl_vector_complex_get( &evec_i.vector, j );	
					eigenvectors[iEvent][i][j] = complex<double>( GSL_REAL(z), GSL_IMAG(z) );
				}
			}
		}
		
		gsl_vector_free(eval);
		gsl_matrix_complex_free(evec);
	
		return;
	}
	


	//-----------------------
	void initialize_vectors()
	{
	
		Npairs                   = vector<vector<long> >(             n_pT, vector<long>(             n_pT, 0.0 ) );
		average_Npairs_per_event = vector<vector<double> >(           n_pT, vector<double>(           n_pT, 0.0 ) );
		V_n_vec                  = vector<vector<complex<double> > >( n_pT, vector<complex<double> >( n_pT, 0.0 ) );
	
		Q_i  = vector<vector<complex<double> > >( n_pT, vector<complex<double> >( n_pT, 0.0 ) );
		Q_j  = vector<vector<complex<double> > >( n_pT, vector<complex<double> >( n_pT, 0.0 ) );
		Q_ij = vector<vector<complex<double> > >( n_pT, vector<complex<double> >( n_pT, 0.0 ) );

		// Needed for error estimation
		Q_i_jackknife_estimates
			= vector<vector<vector<complex<double> > > >( N_total_events,
				vector<vector<complex<double> > >( n_pT,
				vector<complex<double> >( n_pT, 0.0 ) ) );
		Q_j_jackknife_estimates
			= vector<vector<vector<complex<double> > > >( N_total_events,
				vector<vector<complex<double> > >( n_pT,
				vector<complex<double> >( n_pT, 0.0 ) ) );
		Q_ij_jackknife_estimates
			= vector<vector<vector<complex<double> > > >( N_total_events,
				vector<vector<complex<double> > >( n_pT,
				vector<complex<double> >( n_pT, 0.0 ) ) );
		jackknife_V_n_estimates
			= vector<vector<vector<complex<double> > > >( N_total_events,
				vector<vector<complex<double> > >( n_pT,
				vector<complex<double> >( n_pT, 0.0 ) ) );
		jackknife_Npairs
			= vector<vector<vector<double> > >( N_total_events,
				vector<vector<double> > ( n_pT,
				vector<double>( n_pT, 0.0 ) ) );

		return;
	}


	//-----------------------
	void update_Q_n()
	{
		const int current_n_events = all_events.size();
		vector<complex<double> > Q_i_jackknife_estimates_flat( N_total_events * n_pT * n_pT, 0.0 );
		vector<complex<double> > Q_j_jackknife_estimates_flat( N_total_events * n_pT * n_pT, 0.0 );
		vector<complex<double> > Q_ij_jackknife_estimates_flat( N_total_events * n_pT * n_pT, 0.0 );
		vector<double> jackknife_Npairs_flat( N_total_events * n_pT * n_pT, 0.0 );

		vector<complex<double> > Q_i_these_events( current_n_events * n_pT * n_pT, 0.0 );
		vector<complex<double> > Q_j_these_events( current_n_events * n_pT * n_pT, 0.0 );
		vector<complex<double> > Q_ij_these_events( current_n_events * n_pT * n_pT, 0.0 );
		vector<double> Npairs_these_events( current_n_events * n_pT * n_pT, 0.0 );

		auto indexer = [](int iEvent, int pT_i, int pT_j, int n_pT)
						{ return (iEvent * n_pT + pT_i ) * n_pT + pT_j; };

		// loop over all events
		for ( int iEvent = 0; iEvent < current_n_events; iEvent++ )
		{
			if ( iEvent > 0 and iEvent % 100 == 0 )
				cout << "  --> " << iEvent << "/" << current_n_events << " events completed." << endl;
			auto & this_event = all_events[iEvent];

			const int eventID_to_skip = this_event.eventID;

			vector<vector<complex<double> > > Q_i_this_event( n_pT, vector<complex<double> >( n_pT, 0.0 ) );
			vector<vector<complex<double> > > Q_j_this_event( n_pT, vector<complex<double> >( n_pT, 0.0 ) );
			vector<vector<complex<double> > > Q_ij_this_event( n_pT, vector<complex<double> >( n_pT, 0.0 ) );
			vector<vector<double> > Npairs_this_event( n_pT, vector<double>( n_pT, 0.0 ) );
	
			for (int pT_i = 0; pT_i < n_pT; pT_i++)
			for (int pT_j = 0; pT_j < n_pT; pT_j++)
			{
				//cout << "Loop: pT_i = " << pT_i << ", pT_j = " << pT_j << endl;
				
				vector<double> eta_bin_i = this_event.eta_vec[pT_i];
				vector<double> eta_bin_j = this_event.eta_vec[pT_j];
				vector<double> phi_bin_i = this_event.pphi_vec[pT_i];
				vector<double> phi_bin_j = this_event.pphi_vec[pT_j];
				
				const int n_part_i = eta_bin_i.size();
				const int n_part_j = eta_bin_j.size();

				for (int i = 0; i < n_part_i; i++)
					Q_i_this_event[pT_i][pT_j] += exp( complex<double> ( 0.0, order * phi_bin_i[i] ) );
				for (int j = 0; j < n_part_j; j++)
					Q_j_this_event[pT_i][pT_j] += exp( complex<double> ( 0.0, -order * phi_bin_j[j] ) );
				
				for (int i = 0; i < n_part_i; i++)
				for (int j = 0; j < n_part_j; j++)
				{
					if ( pT_i==pT_j and i==j ) continue;
				    double delta_eta = eta_bin_i[i] - eta_bin_j[j];
					if ( delta_eta*delta_eta >= eta_low*eta_low)
					{
						Q_ij_this_event[pT_i][pT_j] += exp( complex<double> ( 0.0, order * ( phi_bin_i[i] - phi_bin_j[j] ) ) );
						Npairs[pT_i][pT_j]++;
						Npairs_this_event[pT_i][pT_j]++;
					}
				}

				// Finally, add in contributions from this event
				// (with normalizations as appropriate)
				complex<double> this_event_Q_i = 0.0, this_event_Q_j = 0.0, this_event_Q_ij = 0.0;
				if ( Vn_mode == 1 )
				{
					if ( n_part_i > 0 )
						this_event_Q_i  = Q_i_this_event[pT_i][pT_j] / double(n_part_i);
					if ( n_part_j > 0 )
						this_event_Q_j  = Q_j_this_event[pT_i][pT_j] / double(n_part_j);
					if ( Npairs_this_event[pT_i][pT_j] > 0 )
						this_event_Q_ij = Q_ij_this_event[pT_i][pT_j] / Npairs_this_event[pT_i][pT_j];
				}
				else
				{
					this_event_Q_i  = Q_i_this_event[pT_i][pT_j];
					this_event_Q_j  = Q_j_this_event[pT_i][pT_j];
					this_event_Q_ij = Q_ij_this_event[pT_i][pT_j];
				}

				// update average over all events
				Q_i[pT_i][pT_j]  += this_event_Q_i;
				Q_j[pT_i][pT_j]  += this_event_Q_j;
				Q_ij[pT_i][pT_j] += this_event_Q_ij;

				const double npairs_this_event_this_bin = Npairs_this_event[pT_i][pT_j];

				int index = indexer(iEvent,pT_i,pT_j,n_pT);
				Q_i_these_events[index]    = this_event_Q_i;
				Q_j_these_events[index]    = this_event_Q_j;
				Q_ij_these_events[index]   = this_event_Q_ij;
				Npairs_these_events[index] = npairs_this_event_this_bin;

			}

		}

		/*
		for ( int iEvent = 0; iEvent < current_n_events; iEvent++ )
		{	
			for (int pT_i = 0; pT_i < n_pT; pT_i++)
			for (int pT_j = 0; pT_j < n_pT; pT_j++)
			{
				

				// update jackknife averages (skip only current event)
				for (int jEvent = 0; jEvent < eventID_to_skip; jEvent++)
				{
					int index = indexer(jEvent,pT_i,pT_j,n_pT);
					Q_i_jackknife_estimates_flat[index] += this_event_Q_i;
					Q_j_jackknife_estimates_flat[index] += this_event_Q_j;
					Q_ij_jackknife_estimates_flat[index] += this_event_Q_ij;
					jackknife_Npairs_flat[index] += npairs_this_event_this_bin;
				}
				for (int jEvent = eventID_to_skip + 1; jEvent < N_total_events; jEvent++)
				{
					int index = indexer(jEvent,pT_i,pT_j,n_pT);
					Q_i_jackknife_estimates_flat[index] += this_event_Q_i;
					Q_j_jackknife_estimates_flat[index] += this_event_Q_j;
					Q_ij_jackknife_estimates_flat[index] += this_event_Q_ij;
					jackknife_Npairs_flat[index] += npairs_this_event_this_bin;
				}

			}		// end loop over pT bins

		}			// end events loop
		*/

cout << "Made it before" << endl;

		int n_completed_events = 0;
		#pragma omp parallel for
		for (int jEvent = 0; jEvent < N_total_events; jEvent++)
		{
			for (int pT_i = 0; pT_i < n_pT; pT_i++)
			for (int pT_j = 0; pT_j < n_pT; pT_j++)
			{
				int jindex = indexer(jEvent,pT_i,pT_j,n_pT);
				for (int iEvent = 0; iEvent < current_n_events; iEvent++)
				{
					int iindex = indexer(iEvent,pT_i,pT_j,n_pT);
					Q_i_jackknife_estimates_flat[jindex]  += Q_i_these_events[iindex];
					Q_j_jackknife_estimates_flat[jindex]  += Q_j_these_events[iindex];
					Q_ij_jackknife_estimates_flat[jindex] += Q_ij_these_events[iindex];
					jackknife_Npairs_flat[jindex]         += Npairs_these_events[iindex];
				}
			}
			#pragma omp critical
			{
				n_completed_events++;
				if ( n_completed_events%10000 == 0 )
					cout << "Completed " << n_completed_events << "/" << N_total_events << " events." << endl;
			}
		}			// end events loop

		// Then subtract off the one event I should have skipped
		for (int iEvent = 0; iEvent < current_n_events; iEvent++)
		for (int pT_i = 0; pT_i < n_pT; pT_i++)
		for (int pT_j = 0; pT_j < n_pT; pT_j++)
		{
			const int eventID_to_skip = all_events[iEvent].eventID;
			int iindex = indexer(iEvent,pT_i,pT_j,n_pT);
			int jindex = indexer(eventID_to_skip,pT_i,pT_j,n_pT);
			Q_i_jackknife_estimates_flat[jindex]  -= Q_i_these_events[iindex];
			Q_j_jackknife_estimates_flat[jindex]  -= Q_j_these_events[iindex];
			Q_ij_jackknife_estimates_flat[jindex] -= Q_ij_these_events[iindex];
			jackknife_Npairs_flat[jindex]         -= Npairs_these_events[iindex];
		}

cout << "Made it after" << endl;

		int index = 0;
		for (int iEvent = 0; iEvent < N_total_events; iEvent++)
		for (int pT_i = 0; pT_i < n_pT; pT_i++)
		for (int pT_j = 0; pT_j < n_pT; pT_j++)
		{
			Q_i_jackknife_estimates[iEvent][pT_i][pT_j]  += Q_i_jackknife_estimates_flat[index];
			Q_j_jackknife_estimates[iEvent][pT_i][pT_j]  += Q_j_jackknife_estimates_flat[index];
			Q_ij_jackknife_estimates[iEvent][pT_i][pT_j] += Q_ij_jackknife_estimates_flat[index];
			jackknife_Npairs[iEvent][pT_i][pT_j]         += jackknife_Npairs_flat[index];
			index++;
		}

		return;
	}


	//-----------------------
	void get_V_n_with_errors_and_eigendecomposition()
	{
		
		cout << endl << endl << endl;
		
		// Going through all possible bin-pairs.
		// Defining new z=0 and N_pairs=0.
		for (int pT_i = 0; pT_i < n_pT; pT_i++)
		for (int pT_j = 0; pT_j < n_pT; pT_j++)
		{

			Q_i[pT_i][pT_j] /= static_cast<double>(N_total_events);
			Q_j[pT_i][pT_j] /= static_cast<double>(N_total_events);
			Q_ij[pT_i][pT_j] /= static_cast<double>(N_total_events);

			// Normalize jackknife estimates too
			for (int iEvent = 0; iEvent < N_total_events; iEvent++)
			{
				Q_i_jackknife_estimates[iEvent][pT_i][pT_j] /= static_cast<double>( N_total_events - 1 );
				Q_j_jackknife_estimates[iEvent][pT_i][pT_j] /= static_cast<double>( N_total_events - 1 );
				Q_ij_jackknife_estimates[iEvent][pT_i][pT_j] /= static_cast<double>( N_total_events - 1 );
				jackknife_Npairs[iEvent][pT_i][pT_j] /= static_cast<double>( N_total_events - 1 );
			}


			average_Npairs_per_event[pT_i][pT_j] = Npairs[pT_i][pT_j] / static_cast<double>(N_total_events);

			V_n_vec[pT_i][pT_j] = Q_ij[pT_i][pT_j] - Q_i[pT_i][pT_j]*Q_j[pT_i][pT_j];

			/*cout << " --> " << Q_ij[pT_i][pT_j] << "   "
							<< Q_i[pT_i][pT_j] << "   "
							<< Q_j[pT_i][pT_j] << "   "
							<< V_n_vec[pT_i][pT_j] << "   "
							<< static_cast<double>(N_total_events) << endl;*/

			for (int iEvent = 0; iEvent < N_total_events; iEvent++)
				jackknife_V_n_estimates[iEvent][pT_i][pT_j]
					= Q_ij_jackknife_estimates[iEvent][pT_i][pT_j]
						- Q_i_jackknife_estimates[iEvent][pT_i][pT_j]
						* Q_j_jackknife_estimates[iEvent][pT_i][pT_j];


		}

//for (int iEvent = 0; iEvent < N_total_events; iEvent++)
//{
//	cout << "JK " << iEvent << " = {";
//	for (int pT_i = 0; pT_i < n_pT; pT_i++)
//	{
//		cout << "{";
//		for (int pT_j = 0; pT_j < n_pT; pT_j++)
//		{
//			cout << real(V_n_vec[pT_i][pT_j]) << "+(" << imag(V_n_vec[pT_i][pT_j]) << ") I";
//			if ( pT_j + 1 == n_pT ) cout << "}";
//			else cout << ", ";
//		}
//		cout << endl;
//	}
//	cout << "}" << endl << endl;
//}

		////////////////////////////////////////////////////
		// Diagonalize Vn matrix averaged over all events
		////////////////////////////////////////////////////

		// Makes it easier to pass this to get eigen
		vector<vector<vector<complex<double> > > > V_n_vec_EnsembleAverage;
		vector<vector<vector<complex<double> > > >
			eigenvectors_EnsembleAverage( 1, vector<vector<complex<double> > > ( n_pT, vector<complex<double> >( n_pT, 0.0 ) ) );
		vector<vector<double> >
			eigenvalues_EnsembleAverage( 1, vector<double>( n_pT, 0.0 ) );
		V_n_vec_EnsembleAverage.push_back( V_n_vec );
		
		// diagonalize and get eigensystem for the average V_n over full ensemble
		get_eigen( V_n_vec_EnsembleAverage, eigenvectors_EnsembleAverage, eigenvalues_EnsembleAverage );


		///////////////////////////////////////////////////////
		// Diagonalize Vn matrix from each jackknife estimate
		///////////////////////////////////////////////////////

		// Next, pass jackknife estimates to eigendecomposition routine
		vector<vector<vector<complex<double> > > >
			jackknife_eigenvectors( N_total_events,  vector<vector<complex<double> > > ( n_pT, vector<complex<double> >( n_pT, 0.0 ) ) );
		vector<vector<double> >
			jackknife_eigenvalues( N_total_events, vector<double>( n_pT, 0.0 ) );


		// diagonalize each jackknife matrix and store the resulting eigensystems
		get_eigen( jackknife_V_n_estimates, jackknife_eigenvectors, jackknife_eigenvalues );

		// ensure eigenvector signs are consistent for each estimate and average
		for (int iEvent = 0; iEvent < N_total_events; iEvent++)
		for (int pT_i = 0; pT_i < n_pT; pT_i++)
			if ( real( jackknife_eigenvectors[iEvent][pT_i][0] ) < 0.0 )
				for (int pT_j = 0; pT_j < n_pT; pT_j++)
					jackknife_eigenvectors[iEvent][pT_i][pT_j] *= -1.0;


//cout << "Check jackknife eigenvalue 0: " << endl;
//for (int iEvent = 0; iEvent < N_total_events; iEvent++)
//	cout << iEvent << "   " << jackknife_eigenvalues[iEvent][0] << endl;


		////////////////////////////////////
		// Only output after this point
		////////////////////////////////////

		// ------------------------------
		// Print results for complex V_n
		cout << "<<<=====================================>>>" << endl;
		cout << "Complex V_n: " << endl;
		// Print V_n matrix
		cout << "{";
		for (int pT_i = 0; pT_i < n_pT; pT_i++)
		{
			cout << "{";
			for (int pT_j = 0; pT_j < n_pT; pT_j++)
			{
				cout << real(V_n_vec[pT_i][pT_j]) << "+(" << imag(V_n_vec[pT_i][pT_j]) << ") I";
				if ( pT_j + 1 == n_pT ) cout << "}";
				else cout << ", ";
			}
			cout << endl;
		}
		cout << "}" << endl << endl;



		if ( Vn_mode == 1 )
		{
			cout << endl << "v^2_2{2} = " << endl << "   { ";
			for (int pT_i = 0; pT_i < n_pT; pT_i++)
				cout << real(V_n_vec[pT_i][pT_i]) << "+(" << imag(V_n_vec[pT_i][pT_i]) << ") I"
						<< ( (pT_i + 1 == n_pT) ? "}" : ", " );
			cout << endl << endl;
		}
		else
		{
			cout << endl << "v^2_2{2} = " << endl << "   { ";
			for (int pT_i = 0; pT_i < n_pT; pT_i++)
				cout << real(V_n_vec[pT_i][pT_i])/average_Npairs_per_event[pT_i][pT_i]
						<< "+(" << imag(V_n_vec[pT_i][pT_i])/average_Npairs_per_event[pT_i][pT_i] << ") I"
						<< ( (pT_i + 1 == n_pT) ? "}" : ", " );
			cout << endl << endl;
		}



//if (1) exit(8);

		for (int pT_i = 0; pT_i < n_pT; pT_i++)
			cout << "Eigenvalue #" << pT_i << " = " << eigenvalues_EnsembleAverage[0][pT_i] << endl;
		cout << endl;
		for (int pT_i = 0; pT_i < n_pT; pT_i++)
		{
			cout << "Eigenvector #" << pT_i << " =" << endl << "{ ";
			for (int pT_j = 0; pT_j < n_pT; pT_j++)
				cout << eigenvectors_EnsembleAverage[0][pT_i][pT_j]
						<< ( (pT_j + 1 == n_pT) ? "}" : ", " );
			cout << endl;
		}
		cout << "<<<=====================================>>>" << endl;


		cout << "Now check jackknife error estimates:" << endl;
		// Same thing with errors using jackknife estimates
		const double prefactor = (N_total_events - 1.0) / static_cast<double>(N_total_events);
		for (int pT_i = 0; pT_i < n_pT; pT_i++)
		{
			complex<double> mean = 0.0;
			for (int iEvent = 0; iEvent < N_total_events; iEvent++)
				mean += jackknife_eigenvalues[iEvent][pT_i];
			mean /= static_cast<double>(N_total_events);

			double variance = 0.0;
			for (int iEvent = 0; iEvent < N_total_events; iEvent++)
			{
//				cout << "Jackknife check: " << iEvent << "   " << pT_i << "  " << jackknife_eigenvalues[iEvent][pT_i] << endl;
				variance += abs( ( jackknife_eigenvalues[iEvent][pT_i] - mean )
							* ( jackknife_eigenvalues[iEvent][pT_i] - mean ) );
			}
			variance *= prefactor;
			cout << "Eigenvalue #" << pT_i << " = "
					<< eigenvalues_EnsembleAverage[0][pT_i] << "; "
					<< mean << " +/- " << sqrt(variance) << endl;
		}
		cout << endl;

		bool store_output = true;
		//string output_name = "./datasets/pp-V1-rapgap-" + to_string(N_total_events) + "-" + to_string(n_pT) + "-" + to_string(Vn_mode) + ".dat";
		//string output_name = resultsDirectory + resultsFilename;
		ofstream output(output_name);

		for (int pT_i = 0; pT_i < n_pT; pT_i++)
		for (int pT_j = 0; pT_j < n_pT; pT_j++)
		{

			//if ( /*pT_i+1==n_pT and */pT_j==0 )
			//	for (int iEvent = 0; iEvent < N_total_events; iEvent++)
			//		cout << "Check eigenvector elements: " << pT_i << "   " << iEvent << "   "
			//				<< real(jackknife_eigenvectors[iEvent][pT_i][pT_j]) << "   "
			//				<< imag(jackknife_eigenvectors[iEvent][pT_i][pT_j]) << endl;

			complex<double> mean = 0.0;
			for (int iEvent = 0; iEvent < N_total_events; iEvent++)
				mean += jackknife_eigenvectors[iEvent][pT_i][pT_j];
			mean /= static_cast<double>(N_total_events);

			double variance = 0.0;
			double variance_re = 0.0, variance_im = 0.0;
			for (int iEvent = 0; iEvent < N_total_events; iEvent++)
			{
				complex<double> diff = jackknife_eigenvectors[iEvent][pT_i][pT_j] - mean;
				variance_re += real(diff) * real(diff);
				variance_im += imag(diff) * imag(diff);
				variance    += abs(diff) * abs(diff);
			}
			variance_re *= prefactor;
			variance_im *= prefactor;
			variance    *= prefactor;
			cout << "Eigenvector #" << pT_i << ", element #" << pT_j << " = "
					<< eigenvectors_EnsembleAverage[0][pT_i][pT_j] << "; "
					<< mean << " +/- " << sqrt(variance_re) << " +/- " << sqrt(variance_im) << " I" << endl;
			if ( store_output )
				output << pT_i << "   " << pT_j << "   "
					<< real(eigenvectors_EnsembleAverage[0][pT_i][pT_j]) << "   "
					<< imag(eigenvectors_EnsembleAverage[0][pT_i][pT_j]) << "   "
					<< real(mean) << "   " << sqrt(variance_re) << "   "
					<< sqrt(variance_im) << "   " << sqrt(variance) << endl;


		}
		if ( store_output ) { output.close(); }
	
		return;
	}




}


#endif
