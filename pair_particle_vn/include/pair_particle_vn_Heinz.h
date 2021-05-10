#ifndef PAIR_PARTICLE_VN_HEINZ_H
#define PAIR_PARTICLE_VN_HEINZ_H

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
#include "sampled_distribution.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

using namespace std;

namespace heinz
{
	bool use_seed = true;
	int seed_index = -1;

	const int n_pT = 6;
	const double max_pT = 3.0;
	const int order = 2;
	const double eta_low = 0.0;
	const double eta_high = 2.0;

	long Nevents;
	const long chosen_multiplicity_per_event = 10;
	const long chosen_number_of_events = 100;
	
	typedef struct
	{
		int eventID, pTbin;
		double pT, pphi, eta;
		
	} particle_data;
	
	
	typedef struct
	{
		int eventID;
		int Nparticles;
	
		//vector<particle_data> particles_in_this_event;
		vector<vector<double> > eta_vec;
		vector<vector<double> > pphi_vec;
		vector<vector<double> > pT_vec;
	} event_data;
	

	vector<event_data> all_events;
	
	vector<vector<double> > pT_vec;
	vector<vector<double> > eta_vec;
	vector<vector<double> > phi_vec;
	vector<vector<int> > event_vec;

	vector<vector<long> >
		Npairs( n_pT, vector<long>( n_pT, 0.0 ) );
	vector<vector<double> >
		average_Npairs_per_event( n_pT, vector<double>( n_pT, 0.0 ) );

	vector<vector<complex<double> > >
		V_n_vec( n_pT, vector<complex<double> >( n_pT, 0.0 ) );

	vector<vector<complex<double> > >
		Q_i( n_pT, vector<complex<double> >( n_pT, 0.0 ) );
	vector<vector<complex<double> > >
		Q_j( n_pT, vector<complex<double> >( n_pT, 0.0 ) );
	vector<vector<complex<double> > >
		Q_ij( n_pT, vector<complex<double> >( n_pT, 0.0 ) );

	void generate_random_data()
	{
		double delta_pT = max_pT/(n_pT);

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


		for (int iEvent = 0; iEvent < chosen_number_of_events; iEvent++)
		{
	
			//double v2 = v2_distribution(generator);
			double v2 = 0.25;
			double psi2 = psi2_distribution(generator);
			double psi2_pTslope = psi2_pTslope_distribution(generator);
			
			for (int iParticle = 0; iParticle < chosen_multiplicity_per_event; iParticle++)
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
				else
				{
					int bin = int(pT/delta_pT);
					eta_vec[bin].push_back(eta);
					phi_vec[bin].push_back(phi);
					pT_vec[bin].push_back(pT);
					event_vec[bin].push_back(event);
					//cout << "Check: " << pT << "   " << eta << "   " << phi << "   " << event << "   " << bin << endl;
				}
			}
		}
						
		//----------------------------------------------
		all_events.clear();
		all_events = vector<event_data>(chosen_number_of_events);
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
/*cout << "CHECK: " << bin << "   " << event_ID_of_this_particle << "   "
					<< pT_vec[bin][iParticle] << "   "
					<< eta_vec[bin][iParticle] << "   "
					<< phi_vec[bin][iParticle] << endl;*/
			}
		}
		//----------------------------------------------
		
		return;
	}
	
	
	//-----------------------
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
	void update_Q_n()
	{
		Nevents += chosen_number_of_events;

		// loop over all events
		for ( auto & this_event : all_events )
		{

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
					Q_ij_this_event[pT_i][pT_j] += exp( complex<double> ( 0.0, order * ( phi_bin_i[i] - phi_bin_j[j] ) ) );
					Npairs[pT_i][pT_j]++;
					Npairs_this_event[pT_i][pT_j]++;
				}

				// Normalize averages
				if ( n_part_i > 0 )
					Q_i[pT_i][pT_j] += Q_i_this_event[pT_i][pT_j] / double(n_part_i);
				if ( n_part_j > 0 )
					Q_j[pT_i][pT_j] += Q_j_this_event[pT_i][pT_j] / double(n_part_j);
				if ( Npairs_this_event[pT_i][pT_j] > 0 )
					Q_ij[pT_i][pT_j] += Q_ij_this_event[pT_i][pT_j] / double(Npairs_this_event[pT_i][pT_j]);

			}

		}

		return;
	}



	//-----------------------
	void get_V_n_and_eigendecomposition()
	{
		
		cout << endl << endl << endl;
		
		// Going through all possible bin-pairs.
		// Defining new z=0 and N_pairs=0.
		for (int pT_i = 0; pT_i < n_pT; pT_i++)
		for (int pT_j = 0; pT_j < n_pT; pT_j++)
		{
			//cout << "Loop: pT_i = " << pT_i << ", pT_j = " << pT_j << endl;
			//cout << " --> " << Q_ij[pT_i][pT_j] << "   " << Q_i[pT_i][pT_j] << "   " << Q_j[pT_i][pT_j] << endl;

			Q_i[pT_i][pT_j] /= static_cast<double>(Nevents);
			Q_j[pT_i][pT_j] /= static_cast<double>(Nevents);
			Q_ij[pT_i][pT_j] /= static_cast<double>(Nevents);

			average_Npairs_per_event[pT_i][pT_j] = Npairs[pT_i][pT_j] / static_cast<double>(Nevents);

			V_n_vec[pT_i][pT_j] = Q_ij[pT_i][pT_j] - Q_i[pT_i][pT_j]*Q_j[pT_i][pT_j];
			//cout << " --> " << Q_ij[pT_i][pT_j] << "   " << Q_i[pT_i][pT_j] << "   " << Q_j[pT_i][pT_j] << "   "
			//		<< V_n_vec[pT_i][pT_j] << "   " << static_cast<double>(Nevents) << endl;
		}	

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
		cout << "}" << endl;

//if (1) exit(8);

		///*
		// Print flow estimates from V_n's diagonal
		cout << endl << "v^2_2{2} = " << endl << "   { ";
		for (int pT_i = 0; pT_i < n_pT; pT_i++)
			cout << real(V_n_vec[pT_i][pT_i]) << "+(" << imag(V_n_vec[pT_i][pT_i]) << ") I"
					<< ( (pT_i + 1 == n_pT) ? "}" : ", " );
		cout << endl << endl;
		//*/

		// Makes it easier to pass this to get eigen
		vector<vector<vector<complex<double> > > > V_n_vec_EnsembleAverage;
		vector<vector<vector<complex<double> > > >
			eigenvectors_EnsembleAverage( 1, vector<vector<complex<double> > > ( n_pT, vector<complex<double> >( n_pT, 0.0 ) ) );
		vector<vector<double> >
			eigenvalues_EnsembleAverage( 1, vector<double>( n_pT, 0.0 ) );
		V_n_vec_EnsembleAverage.push_back( V_n_vec );
		
		// diagonalize and get eigensystem for the average V_n over full ensemble
		get_eigen( V_n_vec_EnsembleAverage, eigenvectors_EnsembleAverage, eigenvalues_EnsembleAverage );

		// Print results for complex V_n
		cout << "<<<=====================================>>>" << endl;
		cout << "Complex V_n: " << endl;
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
		cout << "<<<=====================================>>>" << endl << endl;


		// Print V_n matrix
		cout << "r_n(pT_i, pT_j):" << endl;
		cout << "{";
		for (int pT_i = 0; pT_i < n_pT; pT_i++)
		{
			cout << "{";
			for (int pT_j = 0; pT_j < n_pT; pT_j++)
			{
				cout << V_n_vec[pT_i][pT_j] / sqrt(abs(V_n_vec[pT_i][pT_i]*V_n_vec[pT_j][pT_j]));
				if ( pT_j + 1 == n_pT ) cout << "}";
				else cout << ", ";
			}
			cout << endl;
		}
		cout << "}" << endl;



		/*
		// Re-do for real V_n
		for (int pT_i = 0; pT_i < n_pT; pT_i++)
		for (int pT_j = 0; pT_j < n_pT; pT_j++)
			V_n_vec[pT_i][pT_j] = real( V_n_vec[pT_i][pT_j] );

		V_n_vec_EnsembleAverage.clear();
		eigenvectors_EnsembleAverage.clear();
		eigenvalues_EnsembleAverage.clear();

		
		eigenvectors_EnsembleAverage = vector<vector<vector<complex<double> > > >( 1,
										vector<vector<complex<double> > > ( n_pT,
										vector<complex<double> >( n_pT, 0.0 ) ) );
		
		eigenvalues_EnsembleAverage = vector<vector<double> >( 1,
										vector<double>( n_pT, 0.0 ) );
		V_n_vec_EnsembleAverage.push_back( V_n_vec );


		// diagonalize and get eigensystem for the average V_n over full ensemble
		get_eigen( V_n_vec_EnsembleAverage, eigenvectors_EnsembleAverage, eigenvalues_EnsembleAverage );

		// Print results for real V_n
		cout << "<<<=====================================>>>" << endl;
		cout << "Real V_n: " << endl;
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
		*/
	
		return;
	}




}


#endif
