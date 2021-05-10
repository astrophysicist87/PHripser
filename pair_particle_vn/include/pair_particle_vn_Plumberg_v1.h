#ifndef PAIR_PARTICLE_VN_PLUMBERG_H
#define PAIR_PARTICLE_VN_PLUMBERG_H

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
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

using namespace std;

namespace plumberg
{
	
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

	vector<vector<complex<double> > >
	V_n_vec( n_pT, vector<complex<double> >( n_pT, 0.0 ) );

	
	
	void arrange_data(
	  vector<vector<double> > &eta_vec,
	  vector<vector<double> > &phi_vec,
	  vector<vector<int> > &event_vec,
	  int &N_events,
	  int n_pT, 
	  double max_pT,
	  string filename )
	{
	  double delta_pT = max_pT/(n_pT);
	
	  // Chris added this line
	  pT_vec = vector<vector<double> >(n_pT, vector<double> (0));
	
	  ifstream input( filename );
	  string line;
	  int N_particles = 0;
	  int event;
	  while ( getline( input, line ) )
	  {
	    N_particles++;
	    string value;
	    vector<double> value_vector(0);
	    for( istringstream is (line); is >> value; )
	      value_vector.push_back( stod( value ) );
	    
	    // Sort into bins.
	    double pT = value_vector[0];
	    double eta = value_vector[1];
	    double phi = value_vector[2];
	    event = value_vector[3];
	
	    if(pT >= max_pT)
	    {
	      continue;
	    }
	    else
	    {
	    int bin = int((n_pT/max_pT)*pT);
	    eta_vec[bin].push_back(eta);
	    phi_vec[bin].push_back(phi);
	    pT_vec[bin].push_back(pT);
	    event_vec[bin].push_back(event);
	    }
	  }
	
	  int k;
	  for(k = 0; k < n_pT; k++)
	  {
	    if(eta_vec[k].size() <= 1000){cout << "Bin " << k << " contains <1000 datapoints" << endl;} 
	  }
	
	  input.close();
	
	  N_events = event + 1;
	
	  cout << N_particles << " " << N_events << endl;
	
	  //----------------------------------------------
	  // Chris added this part
	  all_events = vector<event_data>(N_events);
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
		}
	  }
	  //----------------------------------------------
	
	  return;
	}



	void generate_random_data(
	  vector<vector<double> > &eta_vec,
	  vector<vector<double> > &phi_vec,
	  vector<vector<int> > &event_vec,
	  int &N_events,
	  int n_pT, 
	  double max_pT )
	{
	  double delta_pT = max_pT/(n_pT);

		long chosen_multiplicity_per_event = 10;
		long chosen_number_of_events = 500;

	  // Chris added this line
	  pT_vec = vector<vector<double> >(n_pT, vector<double> (0));

	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	default_random_engine generator(seed);
	exponential_distribution<double> pT_distribution(1.0);
	normal_distribution<double> eta_distribution(0.0, 2.0);
	uniform_real_distribution<double> v2_distribution(0.2, 0.3);
	uniform_real_distribution<double> psi2_distribution(0.0, 2.0*M_PI);
	uniform_real_distribution<double> pphi_distribution(0.0, 2.0*M_PI);

	
	  //ifstream input( filename );
	  string line;
	  long N_particles = 0;
	  int event;
vector<double> intervals, weights(36);
for (int interval = 0; interval<=36; interval++)
	intervals.push_back( 2.0*M_PI*interval/36.0 );
	  //while ( getline( input, line ) )
      //while ( N_particles < chosen_number_of_events * chosen_multiplicity_per_event )
	  for (int iEvent = 0; iEvent < chosen_number_of_events; iEvent++)
		{
			//double v2 = v2_distribution(generator);
			double v2 = 0.25;
			double psi2 = psi2_distribution(generator);
		
		  for (int iParticle = 0; iParticle < chosen_multiplicity_per_event; iParticle++)
		  {
		    N_particles++;
		    string value;
		    vector<double> value_vector(0);
		    //for( istringstream is (line); is >> value; )
		    //  value_vector.push_back( stod( value ) );

			//const double pphi_value = pphi_distribution(generator);
			//const double pT_value = (1 + 2.0*v2*cos(2.0*(pphi_value - psi2)))*pT_distribution(generator);
			const double pT_value = pT_distribution(generator);
			const double v2_at_pT = v2*pT_value*pT_value / (1.0+pT_value*pT_value);
			const double psi2_at_pT = psi2+0.5*M_PI*pT_value*pT_value / (1.0+pT_value*pT_value);;
			for (int iWeight = 0; iWeight < 36; iWeight++)
				weights[iWeight] = 1.0 + 2.0*v2_at_pT*cos(2.0*( 0.5*(intervals[iWeight]+intervals[iWeight+1]) - psi2_at_pT));
			std::piecewise_constant_distribution<double> new_pphi_distribution ( intervals.begin(),
																					intervals.end(),
																					weights.begin());
			const double pphi_value = new_pphi_distribution(generator);

			value_vector.push_back( pT_value );
			value_vector.push_back( eta_distribution(generator) );
			value_vector.push_back( pphi_value );
			value_vector.push_back( iEvent );
		    
		    // Sort into bins.
		    double pT = value_vector[0];
		    double eta = value_vector[1];
		    double phi = value_vector[2];
		    event = value_vector[3];
		
		    if(pT >= max_pT)
		    {
		      continue;
		    }
		    else
		    {
		    int bin = int((n_pT/max_pT)*pT);
		    eta_vec[bin].push_back(eta);
		    phi_vec[bin].push_back(phi);
		    pT_vec[bin].push_back(pT);
		    event_vec[bin].push_back(event);
		    }
		  }
	}
	
	  int k;
	  for(k = 0; k < n_pT; k++)
	  {
	    if(eta_vec[k].size() <= 1000){cout << "Bin " << k << " contains <1000 datapoints" << endl;} 
	  }
	
	  //input.close();
	
	  N_events = event + 1;
	
	  cout << N_particles << " " << N_events << endl;
	
	  //----------------------------------------------
	  // Chris added this part
	  all_events = vector<event_data>(N_events);
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
		}
	  }
	  //----------------------------------------------
	
	  return;
	}
	
	
	//-----------------------
	void get_eigen( 
	  vector<vector<vector<double> > > & V_n_vec,
	  vector<vector<vector<double> > > & eigenvectors,
	  vector<vector<double> > & eigenvalues )
	{
		
		int Nevents = V_n_vec.size();
		int dim = V_n_vec[0][0].size();
		gsl_vector_complex *eval = gsl_vector_complex_alloc(dim);
		gsl_matrix_complex *evec = gsl_matrix_complex_alloc(dim,dim);
		double V_n_1D[dim*dim];

		for (int iEvent = 0; iEvent < Nevents; iEvent++)
		{
			int V_n_1D_index = 0;
			
			// Converts nxn V_n matrix to a 1D dataset for GSL functions
			for (int i=0; i<dim; i++)
			for (int j=0; j<dim; j++)
			{
				V_n_1D[V_n_1D_index] = V_n_vec[iEvent][i][j];
				V_n_1D_index++;
			 }
			
			// Solves eigensystem
			gsl_matrix_view m = gsl_matrix_view_array(V_n_1D, dim, dim);
			gsl_eigen_nonsymmv_workspace *w = gsl_eigen_nonsymmv_alloc(dim);
			gsl_eigen_nonsymmv (&m.matrix, eval, evec, w);
			gsl_eigen_nonsymmv_free(w);
			gsl_eigen_nonsymmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC);
			
			// Store results
			for (int i=0; i<dim; i++)
			{
				gsl_complex eval_i = gsl_vector_complex_get(eval, i);
				gsl_vector_complex_view evec_i = gsl_matrix_complex_column(evec, i);
		
			    double lambda = GSL_REAL(eval_i);
			    eigenvalues[iEvent][i] = lambda;
			
				for (int j = 0; j < dim; ++j)
			    {
					gsl_complex z = gsl_vector_complex_get( &evec_i.vector, j );	
					eigenvectors[iEvent][i][j] = GSL_REAL(z);
				}
			}
		}
		
		gsl_vector_complex_free(eval);
		gsl_matrix_complex_free(evec);
	
		return;
	}
	void get_eigen( 
	  vector<vector<vector<complex<double> > > > & V_n_vec,
	  vector<vector<vector<complex<double> > > > & eigenvectors,
	  vector<vector<double> > & eigenvalues )
	{
		
		int Nevents = V_n_vec.size();
		int dim = V_n_vec[0][0].size();
		gsl_vector *eval = gsl_vector_alloc(dim);
		gsl_matrix_complex *evec = gsl_matrix_complex_alloc(dim,dim);
		//complex<double>  V_n_1D[dim*dim];

		for (int iEvent = 0; iEvent < Nevents; iEvent++)
		{
			// Solves eigensystem
			//gsl_matrix_view m = gsl_matrix_view_array(V_n_1D, dim, dim);
			gsl_matrix_complex * m = gsl_matrix_complex_alloc (dim, dim);
			for (int i=0; i<dim; i++)
			for (int j=0; j<dim; j++)
			{
				complex<double> value = V_n_vec[iEvent][i][j];
				//complex<double> value = pow(i+j+2.0,5.0);
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
	void get_V_n_with_errors(
		int order, double eta_low, double eta_high, int Nevents, int n_pT,
		bool print_process = true )
	{
		//vector<vector<complex<double> > >
		//	V_n_vec( n_pT,
		//		vector<complex<double> >( n_pT, 0.0 ) );
	
		// This is a matrix with dimensions Nevents x n_pT x n_pT
		// to hold one jackknife estimate for each bin in each event
		vector<vector<vector<double> > >
			jackknife_V_n_estimates( Nevents, 
				vector<vector<double> > ( n_pT,
					vector<double>( n_pT, 0.0 ) ) );
		
		vector<vector<vector<double> > >
			jackknife_Npairs( Nevents,
				vector<vector<double> > ( n_pT,
					vector<double>( n_pT, 0.0 ) ) );
		vector<vector<double> >
			Npairs( n_pT,
				vector<double>( n_pT, 0.0 ) );
	
		//first, loop over all events
		for (int iEvent = 0; iEvent < Nevents; iEvent++)
		{
	
			event_data & this_event = all_events[iEvent];
	
			// Going through all possible bin-pairs.
			// Defining new z=0 and N_pairs=0.
			for (int pT_i = 0; pT_i < n_pT; pT_i++)
			for (int pT_j = 0; pT_j < n_pT; pT_j++)
			{
				double Npairs_this_bin = 0;
				complex<double> z(0.0, 0.0);
				
				vector<double> eta_bin_i = this_event.eta_vec[pT_i];
				vector<double> eta_bin_j = this_event.eta_vec[pT_j];
				vector<double> phi_bin_i = this_event.pphi_vec[pT_i];
				vector<double> phi_bin_j = this_event.pphi_vec[pT_j];
				
				const int n_part_i = eta_bin_i.size();
				const int n_part_j = eta_bin_j.size();
				
				// Going through all pairs of particles between bins.
				// Discarding if delta_eta is outside limits.
				// Counts N_pairs and adds exponential to z.
				if ( pT_i == pT_j )
				  for (int i = 0; i < n_part_i; i++)
				  for (int j = 0; j < n_part_j; j++)
				  {
					if (i==j) continue;
				    double delta_eta = eta_bin_i[i] - eta_bin_j[j];
				    //cout << delta_eta << endl;
				    if (     delta_eta*delta_eta >= eta_low*eta_low
						       and delta_eta*delta_eta <= eta_high*eta_high )
					{
				      //z += cos( order * ( phi_bin_i[i] - phi_bin_j[j] ) );
				      z += exp( complex<double> ( 0.0, order * ( phi_bin_i[i] - phi_bin_j[j] ) ) );
					  Npairs_this_bin++;
					}
				  } // Finished all pairs of particles within bin-pair
				else
				  for (int i = 0; i < n_part_i; i++)
				  for (int j = 0; j < n_part_j; j++)
				  {
				    double delta_eta = eta_bin_i[i] - eta_bin_j[j];
				    if (     delta_eta*delta_eta >= eta_low*eta_low
						       and delta_eta*delta_eta <= eta_high*eta_high )
					{
				      z += exp( complex<double> ( 0.0, order * ( phi_bin_i[i] - phi_bin_j[j] ) ) );
					  Npairs_this_bin++;
					}
				  } // Finished all pairs of particles within bin-pair
	
				// After computing the contribution to this bin from this event,
				// add it to all but the current event's V_n estimates
				for (int jEvent = 0; jEvent < iEvent; jEvent++)
				{
					jackknife_V_n_estimates[jEvent][pT_i][pT_j] += real(z);
					jackknife_Npairs[jEvent][pT_i][pT_j] += Npairs_this_bin;
				}
				for (int jEvent = iEvent + 1; jEvent < Nevents; jEvent++)
				{
					jackknife_V_n_estimates[jEvent][pT_i][pT_j] += real(z);
					jackknife_Npairs[jEvent][pT_i][pT_j] += Npairs_this_bin;
				}
	
				V_n_vec[pT_i][pT_j] += z;
				Npairs[pT_i][pT_j] += Npairs_this_bin;
		
			} // Finished all bin-pairs for this event			
		}
	
		// Normalize jackknife estimates
		for (int iEvent = 0; iEvent < Nevents; iEvent++)
		for (int pT_i = 0; pT_i < n_pT; pT_i++)
		for (int pT_j = 0; pT_j < n_pT; pT_j++)
		{
			jackknife_V_n_estimates[iEvent][pT_i][pT_j] /= ( Nevents - 1 );
			jackknife_Npairs[iEvent][pT_i][pT_j] /= ( Nevents - 1 );
			//jackknife_V_n_estimates[iEvent][pT_i][pT_j] /= jackknife_Npairs[iEvent][pT_i][pT_j];
		}
			
		for (int pT_i = 0; pT_i < n_pT; pT_i++)
		for (int pT_j = 0; pT_j < n_pT; pT_j++)
		{
			V_n_vec[pT_i][pT_j] /= Nevents;
			Npairs[pT_i][pT_j] /= Nevents;
		}
	
	
		// Variance of jackknife estimates in each bin
		// gives error estimate for that bin
		//cout << "Vn = " << endl;
		const double prefactor = (Nevents - 1.0) / static_cast<double>(Nevents);
		for (int pT_i = 0; pT_i < n_pT; pT_i++)
		for (int pT_j = 0; pT_j < n_pT; pT_j++)
		{
			double mean = 0.0;
			for (int iEvent = 0; iEvent < Nevents; iEvent++)
				mean += jackknife_V_n_estimates[iEvent][pT_i][pT_j];
			mean /= Nevents;

			double variance = 0.0;
			for (int iEvent = 0; iEvent < Nevents; iEvent++)
				variance += abs(( jackknife_V_n_estimates[iEvent][pT_i][pT_j] - mean )
							* ( jackknife_V_n_estimates[iEvent][pT_i][pT_j] - mean ));
//				variance += abs(( jackknife_V_n_estimates[iEvent][pT_i][pT_j] - V_n_vec[pT_i][pT_j] )
//							* ( jackknife_V_n_estimates[iEvent][pT_i][pT_j] - V_n_vec[pT_i][pT_j] ));
			variance *= prefactor;
			//cout << "Vn in bin (" << pT_i << ", " << pT_j << ") = "
			//		<< V_n_vec[pT_i][pT_j] << "; "
			//		<< mean << " +/- " << sqrt(variance) << endl;
			//cout << real(V_n_vec[pT_i][pT_j]) << " + " << imag(V_n_vec[pT_i][pT_j]) << " I,";
			
		}

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
	
		/*for (int pT_i = 0; pT_i < n_pT; pT_i++)
		for (int pT_j = 0; pT_j < n_pT; pT_j++)
		{
			if (pT_i!=pT_j) continue;
	
			double mean = 0.0;
			for (int iEvent = 0; iEvent < Nevents; iEvent++)
				mean += jackknife_V_n_estimates[iEvent][pT_i][pT_j]
						/ double(jackknife_Npairs[iEvent][pT_i][pT_j]);
			mean /= Nevents;

			double variance = 0.0;
			for (int iEvent = 0; iEvent < Nevents; iEvent++)
			{
				const complex<double> diff = ( jackknife_V_n_estimates[iEvent][pT_i][pT_j]
												/double(jackknife_Npairs[iEvent][pT_i][pT_j]) - mean );
				variance += abs(diff*diff);
			}
			variance *= prefactor;
			cout << "v^2_2{2} in bin (" << pT_i << ", " << pT_j << ") = "
					<< V_n_vec[pT_i][pT_j]/double(Npairs[pT_i][pT_j]) << "; "
					<< mean << " +/- " << sqrt(variance) << endl;
		}*/
	



		/*
		//--------------------------------------------------------------
		//--------------------------------------------------------------
		//--------------------------------------------------------------
		// Next, pass jackknife estimates to eigendecomposition routine
		vector<vector<vector<double> > >
			jackknife_eigenvectors( Nevents,  vector<vector<double> > ( n_pT, vector<double>( n_pT, 0.0 ) ) );
		vector<vector<double> >
			jackknife_eigenvalues( Nevents, vector<double>( n_pT, 0.0 ) );


		// diagonalize each jackknife matrix and store the resulting eigensystems
		get_eigen( jackknife_V_n_estimates, jackknife_eigenvectors, jackknife_eigenvalues );

		// ensure eigenvector signs are consistent for each estimate and average
		for (int iEvent = 0; iEvent < Nevents; iEvent++)
		for (int pT_i = 0; pT_i < n_pT; pT_i++)
			if ( jackknife_eigenvectors[iEvent][pT_i][0] < 0.0 )
				for (int pT_j = 0; pT_j < n_pT; pT_j++)
					jackknife_eigenvectors[iEvent][pT_i][pT_j] *= -1.0;
		*/


		/*
		// Makes it easier to pass this to get eigen
		vector<vector<vector<complex<double> > > > V_n_vec_EnsembleAverage;
		vector<vector<vector<double> > >
			eigenvectors_EnsembleAverage( 1, vector<vector<double> > ( n_pT, vector<double>( n_pT, 0.0 ) ) );
		vector<vector<double> >
			eigenvalues_EnsembleAverage( 1, vector<double>( n_pT, 0.0 ) );
		V_n_vec_EnsembleAverage.push_back( V_n_vec );
		
		// diagonalize and get eigensystem for the average V_n over full ensemble
		get_eigen( V_n_vec_EnsembleAverage, eigenvectors_EnsembleAverage, eigenvalues_EnsembleAverage );

		// ensure eigenvector signs are consistent for each estimate and average
		for (int pT_i = 0; pT_i < n_pT; pT_i++)
			if ( eigenvectors_EnsembleAverage[0][pT_i][0] < 0.0 )
				for (int pT_j = 0; pT_j < n_pT; pT_j++)
					eigenvectors_EnsembleAverage[0][pT_i][pT_j] *= -1.0;

		for (int pT_i = 0; pT_i < n_pT; pT_i++)
			cout << "Eigenvalue #" << pT_i << " = " << eigenvalues_EnsembleAverage[0][pT_i] << endl;
		cout << endl;
		for (int pT_i = 0; pT_i < n_pT; pT_i++)
		for (int pT_j = 0; pT_j < n_pT; pT_j++)
			cout << "Eigenvector #" << pT_i << ", element #" << pT_j << " = "
					<< eigenvectors_EnsembleAverage[0][pT_i][pT_j] << endl;

		*/



		/*
		// Print out results
		for (int pT_i = 0; pT_i < n_pT; pT_i++)
		{
			double mean = 0.0;
			for (int iEvent = 0; iEvent < Nevents; iEvent++)
				mean += jackknife_eigenvalues[iEvent][pT_i];
			mean /= Nevents;

			double variance = 0.0;
			for (int iEvent = 0; iEvent < Nevents; iEvent++)
				variance += abs( ( jackknife_eigenvalues[iEvent][pT_i] - mean )
							* ( jackknife_eigenvalues[iEvent][pT_i] - mean ) );
			variance *= prefactor;
			cout << "Eigenvalue #" << pT_i << " = "
					<< eigenvalues_EnsembleAverage[0][pT_i] << "; "
					<< mean << " +/- " << sqrt(variance) << endl;
		}
		cout << endl;
		for (int pT_i = 0; pT_i < n_pT; pT_i++)
		for (int pT_j = 0; pT_j < n_pT; pT_j++)
		{

			double mean = 0.0;
			for (int iEvent = 0; iEvent < Nevents; iEvent++)
				mean += jackknife_eigenvectors[iEvent][pT_i][pT_j];
			mean /= Nevents;

			double variance = 0.0;
			for (int iEvent = 0; iEvent < Nevents; iEvent++)
				variance += abs( ( jackknife_eigenvectors[iEvent][pT_i][pT_j] - mean )
							* ( jackknife_eigenvectors[iEvent][pT_i][pT_j] - mean ) );
			variance *= prefactor;
			cout << "Eigenvector #" << pT_i << ", element #" << pT_j << " = "
					<< eigenvectors_EnsembleAverage[0][pT_i][pT_j] << "; "
					<< mean << " +/- " << sqrt(variance) << endl;
		}
		*/


		// check jackknife estimates directly
		/*
		cout << "Check jackknife estimates directly: " << endl;
		for (int iEvent = 0; iEvent < Nevents; iEvent++)
		{
			for (int pT_i = 0; pT_i < n_pT; pT_i++)
			{
				for (int pT_j = 0; pT_j < n_pT; pT_j++)
					cout << jackknife_V_n_estimates[iEvent][pT_i][pT_j] << "   ";
				cout << endl;
			}
			//cout << endl;
			for (int pT_i = 0; pT_i < n_pT; pT_i++)
			{
				cout << " --> ";
				for (int pT_j = 0; pT_j < n_pT; pT_j++)
					cout << jackknife_eigenvectors[iEvent][pT_i][pT_j] << "   ";
				cout << endl;
			}
			cout << endl;
		}
		*/


	
		return;
	}
	
	





















	//-----------------------
	void get_V_n(
		int order, double eta_low, double eta_high, int Nevents, int n_pT,
		bool print_process = true )
	{
		//vector<vector<complex<double> > >
		//	V_n_vec( n_pT, vector<complex<double> >( n_pT, 0.0 ) );
		
		// Going through all possible bin-pairs.
		// Defining new z=0 and N_pairs=0.
		for (int pT_i = 0; pT_i < n_pT; pT_i++)
		for (int pT_j = 0; pT_j < n_pT; pT_j++)
		{
			cout << "Loop: pT_i = " << pT_i << ", pT_j = " << pT_j << endl;

			complex<double> z_i(0.0, 0.0);
			complex<double> z_j(0.0, 0.0);
			complex<double> z_ij(0.0, 0.0);
				
			// loop over all events
			for (int iEvent = 0; iEvent < Nevents; iEvent++)
			{
		
				event_data & this_event = all_events[iEvent];
	
				vector<double> eta_bin_i = this_event.eta_vec[pT_i];
				vector<double> eta_bin_j = this_event.eta_vec[pT_j];
				vector<double> phi_bin_i = this_event.pphi_vec[pT_i];
				vector<double> phi_bin_j = this_event.pphi_vec[pT_j];
				
				const int n_part_i = eta_bin_i.size();
				const int n_part_j = eta_bin_j.size();

				for (int i = 0; i < n_part_i; i++)
					z_i += exp( complex<double> ( 0.0, order * phi_bin_i[i] ) );
				for (int j = 0; j < n_part_j; j++)
					z_j += exp( complex<double> ( 0.0, -order * phi_bin_j[j] ) );
				
				// Going through all pairs of particles between bins.
				// Discarding if delta_eta is outside limits.
				// Counts N_pairs and adds exponential to z.
				for (int i = 0; i < n_part_i; i++)
				for (int j = 0; j < n_part_j; j++)
				{
					if ( pT_i==pT_j and i==j ) continue;
				    //double delta_eta = eta_bin_i[i] - eta_bin_j[j];
					//if ( delta_eta*delta_eta >= eta_low*eta_low
					//	       and delta_eta*delta_eta <= eta_high*eta_high )
						z_ij += exp( complex<double> ( 0.0, order * ( phi_bin_i[i] - phi_bin_j[j] ) ) );
				} // Finished all pairs of particles within bin-pair
			} // Finished all bin-pairs for this event
		
			z_i /= Nevents;
			z_j /= Nevents;
			z_ij /= Nevents;

			V_n_vec[pT_i][pT_j] = z_ij - z_i*z_j;
			cout << "CHECK: " << z_ij << "   " << z_i << "   " << z_j << "   " << V_n_vec[pT_i][pT_j] << endl;
		}
		
		/*for (int pT_i = 0; pT_i < n_pT; pT_i++)
		for (int pT_j = 0; pT_j < n_pT; pT_j++)
		{
			V_n_vec[pT_i][pT_j] /= Nevents;
			Npairs[pT_i][pT_j] /= Nevents;
		}*/
	

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
	


		// Makes it easier to pass this to get eigen
		vector<vector<vector<complex<double> > > > V_n_vec_EnsembleAverage;
		vector<vector<vector<complex<double> > > >
			eigenvectors_EnsembleAverage( 1, vector<vector<complex<double> > > ( n_pT, vector<complex<double> >( n_pT, 0.0 ) ) );
		vector<vector<double> >
			eigenvalues_EnsembleAverage( 1, vector<double>( n_pT, 0.0 ) );
		V_n_vec_EnsembleAverage.push_back( V_n_vec );
		
		// diagonalize and get eigensystem for the average V_n over full ensemble
		get_eigen( V_n_vec_EnsembleAverage, eigenvectors_EnsembleAverage, eigenvalues_EnsembleAverage );

		for (int pT_i = 0; pT_i < n_pT; pT_i++)
			cout << "Eigenvalue #" << pT_i << " = " << eigenvalues_EnsembleAverage[0][pT_i] << endl;
		/*cout << endl;
		for (int pT_i = 0; pT_i < n_pT; pT_i++)
		for (int pT_j = 0; pT_j < n_pT; pT_j++)
			cout << "Eigenvector #" << pT_i << ", element #" << pT_j << " = "
					<< eigenvectors_EnsembleAverage[0][pT_i][pT_j] << endl;*/


	
		return;
	}




}


#endif
