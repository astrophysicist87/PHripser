#ifndef PAIR_PARTICLE_VN_H
#define PAIR_PARTICLE_VN_H

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <complex>
#include <vector>
#include <iomanip>
#include <stdio.h>
#include <algorithm>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

using namespace std;

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



void get_V_n(
	vector<vector<double> > &V_n_vec,
  vector<vector<double> > &eta_vec,
  vector<vector<double> > &phi_vec,
  int order, double eta_low, double eta_high, bool print_process = true)
{
  int n_pT = V_n_vec.size();

  // Going through all possible bin-pairs.
  // Defining new z=0 and N_pairs=0.
  for (int pT_i = 0; pT_i < n_pT; pT_i++){
    if(print_process == true){cout << "Covariance Matrix Process : " << pT_i << "/" << n_pT << endl;}
  for (int pT_j = 0; pT_j < n_pT; pT_j++)
  {
    //cout << pT_i << " " << pT_j << endl;
    complex<double> z(0.0, 0.0);
    long N_pairs = 0;

    vector<double> eta_bin_i = eta_vec[pT_i];
    vector<double> eta_bin_j = eta_vec[pT_j];
    vector<double> phi_bin_i = phi_vec[pT_i];
    vector<double> phi_bin_j = phi_vec[pT_j];

    const int n_part_i = eta_bin_i.size();
    const int n_part_j = eta_bin_j.size();

    // Going through all pairs of particles between bins.
    // Discarding if delta_eta is outside limits.
    // Counts N_pairs and adds exponential to z.
    if ( pT_i == pT_j )
      for (int i = 0  ; i < n_part_i; i++)
      for (int j = i+1; j < n_part_j; j++)
      {
        double delta_eta = eta_bin_i[i] - eta_bin_j[j];
        //cout << delta_eta << endl;
        if (     delta_eta*delta_eta >= eta_low*eta_low
			       and delta_eta*delta_eta <= eta_high*eta_high )
        {
          z += cos( order * ( phi_bin_i[i] - phi_bin_j[j] ) );
          N_pairs++;
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
          N_pairs++;
        }
      } // Finished all pairs of particles within bin-pair

    // Normalize and store.
    z /= N_pairs;
    double V_n = abs(z);
    //if(isnan(V_n)==true){V_n_vec[pT_i][pT_j] = 0;}
    //else{
	V_n_vec[pT_i][pT_j] = V_n;
	//}

  }} // Finished all bin-pairs
}



void get_eigen( 
  vector<vector<double> > &V_n_vec,
  vector<vector<double> > &eigenvectors,
  vector<double> &eigenvalues)
{
	
	int dim = V_n_vec[0].size();
	double V_n_1D[dim*dim];
	int V_n_1D_index = 0;
	
	// Converts nxn mV_n matrix to a 1D dataset for GSL functions
	for (int i=0; i<dim; i++)
	for (int j=0; j<dim; j++)
	{
		V_n_1D[V_n_1D_index]=V_n_vec[i][j];
		V_n_1D_index++;
	 }
	
	// Solves eigensystem
	gsl_matrix_view m = gsl_matrix_view_array(V_n_1D, dim, dim);
	gsl_vector_complex *eval = gsl_vector_complex_alloc(dim);
	gsl_matrix_complex *evec = gsl_matrix_complex_alloc(dim,dim);
	gsl_eigen_nonsymmv_workspace *w = gsl_eigen_nonsymmv_alloc(dim);
	gsl_eigen_nonsymmv (&m.matrix, eval, evec, w);
	gsl_eigen_nonsymmv_free(w);
	gsl_eigen_nonsymmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC);
	
	// Prints results out neatly
	for (int i=0; i<dim; i++)
	{
		gsl_complex eval_i = gsl_vector_complex_get(eval, i);
		gsl_vector_complex_view evec_i = gsl_matrix_complex_column(evec, i);

    double lambda = GSL_REAL(eval_i);
    eigenvalues[i] = lambda;

    // printf("eigenvalue = %g + %gi\n", GSL_REAL(eval_i), GSL_IMAG(eval_i));
    // printf("eigenvector = \n");
		
		for (int j = 0; j < dim; ++j)
    {
			gsl_complex z = gsl_vector_complex_get( &evec_i.vector, j );
      double z_real = GSL_REAL(z);

      eigenvectors[i][j] = z_real;
      // printf("%g + %gi\n", GSL_REAL(z), GSL_IMAG(z));
	}}
	
	gsl_vector_complex_free(eval);
	gsl_matrix_complex_free(evec);

	return;
}


void get_eigenmode(
  vector<vector<double> > eigenvectors,
  vector<double> eigenvalues,
  int mode,
  string data_title,
  double &eigenvalue_mode,
  vector<double> &eigenvector_mode,
  bool write_result = true)
{
  vector<double> sorted_eigen = eigenvalues;
  sort(sorted_eigen.begin(), sorted_eigen.end(), greater<double>());
  double mode_value = sorted_eigen[mode];

  int index, mode_index;
  for(index = 0; index < eigenvalues.size(); index++){
    if(mode_value == eigenvalues[index]){
      mode_index = index;
      break;
    }
  }

  eigenvalue_mode = eigenvalues[mode_index];
  eigenvector_mode = eigenvectors[mode_index];

  int i;
  // Multiplying each element of the eigenvector with the eigenvalue of the mode
  for(i = 0; i < eigenvector_mode.size(); i++){
    eigenvector_mode[i] *= eigenvalue_mode;}

  // cout << eigenvalues[mode_index] << endl;

  if(write_result == true){
    ofstream output(data_title, ios::out);
    const int w = 12;
    for(i = 0; i < eigenvector_mode.size(); i++)
    {
      output << setprecision(8);
      output << fixed
      << setw(w) << eigenvector_mode[i] << endl;
    }
    output.close();
  }

  return;
}



void error(
  vector<vector<double> > eta_vec,
  vector<vector<double> > phi_vec,
  vector<vector<int> > event_vec,
  vector<double> eigenvector_average,
  vector<vector<double> > V_n_vec,
  int N_events, int n_pT,
  int order, double eta_low, double eta_high)
{
  vector<vector<vector<double> > > V_n_single_event(N_events, vector<vector<double> > (n_pT,vector<double> (n_pT)));
  vector<vector<vector<double> > > V_n_excl_event(N_events, vector<vector<double> > (n_pT,vector<double> (n_pT)));
  vector<int> event_counts(n_pT);
  int ev = 0;
  int i, j, pT_i, current_event;

  // Goes through each event, N out of N_events
  for(current_event = 0; current_event < N_events; current_event++){
    // if(current_event%100 == 0){cout <<  "Error process - loop 1 : " << current_event << "/" << N_events << endl;}
    // cout << endl << "current_event = " << current_event << endl;

    // Copies eta- and phi-values for current event
    vector<vector<double> > eta_copy(n_pT, vector<double> (0));
    vector<vector<double> > phi_copy(n_pT, vector<double> (0));
    for(pT_i = 0; pT_i < n_pT; pT_i++){
      ev = event_counts[pT_i];
      while(event_vec[pT_i][ev] == current_event){
        eta_copy[pT_i].push_back(eta_vec[pT_i][ev]);
        phi_copy[pT_i].push_back(phi_vec[pT_i][ev]);
        ev++;
      }
      event_counts[pT_i] = ev;
    }

    get_V_n(V_n_single_event[current_event], eta_copy, phi_copy, order, eta_low, eta_high, false);

    // Adds contribution to flow matrix to all matrices but current_event
    int matrix;
    for(matrix = 0; matrix < N_events; matrix++){
      if(matrix==current_event){continue;}
      for(i = 0; i < n_pT; i++)
      for(j = 0; j < n_pT; j++)
      {
        V_n_excl_event[current_event][i][j] += V_n_single_event[current_event][i][j];
      }
    }

    // int w = 14;
    // cout << endl << endl;
    // cout << setprecision(8);
    // for(int i = 0 ; i < int(V_n_excl_event[current_event].size()) ; i++)
    // {
    //   for(int j = 0 ; j < int(V_n_excl_event[current_event][i].size()) ; j++)
    //     cout << setw(w) << V_n_excl_event[current_event][i][j] << " ";
    //   cout << endl;
    // }

  } // Finished 1st loop


  vector<double> variance(n_pT);

  // 2nd loop over N_events : Will gdo eigendecomposition and find leading mode, has a running sum of variance for each eigen-dimension
  for(current_event = 0; current_event < N_events; current_event++){
    // if(current_event%100 == 0){cout <<  "Error process - loop 2 : " << current_event << "/" << N_events << endl;}

    vector<vector<double> > eigenvectors(n_pT, vector<double> (n_pT));
    vector<double> eigenvalues(n_pT);
    get_eigen(V_n_excl_event[current_event], eigenvectors, eigenvalues);

    int mode = 0;
    string data_title = "./datasets/eigenvector_error_500N_6pT.dat";
    double eigenvalue_mode;
    vector<double> eigenvector_mode;
    get_eigenmode(eigenvectors, eigenvalues, mode, data_title, eigenvalue_mode, eigenvector_mode, true);

    for(i = 0; i < n_pT; i++)
    {
      variance[i] += (eigenvector_mode[i]-eigenvector_average[i])*(eigenvector_mode[i]-eigenvector_average[i]);
    }


  } // Finished 2nd loop


  // Divides variance by N_events
  for(i = 0; i < n_pT; i++)
  {
    variance[i] /= N_events;
    variance[i] = sqrt(variance[i]);
    cout << "Eigenvector variance element " << i << " = " << variance[i] << endl;
  }


  // Writes the variance to a file
  ofstream output("./datasets/eigen_error_500N_6pT.dat", ios::out);
  const int w = 12;
  for(i = 0; i < n_pT; i++)
  {
    output << setprecision(8);
    output << fixed
    << setw(w) << variance[i] << endl;
  }
  output.close();


  return;

}








#endif
