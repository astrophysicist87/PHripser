#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

//#include <gsl/gsl_sf_gamma.h>

using namespace std;

vector<long double> log_cache, log_nFact_kFact_cache;

long double log_nFact_kFact( long int n, long int k )
{
	double log_k = 0.0;
	//for (int i = 1; i <= k; i++) log_k += log( i );
	for (int i = 1; i <= k; i++) log_k += log_cache[i];
	double log_n = log_k;
	//for (int i = k + 1; i <= n; i++) log_n += log( i );
	for (int i = k + 1; i <= n; i++) log_n += log_cache[i];
	return ( log_k / log_n );
}

long double binomial_coefficient( long int n, long int k )
{
	long int n0 = n, k0 = k, d = n-k, d0 = n-k;
	long double result = 1.0;
	for (int i = 1; i <= k0; i++) result *= static_cast<long double>(n--)
											/static_cast<long double>(k--);
	for (int i = 1; i <= d0; i++) result *= static_cast<long double>(n--)
											/static_cast<long double>(d--);
	return ( result );
}

/*long double binomial_coefficient( long int n, long int k )
{
	long int n0 = n, k0 = k, d = n-k, d0 = n-k;
	long double result = 0.0;
	for (int i = 1; i <= n0; i++) result += log_cache[i];
	for (int i = 1; i <= d0; i++) result -= log_cache[i];
	for (int i = 1; i <= k0; i++) result -= log_cache[i];

	return ( exp( result ) );
}*/

int main(int argc, char** argv)
{
	cout << setprecision(16);
	const long int n = (argc>1) ? std::stoi(argv[1]) : 100;
	log_cache.resize( n + 1 );
	log_nFact_kFact_cache.resize( n + 1 );
	for ( int i = 1; i <= n; i++ ) log_cache[i] = log( i );
	for ( int k = 1; k <= n; k++ ) log_nFact_kFact_cache[k] = log_nFact_kFact( n, k );

	for ( int MM = n; MM > 1; MM-- )
	{
		const long double prefactor = MM / binomial_coefficient( n-1, MM-1 );
		//const long double prefactor = MM / gsl_sf_choose( n-1, MM-1 );
		long double sum = 0.0;
		for ( int k = 1; k < n; k++ ) sum += prefactor * log_nFact_kFact_cache[k]
												* binomial_coefficient( n-k-1, MM-2 );
		//for ( int k = 1; k < n; k++ ) sum += prefactor * log_nFact_kFact_cache[k]
		//										* gsl_sf_choose( n-k-1, MM-2 );
		cout << (static_cast<double>(n) - MM)/(n-1.0) << "   " << sum << endl;
	}

	cout << 1.0 << "   " << 1.0 << endl;

	return 0;
}
