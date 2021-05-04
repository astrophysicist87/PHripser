#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

void compute_PD( const vector<double & birth_radii, const vector<double & death_radii,
				 const vector<double & pVec, vector<double & PD )
{
	PD.clear();
	for ( double p : pVec )
	{
		double PDresult = 0.0;
		for ( int iPair = 0; iPair < birth_radii.size(); iPair++ )
			PDresult += pow( abs( birth_radii[iPair] - death_radii[iPair] ), p );
		PD.push_back( pow( PDresult, 1.0 / p ) );
	}
	return;
}

double get_SD( const vector<double & birth_radii, const vector<double & death_radii )
{
	double SD = 0.0;
	for ( int iPair = 0; iPair < birth_radii.size(); iPair++ )
		SD += abs( death_radii[iPair] - birth_radii[iPair] );
	return ( SD );
}

double get_ED( const vector<double & birth_radii, const vector<double & death_radii )
{
	double SD = get_SD( birth_radii, death_radii );
	double ED = 0.0;
	for ( int iPair = 0; iPair < birth_radii.size(); iPair++ )
	{
		double argument = abs( death_radii[iPair] - birth_radii[iPair] ) / SD;
		ED += argument * log( argument );
	}

	return ( -ED / log(SD) );	
}

int main(int argc, char** argv)
{
	vector<double> pVec;
	for ( double p = 0.1; p <= 10.0+1e-10; p += 0.1 ) pVec.push_back( p );

	vector<double> PD, birth_radii, death_radii;

	// read path to input file from command line
	string path_to_file = string(argv[1]);

	// then read in file itself
	ifstream infile(path_to_file.c_str());
	if (infile.is_open())
	{
		string line;
		while ( getline (infile, line) )
		{
			double birth = 0.0, death = 0.0;
			istringstream iss(line);
			iss >> birth >> death;
			birth_radii.push_back( birth );
			death_radii.push_back( death );
		}
	}

	// complexity measure #1
	compute_PD( birth_radii, death_radii, pVec, PD );

	// complexity measure #2
	double ED = get_ED( birth_radii, death_radii );

	cout << "ED: " << ED << endl;
	for (int ip = 0; ip < pVec.size(); ip++)
		cout << "PD: " << pVec[ip] << "   " << PD[ip] << endl;

	return 0;
}

