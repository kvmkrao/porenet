
#include <iostream>
#include <cstdlib>
#include <vector>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <string>

#include "network01.h"

using namespace std;
int main(int argc, char *argv[])
{
	cout << endl << "Equivalent Random Network Generator for Berea Network" << endl
		 << "Version: 001 - Work in Progress" << endl;
	srand ((unsigned) time (NULL));
	ifstream in;
	in.open("default.dat");

	if(!in)
	{
		cout << "Error opening input data file." << endl;
		system("pause");
		exit(1);
	}
	
	string net, dim, x, y, z, Pore, Throat;
	double i, j, k;
	int n, tn; 

	in >> net >> dim >> x >> y >> z >> i >> j >> k >> Pore >> Throat >>n >>tn;
		
    Network mynetwork;
	mynetwork.netwsize(i, j, k, n, tn);
	system ("pause");
	return 0;
}
