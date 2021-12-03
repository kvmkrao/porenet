 
#include<cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <sys/stat.h>
#include <sstream>
#include <vector>
#include <iomanip>
#include <stdlib.h>

using namespace std;

class Link
{
	public:
		Link();
		~Link();
		void readdata();

		void DisplayContents();
		std::vector<pair<int,int>> ends;
		std::vector<double> rad;
	private:
		int l1, l2 ;
		int nn, nl;
		double radius;
};

class OrenLink
{
	public:
		OrenLink();
		~OrenLink();
		void readdata();
		void DisplayContents();

		std::vector<pair<int,int>> ends;
		std::vector<double> rad;
		std::vector<int> id;
		std::vector<double> len;
	private:
		int l1, l2, Id ;
		int nn, nl, lid;
		double radius,ld,shfnl,ll;
		double n1l, n2l, volth, clay;
		Link Basecl;
};

