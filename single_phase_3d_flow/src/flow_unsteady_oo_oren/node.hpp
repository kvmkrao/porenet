/*
Author: V Kotteda
Date:   Dec 8, 2021 
*/
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

class Node
{
	public:
		Node();
		~Node();
		void readdata();
		void DisplayContents();

		std::vector<double> rad;
		std::vector<double> xc;
		std::vector<double> yc;
		std::vector<double> zc;

	private:
		double xx, yy, zz, radius;
		int nn;
};

class OrenNode
{
	public:
		OrenNode();
		~OrenNode();
		void readdata();
                
		void DisplayContents();

		std::vector<double> rad;
		std::vector<double> xc;
		std::vector<double> yc;
		std::vector<double> zc;
		std::vector<int> conum;

	private:
		double xx, yy, zz, radius;
		int nn, id, ni;
		double nr, vol, shapeFactor, clay;
		double dd1, dd2, dd3,dd4;
		string line, str1, str2, str3, str4, str5, str6, str7,str8;
		Node Basecn;
};


