
#include<cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <sys/stat.h>
#include <sstream>
#include <vector>
#include <iomanip>      // std::setiosflags
#include <stdlib.h>     /* srand, rand */

using namespace std;

class Node
{
        public:
                Node(); 
                ~Node(); 
                void readdata(); 
                void DisplayVectorContents(); 

        //private:
                std::vector<double> rad;
                std::vector<double> xc;
                std::vector<double> yc;
                std::vector<double> zc;

         private:
                double xx, yy, zz, radius;
                int nn; 
};


class Link
{
        public:
                Link(); 
                ~Link();
                void readdata(); 

                void DisplayVectorContents(); 

        //private:
                std::vector<pair<int,int>> ends;
                std::vector<double> rad;
         private:
                int l1, l2 ;
                int nn, nl; 
                double radius; 
};


