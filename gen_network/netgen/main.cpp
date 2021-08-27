#include <fstream>
#include <iostream>
#include <cstdio>
#include <cmath>
#include <string>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <queue>
#include <list>
#include <algorithm>
using namespace std ;

#include "node.h"
#include "rockElem.h"
#include "problem.h"

void Wait4me();
istream& EatWhiteSpace(istream&);
  
int main(int argc, char *argv[])        //Allow file name as argument in program launch
{
    //atexit(Wait4me);                  //Useful little function while debugging in MSVS

    srand((unsigned)time( NULL ));      // Seed to random number generator

    string inputFile;
    
    cout<<"Network Generation Code, version A01" << endl;
    
    if (argc > 1)
        inputFile = argv[1];
    else
    {
        cout << "Please input data file : "; 
        cin >> inputFile;
    }
    
    Problem network(inputFile);         // All initialixation takes place in the constructor
    network.writeData();                // The network data is written to file  
   
    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////
// This thing is useful during debugging in MS Visual studio since the terminal window wont
// dissapear immediatly after the program finishes.
////////////////////////////////////////////////////////////////////////////////////////////
void Wait4me()
{
    cout << "Press any key to exit" << endl;
    getchar();
}

///////////////////////////////////////////////////////////////////////////////////////
// This little function removes white space from the input file. Made it a global such 
// that any class can use it
///////////////////////////////////////////////////////////////////////////////////////
istream& EatWhiteSpace(istream& is)
{
    char c;
    while(is.get(c))
    {
        if(!isspace(c))
        {
            is.putback(c);
            break;
        }
    }
    return is;
}
