
/*
 * Author: V Kotteda
 * Date  : November 19, 2021
 */

#include "link.hpp"
#include "constants.h"

using namespace std;

Link::Link() {}
Link::~Link() {}
void Link::readdata()
{
    ifstream linkinfo(constants::flink1); //("link.dat");
    linkinfo >> nn;
    ifstream linkdia(constants::flink2);  //("link_dia.dat");
    linkdia >> nl;
    for(int i=0; i<nn; i++){
    	linkinfo >> l1 >> l2;
        linkdia >> radius;
        ends.push_back(std::make_pair(l1, l2));
        rad.push_back(radius);
    }
    linkinfo.close();
    linkdia.close();
}

void Link::DisplayContents()
{
    for( unsigned int i = 0; i < ends.size(); i++ )
    {
    	std::cout << "Element[" << i << "] = " << ends.at(i).first <<"\t" << ends.at(i).second << std::endl;
    }
    std::cout << std::endl;
}

OrenLink::OrenLink() {}
OrenLink::~OrenLink() {}

void OrenLink::readdata()
{
	ifstream linkn1;
	linkn1.open(constants::flink1); 
        ifstream linkn2; 
        linkn2.open(constants::flink2); 
	linkn1 >> nl;
        for(int i=0; i<nl ; i++) {     //( id = 0(inside) 1(inlet) 2(outlet) 
        	linkn1 >>  Id >> l1 >> l2  >>  ld  >> shfnl >>  ll ;
        	linkn2 >>  Id >> l1 >> l2  >>  n1l >>  n2l >> ll >> volth >> clay;
        	lid = 0 ;
        	if( l1*l2 < 0.0)  lid = 1;
        	if( l1*l2 == 0.0) lid = 2;
        	ends.push_back(std::make_pair(l1,l2));
        	rad.push_back(ld*1e6);
        	id.push_back(lid);
        	len.push_back(ll*1e6);
	}

	linkn1.close();
	linkn2.close();
}

void OrenLink::DisplayContents(){Basecl.DisplayContents();}

