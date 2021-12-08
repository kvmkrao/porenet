
/*
 * Author: V Kotteda
 * Date  : November 19, 2021
 */


#include "node.hpp"
#include "constants.h"

Node::Node()  {}
Node::~Node() {}

void Node::readdata()
{
    ifstream  nodeloc(constants::fnode1);
    nodeloc >> nn;
    ifstream noderad(constants::fnode2);

    for(int i=0; i<nn; i++){
    	nodeloc >> xx >> yy >> zz;
        noderad >> radius;
        rad.push_back(radius);
        xc.push_back(xx);
        yc.push_back(yy);
        zc.push_back(zz);
    }
    nodeloc.close();
    noderad.close();
}
    
void Node::DisplayContents() 
{
    for( unsigned int i = 0; i < rad.size(); i++){
	std::cout << "Node[" << i << "] = " << rad[i] << std::endl;

    }
    std::cout << std::endl;
}

OrenNode::OrenNode()  {}
OrenNode::~OrenNode() {}

void OrenNode::readdata()
{
	ifstream nodeloc1, nodeloc2;
	nodeloc1.open(constants::fnode1);
	nodeloc2.open(constants::fnode2);
    	nodeloc1 >> nn >>  str1 >>  str2 >> str3 ;
	for(int i=0; i<nn ; i++) {  
		nodeloc2 >> id >> str5 >> str6 >> str7 >> str8; 
		radius   = std::stod(str6);
                //cout << "radius" << radius << endl ; 
		rad.push_back(radius*1.0e6);
	}

    	while ( getline(nodeloc1, line) )
    	{
    	if (line.size() > 54)
        {
    		str1  = line.substr(0, 7);   // id
           	str2  = line.substr(8, 20);  // xx
           	str3  = line.substr(22,34);  // yy
           	str4  = line.substr(36, 48);  // zz
           	str5  = line.substr(50, 52);  // ni
           	xx = std::stod(str2);
           	yy = std::stod(str3);
           	zz = std::stod(str4);
           	ni = std::stoi(str5);
           	xc.push_back(xx*1.0e6);
           	yc.push_back(yy*1.0e6);
           	zc.push_back(zz*1.0e6);
           	conum.push_back(ni);
       	}
   	}
    	nodeloc1.close();
    	nodeloc2.close();
}

void OrenNode::DisplayContents() {
	Basecn.DisplayContents();
}

