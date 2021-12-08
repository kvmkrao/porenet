/*
Author: V Kotteda
Date: Dec 2021
*/

#ifndef CONSTANTS_H
#define CONSTANTS_H

// define your own namespace to hold constants

#include <string>
#include <iostream>
using namespace std;
namespace constants
{

     inline constexpr double PI  {3.141592653589793238463};
     inline constexpr double nonle = 1e-6;
     inline constexpr double m2nm  = 1.0/nonle;
     inline constexpr double pnond = 1e6;
     inline constexpr const int maxrand = 50000;
/*
    inline constexpr char fnode1[] =  "node.dat";
    inline constexpr char fnode2[] =  "node_dia.dat";
    inline constexpr char flink1[] =  "link.dat"; 
    inline constexpr char flink2[] =  "link_dia.dat";

    inline constexpr char fnode1[]  ="veldsteen_poreus_519x557x861_node1.dat"; 
    inline constexpr char fnode2[]  ="veldsteen_poreus_519x557x861_node2.dat"; 
    inline constexpr char flink1[]  ="veldsteen_poreus_519x557x861_link1.dat";  
    inline constexpr char flink2[]  ="veldsteen_poreus_519x557x861_link2.dat"; 
*/ 
    inline constexpr char fnode1[]  ="cubic_node1.dat"; 
    inline constexpr char fnode2[]  ="cubic_node2.dat"; 
    inline constexpr char flink1[]  ="cubic_link1.dat";  
    inline constexpr char flink2[]  ="cubic_link2.dat"; 
 
    // ... other related constants
}
#endif

