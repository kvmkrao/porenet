#ifdef WIN32        
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#endif

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <algorithm>
using namespace std;

#include "node.h"
#include "rockElem.h"

/////////////////////////////////////////////////////////////////////////////////
// Most of the poreinitialisation is delayed until connecting throats are 
// created. This is because their radius is dependant on all connecting throats.
/////////////////////////////////////////////////////////////////////////////////
Pore::Pore(Node *node, double shapeFactor, int connectionNum) : RockElem(shapeFactor), m_node(node)  
{
    m_throats.resize(connectionNum);
    m_connectionNumber = connectionNum;
}

/////////////////////////////////////////////////////////////////////////////////
// Setting the index (required when writing to file) is delayed until connection
// number has been reduced. 
/////////////////////////////////////////////////////////////////////////////////
Throat::Throat(Pore *poreOne, Pore *poreTwo, double shapeFactor, double radius, double length, 
               double clayProp, bool pbcThroat) : RockElem(shapeFactor), m_pbcThroat(pbcThroat)
{
    m_radius = radius;
    m_length = length;
    m_pores.push_back(poreOne);
    m_pores.push_back(poreTwo);
    m_index = 0;

    double area = pow(m_radius, 2.0) / (4.0 * m_shapeFactor);
    m_volume = area * m_length;
    
    if(poreOne->node()->isInOrOutlet() || poreTwo->node()->isInOrOutlet())
        m_volume = 0.0;

    m_clayVolume = (clayProp * m_volume) / (1.0 - clayProp);  
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The pore data is written to file in following format:
//
// *_node1.dat (outOne):
// index, x_pos, y_pos, z_pos, connection num, connecting nodes..., at inlet?, at outlet?, connecting links...
//
// *_node2.dat (outTwo):
// index, volume, radius, shape factor, clay volume
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Pore::writeData(ostream& outOne, ostream& outTwo) const
{
    size_t i;
    outOne.flags(ios::showpoint);
    outOne.flags(ios::scientific);
    outTwo.flags(ios::showpoint);
    outTwo.flags(ios::scientific);

    outOne << *m_node << setw(5) << m_connectionNumber;

    for(i = 0; i < m_throats.size(); ++i)
        outOne << setw(7) << m_throats[i]->nextPore(this)->node()->indexOren();     // Connecting nodes

    outOne << setw(7) << m_node->isAtInlet() << setw(7) << m_node->isAtOutlet();    // In and outlet?

    for(i = 0; i < m_throats.size(); ++i)
        outOne << setw(7) << m_throats[i]->index();                                 // Connecting throats

    outOne << endl;
    
    outTwo << setw(7) << m_node->index() 
        << setw(15) << m_volume 
        << setw(15) << m_radius
        << setw(15) << m_shapeFactor
        << setw(15) << m_clayVolume
        << endl;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
// The throat data is written to file in following format:
//
// *_link1.dat (outOne):
// index, pore 1 index, pore 2 index, radius, shape factor, total length (pore center to pore center)
//
// *_link2.dat (outTwo):
// index, pore 1 index, pore 2 index, length pore 1, length pore 2, length throat, volume, clay volume
/////////////////////////////////////////////////////////////////////////////////////////////////////
void Throat::writeData(ostream& outOne, ostream& outTwo) const
{
    outOne.flags(ios::showpoint);
    outOne.flags(ios::scientific);
    outTwo.flags(ios::showpoint);
    outTwo.flags(ios::scientific);

    double lenPore1(m_pores[0]->halfLength()), lenPore2(m_pores[1]->halfLength());
    double lenTot(lenPore1+lenPore2+m_length);

    outOne << setw(7)   << m_index 
        << setw(7)      << m_pores[0]->node()->indexOren()
        << setw(7)      << m_pores[1]->node()->indexOren()
        << setw(15)     << m_radius
        << setw(15)     << m_shapeFactor
        << setw(15)     << lenTot
        << endl;

    outTwo << setw(7)   << m_index
        << setw(7)      << m_pores[0]->node()->indexOren()
        << setw(7)      << m_pores[1]->node()->indexOren()
        << setw(15)     << lenPore1
        << setw(15)     << lenPore2
        << setw(15)     << m_length
        << setw(15)     << m_volume
        << setw(15)     << m_clayVolume
        << endl;
}

////////////////////////////////////////////////////////////////////////////////
// The pores are finally initialized when the connecting throats are known. The 
// radius of the pore is calculated using the same procedure that Daryl used, 
// ie. the average throat radius multiplied by an aspect ratio (as long as this
// value is larger than the max throat radius). Pore length is taken to be twice
// the radius.
////////////////////////////////////////////////////////////////////////////////
void Pore::finaliseInit(double aspectRatio, double clayProp)
{
    double throatRadiusSum(0.0), maxThroatRadius(0.0);
    int numConn(0);
    
    for(size_t i = 0; i < m_throats.size(); ++i)
    {
        if(m_throats[i] != NULL)                                                
        {                                                                       
            throatRadiusSum += m_throats[i]->radius();   
            maxThroatRadius = max(maxThroatRadius, m_throats[i]->radius());
            ++numConn; 
        }
    }
    
    m_radius = max(aspectRatio * (throatRadiusSum / numConn), maxThroatRadius);
    m_length = 2.0 * m_radius;

    double area = pow(m_radius, 2.0) / (4.0 * m_shapeFactor);
    m_volume = area * m_length;
    m_clayVolume = (clayProp * m_volume) / (1.0 - clayProp);  
}

////////////////////////////////////////////////////////////////////////////////
// When reducing the connection number it becomes necissary to remove the 
// references to these throats in the pores.
////////////////////////////////////////////////////////////////////////////////
void Pore::removeThroat(Throat *throat)
{
    vector< Throat* >::iterator iter = 
        find(m_throats.begin(), m_throats.end(), throat);

    if(iter != m_throats.end()) 
        m_throats.erase(iter);
    else
    {
        cout << "whoa..." << endl;
        exit(-1);
    }

    m_connectionNumber = static_cast< int >(m_throats.size());
}

////////////////////////////////////////////////////////////////////////////////////
// When calculating average connection number we do not want to include connections
// to in/outlet.
////////////////////////////////////////////////////////////////////////////////////
int Pore::connectionNum() const
{
    if(m_node->isInOrOutlet())
        return 0;
    else
        return m_connectionNumber;
}

bool Pore::containsThroat(const Throat* throat) const
{
    for(size_t i = 0; i < m_throats.size(); ++i)
    {
        Throat* throat = m_throats[i];
        if(m_throats[i] == throat) return true;
    }
    
    return false;
}

const Pore* Throat::nextPore(const Pore *callingPore) const
{
    if(m_pores[0] == callingPore)
        return m_pores[1];
    else
        return m_pores[0];
}

////////////////////////////////////////////////////////////////////////////////////
// Rigerous error checking is done to verify that the state of the network is correct
////////////////////////////////////////////////////////////////////////////////////
void Pore::checkState() const
{
    if(!m_node->isInOrOutlet() && (m_throats.size() > 6 || m_throats.size() < 1))
    {
        cerr <<"================================================"   << endl
            << "Error: Throats are not initialized in pore: " 
            << m_node->index()                                      << endl
            << "================================================"   << endl;
        exit(-1);        
    }
    
    if(m_connectionNumber != static_cast< int>(m_throats.size()))
    {
        cerr <<"================================================"   << endl
            << "Error: Throats are not initialized in pore: " 
            << m_node->index()                                      << endl
            << "================================================"   << endl;
        exit(-1);
    }
    
    for(size_t i = 0; i < m_throats.size(); ++i)
    {
        if(m_throats[i] == NULL)
        {
            cerr <<"================================================"   << endl
                << "Error: Throats are not initialized in pore: " 
                << m_node->index()                                      << endl
                << "================================================"   << endl;
            exit(-1);
        }

        if(!m_node->isInOrOutlet() && m_throats[i]->connectionNum() != 2)
        {
            cerr <<"================================================"   << endl
                << "Error: Throats are not initialized in pore: " 
                << m_node->index()                                      << endl
                << "================================================"   << endl;
            exit(-1);
        }

        if(m_throats[i]->radius() > m_radius)
        {
            cerr <<"======================================================" << endl
                << "Error: Throat radius gerater than pore radius in: " 
                << m_node->index()                                          << endl
                << "======================================================" << endl;
            exit(-1);
        }
    }
}

void Throat::checkState() const
{
    for(size_t i = 0; i < m_pores.size(); ++i)
    {
        if(m_pores[i] == NULL)
        {
            cerr <<"================================================"   << endl
                << "Error: Pores are not initialized in throat: " 
                << m_index                                              << endl
                << "================================================"   << endl;
            exit(-1);
        }

        if(!m_pores[i]->containsThroat(this))
        {
            cerr <<"================================================"   << endl
                << "Error: Pores are not initialized in throat: " 
                << m_index                                              << endl
                << "================================================"   << endl;
            exit(-1);
        }      
    }

    if(m_pores.size() != 2)
    {
        cerr <<"================================================"   << endl
            << "Error: Porees are not initialized in throat: " 
            << m_index                                              << endl
            << "================================================"   << endl;
        exit(-1);
    }
}

