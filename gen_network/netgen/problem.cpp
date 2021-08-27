#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <vector>
#include <cassert>
#include <list>
#include <cmath>
using namespace std ;

#include "node.h"
#include "rockElem.h"
#include "problem.h"

const int       Problem::MAX_CONN_NUM   = 6;
const double    Problem::PI             = acos(-1.0);

Problem::Problem(const string fileName)
{  
    ifstream in(fileName.c_str());
    
    if (!in)
    {
        cerr << "==============================================" << '\n' 
             << "Error: Unable to open input file " << fileName  << '\n'
             << "==============================================" << '\n';
        exit( -1 );
    }

    in >> *this;
    in.close();

    initNetworkModel();
}

istream& operator>> (istream& in, Problem& prob) 
{
    char dummy[256];

    in >> prob.m_outFileNameBase;  
    in.getline(dummy,256,'\n');

    in >> prob.m_nX >> prob.m_nY >> prob.m_nZ;  
    in.getline(dummy,256,'\n');
    
    in >> prob.m_throatRadWeibull[Min] 
        >> prob.m_throatRadWeibull[Max] 
        >> prob.m_throatRadWeibull[Delta] 
        >> prob.m_throatRadWeibull[Eta];  
    in.getline(dummy,256,'\n');

    prob.m_throatRadWeibull[0] *= 1.0E-6;       // Get it into m from micro.m
    prob.m_throatRadWeibull[1] *= 1.0E-6;

    in >> prob.m_throatLenWeibull[Min] 
        >> prob.m_throatLenWeibull[Max] 
        >> prob.m_throatLenWeibull[Delta] 
        >> prob.m_throatLenWeibull[Eta];  
    in.getline(dummy,256,'\n');

    prob.m_throatLenWeibull[0] *= 1.0E-6;       // Get it into m from micro.m
    prob.m_throatLenWeibull[1] *= 1.0E-6;

    in >> prob.m_aspectRatioWeibull[Min] 
        >> prob.m_aspectRatioWeibull[Max] 
        >> prob.m_aspectRatioWeibull[Delta] 
        >> prob.m_aspectRatioWeibull[Eta];  
    in.getline(dummy,256,'\n');

    in >> prob.m_triangleGWeibull[Min] 
        >> prob.m_triangleGWeibull[Max] 
        >> prob.m_triangleGWeibull[Delta] 
        >> prob.m_triangleGWeibull[Eta];  
    in.getline(dummy,256,'\n');
    
    if(prob.m_triangleGWeibull[Min] < 0.0 || prob.m_triangleGWeibull[Max] > sqrt(3.0)/36.0) 
    {
        cerr <<"======================================="    << '\n'
            << "Error: Bound on triangular shape factor"    << '\n'
            << "is incorrect"                               << '\n'
            << "======================================="    << '\n';
        exit(-1);
    }
       
    in >> prob.m_throatShapeProp[Square] >> prob.m_throatShapeProp[Circ];
    prob.m_throatShapeProp[Triang] = 1.0 - prob.m_throatShapeProp[Square] - prob.m_throatShapeProp[Circ];
    in.getline(dummy,256,'\n');

    in >> prob.m_poreShapeProp[Square] >> prob.m_poreShapeProp[Circ];
    prob.m_poreShapeProp[Triang] = 1.0 - prob.m_poreShapeProp[Square] - prob.m_poreShapeProp[Circ];
    in.getline(dummy,256,'\n');

    for(int i = 0; i < 3; ++i)
    {
        if(prob.m_throatShapeProp[i] < 0.0 || prob.m_throatShapeProp[i] > 1.0 || 
            prob.m_poreShapeProp[i] < 0.0 || prob.m_poreShapeProp[i] > 1.0)
        {
            cerr << "=======================================" << '\n'
                << "Error: Proportion of square or circular"  << '\n'
                << "elements are incorrect"                   << '\n'
                << "======================================="  << '\n';
            exit(-1);
        }
    }

    in >> prob.m_clayProportion;
    in.getline(dummy,256,'\n');

    if(prob.m_clayProportion > 1.0 || prob.m_clayProportion < 0.0) 
    {
        cerr <<"==============================================" << '\n'
            << "Error: Clay proportion must be between 0 and 1" << '\n'
            << "==============================================" << '\n';
        exit(-1);
    }

    in >> prob.m_averConnectionNum;
    in.getline(dummy,256,'\n');

    if(prob.m_averConnectionNum > 6.0 || prob.m_averConnectionNum <= 0.0) 
    {
        cerr <<"======================================="    << '\n'
            << "Error: Average connection number must"      << '\n'
            << "be between 0 and 6"                         << '\n'
            << "======================================="    << '\n';
        exit(-1);
    }

    char periodicBC;
    in >> periodicBC;
    in.getline(dummy,256,'\n');

    prob.m_periodicBC = (periodicBC == 't' || periodicBC == 'T');


    return in;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
// The network data is written to file in the same format as used by Paal-Eric Oeren. The data is contained
// in four files: *_node1.dat, *_node2.dat, *_link1.dat and *_link2.dat. 
///////////////////////////////////////////////////////////////////////////////////////////////////////////
void Problem::writeData() const
{
    string pOut1FileName(m_outFileNameBase + "_node1.dat"), pOut2FileName(m_outFileNameBase + "_node2.dat");
    
    ofstream pOut1, pOut2;
    pOut1.open(pOut1FileName.c_str());
    pOut2.open(pOut2FileName.c_str());
    pOut1.flags(ios::showpoint);
    pOut1.flags(ios::scientific);

    pOut1 << m_numPores << "   " << m_xDim << "   " << m_yDim << "   " << m_zDim << '\n';

    for(int i = 1; i <= m_numPores; ++i)
        m_pores[i]->writeData(pOut1, pOut2);

    pOut1.close();
    pOut2.close();

    string tOut1FileName(m_outFileNameBase + "_link1.dat"), tOut2FileName(m_outFileNameBase + "_link2.dat");
    ofstream tOut1, tOut2;
    tOut1.open(tOut1FileName.c_str());
    tOut2.open(tOut2FileName.c_str());

    tOut1 << static_cast< int >(m_throats.size()) << "   " << '\n';
    
    list< Throat * >::const_iterator iter;
    for(iter = m_throats.begin(); iter != m_throats.end(); ++iter)
    {
        Throat *throat = *iter;
        throat->writeData(tOut1, tOut2);
    }

    tOut1.close();
    tOut2.close();
}

/////////////////////////////////////////////////////////////////////////////////////////////
// First pores are created. These are then connected by throats in a cubic lattice. The 
// initialisation of the pores is then completed. 
//
// Then the size of the network model is calculated, the connection number reduced and the
// model is checked for possible errors.
/////////////////////////////////////////////////////////////////////////////////////////////
void Problem::initNetworkModel()
{
    m_numPores = m_nX * m_nY * m_nZ;

    createPores();
    connectPoresWithThroats();

    for(size_t elem = 0; elem < m_pores.size(); ++elem)
    {
        assert(m_pores[elem] != NULL);
        double aspectRatio = weibull(m_aspectRatioWeibull[Min], m_aspectRatioWeibull[Max], 
            m_aspectRatioWeibull[Delta], m_aspectRatioWeibull[Eta]);
        m_pores[elem]->finaliseInit(aspectRatio, m_clayProportion);
    }

    findNetworkModelSize();
    reduceConnectionNumber();
    assignThroatIndex();
    setPoreLocation();
    checkNetworkIntegrity();

    cout << '\n'
        << "== Finished generating cubic network =="                            << '\n'
        << "Number of pores:                   " << m_numPores                  << '\n'
        << "Average connection number:         " << m_actualConnectionNumber    << '\n'
        << "Porosity:                          " << m_porosity                  << '\n';
}
////////////////////////////////////////////////////
// The state of both pores and throats are checked.
////////////////////////////////////////////////////
void Problem::checkNetworkIntegrity() const
{
    list< Throat * >::const_iterator iter;
    for(iter = m_throats.begin(); iter != m_throats.end(); ++iter)
    {
        Throat *throat = *iter;
        throat->checkState();
    }

    for(size_t i = 0; i < m_pores.size(); ++i)
        m_pores[i]->checkState();
}

//////////////////////////////////////////////
// The pore location is set on a regular grid.
//////////////////////////////////////////////
void Problem::setPoreLocation()
{
    int numThroats(0);
    for(size_t i = 0; i < m_pores.size(); ++i)
    {
        m_pores[i]->node()->setLocation(m_xDim, m_yDim, m_zDim, m_averageThroatLength);
        numThroats += m_pores[i]->connectionNum();
    }

    m_actualConnectionNumber = ((double) numThroats / m_numPores);
}

//////////////////////////////////////////////////////////////////////////////
// When writing the pore/throat data to file we need to reference the throats
// by an index. The index is set after the connection number has been reduced
// so that they become consecutive.
//////////////////////////////////////////////////////////////////////////////
void Problem::assignThroatIndex()
{
    int index(0);
    double throatVol(0.0);

    list< Throat * >::iterator iter;
    for(iter = m_throats.begin(); iter != m_throats.end(); ++iter)
    {
        Throat *throat = *iter;
        throat->setIndex(++index);
        throatVol += throat->totVolume();
    }

    double poreVolume(0.0);
    for(int i = 1; i <= m_numPores; ++i)
        poreVolume += m_pores[i]->totVolume();

    m_porosity = (throatVol + poreVolume) / (m_xDim * m_yDim * m_zDim);
}

///////////////////////////////////////////////////////////////////////////////////
// A random number between 0 and 1 is drawn. If the number is less than the
// probability for deletion, the throat is deleted. Hence we will only achive 
// the requested connection number if the number of pores is very large.
//////////////////////////////////////////////////////////////////////////////////
void Problem::reduceConnectionNumber()
{
    int numThroats(static_cast< int >(m_throats.size()));
    int numDeletions = static_cast< int >(numThroats * ((MAX_CONN_NUM-m_averConnectionNum)/MAX_CONN_NUM));

    double deletionProb((double)numDeletions/numThroats);

    list< Throat * >::iterator iter;
    for(iter = m_throats.begin(); iter != m_throats.end();)
    {
        Throat *throat = *iter;
        Pore* poreOne = throat->connectingPore(0);
        Pore* poreTwo = throat->connectingPore(1);
                
        if(!m_periodicBC && throat->pbcThroat()) 
        {    
            poreOne->removeThroat(throat);
            poreTwo->removeThroat(throat);            
            iter = m_throats.erase(iter);
            delete throat;
        }
        else
            ++iter;
    }

    for(iter = m_throats.begin(); iter != m_throats.end();)
    {
        double randNum = (double) rand() / RAND_MAX;
        Throat *throat = *iter;
        Pore* poreOne = throat->connectingPore(0);
        Pore* poreTwo = throat->connectingPore(1);
                
        if(poreOne->connectionNum() > 1 && poreTwo->connectionNum() > 1 && randNum < deletionProb) 
        {    
            poreOne->removeThroat(throat);
            poreTwo->removeThroat(throat);            
            iter = m_throats.erase(iter);
            delete throat;
        }
        else
            ++iter;
    }
}

//////////////////////////////////////////////////////////////////////////////
// The size of the model is taken to be the average distance betwen pores at 
// the boundaries.
//////////////////////////////////////////////////////////////////////////////
void Problem::findNetworkModelSize()
{
    int numXFace(m_nY*m_nZ), numYFace(m_nX*m_nZ), numZFace(m_nY*m_nX);
    size_t i;
    double lengthSum;

    lengthSum = 0.0;
    for(i = 1; i < m_pores.size() - 1; ++i)
    {
        lengthSum += m_pores[i]->length();
        lengthSum += m_pores[i]->connectingThroat(0)->length();
    }
    m_xDim = lengthSum / static_cast<double>(numXFace);

    lengthSum = 0.0;
    for(i = 1; i < m_pores.size() - 1; ++i)
    {
        lengthSum += m_pores[i]->length();
        lengthSum += m_pores[i]->connectingThroat(1)->length();
    }
    m_yDim = lengthSum / static_cast<double>(numYFace);

    lengthSum = 0.0;
    for(i = 1; i < m_pores.size() - 1; ++i)
    {
        lengthSum += m_pores[i]->length();
        lengthSum += m_pores[i]->connectingThroat(2)->length();
    }
    m_zDim = lengthSum / static_cast<double>(numZFace);
}


void Problem::connectPoresWithThroats()
{
    double throatLenSum(0.0);
    for(int poreIdx = 1; poreIdx <= m_numPores; ++poreIdx)
    {
        for(int conn = 0; conn < MAX_CONN_NUM; ++conn)
        {        
            const Node* currNode = m_pores[poreIdx]->node();
            bool pbcThroat(false);
            int nextPoreIdx = currNode->nextIndex(conn, pbcThroat);

            if(nextPoreIdx >= 0 && m_pores[poreIdx]->connectingThroat(conn) == NULL)
            {
                assert(m_pores[nextPoreIdx] != NULL && nextPoreIdx < m_numPores + 2);
                
                Pore *nextPore = m_pores[nextPoreIdx];
                
                double radius = weibull(m_throatRadWeibull[Min], m_throatRadWeibull[Max], 
                    m_throatRadWeibull[Delta], m_throatRadWeibull[Eta]);
                double length = weibull(m_throatLenWeibull[Min], m_throatLenWeibull[Max], 
                    m_throatLenWeibull[Delta], m_throatLenWeibull[Eta]);
                double shapeFactor = evalShapeFactor(m_triangleGWeibull, m_throatShapeProp);
                throatLenSum += length;
                                
                Throat *throat = new Throat(nextPore, m_pores[poreIdx], shapeFactor, radius, length, 
                    m_clayProportion, pbcThroat);                
                m_throats.push_back(throat);
                
                m_pores[poreIdx]->addThroat(conn, throat);
                
                int nextConn;                
                if(nextPore->node()->isInOrOutlet())
                    nextConn = (currNode->k()-1) * m_nY + currNode->j() - 1;
                else
                    nextConn = (conn + (MAX_CONN_NUM/2)) % MAX_CONN_NUM;
                
                nextPore->addThroat(nextConn, throat);                
            }
        }
    }
    m_averageThroatLength = throatLenSum/static_cast<double>(m_throats.size());
}

////////////////////////////////////////////////////////////////////////////////////////////////
// wheter a pore/throat is square, circle or triangular is again determined by drawing a 
// random number between 0 and 1
////////////////////////////////////////////////////////////////////////////////////////////////
double Problem::evalShapeFactor(const double propG[], const double proportion[]) const
{
    double randNum = (double) rand() / RAND_MAX;            //Creates a pseudo random number in range [0, 1]
    double shapeFactor;
    
    if(randNum < proportion[Square])
        shapeFactor = 1.0/16.0;                                                     // Square
    else if(randNum < proportion[Square] + proportion[Circ])
        shapeFactor = 1.0/(4.0*PI);                                                 // Circular    
    else
        shapeFactor = weibull(propG[Min], propG[Max], propG[Delta], propG[Eta]);    // Triangular    

    return shapeFactor;
}

//////////////////////////////////////////////////////////////////////////////////////////
// In and outlet pores have been given indicies 0 and numPores + 1. This is because the 
// pores are stored in a vector. Only when writing the data to file their indicies are
// changed to -1 and 0.
//////////////////////////////////////////////////////////////////////////////////////////
void Problem::createPores()
{
    m_pores.resize(m_numPores + 2);

    for(int k = 1; k <= m_nZ; ++k)
    {
        for(int j = 1; j <= m_nY; ++j)
        {
            for(int i = 1; i <= m_nX; ++i)
            {                
                double shapeFactor = evalShapeFactor(m_triangleGWeibull, m_poreShapeProp);
                
                Node *currNode = new Node(i, j, k, m_nX, m_nY, m_nZ);
                Pore *pore = new Pore(currNode, shapeFactor, MAX_CONN_NUM);

                m_pores[currNode->index()] = pore;
            }
        }
    }

    Node *inletNode = new Node(0, 1, 1, m_nX, m_nY, m_nZ);
    Pore *inlet = new Pore(inletNode, 0.04811, m_nY*m_nZ);
    m_pores[inletNode->index()] = inlet;

    Node *outletNode = new Node(m_nX+1, 1, 1, m_nX, m_nY, m_nZ);
    Pore *outlet = new Pore(outletNode, 0.04811, m_nY*m_nZ);
    m_pores[outletNode->index()] = outlet;
}

///////////////////////////////////////////////////////////////////////////////////////////
// Returns a value within the truncated weibull distribution as used by Robin and Fenwick
///////////////////////////////////////////////////////////////////////////////////////////
double Problem::weibull(double min, double max, double delta, double eta) const
{    
    double randNum = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);            //Creates a pseudo random number in range [0, 1]                   
    return (max - min) * pow(-delta*log(randNum*(1.0-exp(-1.0/delta))+exp(-1.0/delta)), 1.0/eta) + min;
}
