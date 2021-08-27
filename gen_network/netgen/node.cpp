#include <iostream>
#include <iomanip>
using namespace std;

#include "node.h"


ostream& operator<< (ostream& out, Node& node) 
{  
    out.flags(ios::showpoint);
    out.flags(ios::scientific);

    out << setprecision(6)
        << setw(7)  << node.m_index
        << setw(15) << node.m_xLoc
        << setw(15) << node.m_yLoc
        << setw(15) << node.m_zLoc;

    return out;
}

/////////////////////////////////////////////////////////////////////////////
// Various constructors for the node class. The node class does make an 
// assumption about the cubic structure of the network when determening
// if a node is within the network or not.
/////////////////////////////////////////////////////////////////////////////
Node::Node(int i, int j, int k, int nX, int nY, int nZ)
    : m_i(i), m_j(j), m_k(k), m_nX(nX), m_nY(nY), m_nZ(nZ)
{
    initNode();
}

///////////////////////////////////////////////////////////////////////////
// Return the next node index given a node and a direction. If the next
// node is outside the lattice an index of -1 is returned. If a pores is
// on a boundary (not in/outlet) the index on the other side is returned
// ie periodic boundary conditions.
//////////////////////////////////////////////////////////////////////////
int Node::nextIndex(int conn, bool& pbcConn) const
{
    int nextIndex(-1);
    pbcConn = false;    

    switch(conn)
    {
    case 0:                                         // iPluss
        if(m_i == m_nX) 
            nextIndex = m_nX * m_nY * m_nZ + 1;     // Outlet
        else if(m_i >= 1 && m_i < m_nX)
            nextIndex = m_index + 1;
        break;
    case 1:                                         // jPluss
        if(m_j == m_nY)
        {
            pbcConn = true;  
            nextIndex = m_index - m_nX * (m_nY-1);
        }
        else if(m_j >= 1 && m_j < m_nY)
            nextIndex = m_index + m_nX;
        break;
    case 2:                                         // kPluss
        if(m_k == m_nZ) 
        {
            pbcConn = true;  
            nextIndex = m_index - m_nX * m_nY * (m_nZ-1); 
        }
        else if(m_k >= 1 && m_k < m_nZ)
            nextIndex = m_index + m_nX * m_nY;
        break;
    case 3:                                         // iMinus
        if(m_i == 1) 
            nextIndex = 0; 
        else if(m_i > 1 && m_i <= m_nX)
            nextIndex = m_index - 1;
        break;
    case 4:                                         // jMinus
        if(m_j == 1) 
        {
            pbcConn = true;  
            nextIndex = m_index + m_nX * (m_nY-1); 
        }
        else if(m_j > 1 && m_j <= m_nY)
            nextIndex = m_index - m_nX;
        break;
    case 5:                                         // kMinus
        if(m_k == 1) 
        {
            pbcConn = true;  
            nextIndex = m_index + m_nX * m_nY * (m_nZ-1); 
        }
        else if(m_k > 1 && m_k <= m_nZ)
            nextIndex = m_index - m_nX * m_nY;
        break;
    default:
        cerr << "=========================" << endl
             << "Error: Crap programmer..." << endl
             << "=========================" << endl;
        exit(-1);
    }
    
    return nextIndex;
}

Node::Node(int index, int nX, int nY, int nZ)
{
    m_nX = nX;
    m_nY = nY;
    m_nZ = nZ;

    if(index == 0)
    {
        m_i = 0;
        m_j = 1;
        m_k = 1;
    }
    else if(index == m_nX*m_nY*m_nZ + 1)
    {
        m_i = m_nX + 1;
        m_j = 1;
        m_k = 1;
    }
    else
    {   
        m_k = (index-1) / (m_nX*m_nY) + 1;
        m_j = ((index-1) - (m_k-1)*(m_nX*m_nY)) / m_nX + 1;
        m_i = index - (m_k-1)*(m_nX*m_nY) - (m_j-1)*m_nX;
    }

    initNode();
}

/////////////////////////////////////////////////////////////////////////////
// The node is initialized by determening if it is inside the network,
// outside or an in/outlet.
/////////////////////////////////////////////////////////////////////////////
void Node::initNode()
{
    if(m_i > 0 && m_i <= m_nX && 
       m_j > 0 && m_j <= m_nY && 
       m_k > 0 && m_k <= m_nZ) 
    {
        m_isOutsideLattice = false;             // We are inside lattice
        m_isInlet = false;
        m_isOutlet = false;
        m_isAtInlet = (m_i == 1);
        m_isAtOutlet = (m_i == m_nX);

        m_index = (m_k-1)*(m_nX*m_nY) + (m_j-1)*m_nX + m_i;
    }
    else if(m_i == 0 && 
            m_j > 0 && m_j <= m_nY && 
            m_k > 0 && m_k <= m_nZ) 
    {
        m_isOutsideLattice = false;             // Inlet node
        m_isInlet = true;
        m_isOutlet = false;
        m_isAtInlet = false;
        m_isAtOutlet = false;

        m_index = 0;
    }
    else if(m_i == (m_nX + 1) && 
            m_j > 0 && m_j <= m_nY && 
            m_k > 0 && m_k <= m_nZ) 
    {
        m_isOutsideLattice = false;             // Outlet node
        m_isInlet = false;
        m_isOutlet = true;
        m_isAtInlet = false;
        m_isAtOutlet = false;

        m_index = m_nX * m_nY * m_nZ + 1;
    }
    else
    {
        m_isOutsideLattice = true;             // We are outside lattice
        m_isInlet = false;
        m_isOutlet = false;
        m_isAtInlet = false;
        m_isAtOutlet = false;

        m_index = -1;
    }
}

void Node::setLocation(double xDim, double yDim, double zDim, double avrThroatLen)
{
    double intXDim(xDim-avrThroatLen), intYDim(yDim-avrThroatLen), intZDim(zDim-avrThroatLen);
    m_xLoc = m_nX > 1 ? (intXDim/(m_nX-1))*(m_i-1)+avrThroatLen/2.0: xDim/2.0;
    m_yLoc = m_nY > 1 ? (intYDim/(m_nY-1))*(m_j-1)+avrThroatLen/2.0: yDim/2.0;
    m_zLoc = m_nZ > 1 ? (intZDim/(m_nZ-1))*(m_k-1)+avrThroatLen/2.0: zDim/2.0;
}

///////////////////////////////////////////////////////////////////
// In/outlet internally has indecies 0 and numPores+1. However when
// writing the data to file this has to be changed to -1 and 0 to
// be compatible with Oren's data format
///////////////////////////////////////////////////////////////////
int Node::indexOren() const
{
    if(m_isInlet)
        return -1;
    else if(m_isOutlet)
        return 0;
    else
        return m_index;
}


