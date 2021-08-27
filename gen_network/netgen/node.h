#ifndef NODE_H
#define NODE_H


//////////////////////////////////////////////////////////////////////
// All connection are aligned in an right handed coord syst, with 
// connections going from 0 to 5 (0:i+, 1:j+, 2:k+, 3:i-, 4:j-, 5:k- 
//
//                    2   1
//                    | /
//                3 --|--  0
//                   /|
//                 4  5
//
/////////////////////////////////////////////////////////////////////
class Node
{

    friend ostream& operator<< (ostream&, Node&);

public:
    
    Node(int i, int j, int k, int nX, int nY, int nZ);
    Node(int index, int nX, int nY, int nZ);

    bool operator==(const Node& rhs) const {return m_index == rhs.index();}
    bool operator>=(const Node& rhs) const {return m_index >= rhs.index();}
    bool operator<=(const Node& rhs) const {return m_index <= rhs.index();}
    bool operator>(const Node& rhs) const {return m_index > rhs.index();}
    bool operator<(const Node& rhs) const {return m_index < rhs.index();}

    int index() const {return m_index;}
    int indexOren() const;
    int nextIndex(int conn, bool& pbcConn) const;
    int i() const {return m_i;}
    int j() const {return m_j;}
    int k() const {return m_k;}
    void setLocation(double xDim, double yDim, double zDim, double avrThroatLen);

    bool isOutsideLattice() const {return m_isOutsideLattice;}
    bool isInlet() const {return m_isInlet;}
    bool isOutlet() const {return m_isOutlet;}
    bool isInOrOutlet() const {return m_isInlet || m_isOutlet;}
    bool isAtInlet() const {return m_isAtInlet;}
    bool isAtOutlet() const {return m_isAtOutlet;}

private:

    void initNode();

    double                      m_xLoc, m_yLoc, m_zLoc;     //The node position in physical space
    int                         m_i, m_j, m_k;              //The node coordinate
    int                         m_nX, m_nY, m_nZ;           //The size of the lattice
    bool                        m_isInlet, m_isOutlet;      //Is the node inlet or outlet node
    bool                        m_isAtInlet, m_isAtOutlet;  //Is the node at inlet or outlet
    bool                        m_isOutsideLattice;         //Is node outside lattice
    int                         m_index;                    //Single consecutive index

};

#endif

