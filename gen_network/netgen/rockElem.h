#ifndef ROCKELEM_H
#define ROCKELEM_H

class Node;
class Throat;
class Pore;

//======================= Base Class ===============================

class RockElem   
{
public:

    RockElem(double shapeFactor) : m_shapeFactor(shapeFactor) {}        
    virtual ~RockElem(){}

    double radius() const {return m_radius;}
    double length() const {return m_length;}
    double totVolume() const {return m_volume + m_clayVolume;}
 
    
protected:

    const double                m_shapeFactor;
    double                      m_radius;
    double                      m_volume;
    double                      m_clayVolume;
    double                      m_length;

};

//==================================== Pores =====================================

class Pore : public RockElem
{
public:

    Pore(Node *, double, int);
    virtual ~Pore(){delete m_node;}

    void finaliseInit(double, double);
    double halfLength() const {return m_length/2.0;}
    const Node* node() const {return m_node;}
    Node* node() {return m_node;}
    void writeData(std::ostream&, std::ostream&) const;

    Throat* connectingThroat(int conn) {return m_throats[conn];}
    void addThroat(int conn, Throat *throat) {m_throats[conn] = throat;}
    void removeThroat(Throat *);
    void checkState() const;
    int connectionNum() const;
    bool containsThroat(const Throat*) const;

private:

    Node*                       m_node;
    int                         m_connectionNumber;
    std::vector< Throat* >      m_throats;

};

//===================================== Throat ===================================

class Throat : public RockElem
{
public:

    Throat(Pore*, Pore*, double, double, double, double, bool);
    virtual ~Throat(){}

    Pore* connectingPore(int conn) {return m_pores[conn];}
    const Pore* nextPore(const Pore*) const;
    void setIndex(int index) {m_index = index;}
    int index() const {return m_index;}
    void writeData(std::ostream&, std::ostream&) const;
    void checkState() const;
    int connectionNum() const {return static_cast< int >(m_pores.size());}
    bool pbcThroat() const {return m_pbcThroat;}

private:

    std::vector< Pore* >        m_pores;
    int                         m_index;
    const bool                  m_pbcThroat;

};

#endif
