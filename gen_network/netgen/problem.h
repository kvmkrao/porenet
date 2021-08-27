
#ifndef PROBLEM_H
#define PROBLEM_H

enum Proportion {Square = 0, Circ, Triang};
enum Property {Min = 0, Max, Delta, Eta};

class Problem
{
    friend istream& operator>> (istream&, Problem&);
    
public:
    
    Problem(const string inputFile);  //Constructor takes care of loading data

    void writeData() const;
        
private:
  
    void initNetworkModel();
    double evalShapeFactor(const double[], const double[]) const;
    void findNetworkModelSize();
    void reduceConnectionNumber();
    void assignThroatIndex();
    void setPoreLocation();
    void createPores();
    void connectPoresWithThroats();
    double weibull(double, double, double, double) const;
    void checkNetworkIntegrity() const;
   
    static const int                MAX_CONN_NUM;
    static const double             PI;

    double                          m_porosity;
    double                          m_xDim, m_yDim, m_zDim; 
    string                          m_outFileNameBase;
    int                             m_nX, m_nY, m_nZ;
    int                             m_numPores;
    vector< Pore * >                m_pores;
    list< Throat * >                m_throats;
    bool                            m_periodicBC;
    
    double                          m_averageThroatLength;
    double                          m_averConnectionNum;
    double                          m_actualConnectionNumber;
    double                          m_clayProportion;
    double                          m_throatRadWeibull[4];                        
    double                          m_throatLenWeibull[4];      // All these are ordered as follows:                        
    double                          m_aspectRatioWeibull[4];    // min, max, delta, eta
    double                          m_triangleGWeibull[4];
    double                          m_throatShapeProp[3]; // Square, Circular, Triangular
    double                          m_poreShapeProp[3];
    

};

#endif
