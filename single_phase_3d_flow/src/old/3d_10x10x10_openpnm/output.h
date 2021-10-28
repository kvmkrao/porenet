/////////////////////////////////////////////////////////////////////////////
/// Author:      Mohammad Sedghi <mhdsedghi@gmail.com>
/// Created:     2018
/// Copyright:   (c) 2018 Mohammad Sedghi
/// Licence:     Commercial-sedghi
/////////////////////////////////////////////////////////////////////////////



#ifndef OUTPUT_VTP_H
#define OUTPUT_VTP_H

#include <vector>

int output_vtp(const std::vector<double> &val,
                        const std::vector<double> &rhs,
                        std::vector<double> &pressures,
                        const std::vector<int> &row,
                        const std::vector<int> &col,
                        int elementCount,
                        int maxIters,
                        double &resid,
                        int &iters);


#endif // OUTPUT_VTP_H
