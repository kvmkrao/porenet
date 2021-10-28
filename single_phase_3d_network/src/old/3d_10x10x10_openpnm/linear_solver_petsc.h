/////////////////////////////////////////////////////////////////////////////
/// Author:      Mohammad Sedghi <mhdsedghi@gmail.com>
/// Created:     2018
/// Copyright:   (c) 2018 Mohammad Sedghi
/// Licence:     Commercial-sedghi
/////////////////////////////////////////////////////////////////////////////



#ifndef LINEAR_SOLVER_PETSC_H
#define LINEAR_SOLVER_PETSC_H

#include <vector>

int linear_solver_petsc(const std::vector<double> &val,
                        const std::vector<double> &rhs,
                        std::vector<double> &pressures,
                        const std::vector<int> &row,
                        const std::vector<int> &col,
                        int elementCount,
                        int maxIters,
                        double &resid,
                        int &iters);


#endif // LINEAR_SOLVER_PETSC_H
