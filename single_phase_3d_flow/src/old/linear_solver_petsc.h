
#ifndef LINEAR_SOLVER_PETSC_H
#define LINEAR_SOLVER_PETSC_H

#include <vector>

//#include "data_type.h"


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
