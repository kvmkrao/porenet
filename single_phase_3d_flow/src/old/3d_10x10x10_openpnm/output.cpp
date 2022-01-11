
#include "linear_solver_petsc.h"
#include <petscksp.h>
#include <iostream>

#include <petscsys.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscpc.h>
#include <petscviewer.h>

int linear_solver_petsc(const std::vector<double> &val,
                        const std::vector<double> &rhs,
                        std::vector<double> &pressures,
                        const std::vector<int> &row,
                        const std::vector<int> &col,
                        int elementCount,
                        int maxIters,
                        double  &resid,
                        int &iters)
{

   Vec        x_pet, rhs_pet;    /* approx solution, RHS, */
   Mat        A;        /* linear system matrix */
   KSP        ksp;      /* linear solver context */

   PetscReal  norm = 0, tol = 1.0e-10;     /* norm of solution error */
   PetscInt   RHS_NDX, ROW_NDX, COL_NDX,  X_NDX,
              BS = 1, max_iter = maxIters,
              Istart, Iend, its;

   PetscErrorCode ierr;
   PetscBool      flg = PETSC_FALSE;
   PetscScalar    VAL;

#if defined(PETSC_USE_LOG)
   PetscLogStage stage;
#endif

   PetscInt num_row = rhs.size();
   PetscInt num_col = num_row;

   printf("num row col %d %d %d\n", num_row, num_col, row.size()); 
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //     Compute the matrix and right-hand-side vector that define
    //     the linear system, Ax = b.
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // Create parallel matrix, specifying only its global dimensions.
    // When using MatCreate(), the matrix format can be specified at
    // runtime. Also, the parallel partitioning of the matrix is
    // determined by PETSc at runtime.

    // Performance tuning note:  For problems of substantial size,
    // preallocation of matrix memory is crucial for attaining good
    // performance. See the matrix chapter of the users manual for details.

//    PetscInt dia_nnz_pet[dia_nnz.size()];
//    PetscInt off_nnz_pet[off_nnz.size()];

   ierr = MatCreate(MPI_COMM_WORLD,&A); CHKERRQ(ierr);
   ierr = MatSetType(A,MATMPIAIJ); CHKERRQ(ierr);
   ierr = MatSetSizes(A,num_row,num_col,elementCount,elementCount); CHKERRQ(ierr);
//    ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,elementCount,elementCount); CHKERRQ(ierr);
   ierr = MatSetFromOptions(A);
   ierr = MatMPIAIJSetPreallocation(A,NULL,NULL,NULL,NULL); CHKERRQ(ierr);

   ierr = MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE); CHKERRQ(ierr);
   ierr = MatGetOwnershipRange(A,&Istart,&Iend); CHKERRQ(ierr);
//    std::cout<<"\t\t\t Istart: " << Istart << "\t Iend: "<<Iend<<std::endl;
   if ( (Iend-Istart) != num_row)
        std::cerr<<"Iend-Istart: " << (Iend-Istart) << " is not the same as num_row: " << num_row <<std::endl;

   for (auto i=0;i<val.size();i++){
        ROW_NDX = row[i];
        COL_NDX = col[i];
        VAL = val[i];
//	if(COL_NDX <= 12348) {
            ierr = MatSetValues(A,1,&ROW_NDX,1,&COL_NDX,&VAL,INSERT_VALUES); CHKERRQ(ierr);
	    //std::cout << ROW_NDX<<" "<<COL_NDX<<" "<<VAL << std::endl; 
//	}

    }
//return 0; 

    ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

    printf("assembled matrix \n");
//    ierr = MatSetOption(A, MAT_SYMMETRIC, PETSC_TRUE); CHKERRQ(ierr);

//    ierr = VecCreateMPI(PETSC_COMM_WORLD,num_row,elementCount,&rhs_pet); CHKERRQ(ierr);
    ierr = VecCreate(MPI_COMM_WORLD,&rhs_pet); CHKERRQ(ierr);
    ierr = VecSetSizes(rhs_pet,num_row,elementCount);CHKERRQ(ierr);
    ierr = VecSetFromOptions(rhs_pet);CHKERRQ(ierr);

    ierr = VecZeroEntries(rhs_pet); CHKERRQ(ierr);
    
    for (auto i=0;i<rhs.size();i++){
	  RHS_NDX = i; 
//        RHS_NDX = rhs_ndx[i];
        VAL = rhs[i];
        ierr = VecSetValues(rhs_pet,1,&RHS_NDX,&VAL,INSERT_VALUES); CHKERRQ(ierr);
    }
    
    ierr = VecAssemblyBegin(rhs_pet); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(rhs_pet); CHKERRQ(ierr);

    printf("assembled vector \n"); 
//    PetscScalar x_sol[pressures.size()];
    std::vector<PetscScalar> x_sol(pressures.size(),0.0);

    ierr = VecCreateMPIWithArray(MPI_COMM_WORLD,BS,num_row,PETSC_DECIDE,&x_sol[0],&x_pet); CHKERRQ(ierr);
    ierr = VecZeroEntries(x_pet); CHKERRQ(ierr);


    for (auto i=0;i<pressures.size();i++){
        PetscScalar x_guess = pressures[i];
        X_NDX = i; //rhs_ndx[i];
        ierr = VecSetValues(x_pet, 1, &X_NDX, &x_guess, INSERT_VALUES);
    }

    printf("set vector x_pet \n"); 
    ierr = VecAssemblyBegin(x_pet); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(x_pet); CHKERRQ(ierr);

    //    !!!!!!!!!!!!!!!!!
    //    !SOLVE WITH PETSC
    //    !!!!!!!!!!!!!!!!!

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Create the linear solver and set various options
         - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

    //     Set operators. Here the matrix that defines the linear system
    //     also serves as the preconditioning matrix.
    ierr = KSPCreate(MPI_COMM_WORLD,&ksp); CHKERRQ(ierr);
    ierr = KSPSetOperators(ksp,A,A); CHKERRQ(ierr);

    ierr = KSPSetInitialGuessNonzero(ksp,PETSC_TRUE); CHKERRQ(ierr);

    ierr = KSPSetTolerances(ksp,tol,PETSC_DEFAULT, PETSC_DEFAULT,max_iter); CHKERRQ(ierr);
    ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);

    /*
         Set linear solver defaults for this problem (optional).
         - By extracting the KSP and PC contexts from the KSP context,
           we can then directly call any KSP and PC routines to set
           various options.
         - The following two statements are optional; all of these
           parameters could alternatively be specified at runtime via
           KSPSetFromOptions().  All of these defaults can be
           overridden at runtime, as indicated below.
      */

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                          Solve the linear system
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    ierr = KSPSetUp(ksp);
    ierr = KSPSolve(ksp,rhs_pet,x_pet); CHKERRQ(ierr);

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                          Check solution and clean up
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    /*
         Check the error
      */
//        ierr = VecAXPY(x_pet,-1.0,u);
    ierr = KSPGetResidualNorm(ksp,&norm); CHKERRQ(ierr);
    ierr = KSPGetIterationNumber(ksp,&its); CHKERRQ(ierr);

    iters = its;
    resid = (double)norm;
    /*
         Print convergence information.  PetscPrintf() produces a single
         print statement from all processes that share a communicator.
         An alternative is PetscFPrintf(), which prints to a file.
      */

//    ierr = PetscPrintf(PETSC_COMM_WORLD,"Residual Norm %g iterations %D\n",(real_tp)norm,its);

    /*
         Free work space.  All PETSc objects should be destroyed when they
         are no longer needed.
      */
    for (auto i=0;i<pressures.size();i++) {	    
        pressures[i] = x_sol[i];
//	printf("pressure at i is %i %d\n", i, pressures[i]); 
    }

//    return 0; 
    ierr = VecDestroy(&rhs_pet);
    ierr = VecDestroy(&x_pet);
    ierr = MatDestroy(&A);
    ierr = KSPDestroy(&ksp);

    /*
         Always call PetscFinalize() before exiting a program.  This routine
           - finalizes the PETSc libraries as well as MPI
           - provides summary and diagnostic information if certain runtime
             options are chosen (e.g., -log_view).
      */

  return 0; //ierr;
}






