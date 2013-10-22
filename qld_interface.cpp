//DD// Matlab interface for QLD using only bounds 2010/05/07

#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;

#include "qld.h"

#include "mex.h"



#define EPS 1e-12


enum QLDMode
{
    QLD_MODE_CHOLESKY_HESSIAN = 0,
    QLD_MODE_NORMAL_HESSIAN = 1,
};


enum QLDStatus
{
    QLD_OK = 0,
    QLD_MAX_ITER_REACHED = 1,
    QLD_INSUFFICIENT_ACCURACY = 2,
    QLD_INTERNAL_ERROR = 3,
    QLD_INSUFFICIENT_MEMORY = 4,
    QLD_INCONSISTENT_CONSTRAINTS = 5,
    QLD_UNKNOWN_ERROR = 6,
};


class QLDSolver
{
    private:
        QLDMode mode_;

        double* H_; // Hessian
        double* g_; // 

        double* lb_;
        double* ub_;

        double* A_;
        double* b_;


        int num_var_;
        int num_eq_;
        int num_ctr_;

        double eps_;
    

    public:
        QLDSolver()
        {
            mode_ = QLD_MODE_NORMAL_HESSIAN;
            
            num_var_ = 0;
            num_eq_ = 0;
        }

        void setOptions(const double eps, const QLDMode mode)
        {
            eps_ = eps;
            mode_ = mode;
        }

        void setObjective (mxArray * H, mxArray * g)
        {
            num_var_ = mxGetM(H); // row number of Hessian matrix

            H_ = (double*) mxGetPr(H);
            g_ = (double*) mxGetPr(g);
        }

        void setBounds (mxArray * lb, mxArray * ub)
        {
            lb_ = (double*) mxGetPr(lb);
            ub_ = (double*) mxGetPr(ub);
        }


        void setConstraints(mxArray * A, mxArray * b, int num_eq)
        {
            num_eq_ = num_eq;

            if (mxIsEmpty(A) || mxIsEmpty(b))
            {
                A_ = NULL;
                b_ = NULL;
                num_ctr_ = 0;
            }
            else
            {
                A_ = (double*) mxGetPr(A);
                b_ = (double*) mxGetPr(b);
                num_ctr_ = mxGetM(b);
            }
        }


        QLDStatus solve (mxArray *x, mxArray *lambda, int * inconsistent_ctr_idx)
        {
            QLDStatus status;


            int MNN = num_ctr_ + num_var_ + num_var_;
            int IOUT = 0;                           // DESIRED OUTPUT UNIT NUMBER (INTEGER).
            int IFAIL;                              // 
            int IPRINT = 0;                         // OUTPUT CONTROL.
            int LWAR = 3*num_var_*num_var_/2 + 10*num_var_ + 2*num_var_ + 100;
            int LIWAR = num_var_;
            int MODE = static_cast <int> (mode_);

            double *WAR  = new double[LWAR];        // REAL WORKING ARRAY. 
            int    *IWAR = new int   [LIWAR];       // INTEGER WORKING ARRAY. 


            ql_(    &num_ctr_, // Number of constraints. 
                    &num_eq_, // Number of equality constraints. XXX int
                    &num_ctr_, // 1 <= MMAX <= (number of constraints), Row dimension of array A containing linear constraints. 
                    &num_var_, // Number of optimization variables. XXX int
                    &num_var_, // 1 <= NMAX <= (number of variables), Row dimension of C.
                    &MNN, // Must be equal to M+N+N when calling QL, dimension of U.
                    H_, // Hessian or its upper triangular factor of a Cholesky decomposition.
                    g_, // The constant vector of the quadratic objective function.
                    A_, // Matrix of the linear constraints, first ME rows for equalities, others inequalities.
                    b_, // Constant values of linear constraints.
                    lb_, // lower bounds on variables 
                    ub_, // upper bunds on variables
                    (double*) mxGetPr(x), // output: solution
                    (double*) mxGetPr(lambda),  // output:
                        // Lagrange multipliers subject to the 
                        // linear constraints and bounds. The first M locations 
                        // contain the multipliers of the M linear constraints, the 
                        // subsequent N locations the multipliers of the lower 
                        // bounds, and the final N locations the multipliers of the
                        // upper bounds. At the optimal solution, all multipliers 
                        // with respect to inequality constraints should be 
                        // nonnegative.
                    &eps_, // the desired final accuracy
                    &MODE, // 
                        // MODE=0 - The user provides an initial Cholesky factorization
                        //          of C, stored in the upper triangular part of the
                        //          array C.
                        // MODE=1 - A Cholesky decomposition to get the first 
                        //          unconstrained minimizer, is computed internally.
                    &IOUT, // Integer preceding messages.
                    &IFAIL, // output: status
                        // IFAIL=0  : The optimality conditions are satisfied.
                        // IFAIL=1  : The algorithm has been stopped after too many
                        //            MAXIT iterations (40*(N+M)).
                        // IFAIL=2  : Termination accuracy insufficient to satisfy 
                        //            convergence criterion. 
                        // IFAIL=3  : Internal inconsistency of QL, division by zero.
                        // IFAIL=5  : Length of a working array is too short. 
                        // IFAIL>100: Constraints are inconsistent and IFAIL=100+ICON,
                        //            where ICON denotes a constraint causing the conflict.
                    &IPRINT, // 0 -- no output, 1 -- final message
                    WAR, // working array of size LWAR
                    &LWAR, // size of WAR >= 3*NMAX*NMAX/2 + 10*NMAX + MMAX + M + 1
                    IWAR, // working integer array of size LIWAR
                    &LIWAR // size of LIWAR >= N
                );


            *inconsistent_ctr_idx = 0;

            switch (IFAIL)
            {
                case 0:
                    status = QLD_OK;
                    break;
                case 1:
                    status = QLD_MAX_ITER_REACHED;
                    break;
                case 2:
                    status = QLD_INSUFFICIENT_ACCURACY;
                    break;
                case 3:
                    status = QLD_INTERNAL_ERROR;
                    break;
                case 5:
                    status = QLD_INSUFFICIENT_MEMORY;
                    break;
                default:
                    if (IFAIL > 100)
                    {
                        status = QLD_INCONSISTENT_CONSTRAINTS;
                        *inconsistent_ctr_idx = IFAIL - 100;
                    }
                    else
                    {
                        status = QLD_UNKNOWN_ERROR;
                    }
                    break;
            }

            return (status);
        }
};


void mexFunction( int num_output, mxArray *output[], int num_input, const mxArray *input[] )
{
    mxArray *x = mxDuplicateArray(input[0]);
    mxArray *H = mxDuplicateArray(input[1]);
    mxArray *g = mxDuplicateArray(input[2]);
    mxArray *lb = mxDuplicateArray(input[3]);
    mxArray *ub = mxDuplicateArray(input[4]);
    mxArray *A = mxDuplicateArray(input[5]);
    mxArray *b = mxDuplicateArray(input[6]);
    int num_eq = static_cast <int> (round(*mxGetPr(input[7])));
    mxArray *options = mxDuplicateArray(input[8]);

    
    mxArray *lambda = mxCreateDoubleMatrix(mxGetM(A) + mxGetM(lb) + mxGetM(ub), 1, mxREAL);
    mxArray *info = NULL;



// parse options
    double eps = EPS;
    QLDMode hess_mode = QLD_MODE_NORMAL_HESSIAN;    

    if (!mxIsEmpty(options))
    {
        mxArray * hess_mode_option = mxGetField (options, 0, "hessian_mode");
        if (hess_mode_option != NULL)
        {
            int hess_mode_int = static_cast <int> (round(*mxGetPr(hess_mode_option)));
            if (hess_mode_int == 0)
            {
                hess_mode = QLD_MODE_CHOLESKY_HESSIAN;
            } 
            else if (hess_mode_int == 1)
            {
                hess_mode = QLD_MODE_NORMAL_HESSIAN;
            }
            else
            {
                mexErrMsgTxt( "Incorrect hessian mode.");
            }
        }


        mxArray * eps_option = mxGetField (options, 0, "tolerance");
        if (eps_option != NULL)
        {
            if (mxIsDouble(eps_option))
            {
                eps = *mxGetPr(eps_option);
            }
            else
            {
                mexErrMsgTxt( "Incorrect type of the tolerance.");
            }
        }
    }


// solve the problem
    QLDSolver qld;
    

    qld.setObjective (H, g);
    qld.setBounds (lb, ub);
    qld.setConstraints (A, b, num_eq);

    qld.setOptions(eps, hess_mode);


    int inconsistent_ctr_idx = 0;
    QLDStatus qld_status = qld.solve(x, lambda, &inconsistent_ctr_idx);


// process results
    // solution
    output[0] = mxDuplicateArray(x);


    // info
    int num_info_fields = 2;
    const char *info_field_names[] = {
        "info",
        "inconsistent_constraint_index",
    };

    output[1] = mxCreateStructMatrix(1, 1, num_info_fields, info_field_names);

    mxArray *info_status = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
    ((INT32_T *) mxGetData (info_status))[0] = static_cast <int> (qld_status);
    mxSetField (output[1], 0, "info", info_status); 

    mxArray *info_ctr = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
    ((INT32_T *) mxGetData (info_ctr))[0] = static_cast <int> (inconsistent_ctr_idx);
    mxSetField (output[1], 0, "inconsistent_constraint_index", info_ctr); 

    
    // lagrange multipliers
    output[2] = mxDuplicateArray(lambda);
        
    return;
}
