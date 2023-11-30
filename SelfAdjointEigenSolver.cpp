/* SelfAdjointEigenSolver.c calculate the eigenvalues and eigenvectors of a
 * selfadjoint matrix (only upper triangle of the matrix is needed);
 * Syntax: E = SelfAdjointEigenSolver(A)
 * or [E V] = SelfAdjointEigenSolver(A)
 */
#include <math.h>
#include "mex.h"

// use the Eigen library
#include <Eigen/Eigenvalues>
using Eigen::MatrixXcd;
using Eigen::MatrixXd;
using Eigen::VectorXd;
typedef std::complex<double> ComplexType;

#define IS_SQARE_FULL(P) ( mxIsDouble(A_IN) && mxGetNumberOfDimensions(P) == 2 && mxGetM(P) == mxGetN(P) && !mxIsSparse(P))
#define IS_SQARE_FULL_COMPLEX(P) (mxGetNumberOfDimensions(P) == 2 && mxGetM(P) == mxGetN(P) && !mxIsSparse(P))

#define IS_REAL_2D_FULL_DOUBLE(P) (!mxIsComplex(P) && \
mxGetNumberOfDimensions(P) == 2 && !mxIsSparse(P) && mxIsDouble(P))
#define IS_REAL_SCALAR(P) (IS_REAL_2D_FULL_DOUBLE(P) && mxGetNumberOfElements(P) == 1)
//mxIsComplex(A)
//M = mxGetM(A)
//N = mxGetN(A)

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Macros for the ouput and input arguments */
#define V_OUT plhs[0]
#define A_IN prhs[0]
    double *E, *A, *V;
    int M, N;
    if(nrhs < 1 || nrhs > 1) /* Check the number of arguments */
        mexErrMsgTxt("Wrong number of input arguments.");
    else if(nlhs > 2)
        mexErrMsgTxt("Too many output arguments.");
    if(!IS_SQARE_FULL(A_IN))
        mexErrMsgTxt("A must be a 2D square double or complex array.");
    M = mxGetM(A_IN); /* Get the dimensions of A */
    // eigenvalues should be real, so create real output
    plhs[nlhs-1] = mxCreateDoubleMatrix(M, 1, mxREAL); /* Create the output matrix */
    E = mxGetPr(plhs[nlhs-1]); /* Get the pointer to the data of B */
    if (!mxIsComplex(A_IN))
    {
        // define real double matrix as needed by Eigen
        MatrixXd H(M,M);
        A = mxGetPr(A_IN); /* Get the pointer to the data of A */
        // put in data to matrix H (only upper triangle!)
        for (size_t orb1 = 0; orb1 < M; ++orb1)
            for (size_t orb2 = 0; orb2 < orb1+1; ++orb2)
                H(orb1,orb2)=A[orb1 + M*orb2];
        // define the eigensolver object       
        Eigen::SelfAdjointEigenSolver<MatrixXd> es;
        // diagonalize the matrix
        es.compute(H);
        // write back the result
        for (size_t orb1 = 0; orb1 < M; ++orb1)
            E[orb1]=es.eigenvalues()[orb1];
        if (nlhs > 1)   
            // also write back the eigenvectors if needed
        {
            V_OUT = mxCreateDoubleMatrix(M, M, mxREAL);
            double *V = mxGetPr(V_OUT);
            for (size_t orb1 = 0; orb1 < M; ++orb1)
                for (size_t orb2 = 0; orb2 < M; ++orb2)
                    V[orb1 + M*orb2]=es.eigenvectors()(orb1,orb2);
        }
    }
    else
    {
        // define the matrix as needed by the Eigen library
        MatrixXcd H(M,M);
        // pointers to the real and imaginary part
        double *Ar = mxGetPr(A_IN);
        double *Ai = mxGetPi(A_IN);
        // put the values into the matrix H (only upper triangle!)
        for (size_t orb1 = 0; orb1 < M; ++orb1)
            for (size_t orb2 = 0; orb2 < orb1+1; ++orb2)
                H(orb1,orb2)=ComplexType(Ar[orb1 + M*orb2],Ai[orb1 + M*orb2]);
        // define the eigensolver object
        Eigen::SelfAdjointEigenSolver<MatrixXcd> es;
        // diagonalize the matrix
        es.compute(H);
        // write back the result
        for (size_t orb1 = 0; orb1 < M; ++orb1)
            E[orb1]=es.eigenvalues()[orb1];
        // also write back the eigenvectors if needed
        if (nlhs > 1)
        {
            V_OUT = mxCreateDoubleMatrix(M, M, mxCOMPLEX);
            double *Vr = mxGetPr(V_OUT);
            double *Vi = mxGetPi(V_OUT);
            for (size_t orb1 = 0; orb1 < M; ++orb1)
                for (size_t orb2 = 0; orb2 < M; ++orb2)
                {
                    Vr[orb1 + M*orb2]=real(es.eigenvectors()(orb1,orb2));
                    Vi[orb1 + M*orb2]=imag(es.eigenvectors()(orb1,orb2));
                }
        }
    }
    return;
}
