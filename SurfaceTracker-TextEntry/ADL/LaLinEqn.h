#ifndef LA_LIN_EQN_HDR
#define LA_LIN_EQN_HDR
/*******************************************************************************
This header file contains the C/C++ prototypes for the LAPACK routines.
Routines in this file are specific for solving systems of linear equations.

Abed M. Hammoud
23Aug94-Stealth
Copyright 1993-1995
*******************************************************************************/
// NOMANPAGE
#ifdef __cplusplus
extern "C" {
#endif


/*
*  -- LAPACK driver routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993 
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LDB, M, N, NRHS, RANK
      REAL               RCOND
*     ..
*     .. Array Arguments ..
      INTEGER            JPVT( * )
      REAL               A( LDA, * ), B( LDB, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  SGELSX computes the minimum-norm solution to a real linear least
*  squares problem:
*      minimize || A * X - B ||
*  using a complete orthogonal factorization of A.  A is an M-by-N
*  matrix which may be rank-deficient.
*
*  Several right hand side vectors b and solution vectors x can be 
*  handled in a single call; they are stored as the columns of the
*  M-by-NRHS right hand side matrix B and the N-by-NRHS solution
*  matrix X.
*
*  The routine first computes a QR factorization with column pivoting:
*      A * P = Q * [ R11 R12 ]
*                  [  0  R22 ]
*  with R11 defined as the largest leading submatrix whose estimated
*  condition number is less than 1/RCOND.  The order of R11, RANK,
*  is the effective rank of A.
*
*  Then, R22 is considered to be negligible, and R12 is annihilated
*  by orthogonal transformations from the right, arriving at the
*  complete orthogonal factorization:
*     A * P = Q * [ T11 0 ] * Z
*                 [  0  0 ]
*  The minimum-norm solution is then
*     X = P * Z' [ inv(T11)*Q1'*B ]
*                [        0       ]
*  where Q1 consists of the first RANK columns of Q.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of
*          columns of matrices B and X. NRHS >= 0.
*
*  A       (input/output) REAL array, dimension (LDA,N)
*          On entry, the M-by-N matrix A.
*          On exit, A has been overwritten by details of its
*          complete orthogonal factorization.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  B       (input/output) REAL array, dimension (LDB,NRHS)
*          On entry, the M-by-NRHS right hand side matrix B.
*          On exit, the N-by-NRHS solution matrix X.
*          If m >= n and RANK = n, the residual sum-of-squares for
*          the solution in the i-th column is given by the sum of
*          squares of elements N+1:M in that column.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B. LDB >= max(1,M,N).
*
*  JPVT    (input/output) INTEGER array, dimension (N)
*          On entry, if JPVT(i) .ne. 0, the i-th column of A is an
*          initial column, otherwise it is a free column.  Before
*          the QR factorization of A, all initial columns are
*          permuted to the leading positions; only the remaining
*          free columns are moved as a result of column pivoting
*          during the factorization.
*          On exit, if JPVT(i) = k, then the i-th column of A*P
*          was the k-th column of A.
*
*  RCOND   (input) REAL
*          RCOND is used to determine the effective rank of A, which
*          is defined as the order of the largest leading triangular
*          submatrix R11 in the QR factorization with pivoting of A,
*          whose estimated condition number < 1/RCOND.
*
*  RANK    (output) INTEGER
*          The effective rank of A, i.e., the order of the submatrix
*          R11.  This is the same as the order of the submatrix T11
*          in the complete orthogonal factorization of A.
*
*  WORK    (workspace) REAL array, dimension
*                      (max( min(M,N)+3*N, 2*min(M,N)+NRHS )),
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*/

#ifdef __cplusplus

// Single precision.
   extern void F77CALL(sgelsx)(const int& M, const int& N, const int& NRHS, 
	float A[], const int& LDA, float B[], const int& LDB, int JPVT[],
	const float& RCOND, int *RANK, float WORK[], int *INFO); 

// Double precision.
   extern void F77CALL(dgelsx)(const int& M, const int& N, const int& NRHS, 
	double A[], const int& LDA, double B[], const int& LDB, int JPVT[],
	const double& RCOND, int *RANK, double WORK[], int *INFO); 

#else

/* Single precision */
   void F77CALL(sgelsx)(int const *M, int const *N, int const *NRHS, 
	float A[], int const *LDA, float B[], int const *LDB, int JPVT[], 
	float const *RCOND, int *RANK, float WORK[], int *INFO); 

/* Double precision */
   void F77CALL(dgelsx)(int const *M, int const *N, int const *NRHS, 
	double A[], int const *LDA, double B[], int const *LDB, int JPVT[], 
	double const *RCOND, int *RANK, double WORK[], int *INFO);

#endif

#endif

/*
*  -- LAPACK driver routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993 
*
*     .. Scalar Arguments ..
      CHARACTER          EQUED, FACT, TRANS
      INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS
      REAL               RCOND
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * ), IWORK( * )
      REAL               A( LDA, * ), AF( LDAF, * ), B( LDB, * ),
     $                   BERR( * ), C( * ), FERR( * ), R( * ),
     $                   WORK( * ), X( LDX, * )
*     ..
*
*  Purpose
*  =======
*
*  SGESVX uses the LU factorization to compute the solution to a real
*  system of linear equations
*     A * X = B,
*  where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
*
*  Error bounds on the solution and a condition estimate are also
*  provided.
*
*  Description
*  ===========
*
*  The following steps are performed:
*
*  1. If FACT = 'E', real scaling factors are computed to equilibrate
*     the system:
*        TRANS = 'N':  diag(R)*A*diag(C)     *inv(diag(C))*X = diag(R)*B
*        TRANS = 'T': (diag(R)*A*diag(C))**T *inv(diag(R))*X = diag(C)*B
*        TRANS = 'C': (diag(R)*A*diag(C))**H *inv(diag(R))*X = diag(C)*B
*     Whether or not the system will be equilibrated depends on the
*     scaling of the matrix A, but if equilibration is used, A is
*     overwritten by diag(R)*A*diag(C) and B by diag(R)*B (if TRANS='N')
*     or diag(C)*B (if TRANS = 'T' or 'C').
*
*  2. If FACT = 'N' or 'E', the LU decomposition is used to factor the
*     matrix A (after equilibration if FACT = 'E') as
*        A = P * L * U,
*     where P is a permutation matrix, L is a unit lower triangular
*     matrix, and U is upper triangular.
*
*  3. The factored form of A is used to estimate the condition number
*     of the matrix A.  If the reciprocal of the condition number is
*     less than machine precision, steps 4-6 are skipped.
*
*  4. The system of equations is solved for X using the factored form
*     of A.
*
*  5. Iterative refinement is applied to improve the computed solution
*     matrix and calculate error bounds and backward error estimates
*     for it.
*
*  6. If FACT = 'E' and equilibration was used, the matrix X is 
*     premultiplied by diag(C) (if TRANS = 'N') or diag(R) (if
*     TRANS = 'T' or 'C') so that it solves the original system
*     before equilibration.
*/

#ifdef __cplusplus

// Single precision.
   extern void F77CALL(sgesvx)(char const& FACT, char const& TRANS, 
	int const& N, int const& NRHS, float A[], int const& LDA, 
	float AF[], int const& LDAF, int IPIV[], char &EQUED, float R[], 
	float C[], float B[], int const& LDB, float X[], int const& LDX, 
	float *RCOND, float FERR[], float BERR[], float WORK[], 
	int IWORK[], int *INFO); 

// Double precision.
   extern void F77CALL(dgesvx)(const char& FACT, const char& TRANS, 
	const int& N, const int& NRHS, double A[], const int& LDA, 
	double AF[], const int& LAF, int IPIV[], char &EQUED, double R[], 
	double C[], double B[], const int& LDB, double X[], const int& LDX, 
	double *RCOND, double FERR[], double BERR[], double WORK[], 
	int IWORK[], int *INFO);

#else

/* Single precision */
   void F77CALL(sgesvx)(char const *FACT, char const *TRANS, 
	int const *N, int const *NRHS, float A[], int const *LDA, 
	float AF[], int const *LAF, int IPIV[], char *EQUED, float R[], 
	float C[], float B[], int const *LDB, float X[], int const *LDX, 
	float *RCOND, float FERR[], float BERR[], float WORK[], 
	int IWORK[], int *INFO); 

/* Double precision */
   void F77CALL(dgesvx)(char const *FACT, char const *TRANS, 
	int const *N, int const *NRHS, double A[], int const *LDA, 
	double AF[], int const *LAF, int IPIV[], char *EQUED, double R[], 
	double C[], double B[], int const *LDB, double X[], int const *LDX, 
	double *RCOND, double FERR[], double BERR[], double WORK[], 
	int IWORK[], int *INFO); 

#endif


#ifdef __cplusplus
}
#endif
