#ifndef LAPACK_HDR
#define LAPACK_HDR
// NOMANPAGE

// LAPACK header files.

#define adl_name2(a,b) a##b

#define F77CALL(x) adl_name2(x,_)

// #include <ADL/ADL_globals.h>
#include "LaNorm.h"
#include "LaLU.h"
#include "LaInv.h"
#include "LaSVD.h"
#include "LaCond.h"
#include "LaLinEqn.h"
#include "LaBalance.h"

/*******************************************************************************
* This class encapsulate the fortran LAPACK library subroutines.
*
* Abed M. Hammoud
* 24Aug94-Stealth
* Copyright 1993-1995
*******************************************************************************/
class Lapack {

public:

// LU Factorization.
   inline static void xgetrf(int M, int N, float A[],
	int LDA, int IPIV[], int *INFO);

// Double precision.
   inline static void xgetrf(int M, int N, double A[],
	int LDA, int IPIV[], int *INFO);

// Inverse of a matrix.
   inline static void xgetri(int N, float A[],
	int LDA, const int IPIV[], float WORK[],
	int LWORK, int *INFO);

// Double precision.
   inline static void xgetri(int N, double A[],
	int LDA, const int IPIV[], double WORK[],
	int LWORK, int *INFO);

// Singular Value Decomposition.
   inline static void xgesvd(char const &JOBU, char const &JOBVT,
	int M, int N, float A[], int LDA,
	float S[], float U[], int LDU, float VT[], int LDVT,
	float WORK[], int LWORK, int *INFO);

// Double precision.
   inline static void xgesvd(char const &JOBU, char const &JOBVT,
	int M, int N, double A[], int LDA, double S[], double U[],
	int LDU, double VT[], int LDVT, double WORK[], int LWORK, int *INFO);

// Matrix conditions, (RCOND).
   inline static void xgecon(char const &NORM, int N, const float A[],
	int LDA, float ANORM, float *RCOND, float WORK[],
	int IWORK[], int *INFO);

// Double precision.
   inline static void xgecon(char const &NORM, int N, const double A[],
	int LDA, double ANORM, double *RCOND, double WORK[],
	int IWORK[], int *INFO);

// Matrix norms.
/*
   inline static float xlange(char const &NORM, int M, int N,
    const float A[], int LDA, float WORK[]);

// Double precision.
   inline static double xlange(char const &NORM, int M, int N,
    const double A[], int LDA, double WORK[]);
*/

   inline static float xlange(char const &NORM, int M, int N,
    const float *A, int LDA, float *WORK);

// Double precision.
   inline static double xlange(char const &NORM, int M, int N,
    const double *A, int LDA, double *WORK);

   inline static double ktest(char const &NORM, int M, int N,
    const double *A, int LDA, double *WORK);

// Find a least square solution using complete orthogonal factorization.
	inline static void xgelsx(int M, int N, int NRHS,
	 float A[], int LDA, float B[], int LDB, int JPVT[],
	 float RCOND, int *RANK, float WORK[], int *INFO);

// Double precision.
	inline static void xgelsx(int M, int N, int NRHS,
	 double A[], int LDA, double B[], int LDB, int JPVT[],
	 double RCOND, int *RANK, double WORK[], int *INFO);

// Use LU factorization to compute the solution to A X = B.
   inline static void xgesvx(char const &FACT, char const &TRANS,
    int N, int NRHS, float A[], int LDA,
    float AF[], int LDAF, int IPIV[], char &EQUED, float R[],
    float C[], float B[], int LDB, float X[], int LDX,
    float *RCOND, float FERR[], float BERR[], float WORK[],
    int IWORK[], int *INFO);

// Double precision.
   inline static void xgesvx(char const &FACT, char const &TRANS,
    int N, int NRHS, double A[], int LDA,
    double AF[], int LAF, int IPIV[], char &EQUED, double R[],
    double C[], double B[], int LDB, double X[], int LDX,
    double *RCOND, double FERR[], double BERR[], double WORK[],
    int IWORK[], int *INFO);

// Matrix Equilibrate.
   inline static void xgeequ(int M, int N, const float A[],
       int LDA, float R[], float C[], float *ROWCND,
       float *COLCND, float *AMAX, int *INFO);

// Double precision.
   inline static void xgeequ(int M, int N, const double A[],
       int LDA, double R[], double C[], double *ROWCND,
       double *COLCND, double *AMAX, int *INFO);

// Matrix Equilibrate, (cont).
   inline static void xlaqge(int M, int N, const float A[], int LDA,
    const float R[], const float C[], float ROWCND, float COLCND,
    float AMAX, char *EQUED);

// Double precision.
   inline static void xlaqge(int M, int N, const double A[],
       int LDA, const double R[], const double C[],
       double ROWCND, double COLCND, double AMAX, char *EQUED);

};


/*******************************************************************************
*
*  SGETRF computes an LU factorization of a general M-by-N matrix A
*  using partial pivoting with row interchanges.
*
*  The factorization has the form
*     A = P * L * U
*  where P is a permutation matrix, L is lower triangular with unit
*  diagonal elements (lower trapezoidal if m > n), and U is upper
*  triangular (upper trapezoidal if m < n).
*
*  This is the right-looking Level 3 BLAS version of the algorithm.
*
*******************************************************************************/
inline void Lapack::xgetrf(int M, int N, float A[],
	int LDA, int IPIV[], int *INFO) {

    F77CALL(sgetrf)(M, N, A, LDA, IPIV, INFO);
}

// Double precision.
inline void Lapack::xgetrf(int M, int N, double A[],
	int LDA, int IPIV[], int *INFO) {

    F77CALL(dgetrf)(M, N, A, LDA, IPIV, INFO);
}

/*******************************************************************************
*
*  SGETRI computes the inverse of a matrix using the LU factorization
*  computed by SGETRF.
*
*  This method inverts U and then computes inv(A) by solving the system
*  inv(A)*L = inv(U) for inv(A).
*
*******************************************************************************/
inline void Lapack::xgetri(int N, float A[],
	int LDA, const int IPIV[], float WORK[], int LWORK, int *INFO) {

    F77CALL(sgetri)(N, A, LDA, IPIV, WORK, LWORK, INFO);
}

// Double precision.
inline void Lapack::xgetri(int N, double A[],
	int LDA, const int IPIV[], double WORK[], int LWORK, int *INFO) {

    F77CALL(dgetri)(N, A, LDA, IPIV, WORK, LWORK, INFO);
}

/*******************************************************************************
*
*  SGESVD computes the singular value decomposition (SVD) of a real
*  M-by-N matrix A, optionally computing the left and/or right singular
*  vectors. The SVD is written
*
*       A = U * SIGMA * transpose(V)
*
*  where SIGMA is an M-by-N matrix which is zero except for its
*  min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and
*  V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA
*  are the singular values of A; they are real and non-negative, and
*  are returned in descending order.  The first min(m,n) columns of
*  U and V are the left and right singular vectors of A.
*
*  Note that the routine returns V**T, not V.
*
*******************************************************************************/
inline void Lapack::xgesvd(char const &JOBU, char const &JOBVT,
	int M, int N, float A[], int LDA,
	float S[], float U[], int LDU, float VT[], int LDVT,
	float WORK[], int LWORK, int *INFO) {

    F77CALL(sgesvd)(JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT,
	WORK, LWORK, INFO);
}

// Double precision.
inline void Lapack::xgesvd(char const &JOBU, char const &JOBVT,
	int M, int N, double A[], int LDA,
	double S[], double U[], int LDU, double VT[], int LDVT,
	double WORK[], int LWORK, int *INFO) {

    F77CALL(dgesvd)(JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT,
	WORK, LWORK, INFO);

}

/*******************************************************************************
*
*  SGECON estimates the reciprocal of the condition number of a general
*  real matrix A, in either the 1-norm or the infinity-norm, using
*  the LU factorization computed by SGETRF.
*
*  An estimate is obtained for norm(inv(A)), and the reciprocal of the
*  condition number is computed as
*     RCOND = 1 / ( norm(A) * norm(inv(A)) ).
*
*******************************************************************************/
// Single precision.
inline void Lapack::xgecon(char const &NORM, int N, const float A[], int LDA,
	float ANORM, float *RCOND, float WORK[], int IWORK[], int *INFO) {

    F77CALL(sgecon)(NORM, N, A, LDA, ANORM, RCOND, WORK, IWORK, INFO);
}

// Double precision.
inline void Lapack::xgecon(char const &NORM, int N, const double A[], int LDA,
	double ANORM, double *RCOND, double WORK[], int IWORK[], int *INFO) {

    F77CALL(dgecon)(NORM, N, A, LDA, ANORM, RCOND, WORK, IWORK, INFO);
}

/********************************************************************************
*  SLANGE  returns the value of the one norm,  or the Frobenius norm, or
*  the  infinity norm,  or the  element of  largest absolute value  of a
*  real matrix A.
*
*******************************************************************************/// Single precision.
inline float Lapack::xlange(char const &NORM, int M, int N,
    const float A[], int LDA, float WORK[]) {

    return F77CALL(slange)(NORM, M, N, A, LDA, WORK);
}

// Double precision.
inline double Lapack::xlange(char const &NORM, int M, int N,
    const double A[], int LDA, double WORK[]) {

    return F77CALL(dlange)(NORM, M, N, A, LDA, WORK);
}

/********************************************************************************
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
*******************************************************************************/inline void Lapack::xgelsx(int M, int N, int NRHS,
    float A[], int LDA, float B[], int LDB, int JPVT[],
    float RCOND, int *RANK, float WORK[], int *INFO) {

    F77CALL(sgelsx)(M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK, WORK, INFO);
}

// Double precision.
inline void Lapack::xgelsx(int M, int N, int NRHS,
    double A[], int LDA, double B[], int LDB, int JPVT[],
    double RCOND, int *RANK, double WORK[], int *INFO) {

    F77CALL(dgelsx)(M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK, WORK, INFO);
}

/********************************************************************************
*  SGESVX uses the LU factorization to compute the solution to a real
*  system of linear equations
*     A * X = B,
*  where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
*
*  Error bounds on the solution and a condition estimate are also
*  provided.
*
*******************************************************************************/inline void Lapack::xgesvx(char const &FACT, char const &TRANS,
    int N, int NRHS, float A[], int LDA,
    float AF[], int LDAF, int IPIV[], char &EQUED, float R[],
    float C[], float B[], int LDB, float X[], int LDX,
    float *RCOND, float FERR[], float BERR[], float WORK[],
    int IWORK[], int *INFO) {

    F77CALL(sgesvx)(FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV, EQUED, R,
    C, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, IWORK, INFO);
}

// Double precision.
   inline void Lapack::xgesvx(char const &FACT, char const &TRANS,
    int N, int NRHS, double A[], int LDA,
    double AF[], int LDAF, int IPIV[], char &EQUED, double R[],
    double C[], double B[], int LDB, double X[], int LDX,
    double *RCOND, double FERR[], double BERR[], double WORK[],
    int IWORK[], int *INFO) {

    F77CALL(dgesvx)(FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV, EQUED, R,
    C, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, IWORK, INFO);
}

/********************************************************************************
*  SGEEQU computes row and column scalings intended to equilibrate an
*  M-by-N matrix A and reduce its condition number.  R returns the row
*  scale factors and C the column scale factors, chosen to try to make
*  the largest entry in each row and column of the matrix B with
*  elements B(i,j)=R(i)*A(i,j)*C(j) have absolute value 1.
*
*  R(i) and C(j) are restricted to be between SMLNUM = smallest safe
*  number and BIGNUM = largest safe number.  Use of these scaling
*  factors is not guaranteed to reduce the condition number of A but
*  works well in practice.
*
*******************************************************************************/// Single precision.
inline void Lapack::xgeequ(int M, int N, const float A[],
    int LDA, float R[], float C[], float *ROWCND,
    float *COLCND, float *AMAX, int *INFO) {

    F77CALL(sgeequ)(M, N, A, LDA, R, C, ROWCND, COLCND, AMAX, INFO);
}

// Double precision.
inline void Lapack::xgeequ(int M, int N, const double A[],
    int LDA, double R[], double C[], double *ROWCND,
    double *COLCND, double *AMAX, int *INFO) {

    F77CALL(dgeequ)(M, N, A, LDA, R, C, ROWCND, COLCND, AMAX, INFO);
}

/********************************************************************************
*  SLAQGE equilibrates a general M by N matrix A using the row and
*  scaling factors in the vectors R and C.
*
*******************************************************************************/// Single precision.
inline void Lapack::xlaqge(int M, int N, const float A[], int LDA,
    const float R[], const float C[], float ROWCND, float COLCND,
    float AMAX, char *EQUED) {

    F77CALL(slaqge)(M, N, A, LDA, R, C, ROWCND, COLCND, AMAX, EQUED);
}

// Double precision.
inline void Lapack::xlaqge(int M, int N, const double A[],
    int LDA, const double R[], const double C[],
    double ROWCND, double COLCND, double AMAX, char *EQUED) {

    F77CALL(dlaqge)(M, N, A, LDA, R, C, ROWCND, COLCND, AMAX, EQUED);
}

#endif	// LAPACK_HDR
