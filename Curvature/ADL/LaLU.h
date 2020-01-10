#ifndef LALU_HDR
#define LALU_HDR
/*******************************************************************************
This header file contains the C/C++ prototypes for the LAPACK routines.
Routines in this file are specific for LU factorization.

Abed M. Hammoud
15Nov94-Stealth
Copyright 1993-1995
*******************************************************************************/
// NOMANPAGE
/*
*  -- LAPACK routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993 
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL               A( LDA, * )
*     ..
*
*  Purpose
*  =======
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
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) REAL array, dimension (LDA,N)
*          On entry, the M-by-N matrix to be factored.
*          On exit, the factors L and U from the factorization
*          A = P*L*U; the unit diagonal elements of L are not stored.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  IPIV    (output) INTEGER array, dimension (min(M,N))
*          The pivot indices; for 1 <= i <= min(M,N), row i of the
*          matrix was interchanged with row IPIV(i).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
*                has been completed, but the factor U is exactly
*                singular, and division by zero will occur if it is used
*                to solve a system of equations.
*
*  =====================================================================
*/
#ifdef __cplusplus

extern "C" {

// Single precision.
   extern void F77CALL(sgetrf)(const int& M, const int& N, float A[], 
	const int& LDA, int IPIV[], int *INFO); 

// Double precision.
   extern void F77CALL(dgetrf)(const int& M, const int& N, double A[], 
	const int& LDA, int IPIV[], int *INFO); 

}

#else

/* Single precision */
   void F77CALL(sgetrf)(int const *M, int const *N, float A[], int const *LDA, 
	int IPIV[], int *INFO); 

/* Double precision */
   void F77CALL(dgetrf)(int const *M, int const *N, double A[], int const *LDA, 
	int IPIV[], int *INFO); 

#endif

#endif

