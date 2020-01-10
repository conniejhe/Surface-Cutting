#ifndef LAINV_HDR
#define LAINV_HDR
/*******************************************************************************
This header file contains the C/C++ prototypes for the LAPACK routines.
Routines in this file are specific for inverse of a matrix.

Abed M. Hammoud
17Nov94-Stealth
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
      INTEGER            INFO, LDA, LWORK, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL               A( LDA, * ), WORK( LWORK )
*     ..
*
*  Purpose
*  =======
*
*  SGETRI computes the inverse of a matrix using the LU factorization
*  computed by SGETRF.
*
*  This method inverts U and then computes inv(A) by solving the system
*  inv(A)*L = inv(U) for inv(A).
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) REAL array, dimension (LDA,N)
*          On entry, the factors L and U from the factorization
*          A = P*L*U as computed by SGETRF.
*          On exit, if INFO = 0, the inverse of the original matrix A.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  IPIV    (input) INTEGER array, dimension (N)
*          The pivot indices from SGETRF; for 1<=i<=N, row i of the
*          matrix was interchanged with row IPIV(i).
*
*  WORK    (workspace) REAL array, dimension (LWORK)
*          On exit, if INFO=0, then WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.  LWORK >= max(1,N).
*          For optimal performance LWORK >= N*NB, where NB is
*          the optimal blocksize returned by ILAENV.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, U(i,i) is exactly zero; the matrix is
*                singular and its inverse could not be computed.
*
*  =====================================================================
*/
#ifdef __cplusplus

extern "C" {

// Single precision.
   extern void F77CALL(sgetri)(const int& N, float A[], 
	const int& LDA, const int IPIV[], float WORK[], 
	const int& LWORK, int *INFO); 

// Double precision.
   extern void F77CALL(dgetri)(const int& N, double A[], 
	const int& LDA, const int IPIV[], double WORK[], 
	const int& LWORK, int *INFO); 
}

#else

/* Single precision */
   void F77CALL(sgetri)(int const *N, float A[], int const *LDA, 
	const int IPIV[], float WORK[], int const *LWORK, int *INFO); 

/* Double precision */
   void F77CALL(dgetri)(int const *N, double A[], int const *LDA, 
	const int IPIV[], double WORK[], int const *LWORK, int *INFO); 

#endif

#endif

