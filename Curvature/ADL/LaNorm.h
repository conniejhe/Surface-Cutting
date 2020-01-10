#ifndef LANORM_HDR
#define LANORM_HDR
/*******************************************************************************
This header file contains the C/C++ prototypes for the LAPACK routines.
Routines in this file are specific for Norms of a matrix.

Abed M. Hammoud
17Nov94-Stealth
Copyright 1993-1995
*******************************************************************************/
// NOMANPAGE
/*
*  -- LAPACK auxiliary routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          NORM
      INTEGER            LDA, M, N
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  SLANGE  returns the value of the one norm,  or the Frobenius norm, or
*  the  infinity norm,  or the  element of  largest absolute value  of a
*  real matrix A.
*
*  Description
*  ===========
*
*  SLANGE returns the value
*
*     SLANGE = ( max(abs(A(i,j))), NORM = 'M' or 'm'
*              (
*              ( norm1(A),         NORM = '1', 'O' or 'o'
*              (
*              ( normI(A),         NORM = 'I' or 'i'
*              (
*              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
*
*  where  norm1  denotes the  one norm of a matrix (maximum column sum),
*  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
*  normF  denotes the  Frobenius norm of a matrix (square root of sum of
*  squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm.
*
*  Arguments
*  =========
*
*  NORM    (input) CHARACTER*1
*          Specifies the value to be returned in SLANGE as described
*          above.
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.  When M = 0,
*          SLANGE is set to zero.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.  When N = 0,
*          SLANGE is set to zero.
*
*  A       (input) REAL array, dimension (LDA,N)
*          The m by n matrix A.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(M,1).
*
*  WORK    (workspace) REAL array, dimension (LWORK),
*          where LWORK >= M when NORM = 'I'; otherwise, WORK is not
*          referenced.
*
* =====================================================================
*/
#ifdef __cplusplus

extern "C" {

// Single precision.
     extern float F77CALL(slange)(const char& NORM, const int& M, const int& N, 
  	const float A[], const int& LDA, float WORK[]); 

// Double precision.
     extern double F77CALL(dlange)(const char& NORM, const int& M, const int& N, 
  	const double A[], const int& LDA, double WORK[]); 

}

#else

/* Single precision */
   extern float F77CALL(slange)(char const *NORM, int const *M, int const *N, 
	const float A[], int const *LDA, float WORK[]); 

/* Double precision */
   extern double F77CALL(dlange)(char const *NORM, int const *M, int const *N, 
	const double A[], int const *LDA, double WORK[]); 

#endif

#endif

