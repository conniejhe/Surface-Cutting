#ifndef LACON_HDR
#define LACON_HDR
/*******************************************************************************
This header file contains the C/C++ prototypes for the LAPACK routines.
Routines in this file are specific for RCOND of a matrix.

Abed M. Hammoud
17Nov94-Stealth
Copyright 1993-1995
*******************************************************************************/
// NOMANPAGE
/*
*  -- LAPACK routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          NORM
      INTEGER            INFO, LDA, N
      REAL               ANORM, RCOND
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * )
      REAL               A( LDA, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  SGECON estimates the reciprocal of the condition number of a general
*  real matrix A, in either the 1-norm or the infinity-norm, using
*  the LU factorization computed by SGETRF.
*
*  An estimate is obtained for norm(inv(A)), and the reciprocal of the
*  condition number is computed as
*     RCOND = 1 / ( norm(A) * norm(inv(A)) ).
*
*  Arguments
*  =========
*
*  NORM    (input) CHARACTER*1
*          Specifies whether the 1-norm condition number or the
*          infinity-norm condition number is required:
*          = '1' or 'O':  1-norm;
*          = 'I':         Infinity-norm.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input) REAL array, dimension (LDA,N)
*          The factors L and U from the factorization A = P*L*U
*          as computed by SGETRF.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  ANORM   (input) REAL
*          If NORM = '1' or 'O', the 1-norm of the original matrix A.
*          If NORM = 'I', the infinity-norm of the original matrix A.
*
*  RCOND   (output) REAL
*          The reciprocal of the condition number of the matrix A,
*          computed as RCOND = 1/(norm(A) * norm(inv(A))).
*
*  WORK    (workspace) REAL array, dimension (4*N)
*
*  IWORK   (workspace) INTEGER array, dimension (N)
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*/
#ifdef __cplusplus

extern "C" {

// Single precision.
   extern void F77CALL(sgecon)(const char& NORM, const int& N, 
	const float A[], const int& LDA, const float& ANORM, 
	float *RCOND, float WORK[], int IWORK[], int *INFO); 

// Double precision.
   extern void F77CALL(dgecon)(const char& NORM, const int& N, 
	const double A[], const int& LDA, const double& ANORM, 
	double *RCOND, double WORK[], int IWORK[], int *INFO); 

#if 0
// Complex, single precision.
   extern void F77CALL(cgecon)(const char& NORM, const int& N, 
	const complex<float>  A[], const int& LDA, const float& ANORM, 
	float *RCOND, complex<float>  WORK[], float RWORK[], int *INFO); 

// Complex, double precision.
   extern void F77CALL(zgecon)(const char& NORM, const int& N, 
	const complex<double>  A[], const int& LDA, const double& ANORM, 
	double *RCOND, complex<double>  WORK[], double RWORK[], int *INFO); 
#endif
}

#else

/* Single precision */
   void F77CALL(sgecon)(char const *NORM, int const *N, const float A[], 
	int const *LDA, float const *ANORM, float *RCOND, 
	float WORK[], int IWORK[], int *INFO); 

/* Double precision */
   void F77CALL(dgecon)(char const *NORM, int const *N, const double A[], 
	int const *LDA, double const *ANORM, double *RCOND, 
	double WORK[], int IWORK[], int *INFO); 

#if 0
/* Single precision */
   void F77CALL(cgecon)(char const *NORM, int const *N, 
	const complex<float> A[], int const *LDA, float const *ANORM, 
	float *RCOND, complex<float>  WORK[], float RWORK[], int *INFO); 

/* Double precision */
   void F77CALL(zgecon)(char const *NORM, int const *N, 
	const complex<double> A[], int const *LDA, double const *ANORM, 
	double *RCOND, complex<double>  WORK[], double RWORK[], int *INFO); 
#endif
#endif

#endif

#ifndef LAEQUILIBRATE_HDR
#define LAEQUILIBRATE_HDR
/*******************************************************************************
This header file contains the C/C++ prototypes for the LAPACK routines.
Routines in this file are specific for matrix balancing.

Abed M. Hammoud
23Nov94-Stealth
Copyright 1993-1995
*******************************************************************************/
/*
*  -- LAPACK routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993 
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
      REAL               AMAX, COLCND, ROWCND
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), C( * ), R( * )
*     ..
*
*  Purpose
*  =======
*
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
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input) REAL array, dimension (LDA,N)
*          The M-by-N matrix whose equilibration factors are
*          to be computed.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  R       (output) REAL array, dimension (M)
*          If INFO = 0 or INFO > M, R contains the row scale factors
*          for A.
*
*  C       (output) REAL array, dimension (N)
*          If INFO = 0,  C contains the column scale factors for A.
*
*  ROWCND  (output) REAL
*          If INFO = 0 or INFO > M, ROWCND contains the ratio of the
*          smallest R(i) to the largest R(i).  If ROWCND >= 0.1 and
*          AMAX is neither too large nor too small, it is not worth
*          scaling by R.
*
*  COLCND  (output) REAL
*          If INFO = 0, COLCND contains the ratio of the smallest
*          C(i) to the largest C(i).  If COLCND >= 0.1, it is not
*          worth scaling by C.
*
*  AMAX    (output) REAL
*          Absolute value of largest matrix element.  If AMAX is very
*          close to overflow or very close to underflow, the matrix
*          should be scaled.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i,  and i is
*                <= M:  the i-th row of A is exactly zero
*                >  M:  the (i-M)-th column of A is exactly zero
*
*  =====================================================================
*/

#ifdef __cplusplus

extern "C" {

// Single precision.
   extern void F77CALL(sgeequ)(const int& M, const int& N, const float A[], 
	const int& LDA, float R[], float C[], float *ROWCND, 
	float *COLCND, float *AMAX, int *INFO); 
      
// Double precision.
   extern void F77CALL(dgeequ)(const int& M, const int& N, const double A[], 
	const int& LDA, double R[], double C[], double *ROWCND, 
	double *COLCND, double *AMAX, int *INFO); 

#if 0
// Complex, single precision.
   extern void F77CALL(cgeequ)(const int& M, const int& N, 
	const complex<float>  A[], const int& LDA, float R[], float C[], 
	float *ROWCND, float *COLCND, float *AMAX, int *INFO); 

// Complex, double precision.
   extern void F77CALL(zgeequ)(const int& M, const int& N, 
	const complex<double>  A[], const int& LDA, double R[], double C[], 
	double *ROWCND, double *COLCND, double *AMAX, int *INFO); 
#endif
}

#else

/* Single precision */
   extern void F77CALL(sgeequ)(int const *M, int const *N, const float A[], 
	int const *LDA, float R[], float C[], float *ROWCND, 
	float *COLCND, float *AMAX, int *INFO); 

/* Double precision */
   extern void F77CALL(dgeequ)(int const *M, int const *N, const double A[], 
	int const *LDA, double R[], double C[], double *ROWCND, 
	double *COLCND, double *AMAX, int *INFO); 

#if 0
/* Complex, single precision */
   extern void F77CALL(cgeequ)(int const *M, int const *N, 
	const complex<float>  A[], int const *LDA, float R[], float C[], 
	float *ROWCND, float *COLCND, float *AMAX, int *INFO); 

/* Complex, double precision */
   extern void F77CALL(zgeequ)(int const *M, int const *N, 
	const complex<double>  A[], int const *LDA, double R[], double C[], 
	double *ROWCND, double *COLCND, double *AMAX, int *INFO); 
#endif

#endif

/*
*  -- LAPACK auxiliary routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          EQUED
      INTEGER            LDA, M, N
      REAL               AMAX, COLCND, ROWCND
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), C( * ), R( * )
*     ..
*
*  Purpose
*  =======
*
*  SLAQGE equilibrates a general M by N matrix A using the row and
*  scaling factors in the vectors R and C.
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
*          On entry, the M by N matrix A.
*          On exit, the equilibrated matrix.  See EQUED for the form of
*          the equilibrated matrix.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(M,1).
*
*  R       (input) REAL array, dimension (M)
*          The row scale factors for A.
*
*  C       (input) REAL array, dimension (N)
*          The column scale factors for A.
*
*  ROWCND  (input) REAL
*          Ratio of the smallest R(i) to the largest R(i).
*
*  COLCND  (input) REAL
*          Ratio of the smallest C(i) to the largest C(i).
*
*  AMAX    (input) REAL
*          Absolute value of largest matrix entry.
*
*  EQUED   (output) CHARACTER*1
*          Specifies the form of equilibration that was done.
*          = 'N':  No equilibration
*          = 'R':  Row equilibration, i.e., A has been premultiplied by
*                  diag(R).
*          = 'C':  Column equilibration, i.e., A has been postmultiplied
*                  by diag(C).
*          = 'B':  Both row and column equilibration, i.e., A has been
*                  replaced by diag(R) * A * diag(C).
*
*  Internal Parameters
*  ===================
*
*  THRESH is a threshold value used to decide if row or column scaling
*  should be done based on the ratio of the row or column scaling
*  factors.  If ROWCND < THRESH, row scaling is done, and if
*  COLCND < THRESH, column scaling is done.
*
*  LARGE and SMALL are threshold values used to decide if row scaling
*  should be done based on the absolute size of the largest matrix
*  element.  If AMAX > LARGE or AMAX < SMALL, row scaling is done.
*
*  =====================================================================
*/

#ifdef __cplusplus

extern "C" {

// Single precision.
   extern void F77CALL(slaqge)(const int& M, const int& N, const float A[], 
	const int& LDA, const float R[], const float C[], 
	const float& ROWCND, const float& COLCND, 
        const float& AMAX, char *EQUED); 
      
// Double precision.
   extern void F77CALL(dlaqge)(const int& M, const int& N, const double A[], 
	const int& LDA, const double R[], const double C[], 
	const double& ROWCND, const double& COLCND, 
        const double& AMAX, char *EQUED); 

#if 0
// Complex, single precision.
   extern void F77CALL(claqge)(const int& M, const int& N, 
	const complex<float>  A[], const int& LDA, const float R[], 
	const float C[], const float& ROWCND, const float& COLCND, 
        const float& AMAX, char *EQUED); 

// Complex, double precision.
   extern void F77CALL(zlaqge)(const int& M, const int& N, 
	const complex<double>  A[], const int& LDA, const double R[], 
	const double C[], const double& ROWCND, const double& COLCND, 
        const double& AMAX, char *EQUED); 
#endif
}

#else

/* Single precision */
   extern void F77CALL(slaqge)(int const *M, int const *N, const float A[], 
	int const *LDA, const float R[], const float C[], 
	float const *ROWCND, float const *COLCND, 
        float const *AMAX, char *EQUED); 

/* Double precision */
   extern void F77CALL(dlaqge)(int const *M, int const *N, const double A[], 
	int const *LDA, const double R[], const double C[], 
	double const *ROWCND, double const *COLCND, 
        double const *AMAX, char *EQUED); 


#if 0
/* Complex, single precision */
   extern void F77CALL(claqge)(int const *M, int const *N, 
	const complex<float>  A[], int const *LDA, const float R[], 
	const float C[], float const *ROWCND, float const *COLCND, 
        float const *AMAX, char *EQUED); 

/* Complex, double precision */
   extern void F77CALL(zlaqge)(int const *M, int const *N, 
	const complex<double>  A[], int const *LDA, const double R[], 
	const double C[], double const *ROWCND, double const *COLCND, 
        double const *AMAX, char *EQUED); 
#endif

#endif

#endif
