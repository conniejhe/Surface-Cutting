#ifndef ADL_GLOBALS_HDR
#define ADL_GLOBALS_HDR
/*****************************************************************************
*Surgical Navigation Technologies Inc Proprietary.
*Do not reproduce without permission in writing.
*Copyright (c) 1995  Surgical Navigation Technologies Inc.
*All rights reserved
*****************************************************************************/
// NOMANPAGE

// KWD : for BrainWorks
#include <OS.h>

#ifdef GNU_COMPILER
#define ADL_LITTLE_ENDIAN
#endif

// Use Big Endian By Default.
// --------------------------
#if !defined(ADL_LITTLE_ENDIAN) && !defined(ADL_BIG_ENDIAN)
#define ADL_BIG_ENDIAN
#endif

// Use Posix Serial Ports.
// -----------------------
#ifndef POSIX_SERIALPORT
#define POSIX_SERIALPORT
#endif

// Check For Proper Dimensions Of Arrays, Vectors, Matrices, etc.
// On By Default.
// --------------------------------------------------------------
#ifndef ADL_DEBUG_DIMENSIONS
#define ADL_DEBUG_DIMENSIONS 1
#endif

// Define The Floating Point Type
// Used By Various Non-Template Classes
// ------------------------------------
typedef float Coord_t;

#endif // ADL_GLOBALS_HDR
