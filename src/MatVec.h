/*************************************************************************\

  Copyright 1999 The University of North Carolina at Chapel Hill.
  All Rights Reserved.

  Permission to use, copy, modify and distribute this software and its
  documentation for educational, research and non-profit purposes, without
  fee, and without a written agreement is hereby granted, provided that the
  above copyright notice and the following three paragraphs appear in all
  copies.

  IN NO EVENT SHALL THE UNIVERSITY OF NORTH CAROLINA AT CHAPEL HILL BE
  LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR
  CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE
  USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF THE UNIVERSITY
  OF NORTH CAROLINA HAVE BEEN ADVISED OF THE POSSIBILITY OF SUCH
  DAMAGES.

  THE UNIVERSITY OF NORTH CAROLINA SPECIFICALLY DISCLAIM ANY
  WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.  THE SOFTWARE
  PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF
  NORTH CAROLINA HAS NO OBLIGATIONS TO PROVIDE MAINTENANCE, SUPPORT,
  UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

  The authors may be contacted via:

  US Mail:             S. Gottschalk
                       Department of Computer Science
                       Sitterson Hall, CB #3175
                       University of N. Carolina
                       Chapel Hill, NC 27599-3175

  Phone:               (919)962-1749

  EMail:               geom@cs.unc.edu


\**************************************************************************/

#ifndef CKL_MATVEC_H
#define CKL_MATVEC_H

#include <cmath>
#include <iostream>
#include "CKL_Compile.h"

namespace CKL
{

#ifndef M_PI
static const CKL_REAL M_PI = (CKL_REAL)3.14159265359;
#endif


inline std::ostream& operator << (std::ostream& os, const CKL_REAL M[3][3])
{
  os << M[0][0] << " " << M[0][1] << " " << M[0][2]
     << M[1][0] << " " << M[1][1] << " " << M[1][2]
     << M[2][0] << " " << M[2][1] << " " << M[2][2];

  return os;
}

inline std::ostream& operator << (std::ostream& os, const CKL_REAL V[3])
{
  os << V[0] << " " << V[1] << " " << V[2];
  
  return os;
}

inline void Midentity(CKL_REAL M[3][3])
{
  M[0][0] = M[1][1] = M[2][2] = 1.0;
  M[0][1] = M[1][2] = M[2][0] = 0.0;
  M[0][2] = M[1][0] = M[2][1] = 0.0;
}

inline void Videntity(CKL_REAL T[3])
{
  T[0] = T[1] = T[2] = 0.0;
}

inline void McM(CKL_REAL Mr[3][3], const CKL_REAL M[3][3])
{
  Mr[0][0] = M[0][0];
  Mr[0][1] = M[0][1];
  Mr[0][2] = M[0][2];
  Mr[1][0] = M[1][0];
  Mr[1][1] = M[1][1];
  Mr[1][2] = M[1][2];
  Mr[2][0] = M[2][0];
  Mr[2][1] = M[2][1];
  Mr[2][2] = M[2][2];
}

inline void MTcM(CKL_REAL Mr[3][3], const CKL_REAL M[3][3])
{
  Mr[0][0] = M[0][0];
  Mr[1][0] = M[0][1];
  Mr[2][0] = M[0][2];
  Mr[0][1] = M[1][0];
  Mr[1][1] = M[1][1];
  Mr[2][1] = M[1][2];
  Mr[0][2] = M[2][0];
  Mr[1][2] = M[2][1];
  Mr[2][2] = M[2][2];
}

inline void VcV(CKL_REAL Vr[3], const CKL_REAL V[3])
{
  Vr[0] = V[0];
  Vr[1] = V[1];
  Vr[2] = V[2];
}

inline void McolcV(CKL_REAL Vr[3], const CKL_REAL M[3][3], int c)
{
  Vr[0] = M[0][c];
  Vr[1] = M[1][c];
  Vr[2] = M[2][c];
}

inline void VscMrow(CKL_REAL Mr[3][3], const CKL_REAL V1[3], const CKL_REAL V2[3], const CKL_REAL V3[3])
{
  Mr[0][0] = V1[0]; Mr[0][1] = V1[1]; Mr[0][2] = V1[2];
  Mr[1][0] = V2[0]; Mr[1][1] = V2[1]; Mr[1][2] = V2[2];
  Mr[2][0] = V3[0]; Mr[2][1] = V3[1]; Mr[2][2] = V3[2];
}

inline void McolcMcol(CKL_REAL Mr[3][3], int cr, const CKL_REAL M[3][3], int c)
{
  Mr[0][cr] = M[0][c];
  Mr[1][cr] = M[1][c];
  Mr[2][cr] = M[2][c];
}

inline void MxMpV(CKL_REAL Mr[3][3], const CKL_REAL M1[3][3], const CKL_REAL M2[3][3], const CKL_REAL T[3])
{
  Mr[0][0] = (M1[0][0] * M2[0][0] +
              M1[0][1] * M2[1][0] +
              M1[0][2] * M2[2][0] +
              T[0]);
  Mr[1][0] = (M1[1][0] * M2[0][0] +
              M1[1][1] * M2[1][0] +
              M1[1][2] * M2[2][0] +
              T[1]);
  Mr[2][0] = (M1[2][0] * M2[0][0] +
              M1[2][1] * M2[1][0] +
              M1[2][2] * M2[2][0] +
              T[2]);
  Mr[0][1] = (M1[0][0] * M2[0][1] +
              M1[0][1] * M2[1][1] +
              M1[0][2] * M2[2][1] +
              T[0]);
  Mr[1][1] = (M1[1][0] * M2[0][1] +
              M1[1][1] * M2[1][1] +
              M1[1][2] * M2[2][1] +
              T[1]);
  Mr[2][1] = (M1[2][0] * M2[0][1] +
              M1[2][1] * M2[1][1] +
              M1[2][2] * M2[2][1] +
              T[2]);
  Mr[0][2] = (M1[0][0] * M2[0][2] +
              M1[0][1] * M2[1][2] +
              M1[0][2] * M2[2][2] +
              T[0]);
  Mr[1][2] = (M1[1][0] * M2[0][2] +
              M1[1][1] * M2[1][2] +
              M1[1][2] * M2[2][2] +
              T[1]);
  Mr[2][2] = (M1[2][0] * M2[0][2] +
              M1[2][1] * M2[1][2] +
              M1[2][2] * M2[2][2] +
              T[2]);
}

inline void MxM(CKL_REAL Mr[3][3], const CKL_REAL M1[3][3], const CKL_REAL M2[3][3])
{
  Mr[0][0] = (M1[0][0] * M2[0][0] +
              M1[0][1] * M2[1][0] +
              M1[0][2] * M2[2][0]);
  Mr[1][0] = (M1[1][0] * M2[0][0] +
              M1[1][1] * M2[1][0] +
              M1[1][2] * M2[2][0]);
  Mr[2][0] = (M1[2][0] * M2[0][0] +
              M1[2][1] * M2[1][0] +
              M1[2][2] * M2[2][0]);
  Mr[0][1] = (M1[0][0] * M2[0][1] +
              M1[0][1] * M2[1][1] +
              M1[0][2] * M2[2][1]);
  Mr[1][1] = (M1[1][0] * M2[0][1] +
              M1[1][1] * M2[1][1] +
              M1[1][2] * M2[2][1]);
  Mr[2][1] = (M1[2][0] * M2[0][1] +
              M1[2][1] * M2[1][1] +
              M1[2][2] * M2[2][1]);
  Mr[0][2] = (M1[0][0] * M2[0][2] +
              M1[0][1] * M2[1][2] +
              M1[0][2] * M2[2][2]);
  Mr[1][2] = (M1[1][0] * M2[0][2] +
              M1[1][1] * M2[1][2] +
              M1[1][2] * M2[2][2]);
  Mr[2][2] = (M1[2][0] * M2[0][2] +
              M1[2][1] * M2[1][2] +
              M1[2][2] * M2[2][2]);
}


inline void MxMT(CKL_REAL Mr[3][3], const CKL_REAL M1[3][3], const CKL_REAL M2[3][3])
{
  Mr[0][0] = (M1[0][0] * M2[0][0] +
              M1[0][1] * M2[0][1] +
              M1[0][2] * M2[0][2]);
  Mr[1][0] = (M1[1][0] * M2[0][0] +
              M1[1][1] * M2[0][1] +
              M1[1][2] * M2[0][2]);
  Mr[2][0] = (M1[2][0] * M2[0][0] +
              M1[2][1] * M2[0][1] +
              M1[2][2] * M2[0][2]);
  Mr[0][1] = (M1[0][0] * M2[1][0] +
              M1[0][1] * M2[1][1] +
              M1[0][2] * M2[1][2]);
  Mr[1][1] = (M1[1][0] * M2[1][0] +
              M1[1][1] * M2[1][1] +
              M1[1][2] * M2[1][2]);
  Mr[2][1] = (M1[2][0] * M2[1][0] +
              M1[2][1] * M2[1][1] +
              M1[2][2] * M2[1][2]);
  Mr[0][2] = (M1[0][0] * M2[2][0] +
              M1[0][1] * M2[2][1] +
              M1[0][2] * M2[2][2]);
  Mr[1][2] = (M1[1][0] * M2[2][0] +
              M1[1][1] * M2[2][1] +
              M1[1][2] * M2[2][2]);
  Mr[2][2] = (M1[2][0] * M2[2][0] +
              M1[2][1] * M2[2][1] +
              M1[2][2] * M2[2][2]);
}

inline void MTxM(CKL_REAL Mr[3][3], const CKL_REAL M1[3][3], const CKL_REAL M2[3][3])
{
  Mr[0][0] = (M1[0][0] * M2[0][0] +
              M1[1][0] * M2[1][0] +
              M1[2][0] * M2[2][0]);
  Mr[1][0] = (M1[0][1] * M2[0][0] +
              M1[1][1] * M2[1][0] +
              M1[2][1] * M2[2][0]);
  Mr[2][0] = (M1[0][2] * M2[0][0] +
              M1[1][2] * M2[1][0] +
              M1[2][2] * M2[2][0]);
  Mr[0][1] = (M1[0][0] * M2[0][1] +
              M1[1][0] * M2[1][1] +
              M1[2][0] * M2[2][1]);
  Mr[1][1] = (M1[0][1] * M2[0][1] +
              M1[1][1] * M2[1][1] +
              M1[2][1] * M2[2][1]);
  Mr[2][1] = (M1[0][2] * M2[0][1] +
              M1[1][2] * M2[1][1] +
              M1[2][2] * M2[2][1]);
  Mr[0][2] = (M1[0][0] * M2[0][2] +
              M1[1][0] * M2[1][2] +
              M1[2][0] * M2[2][2]);
  Mr[1][2] = (M1[0][1] * M2[0][2] +
              M1[1][1] * M2[1][2] +
              M1[2][1] * M2[2][2]);
  Mr[2][2] = (M1[0][2] * M2[0][2] +
              M1[1][2] * M2[1][2] +
              M1[2][2] * M2[2][2]);
}

inline void MxV(CKL_REAL Vr[3], const CKL_REAL M1[3][3], const CKL_REAL V1[3])
{
  Vr[0] = (M1[0][0] * V1[0] +
           M1[0][1] * V1[1] +
           M1[0][2] * V1[2]);
  Vr[1] = (M1[1][0] * V1[0] +
           M1[1][1] * V1[1] +
           M1[1][2] * V1[2]);
  Vr[2] = (M1[2][0] * V1[0] +
           M1[2][1] * V1[1] +
           M1[2][2] * V1[2]);
}


inline void MxVpV(CKL_REAL Vr[3], const CKL_REAL M1[3][3], const CKL_REAL V1[3], const CKL_REAL V2[3])
{
  Vr[0] = (M1[0][0] * V1[0] +
           M1[0][1] * V1[1] +
           M1[0][2] * V1[2] +
           V2[0]);
  Vr[1] = (M1[1][0] * V1[0] +
           M1[1][1] * V1[1] +
           M1[1][2] * V1[2] +
           V2[1]);
  Vr[2] = (M1[2][0] * V1[0] +
           M1[2][1] * V1[1] +
           M1[2][2] * V1[2] +
           V2[2]);
}


inline void sMxVpV(CKL_REAL Vr[3], CKL_REAL s1, const CKL_REAL M1[3][3], const CKL_REAL V1[3], const CKL_REAL V2[3])
{
  Vr[0] = s1 * (M1[0][0] * V1[0] +
                M1[0][1] * V1[1] +
                M1[0][2] * V1[2]) +
    V2[0];
  Vr[1] = s1 * (M1[1][0] * V1[0] +
                M1[1][1] * V1[1] +
                M1[1][2] * V1[2]) +
    V2[1];
  Vr[2] = s1 * (M1[2][0] * V1[0] +
                M1[2][1] * V1[1] +
                M1[2][2] * V1[2]) +
    V2[2];
}

inline void MTxV(CKL_REAL Vr[3], const CKL_REAL M1[3][3], const CKL_REAL V1[3])
{
  Vr[0] = (M1[0][0] * V1[0] +
           M1[1][0] * V1[1] +
           M1[2][0] * V1[2]);
  Vr[1] = (M1[0][1] * V1[0] +
           M1[1][1] * V1[1] +
           M1[2][1] * V1[2]);
  Vr[2] = (M1[0][2] * V1[0] +
           M1[1][2] * V1[1] +
           M1[2][2] * V1[2]);
}

inline void sMTxV(CKL_REAL Vr[3], CKL_REAL s1, const CKL_REAL M1[3][3], const CKL_REAL V1[3])
{
  Vr[0] = s1 * (M1[0][0] * V1[0] +
                M1[1][0] * V1[1] +
                M1[2][0] * V1[2]);
  Vr[1] = s1 * (M1[0][1] * V1[0] +
                M1[1][1] * V1[1] +
                M1[2][1] * V1[2]);
  Vr[2] = s1 * (M1[0][2] * V1[0] +
                M1[1][2] * V1[1] +
                M1[2][2] * V1[2]);
}

inline void sMxV(CKL_REAL Vr[3], CKL_REAL s1, const CKL_REAL M1[3][3], const CKL_REAL V1[3])
{
  Vr[0] = s1 * (M1[0][0] * V1[0] +
                M1[0][1] * V1[1] +
                M1[0][2] * V1[2]);
  Vr[1] = s1 * (M1[1][0] * V1[0] +
                M1[1][1] * V1[1] +
                M1[1][2] * V1[2]);
  Vr[2] = s1 * (M1[2][0] * V1[0] +
                M1[2][1] * V1[1] +
                M1[2][2] * V1[2]);
}


inline void VmV(CKL_REAL Vr[3], const CKL_REAL V1[3], const CKL_REAL V2[3])
{
  Vr[0] = V1[0] - V2[0];
  Vr[1] = V1[1] - V2[1];
  Vr[2] = V1[2] - V2[2];
}

inline void VpV(CKL_REAL Vr[3], const CKL_REAL V1[3], const CKL_REAL V2[3])
{
  Vr[0] = V1[0] + V2[0];
  Vr[1] = V1[1] + V2[1];
  Vr[2] = V1[2] + V2[2];
}

inline void VpVxS(CKL_REAL Vr[3], const CKL_REAL V1[3], const CKL_REAL V2[3], CKL_REAL s)
{
  Vr[0] = V1[0] + V2[0] * s;
  Vr[1] = V1[1] + V2[1] * s;
  Vr[2] = V1[2] + V2[2] * s;
}

inline void MskewV(CKL_REAL M[3][3], const CKL_REAL v[3])
{
  M[0][0] = M[1][1] = M[2][2] = 0.0;
  M[1][0] = v[2];
  M[0][1] = -v[2];
  M[0][2] = v[1];
  M[2][0] = -v[1];
  M[1][2] = -v[0];
  M[2][1] = v[0];
}


inline void VcrossV(CKL_REAL Vr[3], const CKL_REAL V1[3], const CKL_REAL V2[3])
{
  Vr[0] = V1[1] * V2[2] - V1[2] * V2[1];
  Vr[1] = V1[2] * V2[0] - V1[0] * V2[2];
  Vr[2] = V1[0] * V2[1] - V1[1] * V2[0];
}

inline CKL_REAL Vlength(CKL_REAL V[3])
{
  return sqrt(V[0] * V[0] + V[1] * V[1] + V[2] * V[2]);
}

inline void Vnormalize(CKL_REAL V[3])
{
  CKL_REAL d = (CKL_REAL)1.0 / sqrt(V[0] * V[0] + V[1] * V[1] + V[2] * V[2]);
  V[0] *= d;
  V[1] *= d;
  V[2] *= d;
}

inline CKL_REAL VdotV(const CKL_REAL V1[3], const CKL_REAL V2[3])
{
  return (V1[0] * V2[0] + V1[1] * V2[1] + V1[2] * V2[2]);
}

inline CKL_REAL VdistV2(const CKL_REAL V1[3], const CKL_REAL V2[3])
{
  return ((V1[0] - V2[0]) * (V1[0] - V2[0]) +
          (V1[1] - V2[1]) * (V1[1] - V2[1]) +
          (V1[2] - V2[2]) * (V1[2] - V2[2]));
}

inline void VxS(CKL_REAL Vr[3], const CKL_REAL V[3], CKL_REAL s)
{
  Vr[0] = V[0] * s;
  Vr[1] = V[1] * s;
  Vr[2] = V[2] * s;
}

inline void MxS(CKL_REAL Mr[3][3], const CKL_REAL M[3][3], CKL_REAL s)
{
  Mr[0][0] = M[0][0] * s; Mr[0][1] = M[0][1] * s; Mr[0][2] = M[0][2] * s;
  Mr[1][0] = M[1][0] * s; Mr[1][1] = M[1][1] * s; Mr[1][2] = M[1][2] * s;
  Mr[2][0] = M[2][0] * s; Mr[2][1] = M[2][1] * s; Mr[2][2] = M[2][2] * s;
}

inline void MmM(CKL_REAL Mr[3][3], const CKL_REAL M1[3][3], const CKL_REAL M2[3][3])
{
  Mr[0][0] = M1[0][0] - M2[0][0];   Mr[0][1] = M1[0][1] - M2[0][1];   Mr[0][2] = M1[0][2] - M2[0][2];
  Mr[1][0] = M1[1][0] - M2[1][0];   Mr[1][1] = M1[1][1] - M2[1][1];   Mr[1][2] = M1[1][2] - M2[1][2];
  Mr[2][0] = M1[2][0] - M2[2][0];   Mr[2][1] = M1[2][1] - M2[2][1];   Mr[2][2] = M1[2][2] - M2[2][2];
}

inline void MpM(CKL_REAL Mr[3][3], const CKL_REAL M1[3][3], const CKL_REAL M2[3][3])
{
  Mr[0][0] = M1[0][0] + M2[0][0];   Mr[0][1] = M1[0][1] + M2[0][1];   Mr[0][2] = M1[0][2] + M2[0][2];
  Mr[1][0] = M1[1][0] + M2[1][0];   Mr[1][1] = M1[1][1] + M2[1][1];   Mr[1][2] = M1[1][2] + M2[1][2];
  Mr[2][0] = M1[2][0] + M2[2][0];   Mr[2][1] = M1[2][1] + M2[2][1];   Mr[2][2] = M1[2][2] + M2[2][2];
}


inline void MRotZ(CKL_REAL Mr[3][3], CKL_REAL t)
{
  Mr[0][0] = cos(t);
  Mr[1][0] = sin(t);
  Mr[0][1] = -Mr[1][0];
  Mr[1][1] = Mr[0][0];
  Mr[2][0] = Mr[2][1] = 0.0;
  Mr[0][2] = Mr[1][2] = 0.0;
  Mr[2][2] = 1.0;
}

inline void MRotX(CKL_REAL Mr[3][3], CKL_REAL t)
{
  Mr[1][1] = cos(t);
  Mr[2][1] = sin(t);
  Mr[1][2] = -Mr[2][1];
  Mr[2][2] = Mr[1][1];
  Mr[0][1] = Mr[0][2] = 0.0;
  Mr[1][0] = Mr[2][0] = 0.0;
  Mr[0][0] = 1.0;
}

inline void MRotY(CKL_REAL Mr[3][3], CKL_REAL t)
{
  Mr[2][2] = cos(t);
  Mr[0][2] = sin(t);
  Mr[2][0] = -Mr[0][2];
  Mr[0][0] = Mr[2][2];
  Mr[1][2] = Mr[1][0] = 0.0;
  Mr[2][1] = Mr[0][1] = 0.0;
  Mr[1][1] = 1.0;
}

inline void MVtoOGL(double oglm[16], const CKL_REAL R[3][3], const CKL_REAL T[3])
{
  oglm[0] = (double)R[0][0];
  oglm[1] = (double)R[1][0];
  oglm[2] = (double)R[2][0];
  oglm[3] = 0.0;
  oglm[4] = (double)R[0][1];
  oglm[5] = (double)R[1][1];
  oglm[6] = (double)R[2][1];
  oglm[7] = 0.0;
  oglm[8] = (double)R[0][2];
  oglm[9] = (double)R[1][2];
  oglm[10] = (double)R[2][2];
  oglm[11] = 0.0;
  oglm[12] = (double)T[0];
  oglm[13] = (double)T[1];
  oglm[14] = (double)T[2];
  oglm[15] = 1.0;
}

inline void OGLtoMV(CKL_REAL R[3][3], CKL_REAL T[3], const double oglm[16])
{
  R[0][0] = (CKL_REAL)oglm[0];
  R[1][0] = (CKL_REAL)oglm[1];
  R[2][0] = (CKL_REAL)oglm[2];
  
  R[0][1] = (CKL_REAL)oglm[4];
  R[1][1] = (CKL_REAL)oglm[5];
  R[2][1] = (CKL_REAL)oglm[6];
  
  R[0][2] = (CKL_REAL)oglm[8];
  R[1][2] = (CKL_REAL)oglm[9];
  R[2][2] = (CKL_REAL)oglm[10];
  
  T[0] = (CKL_REAL)oglm[12];
  T[1] = (CKL_REAL)oglm[13];
  T[2] = (CKL_REAL)oglm[14];
}

// taken from quatlib, written by Richard Holloway
static const int QX = 0;
static const int QY = 1;
static const int QZ = 2;
static const int QW = 3;

inline void MRotQ(CKL_REAL destMatrix[3][3], CKL_REAL srcQuat[4])
{
  CKL_REAL  s;
  CKL_REAL  xs, ys, zs,
    wx, wy, wz,
    xx, xy, xz,
    yy, yz, zz;
            
  /*
   * For unit srcQuat, just set s = 2.0; or set xs = srcQuat[QX] +
   *   srcQuat[QX], etc.
   */
  
  s = (CKL_REAL)2.0 / (srcQuat[QX] * srcQuat[QX] + srcQuat[QY] * srcQuat[QY] +
                       srcQuat[QZ] * srcQuat[QZ] + srcQuat[QW] * srcQuat[QW]);
                       
  xs = srcQuat[QX] * s;
  ys = srcQuat[QY] * s;
  zs = srcQuat[QZ] * s;
  wx = srcQuat[QW] * xs;
  wy = srcQuat[QW] * ys;
  wz = srcQuat[QW] * zs;
  xx = srcQuat[QX] * xs;
  xy = srcQuat[QX] * ys;
  xz = srcQuat[QX] * zs;
  yy = srcQuat[QY] * ys;
  yz = srcQuat[QY] * zs;
  zz = srcQuat[QZ] * zs;
  
  destMatrix[QX][QX] = (CKL_REAL)1.0 - (yy + zz);
  destMatrix[QX][QY] = xy + wz;
  destMatrix[QX][QZ] = xz - wy;
  
  destMatrix[QY][QX] = xy - wz;
  destMatrix[QY][QY] = (CKL_REAL)1.0 - (xx + zz);
  destMatrix[QY][QZ] = yz + wx;
  
  destMatrix[QZ][QX] = xz + wy;
  destMatrix[QZ][QY] = yz - wx;
  destMatrix[QZ][QZ] = (CKL_REAL)1.0 - (xx + yy);
}

inline void Mqinverse(CKL_REAL Mr[3][3], CKL_REAL m[3][3])
{
  int i, j;
  
  for(i = 0; i < 3; i++)
    for(j = 0; j < 3; j++)
    {
      int i1 = (i + 1) % 3;
      int i2 = (i + 2) % 3;
      int j1 = (j + 1) % 3;
      int j2 = (j + 2) % 3;
      Mr[i][j] = (m[j1][i1] * m[j2][i2] - m[j1][i2] * m[j2][i1]);
    }
}


#define ROTATE(a,i,j,k,l) g=a[i][j]; h=a[k][l]; a[i][j]=g-s*(h+g*tau); a[k][l]=h+s*(g-h*tau);

void inline Meigen(CKL_REAL vout[3][3], CKL_REAL dout[3], CKL_REAL a[3][3])
{
  int n = 3;
  int j, iq, ip, i;
  CKL_REAL tresh, theta, tau, t, sm, s, h, g, c;
  int nrot;
  CKL_REAL b[3];
  CKL_REAL z[3];
  CKL_REAL v[3][3];
  CKL_REAL d[3];

  Midentity(v);
  for(ip = 0; ip < n; ip++)
  {
    b[ip] = a[ip][ip];
    d[ip] = a[ip][ip];
    z[ip] = 0.0;
  }

  nrot = 0;

  for(i = 0; i < 50; i++)
  {

    sm = 0.0;
    for(ip = 0; ip < n; ip++) for(iq = ip + 1; iq < n; iq++) sm += fabs(a[ip][iq]);
    if(sm == 0.0)
    {
      McM(vout, v);
      VcV(dout, d);
      return;
    }


    if(i < 3) tresh = (CKL_REAL)0.2 * sm / (n * n);
    else tresh = 0.0;

    for(ip = 0; ip < n; ip++) for(iq = ip + 1; iq < n; iq++)
                              {
                                g = (CKL_REAL)100.0 * fabs(a[ip][iq]);
                                if(i > 3 &&
                                   fabs(d[ip]) + g == fabs(d[ip]) &&
                                   fabs(d[iq]) + g == fabs(d[iq]))
                                  a[ip][iq] = 0.0;
                                else if(fabs(a[ip][iq]) > tresh)
                                {
                                  h = d[iq] - d[ip];
                                  if(fabs(h) + g == fabs(h)) t = (a[ip][iq]) / h;
                                  else
                                  {
                                    theta = (CKL_REAL)0.5 * h / (a[ip][iq]);
                                    t = (CKL_REAL)(1.0 / (fabs(theta) + sqrt(1.0 + theta * theta)));
                                    if(theta < 0.0) t = -t;
                                  }
                                  c = (CKL_REAL)1.0 / sqrt(1 + t * t);
                                  s = t * c;
                                  tau = s / ((CKL_REAL)1.0 + c);
                                  h = t * a[ip][iq];
                                  z[ip] -= h;
                                  z[iq] += h;
                                  d[ip] -= h;
                                  d[iq] += h;
                                  a[ip][iq] = 0.0;
                                  for(j = 0; j < ip; j++)
                                  {
                                    ROTATE(a, j, ip, j, iq);
                                  }
                                  for(j = ip + 1; j < iq; j++)
                                  {
                                    ROTATE(a, ip, j, j, iq);
                                  }
                                  for(j = iq + 1; j < n; j++)
                                  {
                                    ROTATE(a, ip, j, iq, j);
                                  }
                                  for(j = 0; j < n; j++)
                                  {
                                    ROTATE(v, j, ip, j, iq);
                                  }
                                  nrot++;
                                }
                              }
    for(ip = 0; ip < n; ip++)
    {
      b[ip] += z[ip];
      d[ip] = b[ip];
      z[ip] = 0.0;
    }
  }

  std::cerr << "eigen: too many iterations in Jacobi transform.\n";

  return;
}


inline void QRotM(CKL_REAL destQuat[4], CKL_REAL srcMatrix[3][3])
{
  const int next[3] = {1, 2, 0};

  CKL_REAL trace = srcMatrix[0][0] + srcMatrix[1][1] + srcMatrix[2][2];
  CKL_REAL root;

  if(trace > 0.0)
  {
    // |w| > 1/2, may as well choose w > 1/2
    root = sqrt(trace + 1.0);  // 2w
    destQuat[0] = 0.5 * root;
    root = 0.5 / root;  // 1/(4w)
    destQuat[1] = (srcMatrix[2][1] - srcMatrix[1][2])*root;
    destQuat[2] = (srcMatrix[0][2] - srcMatrix[2][0])*root;
    destQuat[3] = (srcMatrix[1][0] - srcMatrix[0][1])*root;
  }
  else
  {
    // |w| <= 1/2
    int i = 0;
    if(srcMatrix[1][1] > srcMatrix[0][0])
    {
      i = 1;
    }
    if(srcMatrix[2][2] > srcMatrix[i][i])
    {
      i = 2;
    }
    int j = next[i];
    int k = next[j];

    root = sqrt(srcMatrix[i][i] - srcMatrix[j][j] - srcMatrix[k][k] + 1.0);
    CKL_REAL* quat[3] = { &destQuat[1], &destQuat[2], &destQuat[3] };
    *quat[i] = 0.5 * root;
    root = 0.5 / root;
    destQuat[0] = (srcMatrix[k][j] - srcMatrix[j][k]) * root;
    *quat[j] = (srcMatrix[j][i] + srcMatrix[i][j]) * root;
    *quat[k] = (srcMatrix[k][i] + srcMatrix[i][k]) * root;
  }
}

inline void Quatinverse(CKL_REAL destInvQuat[4], CKL_REAL srcQuat[4])
{
  CKL_REAL sqr_length = srcQuat[0] * srcQuat[0] + srcQuat[1] * srcQuat[1] + srcQuat[2] * srcQuat[2] + srcQuat[3] * srcQuat[3];
  if(sqr_length > 0)
  {
    CKL_REAL inv_length = 1 / sqrt(sqr_length);
    destInvQuat[0] = srcQuat[0] * inv_length;
    destInvQuat[1] = -srcQuat[1] * inv_length;
    destInvQuat[2] = -srcQuat[2] * inv_length;
    destInvQuat[3] = -srcQuat[3] * inv_length;
  }
  else
  {
    destInvQuat[0] = srcQuat[0];
    destInvQuat[1] = -srcQuat[1];
    destInvQuat[2] = -srcQuat[2];
    destInvQuat[3] = -srcQuat[3];
  }
}

inline void QuatxQuat(CKL_REAL destQuat[4], CKL_REAL quat1[4], CKL_REAL quat2[4])
{

    destQuat[0] = quat1[0] * quat2[0] - quat1[1] * quat2[1] - quat1[2] * quat2[2] - quat1[3] * quat2[3];
    destQuat[1] = quat1[0] * quat2[1] + quat1[1] * quat2[0] + quat1[2] * quat2[3] - quat1[3] * quat2[2];
    destQuat[2] = quat1[0] * quat2[2] - quat1[1] * quat2[3] + quat1[2] * quat2[0] + quat1[3] * quat2[1];
    destQuat[3] = quat1[0] * quat2[3] + quat1[1] * quat2[2] - quat1[2] * quat2[1] + quat1[3] * quat2[0];
}

inline CKL_REAL QuatDot(CKL_REAL quat1[4], CKL_REAL quat2[4])
{
  return quat1[0] * quat2[0] + quat1[1] * quat2[1] + quat1[2] * quat2[2] + quat1[3] * quat2[3];
}

inline void VInterp(CKL_REAL Vr[3], CKL_REAL V1[3], CKL_REAL V2[3], CKL_REAL t)
{
  Vr[0] = V1[0] + t * (V2[0] - V1[0]);
  Vr[1] = V1[1] + t * (V2[1] - V1[1]);
  Vr[2] = V1[2] + t * (V2[1] - V1[2]);
}

inline void VRay(CKL_REAL Vr[3], CKL_REAL V1[3], CKL_REAL ray[3], CKL_REAL t)
{
  Vr[0] = V1[0] + t * ray[0];
  Vr[1] = V1[1] + t * ray[1];
  Vr[2] = V1[2] + t * ray[2];
}

inline void MRotEuler(CKL_REAL destMatrix[3][3], CKL_REAL euler[3])
{
    CKL_REAL ci(cos(euler[0]));
    CKL_REAL cj(cos(euler[1]));
    CKL_REAL ch(cos(euler[2]));
    CKL_REAL si(sin(euler[0]));
    CKL_REAL sj(sin(euler[1]));
    CKL_REAL sh(sin(euler[2]));
    CKL_REAL cc = ci * ch;
    CKL_REAL cs = ci * sh;
    CKL_REAL sc = si * ch;
    CKL_REAL ss = si * sh;

    destMatrix[0][0] = cj * ch;
    destMatrix[0][1] = sj * sc - cs;
    destMatrix[0][2] = sj * cc + ss;
    destMatrix[1][0] = cj * sh;
    destMatrix[1][1] = sj * ss + cc;
    destMatrix[1][2] = sj * cs - sc;
    destMatrix[2][0] = -sj;
    destMatrix[2][1] = cj * si;
    destMatrix[2][2] = cj * ci;
}

inline void EulerRotQuat(CKL_REAL destEuler[3], CKL_REAL quat[4])
{
  CKL_REAL R[3][3];
  MRotQ(R, quat);
  CKL_REAL a = atan2(R[1][0], R[0][0]);
  CKL_REAL b = asin(-R[2][0]);
  CKL_REAL c = atan2(R[2][1], R[2][2]);

  if(b == M_PI * 0.5)
  {
    if(a > 0)
      a -= M_PI;
    else
      a += M_PI;
    if(c > 0)
      c -= M_PI;
    else
      c += M_PI;
  }

  destEuler[0] = a;
  destEuler[1] = b;
  destEuler[2] = c;
}

inline void Slerp(CKL_REAL destQuat[4], CKL_REAL srcQuat1[4], CKL_REAL srcQuat2[4], CKL_REAL t)
{
  CKL_REAL cs = QuatDot(srcQuat1, srcQuat2);
  CKL_REAL angle = acos(cs);

  if(fabs(angle) > 0)
  {
    CKL_REAL sn = sin(angle);
    CKL_REAL invSn = 1.0 / sn;
    CKL_REAL tAngle = t * angle;
    CKL_REAL coeff0 = sin(angle - tAngle) * invSn;
    CKL_REAL coeff1 = sin(tAngle) * invSn;

    destQuat[0] = coeff0 * srcQuat1[0] + coeff1 * srcQuat2[0];
    destQuat[1] = coeff0 * srcQuat1[1] + coeff1 * srcQuat2[1];
    destQuat[2] = coeff0 * srcQuat1[2] + coeff1 * srcQuat2[2];
    destQuat[3] = coeff0 * srcQuat1[3] + coeff1 * srcQuat2[3];
  }
  else
  {
    destQuat[0] = srcQuat1[0];
    destQuat[1] = srcQuat1[1];
    destQuat[2] = srcQuat1[2];
    destQuat[3] = srcQuat1[3];
  }
}


}

#endif
