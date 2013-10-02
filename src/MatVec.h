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

}

#endif
