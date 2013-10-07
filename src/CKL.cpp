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

  US Mail:             S. Gottschalk, E. Larsen
                       Department of Computer Science
                       Sitterson Hall, CB #3175
                       University of N. Carolina
                       Chapel Hill, NC 27599-3175

  Phone:               (919)962-1749

  EMail:               geom@cs.unc.edu


\**************************************************************************/

#include <cstdio>
#include <string.h>
#include <iostream>
#include "CKL.h"
#include "BVTQ.h"
#include "Build.h"
#include "MatVec.h"
#include "GetTime.h"
#include "TriDist.h"
#include "Classifier.h"
#include "NearestNeighbors.h"

namespace CKL
{

CKL_REAL distanceToPlane(const CKL_REAL n[3], CKL_REAL t, const CKL_REAL v[3])
{
  return VdotV(n, v) - t;
}


bool buildTrianglePlane(const CKL_REAL V1[3], const CKL_REAL V2[3], const CKL_REAL V3[3], CKL_REAL n[3], CKL_REAL& t)
{
  CKL_REAL V2mV1[3];
  CKL_REAL V3mV1[3];
  VmV(V2mV1, V2, V1);
  VmV(V3mV1, V3, V1);
  VcrossV(n, V2mV1, V3mV1);
  Vnormalize(n);
  t = VdotV(n, V1);
  return true;
}

void computeDeepestPoints(const CKL_REAL* clipped_points, std::size_t num_clipped_points,
                          const CKL_REAL n[3], CKL_REAL t,
                          CKL_REAL* penetration_depth,
                          CKL_REAL* deepest_points,
                          std::size_t* num_deepest_points)
{
  *num_deepest_points = 0;
  CKL_REAL max_depth = -std::numeric_limits<CKL_REAL>::max();
  std::size_t num_deepest_points_ = 0;
  std::size_t num_neg = 0;
  std::size_t num_pos = 0;
  std::size_t num_zero = 0;

  for(std::size_t i = 0; i < num_clipped_points; ++i)
  {
    CKL_REAL dist = -distanceToPlane(n, t, clipped_points + 3 * i);
    if(dist > 1e-6) num_pos++;
    else if(dist < -1e-6) num_neg++;
    else num_zero++;
    if(dist > max_depth)
    {
      max_depth = dist;
      num_deepest_points_ = 1;
      VcV(deepest_points + 3 * (num_deepest_points_ - 1), clipped_points + 3 * i);
    }
    else if(dist + 1e-6 >= max_depth)
    {
      num_deepest_points_++;
      VcV(deepest_points + 3 * (num_deepest_points_ - 1), clipped_points + 3 * i);
    }
  }

  if(max_depth < -1e-6)
    num_deepest_points_ = 0;

  if(num_zero == 0 && ((num_neg == 0) || (num_pos == 0)))
    num_deepest_points_ = 0;

  *penetration_depth = max_depth;
  *num_deepest_points = num_deepest_points_;
}







                       
enum BUILD_STATE
  {
    CKL_BUILD_STATE_EMPTY,     // empty state, immediately after constructor
    CKL_BUILD_STATE_BEGUN,     // after BeginModel(), state for adding triangles
    CKL_BUILD_STATE_PROCESSED  // after tree has been built, ready to use
  };

CKL_Model::CKL_Model()
{
  // no bounding volume tree yet
  
  b = 0;
  num_bvs_alloced = 0;
  num_bvs = 0;
  
  // no tri list yet
  
  tris = 0;
  num_tris = 0;
  num_tris_alloced = 0;
  
  last_tri = 0;
  
  build_state = CKL_BUILD_STATE_EMPTY;
}

CKL_Model::~CKL_Model()
{
  if(b != NULL)
    delete [] b;
  if(tris != NULL)
    delete [] tris;
}

int CKL_Model::BeginModel(int n)
{
  // reset to initial state if necessary
  
  if(build_state != CKL_BUILD_STATE_EMPTY)
  {
    delete [] b;
    delete [] tris;
    
    num_tris = num_bvs = num_tris_alloced = num_bvs_alloced = 0;
  }
  
  // prepare model for addition of triangles
  
  if(n <= 0) n = 8;
  num_tris_alloced = n;
  tris = new Tri[n];
  if(!tris)
  {
    std::cerr << "CKL Error!  Out of memory for tri array on BeginModel() call!" << std::endl;
    return CKL_ERR_MODEL_OUT_OF_MEMORY;
  }
  
  // give a warning if called out of sequence
  
  if(build_state != CKL_BUILD_STATE_EMPTY)
  {
    std::cerr << "CKL Warning! Called BeginModel() on a CKL_Model that\n"
              << "was not empty. This model was cleared and previous\n"
              << "triangle additions were lost.\n";
    build_state = CKL_BUILD_STATE_BEGUN;
    return CKL_ERR_BUILD_OUT_OF_SEQUENCE;
  }
  
  build_state = CKL_BUILD_STATE_BEGUN;
  return CKL_OK;
}

int CKL_Model::AddTri(const CKL_REAL *p1,
                      const CKL_REAL *p2,
                      const CKL_REAL *p3,
                      int id)
{
  if(build_state == CKL_BUILD_STATE_EMPTY)
  {
    BeginModel();
  }
  else if(build_state == CKL_BUILD_STATE_PROCESSED)
  {
    std::cerr << "CKL Warning! Called AddTri() on CKL_Model \n"
              << "object that was already ended. AddTri() was\n"
              << "ignored.  Must do a BeginModel() to clear the\n"
              << "model for addition of new triangles\n";
    return CKL_ERR_BUILD_OUT_OF_SEQUENCE;
  }
  
  // allocate for new triangles
  
  if(num_tris >= num_tris_alloced)
  {
    Tri *temp;
    temp = new Tri[num_tris_alloced * 2];
    if(!temp)
    {
      std::cerr << "CKL Error!  Out of memory for tri array on"
                << " AddTri() call!\n";
      return CKL_ERR_MODEL_OUT_OF_MEMORY;
    }
    memcpy(temp, tris, sizeof(Tri)*num_tris);
    delete [] tris;
    tris = temp;
    num_tris_alloced = num_tris_alloced * 2;
  }
  
  // initialize the new triangle
  
  tris[num_tris].p1[0] = p1[0];
  tris[num_tris].p1[1] = p1[1];
  tris[num_tris].p1[2] = p1[2];
  
  tris[num_tris].p2[0] = p2[0];
  tris[num_tris].p2[1] = p2[1];
  tris[num_tris].p2[2] = p2[2];
  
  tris[num_tris].p3[0] = p3[0];
  tris[num_tris].p3[1] = p3[1];
  tris[num_tris].p3[2] = p3[2];
  
  tris[num_tris].id = id;
  
  num_tris += 1;
  
  return CKL_OK;
}

int CKL_Model::EndModel()
{
  if(build_state == CKL_BUILD_STATE_PROCESSED)
  {
    std::cerr << "CKL Warning! Called EndModel() on CKL_Model \n"
              << "object that was already ended. EndModel() was\n"
              << "ignored.  Must do a BeginModel() to clear the\n"
              << "model for addition of new triangles\n";
    return CKL_ERR_BUILD_OUT_OF_SEQUENCE;
  }
  
  // report error is no tris
  
  if(num_tris == 0)
  {
    std::cerr << "CKL Error! EndModel() called on model with"
              << " no triangles\n";
    return CKL_ERR_BUILD_EMPTY_MODEL;
  }
  
  // shrink fit tris array
  
  if(num_tris_alloced > num_tris)
  {
    Tri *new_tris = new Tri[num_tris];
    if(!new_tris)
    {
      std::cerr << "CKL Error!  Out of memory for tri array "
                << "in EndModel() call!\n";
      return CKL_ERR_MODEL_OUT_OF_MEMORY;
    }
    memcpy(new_tris, tris, sizeof(Tri)*num_tris);
    delete [] tris;
    tris = new_tris;
    num_tris_alloced = num_tris;
  }
  
  // create an array of BVs for the model
  
  b = new BV[2 * num_tris - 1];
  if(!b)
  {
    std::cerr << "CKL Error! out of memory for BV array "
              << "in EndModel()\n";
    return CKL_ERR_MODEL_OUT_OF_MEMORY;
  }
  num_bvs_alloced = 2 * num_tris - 1;
  num_bvs = 0;
  
  // we should build the model now.
  
  build_model(this);
  build_state = CKL_BUILD_STATE_PROCESSED;
  
  last_tri = tris;
  
  return CKL_OK;
}

int CKL_Model::MemUsage(int msg)
{
  int mem_bv_list = sizeof(BV) * num_bvs;
  int mem_tri_list = sizeof(Tri) * num_tris;
  
  int total_mem = mem_bv_list + mem_tri_list + sizeof(CKL_Model);
  
  if(msg)
  {
    std::cerr << "Total for model " << std::hex << this << ": " << total_mem << "bytes\n";
    std::cerr << "BVs: " << num_bvs << " alloced, take " << sizeof(BV) << " bytes each\n";
    std::cerr << "Tris: " << num_tris << " alloced, take " << sizeof(Tri) << " bytes each\n";
  }
  
  return total_mem;
}

//  COLLIDE STUFF
//
//--------------------------------------------------------------------------

CKL_CollideResult::CKL_CollideResult()
{
  pairs = 0;
  num_pairs = num_pairs_alloced = 0;
  num_bv_tests = 0;
  num_tri_tests = 0;
}

CKL_CollideResult::~CKL_CollideResult()
{
  delete [] pairs;
}

void CKL_CollideResult::FreePairsList()
{
  num_pairs = num_pairs_alloced = 0;
  delete [] pairs;
  pairs = 0;
}

// may increase OR reduce mem usage
void CKL_CollideResult::SizeTo(int n)
{
  CollisionPair *temp;
  
  if(n < num_pairs)
  {
    std::cerr << "CKL Error: Internal error in "
              << "'CKL_CollideResult::SizeTo(int n)'\n";
    std::cerr << "       n = " << n << ", but num_pairs = " << num_pairs << "\n";
    return;
  }
  
  temp = new CollisionPair[n];
  memcpy(temp, pairs, num_pairs * sizeof(CollisionPair));
  delete [] pairs;
  pairs = temp;
  num_pairs_alloced = n;
  return;
}

void CKL_CollideResult::Add(int a, int b)
{
  if(num_pairs >= num_pairs_alloced)
  {
    // allocate more
    
    SizeTo(num_pairs_alloced * 2 + 8);
  }
  
  // now proceed as usual
  
  pairs[num_pairs].id1 = a;
  pairs[num_pairs].id2 = b;
  num_pairs++;
}

void CKL_CollideResult::Add(int a, int b, CKL_REAL contact_point[3], CKL_REAL contact_normal[3])
{
  if(num_pairs >= num_pairs_alloced)
  {
    // allocate more
    
    SizeTo(num_pairs_alloced * 2 + 8);
  }
  
  // now proceed as usual
  
  pairs[num_pairs].id1 = a;
  pairs[num_pairs].id2 = b;
  VcV(pairs[num_pairs].point, contact_point);
  VcV(pairs[num_pairs].normal, contact_normal);
  num_pairs++;
}


// TRIANGLE OVERLAP TEST

inline CKL_REAL max(CKL_REAL a, CKL_REAL b, CKL_REAL c)
{
  CKL_REAL t = a;
  if(b > t) t = b;
  if(c > t) t = c;
  return t;
}

inline CKL_REAL min(CKL_REAL a, CKL_REAL b, CKL_REAL c)
{
  CKL_REAL t = a;
  if(b < t) t = b;
  if(c < t) t = c;
  return t;
}

int project6(CKL_REAL *ax,
             CKL_REAL *p1, CKL_REAL *p2, CKL_REAL *p3,
             CKL_REAL *q1, CKL_REAL *q2, CKL_REAL *q3)
{
  CKL_REAL P1 = VdotV(ax, p1);
  CKL_REAL P2 = VdotV(ax, p2);
  CKL_REAL P3 = VdotV(ax, p3);
  CKL_REAL Q1 = VdotV(ax, q1);
  CKL_REAL Q2 = VdotV(ax, q2);
  CKL_REAL Q3 = VdotV(ax, q3);
  
  CKL_REAL mx1 = max(P1, P2, P3);
  CKL_REAL mn1 = min(P1, P2, P3);
  CKL_REAL mx2 = max(Q1, Q2, Q3);
  CKL_REAL mn2 = min(Q1, Q2, Q3);
  
  if(mn1 > mx2) return 0;
  if(mn2 > mx1) return 0;
  return 1;
}

// very robust triangle intersection test
// uses no divisions
// works on coplanar triangles
int TriContact(CKL_REAL *P1, CKL_REAL *P2, CKL_REAL *P3,
               CKL_REAL *Q1, CKL_REAL *Q2, CKL_REAL *Q3,
               CKL_REAL* contact_points = NULL,
               std::size_t* num_contact_points = NULL,
               CKL_REAL* penetration_depth = NULL,
               CKL_REAL* normal = NULL)
{

  // One triangle is (p1,p2,p3).  Other is (q1,q2,q3).
  // Edges are (e1,e2,e3) and (f1,f2,f3).
  // Normals are n1 and m1
  // Outwards are (g1,g2,g3) and (h1,h2,h3).
  //
  // We assume that the triangle vertices are in the same coordinate system.
  //
  // First thing we do is establish a new c.s. so that p1 is at (0,0,0).
  
  CKL_REAL p1[3], p2[3], p3[3];
  CKL_REAL q1[3], q2[3], q3[3];
  CKL_REAL e1[3], e2[3], e3[3];
  CKL_REAL f1[3], f2[3], f3[3];
  CKL_REAL g1[3], g2[3], g3[3];
  CKL_REAL h1[3], h2[3], h3[3];
  CKL_REAL n1[3], m1[3];
  
  CKL_REAL ef11[3], ef12[3], ef13[3];
  CKL_REAL ef21[3], ef22[3], ef23[3];
  CKL_REAL ef31[3], ef32[3], ef33[3];
  
  p1[0] = P1[0] - P1[0];
  p1[1] = P1[1] - P1[1];
  p1[2] = P1[2] - P1[2];
  p2[0] = P2[0] - P1[0];
  p2[1] = P2[1] - P1[1];
  p2[2] = P2[2] - P1[2];
  p3[0] = P3[0] - P1[0];
  p3[1] = P3[1] - P1[1];
  p3[2] = P3[2] - P1[2];
  
  q1[0] = Q1[0] - P1[0];
  q1[1] = Q1[1] - P1[1];
  q1[2] = Q1[2] - P1[2];
  q2[0] = Q2[0] - P1[0];
  q2[1] = Q2[1] - P1[1];
  q2[2] = Q2[2] - P1[2];
  q3[0] = Q3[0] - P1[0];
  q3[1] = Q3[1] - P1[1];
  q3[2] = Q3[2] - P1[2];
  
  e1[0] = p2[0] - p1[0];
  e1[1] = p2[1] - p1[1];
  e1[2] = p2[2] - p1[2];
  e2[0] = p3[0] - p2[0];
  e2[1] = p3[1] - p2[1];
  e2[2] = p3[2] - p2[2];
  e3[0] = p1[0] - p3[0];
  e3[1] = p1[1] - p3[1];
  e3[2] = p1[2] - p3[2];
  
  f1[0] = q2[0] - q1[0];
  f1[1] = q2[1] - q1[1];
  f1[2] = q2[2] - q1[2];
  f2[0] = q3[0] - q2[0];
  f2[1] = q3[1] - q2[1];
  f2[2] = q3[2] - q2[2];
  f3[0] = q1[0] - q3[0];
  f3[1] = q1[1] - q3[1];
  f3[2] = q1[2] - q3[2];
  
  VcrossV(n1, e1, e2);
  VcrossV(m1, f1, f2);
  
  VcrossV(g1, e1, n1);
  VcrossV(g2, e2, n1);
  VcrossV(g3, e3, n1);
  VcrossV(h1, f1, m1);
  VcrossV(h2, f2, m1);
  VcrossV(h3, f3, m1);
  
  VcrossV(ef11, e1, f1);
  VcrossV(ef12, e1, f2);
  VcrossV(ef13, e1, f3);
  VcrossV(ef21, e2, f1);
  VcrossV(ef22, e2, f2);
  VcrossV(ef23, e2, f3);
  VcrossV(ef31, e3, f1);
  VcrossV(ef32, e3, f2);
  VcrossV(ef33, e3, f3);
  
  // now begin the series of tests
  
  if(!project6(n1, p1, p2, p3, q1, q2, q3)) return 0;
  if(!project6(m1, p1, p2, p3, q1, q2, q3)) return 0;
  
  if(!project6(ef11, p1, p2, p3, q1, q2, q3)) return 0;
  if(!project6(ef12, p1, p2, p3, q1, q2, q3)) return 0;
  if(!project6(ef13, p1, p2, p3, q1, q2, q3)) return 0;
  if(!project6(ef21, p1, p2, p3, q1, q2, q3)) return 0;
  if(!project6(ef22, p1, p2, p3, q1, q2, q3)) return 0;
  if(!project6(ef23, p1, p2, p3, q1, q2, q3)) return 0;
  if(!project6(ef31, p1, p2, p3, q1, q2, q3)) return 0;
  if(!project6(ef32, p1, p2, p3, q1, q2, q3)) return 0;
  if(!project6(ef33, p1, p2, p3, q1, q2, q3)) return 0;
  
  if(!project6(g1, p1, p2, p3, q1, q2, q3)) return 0;
  if(!project6(g2, p1, p2, p3, q1, q2, q3)) return 0;
  if(!project6(g3, p1, p2, p3, q1, q2, q3)) return 0;
  if(!project6(h1, p1, p2, p3, q1, q2, q3)) return 0;
  if(!project6(h2, p1, p2, p3, q1, q2, q3)) return 0;
  if(!project6(h3, p1, p2, p3, q1, q2, q3)) return 0;

  if(contact_points && num_contact_points && penetration_depth && normal)
  {
    CKL_REAL n1[3], n2[3];
    CKL_REAL t1, t2;
    buildTrianglePlane(P1, P2, P3, n1, t1);
    buildTrianglePlane(Q1, Q2, Q3, n2, t2);

    CKL_REAL deepest_points1[9];
    std::size_t num_deepest_points1 = 0;
    CKL_REAL deepest_points2[9];
    std::size_t num_deepest_points2 = 0;
    CKL_REAL penetration_depth1, penetration_depth2;

    CKL_REAL P[9];
    CKL_REAL Q[9];
    VcV(P, P1); VcV(P + 3, P2); VcV(P + 6, P3);
    VcV(Q, Q1); VcV(Q + 3, Q2); VcV(Q + 6, Q3);
    
    computeDeepestPoints(Q, 3, n1, t1, &penetration_depth2, deepest_points2, &num_deepest_points2);
    computeDeepestPoints(P, 3, n2, t2, &penetration_depth1, deepest_points1, &num_deepest_points1);

    if(penetration_depth1 > penetration_depth2)
    {
      *num_contact_points = std::min(num_deepest_points2, (std::size_t)2);
      for(std::size_t i = 0; i < *num_contact_points; ++i)
      {
        VcV(contact_points + 3 * i, deepest_points2 + 3 * i);
      }

      VcV(normal, n1);
      *penetration_depth = penetration_depth2;
    }
    else
    {
      *num_contact_points = std::min(num_deepest_points1, (std::size_t)2);
      for(std::size_t i = 0; i < *num_contact_points; ++i)
      {
        VcV(contact_points + 3 * i, deepest_points1 + 3 * i);
      }

      VxS(normal, n2, -1.0);
      *penetration_depth = penetration_depth1;
    }
  }
  
  return 1;
}

inline CKL_REAL TriDistance(CKL_REAL R[3][3], CKL_REAL T[3], Tri *t1, Tri *t2,
                            CKL_REAL p[3], CKL_REAL q[3])
{
  // transform tri 2 into same space as tri 1
  
  CKL_REAL tri1[3][3], tri2[3][3];
  
  VcV(tri1[0], t1->p1);
  VcV(tri1[1], t1->p2);
  VcV(tri1[2], t1->p3);
  MxVpV(tri2[0], R, t2->p1, T);
  MxVpV(tri2[1], R, t2->p2, T);
  MxVpV(tri2[2], R, t2->p3, T);
  
  return TriDist(p, q, tri1, tri2);
}


void CollideRecurse(CKL_CollideResult *res,
                    CKL_REAL R[3][3], CKL_REAL T[3], // b2 relative to b1
                    CKL_Model *o1, int b1,
                    CKL_Model *o2, int b2, int flag)
{
  // first thing, see if we're overlapping
  
  res->num_bv_tests++;
  
  if(!BV_Overlap(R, T, o1->child(b1), o2->child(b2))) return;
  
  // if we are, see if we test triangles next
  
  int l1 = o1->child(b1)->Leaf();
  int l2 = o2->child(b2)->Leaf();
  
  if(l1 && l2)
  {
    res->num_tri_tests++;
    
    // transform the points in b2 into space of b1, then compare
    
    Tri *t1 = &o1->tris[-o1->child(b1)->first_child - 1];
    Tri *t2 = &o2->tris[-o2->child(b2)->first_child - 1];
    CKL_REAL q1[3], q2[3], q3[3];
    CKL_REAL *p1 = t1->p1;
    CKL_REAL *p2 = t1->p2;
    CKL_REAL *p3 = t1->p3;
    MxVpV(q1, res->R, t2->p1, res->T);
    MxVpV(q2, res->R, t2->p2, res->T);
    MxVpV(q3, res->R, t2->p3, res->T);

    CKL_REAL contact_point[6];
    CKL_REAL contact_normal[6];
    std::size_t num_contact_point;
    CKL_REAL penetration_depth;
    
    if(TriContact(p1, p2, p3, q1, q2, q3, contact_point, &num_contact_point, &penetration_depth, contact_normal))
    {
      // add this to result

      for(int j = 0; j < num_contact_point; ++j)
      {
        res->Add(t1->id, t2->id, contact_point + 3 * j, contact_normal + 3 * j);
      }
    }
    
    return;
  }
  
  // we dont, so decide whose children to visit next
  
  CKL_REAL sz1 = o1->child(b1)->GetSize();
  CKL_REAL sz2 = o2->child(b2)->GetSize();
  
  CKL_REAL Rc[3][3], Tc[3], Ttemp[3];
  
  if(l2 || (!l1 && (sz1 > sz2)))
  {
    int c1 = o1->child(b1)->first_child;
    int c2 = c1 + 1;
    
    MTxM(Rc, o1->child(c1)->R, R);
#if CKL_BV_TYPE & OBB_TYPE
    VmV(Ttemp, T, o1->child(c1)->To);
#else
    VmV(Ttemp, T, o1->child(c1)->Tr);
#endif
    MTxV(Tc, o1->child(c1)->R, Ttemp);
    CollideRecurse(res, Rc, Tc, o1, c1, o2, b2, flag);
    
    if((flag == CKL_FIRST_CONTACT) && (res->num_pairs > 0)) return;
    
    MTxM(Rc, o1->child(c2)->R, R);
#if CKL_BV_TYPE & OBB_TYPE
    VmV(Ttemp, T, o1->child(c2)->To);
#else
    VmV(Ttemp, T, o1->child(c2)->Tr);
#endif
    MTxV(Tc, o1->child(c2)->R, Ttemp);
    CollideRecurse(res, Rc, Tc, o1, c2, o2, b2, flag);
  }
  else
  {
    int c1 = o2->child(b2)->first_child;
    int c2 = c1 + 1;
    
    MxM(Rc, R, o2->child(c1)->R);
#if CKL_BV_TYPE & OBB_TYPE
    MxVpV(Tc, R, o2->child(c1)->To, T);
#else
    MxVpV(Tc, R, o2->child(c1)->Tr, T);
#endif
    CollideRecurse(res, Rc, Tc, o1, b1, o2, c1, flag);
    
    if((flag == CKL_FIRST_CONTACT) && (res->num_pairs > 0)) return;
    
    MxM(Rc, R, o2->child(c2)->R);
#if CKL_BV_TYPE & OBB_TYPE
    MxVpV(Tc, R, o2->child(c2)->To, T);
#else
    MxVpV(Tc, R, o2->child(c2)->Tr, T);
#endif
    CollideRecurse(res, Rc, Tc, o1, b1, o2, c2, flag);
  }
}

int CKL_Collide(CKL_CollideResult *res,
                CKL_REAL R1[3][3], CKL_REAL T1[3], CKL_Model *o1,
                CKL_REAL R2[3][3], CKL_REAL T2[3], CKL_Model *o2,
                int flag)
{
  double t1 = GetTime();
  
  // make sure that the models are built
  
  if(o1->build_state != CKL_BUILD_STATE_PROCESSED)
    return CKL_ERR_UNPROCESSED_MODEL;
  if(o2->build_state != CKL_BUILD_STATE_PROCESSED)
    return CKL_ERR_UNPROCESSED_MODEL;
    
  // clear the stats
  
  res->num_bv_tests = 0;
  res->num_tri_tests = 0;
  
  // don't release the memory, but reset the num_pairs counter
  
  res->num_pairs = 0;
  
  // Okay, compute what transform [R,T] that takes us from cs1 to cs2.
  // [R,T] = [R1,T1]'[R2,T2] = [R1',-R1'T][R2,T2] = [R1'R2, R1'(T2-T1)]
  // First compute the rotation part, then translation part
  
  MTxM(res->R, R1, R2);
  CKL_REAL Ttemp[3];
  VmV(Ttemp, T2, T1);
  MTxV(res->T, R1, Ttemp);
  
  // compute the transform from o1->child(0) to o2->child(0)
  
  CKL_REAL Rtemp[3][3], R[3][3], T[3];
  
  MxM(Rtemp, res->R, o2->child(0)->R);
  MTxM(R, o1->child(0)->R, Rtemp);
  
#if CKL_BV_TYPE & OBB_TYPE
  MxVpV(Ttemp, res->R, o2->child(0)->To, res->T);
  VmV(Ttemp, Ttemp, o1->child(0)->To);
#else
  MxVpV(Ttemp, res->R, o2->child(0)->Tr, res->T);
  VmV(Ttemp, Ttemp, o1->child(0)->Tr);
#endif
  
  MTxV(T, o1->child(0)->R, Ttemp);
  
  // now start with both top level BVs
  
  CollideRecurse(res, R, T, o1, 0, o2, 0, flag);
  
  double t2 = GetTime();
  res->query_time_secs = t2 - t1;
  
  return CKL_OK;
}

#if CKL_BV_TYPE & RSS_TYPE // distance/tolerance only available with RSS
// unless an OBB distance test is supplied in
// BV.cpp

// DISTANCE STUFF
//
//--------------------------------------------------------------------------

void DistanceRecurse(CKL_DistanceResult *res,
                     CKL_REAL R[3][3], CKL_REAL T[3], // b2 relative to b1
                     CKL_Model *o1, int b1,
                     CKL_Model *o2, int b2)
{
  CKL_REAL sz1 = o1->child(b1)->GetSize();
  CKL_REAL sz2 = o2->child(b2)->GetSize();
  int l1 = o1->child(b1)->Leaf();
  int l2 = o2->child(b2)->Leaf();
  
  if(l1 && l2)
  {
    // both leaves.  Test the triangles beneath them.
    
    res->num_tri_tests++;
    
    CKL_REAL p[3], q[3];
    
    Tri *t1 = &o1->tris[-o1->child(b1)->first_child - 1];
    Tri *t2 = &o2->tris[-o2->child(b2)->first_child - 1];
    
    CKL_REAL d = TriDistance(res->R, res->T, t1, t2, p, q);
    
    if(d < res->distance)
    {
      res->distance = d;
      
      VcV(res->p1, p);         // p already in c.s. 1
      VcV(res->p2, q);         // q must be transformed
      // into c.s. 2 later
      o1->last_tri = t1;
      o2->last_tri = t2;
    }
    
    return;
  }
  
  // First, perform distance tests on the children. Then traverse
  // them recursively, but test the closer pair first, the further
  // pair second.
  
  int a1, a2, c1, c2; // new bv tests 'a' and 'c'
  CKL_REAL R1[3][3], T1[3], R2[3][3], T2[3], Ttemp[3];
  
  if(l2 || (!l1 && (sz1 > sz2)))
  {
    // visit the children of b1
    
    a1 = o1->child(b1)->first_child;
    a2 = b2;
    c1 = o1->child(b1)->first_child + 1;
    c2 = b2;
    
    MTxM(R1, o1->child(a1)->R, R);
#if CKL_BV_TYPE & RSS_TYPE
    VmV(Ttemp, T, o1->child(a1)->Tr);
#else
    VmV(Ttemp, T, o1->child(a1)->To);
#endif
    MTxV(T1, o1->child(a1)->R, Ttemp);
    
    MTxM(R2, o1->child(c1)->R, R);
#if CKL_BV_TYPE & RSS_TYPE
    VmV(Ttemp, T, o1->child(c1)->Tr);
#else
    VmV(Ttemp, T, o1->child(c1)->To);
#endif
    MTxV(T2, o1->child(c1)->R, Ttemp);
  }
  else
  {
    // visit the children of b2
    
    a1 = b1;
    a2 = o2->child(b2)->first_child;
    c1 = b1;
    c2 = o2->child(b2)->first_child + 1;
    
    MxM(R1, R, o2->child(a2)->R);
#if CKL_BV_TYPE & RSS_TYPE
    MxVpV(T1, R, o2->child(a2)->Tr, T);
#else
    MxVpV(T1, R, o2->child(a2)->To, T);
#endif
    
    MxM(R2, R, o2->child(c2)->R);
#if CKL_BV_TYPE & RSS_TYPE
    MxVpV(T2, R, o2->child(c2)->Tr, T);
#else
    MxVpV(T2, R, o2->child(c2)->To, T);
#endif
  }
  
  res->num_bv_tests += 2;
  
  CKL_REAL d1 = BV_Distance(R1, T1, o1->child(a1), o2->child(a2));
  CKL_REAL d2 = BV_Distance(R2, T2, o1->child(c1), o2->child(c2));
  
  if(d2 < d1)
  {
    if((d2 < (res->distance - res->abs_err)) ||
       (d2 * (1 + res->rel_err) < res->distance))
    {
      DistanceRecurse(res, R2, T2, o1, c1, o2, c2);
    }
    
    if((d1 < (res->distance - res->abs_err)) ||
       (d1 * (1 + res->rel_err) < res->distance))
    {
      DistanceRecurse(res, R1, T1, o1, a1, o2, a2);
    }
  }
  else
  {
    if((d1 < (res->distance - res->abs_err)) ||
       (d1 * (1 + res->rel_err) < res->distance))
    {
      DistanceRecurse(res, R1, T1, o1, a1, o2, a2);
    }
    
    if((d2 < (res->distance - res->abs_err)) ||
       (d2 * (1 + res->rel_err) < res->distance))
    {
      DistanceRecurse(res, R2, T2, o1, c1, o2, c2);
    }
  }
}

void DistanceQueueRecurse(CKL_DistanceResult *res,
                          CKL_REAL R[3][3], CKL_REAL T[3],
                          CKL_Model *o1, int b1,
                          CKL_Model *o2, int b2)
{
  BVTQ bvtq(res->qsize);
  
  BVT min_test;
  min_test.b1 = b1;
  min_test.b2 = b2;
  McM(min_test.R, R);
  VcV(min_test.T, T);
  
  while(1)
  {
    int l1 = o1->child(min_test.b1)->Leaf();
    int l2 = o2->child(min_test.b2)->Leaf();
    
    if(l1 && l2)
    {
      // both leaves.  Test the triangles beneath them.
      
      res->num_tri_tests++;
      
      CKL_REAL p[3], q[3];
      
      Tri *t1 = &o1->tris[-o1->child(min_test.b1)->first_child - 1];
      Tri *t2 = &o2->tris[-o2->child(min_test.b2)->first_child - 1];
      
      CKL_REAL d = TriDistance(res->R, res->T, t1, t2, p, q);
      
      if(d < res->distance)
      {
        res->distance = d;
        
        VcV(res->p1, p);         // p already in c.s. 1
        VcV(res->p2, q);         // q must be transformed
        // into c.s. 2 later
        o1->last_tri = t1;
        o2->last_tri = t2;
      }
    }
    else if(bvtq.GetNumTests() == bvtq.GetSize() - 1)
    {
      // queue can't get two more tests, recur
      
      DistanceQueueRecurse(res, min_test.R, min_test.T,
                           o1, min_test.b1, o2, min_test.b2);
    }
    else
    {
      // decide how to descend to children
      
      CKL_REAL sz1 = o1->child(min_test.b1)->GetSize();
      CKL_REAL sz2 = o2->child(min_test.b2)->GetSize();
      
      res->num_bv_tests += 2;
      
      BVT bvt1, bvt2;
      CKL_REAL Ttemp[3];
      
      if(l2 || (!l1 && (sz1 > sz2)))
      {
        // put new tests on queue consisting of min_test.b2
        // with children of min_test.b1
        
        int c1 = o1->child(min_test.b1)->first_child;
        int c2 = c1 + 1;
        
        // init bv test 1
        
        bvt1.b1 = c1;
        bvt1.b2 = min_test.b2;
        MTxM(bvt1.R, o1->child(c1)->R, min_test.R);
#if CKL_BV_TYPE & RSS_TYPE
        VmV(Ttemp, min_test.T, o1->child(c1)->Tr);
#else
        VmV(Ttemp, min_test.T, o1->child(c1)->To);
#endif
        MTxV(bvt1.T, o1->child(c1)->R, Ttemp);
        bvt1.d = BV_Distance(bvt1.R, bvt1.T,
                             o1->child(bvt1.b1), o2->child(bvt1.b2));
                             
        // init bv test 2
        
        bvt2.b1 = c2;
        bvt2.b2 = min_test.b2;
        MTxM(bvt2.R, o1->child(c2)->R, min_test.R);
#if CKL_BV_TYPE & RSS_TYPE
        VmV(Ttemp, min_test.T, o1->child(c2)->Tr);
#else
        VmV(Ttemp, min_test.T, o1->child(c2)->To);
#endif
        MTxV(bvt2.T, o1->child(c2)->R, Ttemp);
        bvt2.d = BV_Distance(bvt2.R, bvt2.T,
                             o1->child(bvt2.b1), o2->child(bvt2.b2));
      }
      else
      {
        // put new tests on queue consisting of min_test.b1
        // with children of min_test.b2
        
        int c1 = o2->child(min_test.b2)->first_child;
        int c2 = c1 + 1;
        
        // init bv test 1
        
        bvt1.b1 = min_test.b1;
        bvt1.b2 = c1;
        MxM(bvt1.R, min_test.R, o2->child(c1)->R);
#if CKL_BV_TYPE & RSS_TYPE
        MxVpV(bvt1.T, min_test.R, o2->child(c1)->Tr, min_test.T);
#else
        MxVpV(bvt1.T, min_test.R, o2->child(c1)->To, min_test.T);
#endif
        bvt1.d = BV_Distance(bvt1.R, bvt1.T,
                             o1->child(bvt1.b1), o2->child(bvt1.b2));
                             
        // init bv test 2
        
        bvt2.b1 = min_test.b1;
        bvt2.b2 = c2;
        MxM(bvt2.R, min_test.R, o2->child(c2)->R);
#if CKL_BV_TYPE & RSS_TYPE
        MxVpV(bvt2.T, min_test.R, o2->child(c2)->Tr, min_test.T);
#else
        MxVpV(bvt2.T, min_test.R, o2->child(c2)->To, min_test.T);
#endif
        bvt2.d = BV_Distance(bvt2.R, bvt2.T,
                             o1->child(bvt2.b1), o2->child(bvt2.b2));
      }
      
      bvtq.AddTest(bvt1);
      bvtq.AddTest(bvt2);
    }
    
    if(bvtq.Empty())
    {
      break;
    }
    else
    {
      min_test = bvtq.ExtractMinTest();
      
      if((min_test.d + res->abs_err >= res->distance) &&
         ((min_test.d * (1 + res->rel_err)) >= res->distance))
      {
        break;
      }
    }
  }
}

int CKL_Distance(CKL_DistanceResult *res,
                 CKL_REAL R1[3][3], CKL_REAL T1[3], CKL_Model *o1,
                 CKL_REAL R2[3][3], CKL_REAL T2[3], CKL_Model *o2,
                 CKL_REAL rel_err, CKL_REAL abs_err,
                 int qsize)
{

  double time1 = GetTime();
  
  // make sure that the models are built
  
  if(o1->build_state != CKL_BUILD_STATE_PROCESSED)
    return CKL_ERR_UNPROCESSED_MODEL;
  if(o2->build_state != CKL_BUILD_STATE_PROCESSED)
    return CKL_ERR_UNPROCESSED_MODEL;
    
  // Okay, compute what transform [R,T] that takes us from cs2 to cs1.
  // [R,T] = [R1,T1]'[R2,T2] = [R1',-R1'T][R2,T2] = [R1'R2, R1'(T2-T1)]
  // First compute the rotation part, then translation part
  
  MTxM(res->R, R1, R2);
  CKL_REAL Ttemp[3];
  VmV(Ttemp, T2, T1);
  MTxV(res->T, R1, Ttemp);
  
  // establish initial upper bound using last triangles which
  // provided the minimum distance
  
  CKL_REAL p[3], q[3];
  res->distance = TriDistance(res->R, res->T, o1->last_tri, o2->last_tri, p, q);
  VcV(res->p1, p);
  VcV(res->p2, q);
  
  // initialize error bounds
  
  res->abs_err = abs_err;
  res->rel_err = rel_err;
  
  // clear the stats
  
  res->num_bv_tests = 0;
  res->num_tri_tests = 0;
  
  // compute the transform from o1->child(0) to o2->child(0)
  
  CKL_REAL Rtemp[3][3], R[3][3], T[3];
  
  MxM(Rtemp, res->R, o2->child(0)->R);
  MTxM(R, o1->child(0)->R, Rtemp);
  
#if CKL_BV_TYPE & RSS_TYPE
  MxVpV(Ttemp, res->R, o2->child(0)->Tr, res->T);
  VmV(Ttemp, Ttemp, o1->child(0)->Tr);
#else
  MxVpV(Ttemp, res->R, o2->child(0)->To, res->T);
  VmV(Ttemp, Ttemp, o1->child(0)->To);
#endif
  MTxV(T, o1->child(0)->R, Ttemp);
  
  // choose routine according to queue size
  
  if(qsize <= 2)
  {
    DistanceRecurse(res, R, T, o1, 0, o2, 0);
  }
  else
  {
    res->qsize = qsize;
    
    DistanceQueueRecurse(res, R, T, o1, 0, o2, 0);
  }
  
  // res->p2 is in cs 1 ; transform it to cs 2
  
  CKL_REAL u[3];
  VmV(u, res->p2, res->T);
  MTxV(res->p2, res->R, u);
  
  double time2 = GetTime();
  res->query_time_secs = time2 - time1;
  
  return CKL_OK;
}

// Tolerance Stuff
//
//---------------------------------------------------------------------------
void ToleranceRecurse(CKL_ToleranceResult *res,
                      CKL_REAL R[3][3], CKL_REAL T[3],
                      CKL_Model *o1, int b1, CKL_Model *o2, int b2)
{
  CKL_REAL sz1 = o1->child(b1)->GetSize();
  CKL_REAL sz2 = o2->child(b2)->GetSize();
  int l1 = o1->child(b1)->Leaf();
  int l2 = o2->child(b2)->Leaf();
  
  if(l1 && l2)
  {
    // both leaves - find if tri pair within tolerance
    
    res->num_tri_tests++;
    
    CKL_REAL p[3], q[3];
    
    Tri *t1 = &o1->tris[-o1->child(b1)->first_child - 1];
    Tri *t2 = &o2->tris[-o2->child(b2)->first_child - 1];
    
    CKL_REAL d = TriDistance(res->R, res->T, t1, t2, p, q);
    
    if(d <= res->tolerance)
    {
      // triangle pair distance less than tolerance
      
      res->closer_than_tolerance = 1;
      res->distance = d;
      VcV(res->p1, p);         // p already in c.s. 1
      VcV(res->p2, q);         // q must be transformed
      // into c.s. 2 later
    }
    
    return;
  }
  
  int a1, a2, c1, c2; // new bv tests 'a' and 'c'
  CKL_REAL R1[3][3], T1[3], R2[3][3], T2[3], Ttemp[3];
  
  if(l2 || (!l1 && (sz1 > sz2)))
  {
    // visit the children of b1
    
    a1 = o1->child(b1)->first_child;
    a2 = b2;
    c1 = o1->child(b1)->first_child + 1;
    c2 = b2;
    
    MTxM(R1, o1->child(a1)->R, R);
#if CKL_BV_TYPE & RSS_TYPE
    VmV(Ttemp, T, o1->child(a1)->Tr);
#else
    VmV(Ttemp, T, o1->child(a1)->To);
#endif
    MTxV(T1, o1->child(a1)->R, Ttemp);
    
    MTxM(R2, o1->child(c1)->R, R);
#if CKL_BV_TYPE & RSS_TYPE
    VmV(Ttemp, T, o1->child(c1)->Tr);
#else
    VmV(Ttemp, T, o1->child(c1)->To);
#endif
    MTxV(T2, o1->child(c1)->R, Ttemp);
  }
  else
  {
    // visit the children of b2
    
    a1 = b1;
    a2 = o2->child(b2)->first_child;
    c1 = b1;
    c2 = o2->child(b2)->first_child + 1;
    
    MxM(R1, R, o2->child(a2)->R);
#if CKL_BV_TYPE & RSS_TYPE
    MxVpV(T1, R, o2->child(a2)->Tr, T);
#else
    MxVpV(T1, R, o2->child(a2)->To, T);
#endif
    MxM(R2, R, o2->child(c2)->R);
#if CKL_BV_TYPE & RSS_TYPE
    MxVpV(T2, R, o2->child(c2)->Tr, T);
#else
    MxVpV(T2, R, o2->child(c2)->To, T);
#endif
  }
  
  res->num_bv_tests += 2;
  
  CKL_REAL d1 = BV_Distance(R1, T1, o1->child(a1), o2->child(a2));
  CKL_REAL d2 = BV_Distance(R2, T2, o1->child(c1), o2->child(c2));
  
  if(d2 < d1)
  {
    if(d2 <= res->tolerance) ToleranceRecurse(res, R2, T2, o1, c1, o2, c2);
    if(res->closer_than_tolerance) return;
    if(d1 <= res->tolerance) ToleranceRecurse(res, R1, T1, o1, a1, o2, a2);
  }
  else
  {
    if(d1 <= res->tolerance) ToleranceRecurse(res, R1, T1, o1, a1, o2, a2);
    if(res->closer_than_tolerance) return;
    if(d2 <= res->tolerance) ToleranceRecurse(res, R2, T2, o1, c1, o2, c2);
  }
}

void ToleranceQueueRecurse(CKL_ToleranceResult *res,
                           CKL_REAL R[3][3], CKL_REAL T[3],
                           CKL_Model *o1, int b1,
                           CKL_Model *o2, int b2)
{
  BVTQ bvtq(res->qsize);
  BVT min_test;
  min_test.b1 = b1;
  min_test.b2 = b2;
  McM(min_test.R, R);
  VcV(min_test.T, T);
  
  while(1)
  {
    int l1 = o1->child(min_test.b1)->Leaf();
    int l2 = o2->child(min_test.b2)->Leaf();
    
    if(l1 && l2)
    {
      // both leaves - find if tri pair within tolerance
      
      res->num_tri_tests++;
      
      CKL_REAL p[3], q[3];
      
      Tri *t1 = &o1->tris[-o1->child(min_test.b1)->first_child - 1];
      Tri *t2 = &o2->tris[-o2->child(min_test.b2)->first_child - 1];
      
      CKL_REAL d = TriDistance(res->R, res->T, t1, t2, p, q);
      
      if(d <= res->tolerance)
      {
        // triangle pair distance less than tolerance
        
        res->closer_than_tolerance = 1;
        res->distance = d;
        VcV(res->p1, p);         // p already in c.s. 1
        VcV(res->p2, q);         // q must be transformed
        // into c.s. 2 later
        return;
      }
    }
    else if(bvtq.GetNumTests() == bvtq.GetSize() - 1)
    {
      // queue can't get two more tests, recur
      
      ToleranceQueueRecurse(res, min_test.R, min_test.T,
                            o1, min_test.b1, o2, min_test.b2);
      if(res->closer_than_tolerance == 1) return;
    }
    else
    {
      // decide how to descend to children
      
      CKL_REAL sz1 = o1->child(min_test.b1)->GetSize();
      CKL_REAL sz2 = o2->child(min_test.b2)->GetSize();
      
      res->num_bv_tests += 2;
      
      BVT bvt1, bvt2;
      CKL_REAL Ttemp[3];
      
      if(l2 || (!l1 && (sz1 > sz2)))
      {
        // add two new tests to queue, consisting of min_test.b2
        // with the children of min_test.b1
        
        int c1 = o1->child(min_test.b1)->first_child;
        int c2 = c1 + 1;
        
        // init bv test 1
        
        bvt1.b1 = c1;
        bvt1.b2 = min_test.b2;
        MTxM(bvt1.R, o1->child(c1)->R, min_test.R);
#if CKL_BV_TYPE & RSS_TYPE
        VmV(Ttemp, min_test.T, o1->child(c1)->Tr);
#else
        VmV(Ttemp, min_test.T, o1->child(c1)->To);
#endif
        MTxV(bvt1.T, o1->child(c1)->R, Ttemp);
        bvt1.d = BV_Distance(bvt1.R, bvt1.T,
                             o1->child(bvt1.b1), o2->child(bvt1.b2));
                             
        // init bv test 2
        
        bvt2.b1 = c2;
        bvt2.b2 = min_test.b2;
        MTxM(bvt2.R, o1->child(c2)->R, min_test.R);
#if CKL_BV_TYPE & RSS_TYPE
        VmV(Ttemp, min_test.T, o1->child(c2)->Tr);
#else
        VmV(Ttemp, min_test.T, o1->child(c2)->To);
#endif
        MTxV(bvt2.T, o1->child(c2)->R, Ttemp);
        bvt2.d = BV_Distance(bvt2.R, bvt2.T,
                             o1->child(bvt2.b1), o2->child(bvt2.b2));
      }
      else
      {
        // add two new tests to queue, consisting of min_test.b1
        // with the children of min_test.b2
        
        int c1 = o2->child(min_test.b2)->first_child;
        int c2 = c1 + 1;
        
        // init bv test 1
        
        bvt1.b1 = min_test.b1;
        bvt1.b2 = c1;
        MxM(bvt1.R, min_test.R, o2->child(c1)->R);
#if CKL_BV_TYPE & RSS_TYPE
        MxVpV(bvt1.T, min_test.R, o2->child(c1)->Tr, min_test.T);
#else
        MxVpV(bvt1.T, min_test.R, o2->child(c1)->To, min_test.T);
#endif
        bvt1.d = BV_Distance(bvt1.R, bvt1.T,
                             o1->child(bvt1.b1), o2->child(bvt1.b2));
                             
        // init bv test 2
        
        bvt2.b1 = min_test.b1;
        bvt2.b2 = c2;
        MxM(bvt2.R, min_test.R, o2->child(c2)->R);
#if CKL_BV_TYPE & RSS_TYPE
        MxVpV(bvt2.T, min_test.R, o2->child(c2)->Tr, min_test.T);
#else
        MxVpV(bvt2.T, min_test.R, o2->child(c2)->To, min_test.T);
#endif
        bvt2.d = BV_Distance(bvt2.R, bvt2.T,
                             o1->child(bvt2.b1), o2->child(bvt2.b2));
      }
      
      // put children tests in queue
      
      if(bvt1.d <= res->tolerance) bvtq.AddTest(bvt1);
      if(bvt2.d <= res->tolerance) bvtq.AddTest(bvt2);
    }
    
    if(bvtq.Empty() || (bvtq.MinTest() > res->tolerance))
    {
      res->closer_than_tolerance = 0;
      return;
    }
    else
    {
      min_test = bvtq.ExtractMinTest();
    }
  }
}

int CKL_Tolerance(CKL_ToleranceResult *res,
                  CKL_REAL R1[3][3], CKL_REAL T1[3], CKL_Model *o1,
                  CKL_REAL R2[3][3], CKL_REAL T2[3], CKL_Model *o2,
                  CKL_REAL tolerance,
                  int qsize)
{
  double time1 = GetTime();
  
  // make sure that the models are built
  
  if(o1->build_state != CKL_BUILD_STATE_PROCESSED)
    return CKL_ERR_UNPROCESSED_MODEL;
  if(o2->build_state != CKL_BUILD_STATE_PROCESSED)
    return CKL_ERR_UNPROCESSED_MODEL;
    
  // Compute the transform [R,T] that takes us from cs2 to cs1.
  // [R,T] = [R1,T1]'[R2,T2] = [R1',-R1'T][R2,T2] = [R1'R2, R1'(T2-T1)]
  
  MTxM(res->R, R1, R2);
  CKL_REAL Ttemp[3];
  VmV(Ttemp, T2, T1);
  MTxV(res->T, R1, Ttemp);
  
  // set tolerance, used to prune the search
  
  if(tolerance < 0.0) tolerance = 0.0;
  res->tolerance = tolerance;
  
  // clear the stats
  
  res->num_bv_tests = 0;
  res->num_tri_tests = 0;
  
  // initially assume not closer than tolerance
  
  res->closer_than_tolerance = 0;
  
  // compute the transform from o1->child(0) to o2->child(0)
  
  CKL_REAL Rtemp[3][3], R[3][3], T[3];
  
  MxM(Rtemp, res->R, o2->child(0)->R);
  MTxM(R, o1->child(0)->R, Rtemp);
#if CKL_BV_TYPE & RSS_TYPE
  MxVpV(Ttemp, res->R, o2->child(0)->Tr, res->T);
  VmV(Ttemp, Ttemp, o1->child(0)->Tr);
#else
  MxVpV(Ttemp, res->R, o2->child(0)->To, res->T);
  VmV(Ttemp, Ttemp, o1->child(0)->To);
#endif
  MTxV(T, o1->child(0)->R, Ttemp);
  
  // find a distance lower bound for trivial reject
  
  CKL_REAL d = BV_Distance(R, T, o1->child(0), o2->child(0));
  
  if(d <= res->tolerance)
  {
    // more work needed - choose routine according to queue size
    
    if(qsize <= 2)
    {
      ToleranceRecurse(res, R, T, o1, 0, o2, 0);
    }
    else
    {
      res->qsize = qsize;
      ToleranceQueueRecurse(res, R, T, o1, 0, o2, 0);
    }
  }
  
  // res->p2 is in cs 1 ; transform it to cs 2
  
  CKL_REAL u[3];
  VmV(u, res->p2, res->T);
  MTxV(res->p2, res->R, u);
  
  double time2 = GetTime();
  res->query_time_secs = time2 - time1;
  
  return CKL_OK;
}

#endif


Transform::Transform(const CKL_REAL R_[3][3], const CKL_REAL T_[3])
{
  McM(R, R_);
  VcV(T, T_);
}

Transform::Transform()
{
  Videntity(T);
  Midentity(R);
}


int CKL_ContinuousCollide(CKL_ContinuousCollideResult *result,
                          CKL_REAL R11[3][3], CKL_REAL T11[3], CKL_REAL R12[3][3], CKL_REAL T12[3], CKL_Model *o1,
                          CKL_REAL R21[3][3], CKL_REAL T21[3], CKL_REAL R22[3][3], CKL_REAL T22[3], CKL_Model *o2,
                          int N)
{
  // convert matrix to quaternion
  CKL_REAL quat11[4], quat12[4], quat21[4], quat22[4];
  QRotM(quat11, R11);
  QRotM(quat12, R12);
  QRotM(quat21, R21);
  QRotM(quat22, R22);

  // temporary matrices, vertors and quaternions for interpolation
  CKL_REAL R1[3][3];
  CKL_REAL T1[3];
  CKL_REAL quat1[4];

  CKL_REAL R2[3][3];
  CKL_REAL T2[3];
  CKL_REAL quat2[4];

  // delta translation
  CKL_REAL delta_T1[3];
  CKL_REAL delta_T2[3];
  VmV(delta_T1, T12, T11);
  VmV(delta_T2, T22, T21);

  
  for(std::size_t i = 0; i < N; ++i)
  {
    CKL_REAL t = i / (CKL_REAL)(N - 1);

    // rotation interpolation
    Slerp(quat1, quat11, quat12, t);
    Slerp(quat2, quat21, quat22, t);

    // quaternion to matrix
    MRotQ(R1, quat1);
    MRotQ(R2, quat2);


    // translation interpolation
    VRay(T1, T11, delta_T1, t);
    VRay(T2, T21, delta_T2, t);

    CKL_CollideResult cresult;
    CKL_Collide(&cresult,
                R1, T1, o1,
                R2, T2, o2,
                CKL_FIRST_CONTACT);
    
    if(cresult.NumPairs() > 0)
    {
      result->is_collide = true;
      result->time_of_contact = t;
      McM(result->R1, R1);
      VcV(result->T1, T1);
      McM(result->R2, R2);
      VcV(result->T2, T2);

      return CKL_OK;
    }
  }

  result->is_collide = false;
  result->time_of_contact = 1.001;

  return CKL_OK;
}

void sampleSE3Euler(Vecnf<6>& v, CKL_REAL trans_lower[3], CKL_REAL trans_upper[3])
{
  CKL_REAL randv0 = (CKL_REAL)(rand() / double(RAND_MAX));
  CKL_REAL randv1 = (CKL_REAL)(rand() / double(RAND_MAX));
  CKL_REAL randv2 = (CKL_REAL)(rand() / double(RAND_MAX));

  v[0] = trans_lower[0] + randv0 * (trans_upper[0] - trans_lower[0]);
  v[1] = trans_lower[1] + randv1 * (trans_upper[1] - trans_lower[1]);
  v[2] = trans_lower[2] + randv2 * (trans_upper[2] - trans_lower[2]);

  CKL_REAL quat[4];
  CKL_REAL randq0 = (CKL_REAL)(rand() / double(RAND_MAX));
  CKL_REAL randq1 = (CKL_REAL)(rand() / double(RAND_MAX));
  CKL_REAL randq2 = (CKL_REAL)(rand() / double(RAND_MAX));

  CKL_REAL x0 = randq0;
  
  CKL_REAL r1 = sqrt(1 - x0), r2 = sqrt(x0);
  CKL_REAL t1 = 2 * M_PI * randq1, t2 = 2 * M_PI * randq2;
  CKL_REAL c1 = cos(t1), s1 = sin(t1);
  CKL_REAL c2 = cos(t2), s2 = sin(t2);

  quat[0] = s1 * r1;
  quat[1] = c1 * r1;
  quat[2] = s2 * r2;
  quat[3] = c2 * r2;

  CKL_REAL euler[3];
  EulerRotQuat(euler, quat);

  v[3] = euler[0];
  v[4] = euler[1];
  v[5] = euler[2];
}


void ModelRadiusAndCenter(CKL_REAL& radius, CKL_REAL center[3], const CKL_Model* o)
{
  CKL_REAL bmin[3], bmax[3];
  bmin[0] = bmin[1] = bmin[2] = std::numeric_limits<CKL_REAL>::max();
  bmax[0] = bmax[1] = bmax[2] = -std::numeric_limits<CKL_REAL>::max();

  for(int i = 0; i < o->num_tris; ++i)
  {
    Tri tri = o->tris[i];
    for(int j = 0; j < 3; ++j)
    {
      if(tri.p1[j] < bmin[j]) bmin[j] = tri.p1[j];
      if(tri.p1[j] > bmax[j]) bmax[j] = tri.p1[j];
    }

    for(int j = 0; j < 3; ++j)
    {
      if(tri.p2[j] < bmin[j]) bmin[j] = tri.p2[j];
      if(tri.p2[j] > bmax[j]) bmax[j] = tri.p2[j];
    }

    for(int j = 0; j < 3; ++j)
    {
      if(tri.p3[j] < bmin[j]) bmin[j] = tri.p3[j];
      if(tri.p3[j] > bmax[j]) bmax[j] = tri.p3[j];
    }
  }

  CKL_REAL diag[3];
  VmV(diag, bmax, bmin);
  radius = 0.5 * Vlength(diag);

  VInterp(center, bmin, bmax, 0.5);
}

CKL_REAL ModelVolume(const CKL_Model* o)
{
  CKL_REAL vol = 0;
  for(int i = 0; i < o->num_tris; ++i)
  {
    const Tri& tri = o->tris[i];
    CKL_REAL cross[3];
    VcrossV(cross, tri.p1, tri.p2);
    CKL_REAL d_six_vol = VdotV(cross, tri.p3);
    vol += d_six_vol;
  }

  return vol / 6;
}

void ModelCoM(CKL_REAL com[3], const CKL_Model* o)
{
  CKL_REAL vol = 0;
  Videntity(com);
  for(int i = 0; i < o->num_tris; ++i)
  {
    const Tri& tri = o->tris[i];
    CKL_REAL cross[3];
    VcrossV(cross, tri.p1, tri.p2);
    CKL_REAL d_six_vol = VdotV(cross, tri.p3);
    vol += d_six_vol;
    CKL_REAL ave[3];
    VpV(ave, tri.p1, tri.p2);
    VpV(ave, ave, tri.p3);
    VxS(ave, ave, d_six_vol);
    VpV(com, com, ave);
  }

  VxS(com, com, 1.0/(vol * 4));
}

void computeMomentofInertia(CKL_REAL C[3][3], const CKL_Model* o)
{
  for(int i = 0; i < 3; ++i)
    for(int j = 0; j < 3; ++j)
      C[i][j] = 0;

  CKL_REAL C_canonical[3][3] = {{1/60.0, 1/120.0, 1/120.0},
                                {1/120.0, 1/60.0, 1/120.0},
                                {1/120.0, 1/120.0, 1/60.0}};

  for(int i = 0; i < o->num_tris; ++i)
  {
    const Tri& tri = o->tris[i];
    CKL_REAL cross[3];
    VcrossV(cross, tri.p1, tri.p2);
    CKL_REAL d_six_vol = VdotV(cross, tri.p3);

    CKL_REAL A[3][3];
    VscMrow(A, tri.p1, tri.p2, tri.p3);
    CKL_REAL AT[3][3];
    MTcM(AT, A);
    CKL_REAL ATxC[3][3];
    CKL_REAL ATxCxA[3][3];
    MxM(ATxC, AT, C_canonical);
    MxM(ATxCxA, ATxC, A);
    MxS(ATxCxA, ATxCxA, d_six_vol);
    MpM(C, C, ATxCxA);    
  }

  CKL_REAL trace_C = C[0][0] + C[1][1] + C[2][2];

  C[0][0] = trace_C - C[0][0];
  C[0][1] = - C[0][1];
  C[0][2] = - C[0][2];

  C[1][0] = - C[1][0];
  C[1][1] = trace_C - C[1][1];
  C[1][2] = - C[1][2];

  C[2][0] = - C[2][0];
  C[2][1] = - C[2][1];
  C[2][2] = trace_C - C[2][2];
}

void ModelMomentOfInertiaRelatedToCoM(CKL_REAL C[3][3], const CKL_Model* o)
{
  computeMomentofInertia(C, o);
  CKL_REAL com[3];
  ModelCoM(com, o);
  CKL_REAL V = ModelVolume(o);
  
  C[0][0] -= V * (com[1] * com[1] + com[2] * com[2]);
  C[0][1] += V * com[0] * com[1];
  C[0][2] += V * com[0] * com[2];
  C[1][0] += V * com[1] * com[0];
  C[1][1] -= V * (com[0] * com[0] + com[2] * com[2]);
  C[1][2] += V * com[1] * com[2];
  C[2][0] += V * com[2] * com[0];
  C[2][1] += V * com[2] * com[1];
  C[2][2] -= V * (com[0] * com[0] + com[1] * com[1]);
}

// tf2 = tf1 * tf
Transform relativeTransform(const Transform& tf1, const Transform& tf2)
{
  Transform tf;
  CKL_REAL R1inv[3][3];
  MTcM(R1inv, tf1.R);
  MxM(tf.R, R1inv, tf2.R);
  MxV(tf.T, R1inv, tf2.T);
  VmV(tf.T, tf.T, tf1.T);
  return tf;
}

// tf2 = tf * tf1;
Transform relativeTransform2(const Transform& tf1, const Transform& tf2)
{
  Transform tf;
  MxMT(tf.R, tf2.R, tf1.R);
  CKL_REAL t1rotated[3];
  MxV(t1rotated, tf.R, tf1.T);
  VmV(tf.T, tf2.T, t1rotated);
  return tf;  
}


struct DefaultDistanceFunctor : public DistanceFunctor<Transform>
{
public:
  const CKL_Model* o;

  CKL_REAL rot_x_weight, rot_y_weight, rot_z_weight;

  DefaultDistanceFunctor()
  {
    o = NULL;
    rot_x_weight = rot_y_weight = rot_z_weight = 1;
  }

  DefaultDistanceFunctor(const CKL_Model* o_)
  {
    o = o_;
    init();
  }

  void init()
  {
    rot_x_weight = rot_y_weight = rot_z_weight = 1;

    if(o)
    {
      CKL_REAL V = ModelVolume(o);
      CKL_REAL C[3][3];
      ModelMomentOfInertiaRelatedToCoM(C, o);
      CKL_REAL eigen_values[3];
      CKL_REAL eigen_vectors[3][3];

      Meigen(eigen_vectors, eigen_values, C);
      rot_x_weight = eigen_values[0] * 4 / V;
      rot_y_weight = eigen_values[1] * 4 / V;
      rot_z_weight = eigen_values[2] * 4 / V;

      std::cout << "rotation weight " << rot_x_weight << " " << rot_y_weight << " " << rot_z_weight << std::endl;
    }
  }
  
  virtual double dist(const Transform& tf1, const Transform& tf2) const
  {
    Transform tf = relativeTransform2(tf1, tf2);
    CKL_REAL quat[4];
    QRotM(quat, tf.R);
    CKL_REAL d =
        rot_x_weight * quat[1] * quat[1]
      + rot_y_weight * quat[2] * quat[2]
      + rot_z_weight * quat[3] * quat[3]
      + tf.T[0] * tf.T[0]
      + tf.T[1] * tf.T[1]
      + tf.T[2] * tf.T[2];

    return d;
  }

  
};


std::vector<Transform> CKL_PenetrationDepthModelLearning(CKL_Model* o1, CKL_Model* o2,
                                                         std::size_t n_samples, std::size_t knn_k,
                                                         CKL_REAL margin)
{
  CKL_REAL r1, r2;
  CKL_REAL c1[3], c2[3];
  ModelRadiusAndCenter(r1, c1, o1);
  ModelRadiusAndCenter(r2, c2, o2);
  CKL_REAL r = r1 + r2 + margin;

  CKL_REAL trans_lower[3] = {-r, -r, -r};
  CKL_REAL trans_upper[3] = {r, r, r};

  std::vector<Item<6> > data(n_samples);
  for(std::size_t i = 0; i < n_samples; ++i)
  {
    Vecnf<6> q;
    sampleSE3Euler(q, trans_lower, trans_upper);

    CKL_REAL trans[3] = {q[0], q[1], q[2]};
    CKL_REAL euler[3] = {q[3], q[4], q[5]};

    CKL_REAL R[3][3];
    CKL_REAL T[3];
    CKL_REAL rotated_c2[3];
    
    MRotEuler(R, euler);
    VpV(trans, trans, c1);
    MxV(rotated_c2, R, c2);
    VmV(T, trans, rotated_c2);

    CKL_CollideResult cresult;
    CKL_REAL R1[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    CKL_REAL T1[3] = {0, 0, 0};
    CKL_Collide(&cresult,
                R1, T1, o1,
                R, T, o2,
                CKL_FIRST_CONTACT);

    data[i] = Item<6>(q, (cresult.NumPairs() > 0));
  }

  Classifier<6> classifier;
  classifier.setScaler(computeScaler(data));

  classifier.learn(data);

  classifier.save("model.txt");


  std::vector<Transform> support_transforms_positive;
  std::vector<Transform> support_transforms_negative;
  
  std::vector<Item<6> > svs = classifier.getSupportVectors();
  for(std::size_t i = 0; i < svs.size(); ++i)
  {
    const Vecnf<6>&q = svs[i].q;

    CKL_REAL trans[3] = {q[0], q[1], q[2]};
    CKL_REAL euler[3] = {q[3], q[4], q[5]};
    
    CKL_REAL R[3][3];
    CKL_REAL T[3];
    CKL_REAL rotated_c2[3];
    
    MRotEuler(R, euler);
    VpV(trans, trans, c1);
    MxV(rotated_c2, R, c2);
    VmV(T, trans, rotated_c2);
    
    Transform tf(R, T);    
    if(svs[i].label)
      support_transforms_positive.push_back(tf);
    else
      support_transforms_negative.push_back(tf);
  }


  DefaultDistanceFunctor distance_func(o2);
  NearestNeighbors<Transform> knn_solver_positive;
  NearestNeighbors<Transform> knn_solver_negative;

  knn_solver_positive.setDistanceFunctor(&distance_func);
  knn_solver_negative.setDistanceFunctor(&distance_func);
  
  knn_solver_positive.add(support_transforms_positive);
  knn_solver_negative.add(support_transforms_negative);

  std::vector<Transform> contact_vectors;
  Transform tf_identity;  
  
  for(std::size_t i = 0; i < support_transforms_positive.size(); ++i)
  {
    std::vector<Transform> nbh;
    knn_solver_negative.nearestK(support_transforms_positive[i], knn_k, nbh);

    for(std::size_t j = 0; j < nbh.size(); ++j)
    {
      CKL_ContinuousCollideResult result;
      CKL_ContinuousCollide(&result,
                            tf_identity.R, tf_identity.T, tf_identity.R, tf_identity.T, o1,
                            nbh[j].R, nbh[j].T, support_transforms_positive[i].R, support_transforms_positive[i].T, o2,
                            100);

      // std::cout << result.time_of_contact << std::endl;
      
      contact_vectors.push_back(Transform(result.R2, result.T2));
    }
  }

  for(std::size_t i = 0; i < support_transforms_negative.size(); ++i)
  {
    std::vector<Transform> nbh;
    knn_solver_positive.nearestK(support_transforms_negative[i], knn_k, nbh);

    for(std::size_t j = 0; j < nbh.size(); ++j)
    {
      CKL_ContinuousCollideResult result;
      CKL_ContinuousCollide(&result,
                            tf_identity.R, tf_identity.T, tf_identity.R, tf_identity.T, o1,
                            support_transforms_negative[i].R, support_transforms_negative[i].T, nbh[j].R, nbh[j].T, o2,
                            100);

      // std::cout << result.time_of_contact << std::endl;

      contact_vectors.push_back(Transform(result.R2, result.T2));
    }
  }

  return contact_vectors;
}



int CKL_PenetrationDepth(CKL_PenetrationDepthResult* res,
                         CKL_REAL R1[3][3], CKL_REAL T1[3], CKL_Model* o1,
                         CKL_REAL R2[3][3], CKL_REAL T2[3], CKL_Model* o2,
                         const std::vector<Transform>& contact_space)
{
  DefaultDistanceFunctor distance_func(o2);
  NearestNeighbors<Transform> knn_solver;
  knn_solver.setDistanceFunctor(&distance_func);
  knn_solver.add(contact_space);

  Transform tf1(R1, T1), tf2(R2, T2);
  Transform tf = relativeTransform2(tf1, tf2);
  Transform tf_nearest = knn_solver.nearest(tf);
  
  MxM(res->R, tf1.R, tf_nearest.R);
  MxV(res->T, tf1.R, tf_nearest.T);
  VpV(res->T, res->T, tf1.T);

  res->pd_value = sqrt(distance_func.dist(Transform(res->R, res->T), tf2));

  return CKL_OK;
}






}
