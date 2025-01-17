---------------------------------------------------------------------------


                                CKL v. 1.0

                   UNC - Chapel Hill Computer Science
                           geom@cs.unc.edu

---------------------------------------------------------------------------
Changes:

1.0   - first release of library 
---------------------------------------------------------------------------


I. Introduction

 CKL, which stands for Collision Kernel Package, is a library for three
 types of proximity queries performed on geometric models composed of
 triangles:

   * collision detection - detect whether two models overlap, and
                           optionally, which triangles of the models
                           overlap.

   * distance computation - compute the distance between two models, 
                            i. e., the length of the shortest translation
                            that makes the models overlap

   * tolerance verification - detect whether two models are closer or
                              farther than a tolerance value.

 By default, the library uses "RSS" bounding volumes for distance and
 tolerance queries, and OBBs for collision detection (see CKL_Compile.h).
 Descriptions of the bounding volumes and algorithms used in this package 
 are contained in:

   Eric Larsen, Stefan Gottschalk, Ming Lin, Dinesh Manocha,
   "Fast Proximity Queries with Swept Sphere Volumes", 
   Technical report TR99-018, Department of Computer Science, 
   UNC Chapel Hill 

   S. Gottschalk, M. C. Lin and D. Manocha,
   "OBB-Tree: A Hierarchical Structure for Rapid Interference Detection",
   Technical report TR96-013, Department of Computer Science, University
   of N. Carolina, Chapel Hill.
   Proc. of ACM Siggraph'96.

II. Layout of Files

 CKL_v1.3/
   Makefile           Unix makefile to build CKL library

   src/
     CKL source

   lib/             
     libCKL.a         after Unix compilation
     CKL.lib          after Win32 compilation

   include/
     CKL.h            include this file to use CKL classes and functions.
     CKL_Internal.h   
     CKL_Compile.h    *WARNING* you should only modify CKL_Compile.h in
     Tri.h            the src directory, not here, because these files
     BV.h             are copied here from src when you perform a build
                     
   demos/
     Makefile         Unix makefile for both demos
 
     falling/         source and project files
     sample/          "      "   "       "
     spinning/        "      "   "       "

III. Building the CKL Library

 In the top level directory, there is a Unix Makefile for building the CKL
 library.  Type 'make' to create a 'libCKL.a' in the lib directory.
 The compiler is currently set to g++ with -O2 level optimization. 

 Building on either platform has a side effect of copying the include
 files needed for a client application to the include/ directory.  

IV. Building the Demos

 In the demos directory is a Unix Makefile.  Typing 'make' will perform a
 'make' in the 'sample' and 'spinning' directories.  

   sample

   This demo is adapted from the sample client included with RAPID.  Two 
   tori are created, and proximity queries are performed on them at  
   several configurations

   spinning

   The spinning demo is a GLUT application, so paths to the GLUT & OpenGL
   libraries and includes must be set in spinning/Makefile, or in the 
   VC++ project settings. When run, a bunny and a torus should appear in
   the GLUT window, with a line drawn between their closest points.
   Pressing a key alternately starts and stops them spinning.

   falling

   This demo is also a GLUT application, showing a bent torus
   falling through the center of a knobby torus.  Each of the three 
   proximity query types can be demonstrated.

V. Creating a CKL Client Application

 "CKL.h" contains the most complete information on constructing client
 applications.  Here is a summary of the steps involved.

 1. Include the CKL API header.
 
    #include "CKL.h"

 2. Create two instances of CKL_Model.

    CKL_Model m1, m2;

 3. Specify the triangles of each CKL_Model.

    Note that CKL uses the CKL_REAL type for all its floating point 
    values. This can be set in "CKL_Compile.h", and is "double" by 
    default

    // begin m1

    m1.BeginModel();

    // create some triangles

    CKL_REAL p1[3], p2[3], p3[3];  
    CKL_REAL q1[3], q2[3], q3[3];
    CKL_REAL r1[3], r2[3], r3[3];

    // initialize the points
     .
     . 
     .  

    // add triangles that will belong to m1
 
    m1.AddTri(p1, p2, p3, 0);
    m1.AddTri(q1, q2, q3, 1);
    m1.AddTri(r1, r2, r3, 2);
 
    // end m1, which builds the model

    m1.EndModel();

 4. Specify the orientation and position of each model.

    The position of a model is specified as a 3 vector giving the
    position of its frame in the world, stored in a CKL_REAL [3]. 

    The rotation for a model is specified as a 3x3 matrix, whose columns
    are the model frame's basis vectors, stored in row major order in
    a CKL_REAL [3][3];

    Note that an OpenGL 4x4 matrix has column major storage.

 5. Perform any of the three proximity queries.

    // collision

    CKL_CollideResult cres;
    CKL_Collide(&cres,R1,T1,&m1,R2,T2,&m2);

    // distance

    CKL_DistanceResult dres;
    double rel_err = 0.0, abs_err = 0.0;
    CKL_Distance(&dres,R1,T1,&m1,R2,T2,&m2,rel_err,abs_err);

    // tolerance
 
    CKL_ToleranceResult tres;
    double tolerance = 1.0;
    CKL_Tolerance(&tres,R1,T1,&m1,R2,T2,&m2,tolerance);

    See "CKL.h" for complete information.

 6. Access the result structure passed in the query call.

    int colliding = cres.Colliding();
    double distance = dres.Distance();
    int closer = tres.CloserThanTolerance();

    See "CKL.h" for the complete interface to each result structure.

