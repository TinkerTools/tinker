/*
 * alf_tools_gmp.c   Version 2   3/30/2007   Patrice Koehl
 *
 * This is the C version of alfcx_tools.f, which performs all
 * operations with multi precision arithmetic using the GMP package
 *
 * Copyright (C) 2002 Patrice Koehl
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "uniongmp.h"

/*
 * Naming convention between C and Fortran routines:
 *
 * Let us consider two Fortran subroutines, foo and foo_with_underscore
 *
 * Case 1: the name of the Fortran subroutine (foo) does not contain
 *    an underscore (_)
 *
 * In a C program, we must add one underscore to the Fortran name:
 *    the C program refers to the subroutine as:  foo_(...)
 *    while the Fortran code writes it as:  foo(...)
 * This is independent of the compiler pair, at least for the
 *    gcc/gfortran and Intel icc/ifort compilers
 *
 * Case 2: the name of the Fortran subroutine (foo_with_underscore)
 *    contains at least one underscore (_)
 *
 * Treatment of this case is compiler dependent
 * The Intel compiler treats this case as if it was case 1, i.e., needs
 *    one underscore at the end of the Fortran name in the C program
 * However, GNU gfortran requires two underscores at the end of the
 *    Fortran name in the C program
 *
 * Here this is solved by introducing two functions, FTName1 and FTName2,
 *    where FTName2 is compiler dependent
 */

#define F77Name1(x) x##_    /* Add one underscore to Fortran program */

#if defined intel
#define F77Name2(x) x##_    /* Add one underscore to Fortran program
                                for Intel compiler */
#else
#define F77Name2(x) x##__   /* Add two underscores to Fortran program
                                  for GNU compiler */
#endif

/*------------------------------------------------------------------------*/

void set_edge(double *a, double *b, double ra, double rb, double scale);

void set_triangle(double *a, double *b, double *c, double ra, double rb, 
      double rc, double scale);

void real_to_gmp(double *coord, double scale, int idx, mpz_t val);

void build_weight(mpz_t ax, mpz_t ay, mpz_t az, mpz_t r, mpz_t w);

void scalar_to_gmp(double coord, double scale, mpz_t val);

#ifdef __cplusplus
extern "C"
#endif
void F77Name2(tetra_radius_gmp)(double *a, double *b, double *c,
              double *d, double *ra, double *rb, double *rc, double *rd,
              int *test, double *scale, double *alpha);

#ifdef __cplusplus
extern "C"
#endif
void F77Name2(vertex_attach_gmp)(double *a, double *b, double *ra,
              double *rb, double *scale, int *testa, int *testb);

#ifdef __cplusplus
extern "C"
#endif
void F77Name2(edge_attach_gmp)(double *a, double *b, double *c, double *ra,
              double *rb, double *rc, double *scale, int *test, int *memory);

#ifdef __cplusplus
extern "C"
#endif
void F77Name2(edge_radius_gmp)(double *a, double *b, double *ra, double *rb,
              int *test, double *scale, double *alpha, int *memory);

#ifdef __cplusplus
extern "C"
#endif
void F77Name2(triangle_attach_gmp)(double *a, double *b, double *c,
              double *d, double *ra, double *rb, double *rc, double *rd,
              double *scale, int *test, int *memory);

#ifdef __cplusplus
extern "C"
#endif
void F77Name2(triangle_radius_gmp)(double *a, double *b, double *c,
              double *ra, double *rb, double *rc, int *test, double *scale,
              double *alpha, int *memory);

void F77Name2(set_alf_gmp)();
void F77Name2(clear_alf_gmp)();

/*------------------------------------------------------------------------*/

/* scalar_to_gmp:                                                   */
/* This subroutine converts one coordinate of one point into mpz_t type
 * */

void scalar_to_gmp(double coord, double scale, mpz_t val)
{
        double x;
        int ivalue;

        mpz_set_d(temp3, scale);

        ivalue= (int) coord;
        mpz_set_si(temp1,ivalue);
        mpz_mul(temp1,temp1,temp3);
        x = (coord-ivalue)*(scale);
        ivalue = (int) round(x);
        mpz_set_si(temp2,ivalue);
        mpz_add(val,temp1,temp2);
}

/*------------------------------------------------------------------------*/
/* build_weight:
   This subroutine builds the weight of a point:
   w = x**2 + y**2 + z**2 - ra**2
                           */

void build_weight(mpz_t ax, mpz_t ay, mpz_t az, mpz_t r, mpz_t w)
{
   mpz_mul(temp1,r,r);
   mpz_mul(temp2,ax,ax), mpz_sub(temp1,temp2,temp1);
   mpz_mul(temp2,ay,ay), mpz_add(temp1,temp2,temp1);
   mpz_mul(temp2,az,az), mpz_add(w,temp2,temp1);
}

/*------------------------------------------------------------------------*/
/* set_triangle:
   This subroutine sets all common arrays for gmp calculation over edges
   Input:
      ia,ib: indices of the two points considered
*/
void set_triangle(double *a, double *b, double *c, double ra, double rb, double rc, 
                    double scale)
{
   int i,j;

//       0. define coordinates


   for (i=0; i<3; i++)
   {
      real_to_gmp(a,scale,i,a_mp[i+1]);
      real_to_gmp(b,scale,i,b_mp[i+1]);
      real_to_gmp(c,scale,i,c_mp[i+1]);
   }
   
   scalar_to_gmp(ra,scale,ra_mp);
   scalar_to_gmp(rb,scale,rb_mp);
   scalar_to_gmp(rc,scale,rc_mp);

   build_weight(a_mp[1],a_mp[2],a_mp[3],ra_mp,a_mp[4]);
   build_weight(b_mp[1],b_mp[2],b_mp[3],rb_mp,b_mp[4]);
   build_weight(c_mp[1],c_mp[2],c_mp[3],rc_mp,c_mp[4]);

/*
       1. Computes all Minors Mab(i,j)= M(a,b,i,j)   = Det | a(i)  a(j) |
                                                           | b(i)  b(j) |
*/
   for (i=1;  i<4; i++)
   {
      for (j=i+1; j<5 ; j++)
      {
         mpz_mul(temp1,a_mp[j],b_mp[i]); 
         mpz_mul(temp2,a_mp[i],b_mp[j]);
         mpz_sub(Mab[i][j],temp2,temp1);
         mpz_mul(temp1,a_mp[j],c_mp[i]); 
         mpz_mul(temp2,a_mp[i],c_mp[j]);
         mpz_sub(Mac[i][j],temp2,temp1);
         mpz_mul(temp1,b_mp[j],c_mp[i]); 
         mpz_mul(temp2,b_mp[i],c_mp[j]);
         mpz_sub(Mbc[i][j],temp2,temp1);
      }
   }

/*
       Now compute all Minors
               S(i,j) = M(a,b,c,i,j,0)    = Det | a(i) a(j) 1 |
                                                | b(i) b(j) 1 |
                                                | c(i) c(j) 1 |

       a,b,c are the 3 vertices of the triangle, i and j correspond
       to two of the coordinates of the vertices

       for all i in [1,3] and all j in [i+1,4]
*/

   for (i=1;  i<4; i++)
   {
      for (j=i+1; j<5 ; j++)
      {
         mpz_sub(temp1,Mbc[i][j],Mac[i][j]);
         mpz_add(S[i][j],temp1,Mab[i][j]);
      }
   }

/*
       Now compute all Minors
               T(i,j) = M(a,b,c,i,j,4)    = Det | a(i) a(j) a(4) |
                                                | b(i) b(j) b(4) |
                                                | c(i) c(j) c(4) |

       for all i in [1,2] and all j in [i+1,3]
*/

   for (i=1;  i<3; i++)
   {
      for (j=i+1; j<4 ; j++)
      {
         mpz_mul(temp1,a_mp[4],Mbc[i][j]);
         mpz_mul(temp2,b_mp[4],Mac[i][j]);
         mpz_sub(temp1,temp1,temp2);
         mpz_mul(temp2,c_mp[4],Mab[i][j]);
         mpz_add(T[i][j],temp1,temp2);
      }
   }
/*
       Finally,  need Dabc = M(a,b,c,1,2,3) Det | a(1) a(2) a(3) |
                                                | b(2) b(2) b(3) |
                                                | c(3) c(2) c(3) |
*/
   mpz_mul(temp1,a_mp[1],Mbc[2][3]); mpz_mul(temp2,b_mp[1],Mac[2][3]); mpz_sub(temp1,temp1,temp2);
   mpz_mul(temp2,c_mp[1],Mab[2][3]); mpz_add(Dabc,temp1,temp2);

}
/*------------------------------------------------------------------------*/
/* set_edge:
   This subroutine sets all common arrays for gmp calculation over edges
   Input:
      ia,ib: indices of the two points considered
*/
void set_edge(double *a, double *b, double ra, double rb, double scale)
{
   int i,j,k;

//       0. define coordinates


   for (i=0; i<3; i++)
   {
      real_to_gmp(a,scale,i,a_mp[i+1]);
      real_to_gmp(b,scale,i,b_mp[i+1]);
   }
   
   scalar_to_gmp(ra,scale,ra_mp);
   scalar_to_gmp(rb,scale,rb_mp);

   build_weight(a_mp[1],a_mp[2],a_mp[3],ra_mp,a_mp[4]);
   build_weight(b_mp[1],b_mp[2],b_mp[3],rb_mp,b_mp[4]);

/*
       1. Compute all Minors Dab(i) = M(a,b,i,0) = Det | a(i) 1 |
                                                       | b(i) 1 |
*/
   for (i=1; i<5; i++)
   {
      mpz_sub(Dab[i],a_mp[i],b_mp[i]);
   }

/*
       2. Computes all Minors Sab(i,j)= M(a,b,i,j)   = Det | a(i)  a(j) |
                                                           | b(i)  b(j) |
*/
   for (i=1;  i<3; i++)
   {
      for (j=i+1; j<4 ; j++)
      {
         k=i+j-2;
         mpz_mul(temp1,a_mp[j],b_mp[i]); 
         mpz_mul(temp2,a_mp[i],b_mp[j]);
         mpz_sub(Sab[k],temp2,temp1);
      }
   }

/*
      3. Computes all Minors Tab(i)= M(a,b,i,4)   = Det | a(i)  a(4) |
                                                         | b(i)  b(4) |
*/
   for (i=1; i<4; i++)
   {
      mpz_mul(temp1,a_mp[i],b_mp[4]); 
      mpz_mul(temp2,a_mp[4],b_mp[i]);
      mpz_sub(Tab[i],temp1,temp2);
   }
}

/*------------------------------------------------------------------------*/
/* tetra_radius_gmp   Version 1 11/24/2000   Patrice Koehl
 
    This subroutine computes the radius R of the circumsphere containing
    a tetrahedron [A,B,C,D], as well as check if any fourth point L 
    (in A, B, C, D) of the tetrahedron is "hidden" by its opposite 
    face [I,J,K], i.e. is interior to the circumsphere of [I,J,K]

   The array "coord" contains the coordinates of A,B,C,D in that order
   The array "radius" contains the radii of A,B,C and D
 
    Since we are only interested at how R compares to Alpha, we don't
    output R, rather the result of the comparison
 
    Computation extends to all four faces of the tetrahedron, as
    well as to all six edges of the tetrahedron
 
   This procedure works with Multiple Precision Integer Arithmetics  (MPIA)
    The package GMP is used for MPIA (with a C wrapper)
*/
#ifdef __cplusplus
extern "C" {
#endif
void F77Name2(tetra_radius_gmp)( double *a, double *b, double *c, double *d, double *ra,
      double *rb, double *rc, double *rd, int *testr, double *scale, double *alpha)

{
   int i,j,k,coef;
   int ivalue;
   double value;

/*   Transfer data in multiple precision */

   for (i=0; i<3; i++)
   {
      real_to_gmp(a,*scale,i,a_mp[i+1]);
      real_to_gmp(b,*scale,i,b_mp[i+1]);
      real_to_gmp(c,*scale,i,c_mp[i+1]);
      real_to_gmp(d,*scale,i,d_mp[i+1]);
   }
   
   scalar_to_gmp(*ra,*scale,ra_mp);
   scalar_to_gmp(*rb,*scale,rb_mp);
   scalar_to_gmp(*rc,*scale,rc_mp);
   scalar_to_gmp(*rd,*scale,rd_mp);

   build_weight(a_mp[1],a_mp[2],a_mp[3],ra_mp,wa);
   build_weight(b_mp[1],b_mp[2],b_mp[3],rb_mp,wb);
   build_weight(c_mp[1],c_mp[2],c_mp[3],rc_mp,wc);
   build_weight(d_mp[1],d_mp[2],d_mp[3],rd_mp,wd);

   value = (*alpha)*(*scale); ivalue = (int) floor(value); 
   mpz_set_si(alp,ivalue);

/*   1. Computes all Minors Smn(i+j-2)= M(m,n,i,j) = Det | m(i)  m(j) |
                              | n(i)  n(j) |
   for all i in [1,2] and all j in [i+1,3]                         */

   for (i=1;  i<3; i++)
   {
      for (j=i+1; j<4 ; j++)
      {
         k=i+j-2;
         mpz_mul(temp1,a_mp[j],b_mp[i]); 
         mpz_mul(temp2,a_mp[i],b_mp[j]);
         mpz_sub(Sab[k],temp2,temp1);
         mpz_mul(temp1,a_mp[j],c_mp[i]); 
         mpz_mul(temp2,a_mp[i],c_mp[j]);
         mpz_sub(Sac[k],temp2,temp1);
         mpz_mul(temp1,a_mp[j],d_mp[i]); 
         mpz_mul(temp2,a_mp[i],d_mp[j]);
         mpz_sub(Sad[k],temp2,temp1);
         mpz_mul(temp1,b_mp[j],c_mp[i]); 
         mpz_mul(temp2,b_mp[i],c_mp[j]);
         mpz_sub(Sbc[k],temp2,temp1);
         mpz_mul(temp1,b_mp[j],d_mp[i]); 
         mpz_mul(temp2,b_mp[i],d_mp[j]);
         mpz_sub(Sbd[k],temp2,temp1);
         mpz_mul(temp1,c_mp[j],d_mp[i]); 
         mpz_mul(temp2,c_mp[i],d_mp[j]);
         mpz_sub(Scd[k],temp2,temp1);
      }
   }

/*   Now compute all Minors 
      Sq(i+j-2) = M(m,n,p,i,j,0) = Det | m(i) m(j) 1 |
                               | n(i) n(j) 1 |
                   | p(i) p(j) 1 |

   and all Minors
      Det(i+j-2) = M(m,n,p,q,i,j,4,0) = Det | m(i) m(j) m(4) 1 |
                        | n(i) n(j) n(4) 1 |
                        | p(i) p(j) p(4) 1 |
                        | q(i) q(j) q(4) 1 |

   m,n,p,q are the four vertices of the tetrahedron, i and j correspond
   to two of the coordinates of the vertices, and m(4) refers to the
   "weight" of vertices m                                           */
 
   for (i=1; i<4; i++)
   {
      mpz_sub(temp1,Scd[i],Sbd[i]); mpz_add(Sa[i],temp1,Sbc[i]);
      mpz_mul(temp2,Sa[i],wa);
      mpz_sub(temp1,Scd[i],Sad[i]); mpz_add(Sb[i],temp1,Sac[i]);
      mpz_mul(temp3,Sb[i],wb); mpz_sub(temp2,temp2,temp3);
      mpz_sub(temp1,Sbd[i],Sad[i]); mpz_add(Sc[i],temp1,Sab[i]);
      mpz_mul(temp3,Sc[i],wc); mpz_add(temp2,temp2,temp3);
      mpz_sub(temp1,Sbc[i],Sac[i]); mpz_add(Sd[i],temp1,Sab[i]);
      mpz_mul(temp3,Sd[i],wd); mpz_sub(Deter[i],temp2,temp3);
      mpz_neg(Sam1[i],Sa[i]); mpz_neg(Sbm1[i],Sb[i]);
      mpz_neg(Scm1[i],Sc[i]); mpz_neg(Sdm1[i],Sd[i]);
   }
 
/*
   Now compute the determinant needed to compute the radius of the
   circumsphere of the tetrahedron :

      Det1 = Minor(a,b,c,d,4,2,3,0)
      Det2 = Minor(a,b,c,d,1,3,4,0)
      Det3 = Minor(a,b,c,d,1,2,4,0)
      Det4 = Minor(a,b,c,d,1,2,3,0)
                           */

   mpz_set(det1,Deter[3]);
   mpz_set(det2,Deter[2]);
   mpz_set(det3,Deter[1]);

   mpz_mul(temp1,a_mp[1],Sa[3]);mpz_mul(temp2,b_mp[1],Sb[3]);
   mpz_sub(temp3,temp1,temp2);
   mpz_mul(temp1,c_mp[1],Sc[3]);mpz_mul(temp2,d_mp[1],Sd[3]);
   mpz_sub(temp1,temp1,temp2);
   mpz_add(det4,temp1,temp3);

/*
   Now compute all minors:
      Dmnp = Minor(m,n,p,1,2,3) = Det | m(1) m(2) m(3) |
                  | n(1) n(2) n(3) |
                  | p(1) p(2) p(3) |
                           */

   mpz_mul(temp1,a_mp[1],Sbc[3]); mpz_mul(temp2,b_mp[1],Sac[3]);
   mpz_sub(temp3,temp1,temp2);
   mpz_mul(temp1,c_mp[1],Sab[3]);mpz_add(Dabc,temp3,temp1);

   mpz_mul(temp1,a_mp[1],Sbd[3]); mpz_mul(temp2,b_mp[1],Sad[3]);
   mpz_sub(temp3,temp1,temp2);
   mpz_mul(temp1,d_mp[1],Sab[3]);mpz_add(Dabd,temp3,temp1);

   mpz_mul(temp1,a_mp[1],Scd[3]); mpz_mul(temp2,c_mp[1],Sad[3]);
   mpz_sub(temp3,temp1,temp2);
   mpz_mul(temp1,d_mp[1],Sac[3]);mpz_add(Dacd,temp3,temp1);

   mpz_mul(temp1,b_mp[1],Scd[3]); mpz_mul(temp2,c_mp[1],Sbd[3]);
   mpz_sub(temp3,temp1,temp2);
   mpz_mul(temp1,d_mp[1],Sbc[3]);mpz_add(Dbcd,temp3,temp1);

/*
   We also need :
      Det = Det | m(1) m(2) m(3) m(4) |
           | n(1) n(2) n(3) n(4) |
           | p(1) p(2) p(3) p(4) |
           | q(1) q(2) q(3) q(4) |
                        */

   mpz_mul(temp1,wa,Dbcd); mpz_mul(temp2,wb,Dacd);
   mpz_sub(temp3,temp2,temp1);
   mpz_mul(temp1,wc,Dabd); mpz_mul(temp2,wd,Dabc);
   mpz_sub(temp1,temp2,temp1); mpz_add(Dabcd,temp3,temp1);

/*
   The radius of the circumsphere of the weighted tetrahedron is then:
   r_t = (Det1*Det1 + Det2*Det2 + Det3*Det3 + 4*Det4*Dabcd)/(4*Det4*Det4)
                        */

   mpz_mul(temp1,det4,det4); coef=4; mpz_mul_si(den,temp1,coef);

   mpz_mul(temp1,det1,det1); mpz_mul(temp2,det2,det2);
   mpz_add(temp1,temp1,temp2); mpz_mul(temp2,det3,det3);
   mpz_add(temp1,temp1,temp2); mpz_mul(temp2,det4,Dabcd);
   mpz_mul_si(temp2,temp2,coef); mpz_add(num,temp1,temp2);
   
   mpz_mul(temp1,den,alp); mpz_sub(temp2,num,temp1);
   
/* 
    If tetrahedron is part of the alpha shape, then the 4 triangles,
    the 6 edges and the four vertices are also part of the alpha
    complex
                         */
   if(!(mpz_sgn(temp2) > 0)) 
   {
      (*testr)=1;
   }
   else
   {
      (*testr)=0;
   }
}
#ifdef __cplusplus
}
#endif

/*------------------------------------------------------------------------*/
/* vertex_attach_gmp:
   This subroutine checks if a vertex is attached to another vertex 
   Input:
      a, b   : coordinates of the two points
      ra,rb   : radii of the two points
   Output:
      testa   : flag equal to 1 if a is attached to b
      testb   : flag equal to 1 if b is attached to a
*/
#ifdef __cplusplus
extern "C" {
#endif
void   F77Name2(vertex_attach_gmp)(double *a, double *b, double *ra, double *rb, 
   double *scale, int *testa, int *testb)

{
   int i;

   (*testa = 0);
   (*testb = 0);

   for (i=0; i<4; i++)
   {
      real_to_gmp(a, *scale, i,temp1);
      real_to_gmp(b, *scale, i,temp2);
      mpz_sub(Dab[i],temp1,temp2);
   }

   scalar_to_gmp(*ra, *scale, ra_mp);
   scalar_to_gmp(*rb, *scale, rb_mp);

   mpz_mul(temp1, Dab[0],Dab[0]);
   mpz_mul(temp2, Dab[1],Dab[1]);
   mpz_mul(temp3, Dab[2],Dab[2]);
   mpz_add(temp1, temp1,temp2); mpz_add(dist2,temp1,temp3);

   mpz_mul(ra2, ra_mp, ra_mp); mpz_mul(rb2, rb_mp, rb_mp);

   mpz_add(dtest,dist2,ra2);
   mpz_sub(dtest,dtest,rb2);
   if(mpz_sgn(dtest) < 0) (*testa = 1);

   mpz_sub(dtest,dist2,ra2);
   mpz_add(dtest,dtest,rb2);
   if(mpz_sgn(dtest) < 0) (*testb = 1);

}
#ifdef __cplusplus
}
#endif

/*------------------------------------------------------------------------*/
/* edge_attach_gmp:
   This subroutine checks if an edge of the regular triangulation is "attached"
   to another vertex (i.e.  if the vertex belongs to the smallest circumsphere
   of the edge). 

                           */

#ifdef __cplusplus
extern "C" {
#endif
void F77Name2(edge_attach_gmp)(double *a, double *b, double *c, double *ra, 
         double *rb, double *rc, double *scale, int *testa, int *memory)

{
   int i,j,k,coef;

/* Set up calculation, if not already done) */

   if((*memory) != 1) set_edge(a, b, *ra, *rb, *scale);

   for (i=0; i<3; i++)
   {
      real_to_gmp(c,*scale,i,c_mp[i+1]);
   }

   scalar_to_gmp(*rc,*scale,rc_mp);
   build_weight(c_mp[1],c_mp[2],c_mp[3],rc_mp,c_mp[4]);

/* 
       Need to compute:
       Sc      : minor(a,b,c,i,j,0) for i=1,2 and j = i+1,3
       Tc      : minor(a,b,c,i,4,0) for i = 1,2,3
*/

   for (i=1; i<3 ; i++)
   {
      for (j=i+1; j<4; j++)
      {
         k=i+j-2;
         mpz_mul(temp1,c_mp[i],Dab[j]);
         mpz_mul(temp2,c_mp[j],Dab[i]);
         mpz_sub(temp1,temp1,temp2);
         mpz_add(Sc[k],temp1,Sab[k]);
      }
   }

   for (i=1; i<4; i++)
   {
      mpz_mul(temp1,c_mp[i],Dab[4]);
      mpz_mul(temp2,c_mp[4],Dab[i]);
      mpz_sub(temp1,temp1,temp2);
      mpz_add(Tc[i],temp1,Tab[i]);
   }

/*   This is the "hidden1" part */

   (*testa) = 0;

   if( mpz_cmp(a_mp[1],b_mp[1]) != 0) 
   {
      for (i = 1; i < 4 ; i++)
      {
         mpz_set(res[0][i],Dab[i]);
         mpz_set(res2_c[i][4],Tc[i]);
      }
      mpz_set(res[1][2],Sab[1]); mpz_set(res[1][3],Sab[2]);
      mpz_set(res[2][3],Sab[3]);
      mpz_set(res2_c[1][2],Sc[1]); mpz_set(res2_c[1][3],Sc[2]);
      mpz_set(res2_c[2][3],Sc[3]);
   }
   else if ( mpz_cmp(a_mp[2],b_mp[2]) != 0)
   {
      mpz_set(res[0][1],Dab[2]); mpz_set(res[0][2],Dab[3]);
      mpz_set(res[0][3],Dab[1]);
      mpz_set(res[1][2],Sab[3]);
      mpz_neg(res[1][3],Sab[1]); mpz_neg(res[2][3],Sab[2]);
      mpz_set(res2_c[1][2],Sc[3]);
      mpz_neg(res2_c[1][3],Sc[1]); mpz_neg(res2_c[2][3],Sc[2]);
      mpz_set(res2_c[1][4],Tc[2]); mpz_set(res2_c[2][4],Tc[3]);
      mpz_set(res2_c[3][4],Tc[1]);
   }
   else if (  mpz_cmp(a_mp[3],b_mp[3]) != 0)
   {
      mpz_set(res[0][1],Dab[3]); mpz_set(res[0][2],Dab[1]);
      mpz_set(res[0][3],Dab[2]);
      mpz_neg(res[1][2],Sab[2]);
      mpz_neg(res[1][3],Sab[3]); mpz_set(res[2][3],Sab[1]);
      mpz_neg(res2_c[1][2],Sc[2]);
      mpz_neg(res2_c[1][3],Sc[3]); mpz_set(res2_c[2][3],Sc[1]);
      mpz_set(res2_c[1][4],Tc[3]); mpz_set(res2_c[2][4],Tc[1]);
      mpz_set(res2_c[3][4],Tc[2]);
   }
   else
   {
      exit(1);
   }

   mpz_mul(r_11,res[0][1],res[0][1]);
   mpz_mul(r_22,res[0][2],res[0][2]);
   mpz_mul(r_33,res[0][3],res[0][3]);
   mpz_mul(temp1,res[0][3],res[1][2]); 
   mpz_mul(temp2,res[0][2],res[1][3]);
   mpz_sub(diff,temp1,temp2);

/* Compute det0 */

   mpz_add(temp1,r_22,r_33); mpz_add(temp1,temp1,r_11); 
   mpz_mul(temp1,temp1,res[0][1]); 
   coef = -2;
   mpz_mul_si(det0,temp1,coef);


/* Now check if edge (ab) is attached to c */

   mpz_mul(temp1,res[1][2],res2_c[1][2]);
   mpz_mul(temp2,res[1][3],res2_c[1][3]);
   mpz_add(temp1,temp1,temp2);
   coef = -2;
   mpz_mul_si(temp1,temp1,coef);

   mpz_set_si(temp2,0);
   for (i=1; i<4; i++)
   {
      mpz_mul(temp3,res[0][i],res2_c[i][4]);
      mpz_add(temp2,temp2,temp3);
   }
   mpz_add(temp1,temp2,temp1);
   mpz_mul(temp1,temp1,res[0][1]);

   mpz_mul(temp2,res2_c[2][3],diff);
   mpz_mul_si(temp2,temp2,coef);
   mpz_sub(temp3,temp1,temp2);
   mpz_mul(dtest,temp3,det0);

   if(mpz_sgn(dtest) < 0) (*testa = 1);

}
#ifdef __cplusplus
}
#endif

/*------------------------------------------------------------------------*/
/* edge_radius_gmp:
   This subroutine checks if the radius of the smallest circumsphere of an edge
   of the regular triangulation is smaller than the value of alpha
   For that, it needs:

   Input:
      a,b   : coordinates of the two vertices defining the edge
      ra,rb   : radii of these two vertices
      Dab   : minor(a,b,i,0) for all i=1,2,3
      Sab   : minor(a,b,i,j) for i = 1,2 and j =i+1,3
   Ouput:
      testr   : flag that defines if radius < alpha or not
                             
For comments, see alfcx.f which contains the fortran equivalence of this
routine                                                                           */

#ifdef __cplusplus
extern "C" {
#endif
void F77Name2(edge_radius_gmp)(double *a, double *b, double *ra, double *rb, 
         int *testr, double *scale, double *alpha, int *memory)

{
   int i,coef, ivalue;
   double value;

   value = (*alpha)*(*scale); ivalue = (int) floor(value); 
   mpz_set_si(alp,ivalue);

   if((*memory) != 1) set_edge(a,b,*ra,*rb,*scale);

/*   This is the "hidden1" part */

   (*testr) = 0;
   mpz_set(res[0][4],Dab[4]);


   if( mpz_cmp(a_mp[1],b_mp[1]) != 0) 
   {
      for (i = 1; i < 4 ; i++)
      {
         mpz_set(res[0][i],Dab[i]);
         mpz_mul(temp1,b_mp[i],a_mp[4]);
         mpz_mul(temp2,a_mp[i],b_mp[4]);
         mpz_sub(res[i][4],temp2,temp1);
      }
      mpz_set(res[1][2],Sab[1]); mpz_set(res[1][3],Sab[2]);
      mpz_set(res[2][3],Sab[3]);
   }
   else if ( mpz_cmp(a_mp[2],b_mp[2]) != 0)
   {
      mpz_set(res[0][1],Dab[2]); mpz_set(res[0][2],Dab[3]);
      mpz_set(res[0][3],Dab[1]);
      mpz_set(res[1][2],Sab[3]);
      mpz_neg(res[1][3],Sab[1]); mpz_neg(res[2][3],Sab[2]);
      mpz_mul(temp1,a_mp[2],b_mp[4]); mpz_mul(temp2,b_mp[2],a_mp[4]);
      mpz_sub(res[1][4],temp1,temp2);
      mpz_mul(temp1,a_mp[3],b_mp[4]); mpz_mul(temp2,b_mp[3],a_mp[4]);
      mpz_sub(res[2][4],temp1,temp2);
      mpz_mul(temp1,a_mp[1],b_mp[4]); mpz_mul(temp1,b_mp[1],a_mp[4]);
      mpz_sub(res[3][4],temp1,temp2);
   }
   else if (  mpz_cmp(a_mp[3],b_mp[3]) != 0)
   {
      mpz_set(res[0][1],Dab[3]); mpz_set(res[0][2],Dab[1]);
      mpz_set(res[0][3],Dab[2]);
      mpz_neg(res[1][2],Sab[2]);
      mpz_neg(res[1][3],Sab[3]); mpz_set(res[2][3],Sab[1]);
      mpz_mul(temp1,a_mp[3],b_mp[4]); mpz_mul(temp2,b_mp[3],a_mp[4]);
      mpz_sub(res[1][4],temp1,temp2);
      mpz_mul(temp1,a_mp[1],b_mp[4]); mpz_mul(temp2,b_mp[1],a_mp[4]);
      mpz_sub(res[2][4],temp1,temp2);
      mpz_mul(temp1,a_mp[2],b_mp[4]); mpz_mul(temp1,b_mp[2],a_mp[4]);
      mpz_sub(res[3][4],temp1,temp2);
   }
   else
   {
      exit(1);
   }

   mpz_mul(r_11,res[0][1],res[0][1]);
   mpz_mul(r_22,res[0][2],res[0][2]);
   mpz_mul(r_33,res[0][3],res[0][3]);
   mpz_mul(r_14,res[0][1],res[0][4]);
   mpz_mul(r_313,res[0][3],res[1][3]);
   mpz_mul(r_212,res[0][2],res[1][2]);
   mpz_mul(temp1,res[0][3],res[1][2]); 
   mpz_mul(temp2,res[0][2],res[1][3]);
   mpz_sub(diff,temp1,temp2);

/* Compute det0 */

   mpz_add(temp1,r_22,r_33); mpz_add(temp1,temp1,r_11); 
   mpz_mul(temp1,temp1,res[0][1]); 
   coef = -2;
   mpz_mul_si(det0,temp1,coef);

/* Compute det1 */

   mpz_add(temp1,r_313,r_212);
   coef = 2;
   mpz_mul_si(temp1,temp1,coef);
   mpz_sub(temp1,temp1,r_14);
   mpz_mul(det1,res[0][1],temp1);


/* Compute det2 */

   mpz_add(temp1,r_11,r_33);
   mpz_mul(temp1,temp1,res[1][2]);
   coef = -2;
   mpz_mul_si(temp1,temp1,coef);
   mpz_mul_si(temp2,r_313,coef);
   mpz_add(temp2,temp2,r_14);
   mpz_mul(temp2,temp2,res[0][2]);
   mpz_sub(det2,temp1,temp2);

/* Compute det3 */

   mpz_add(temp1,r_11,r_22);
   mpz_mul(temp1,temp1,res[1][3]);
   mpz_mul_si(temp1,temp1,coef);
   mpz_mul_si(temp2,r_212,coef);
   mpz_add(temp2,temp2,r_14);
   mpz_mul(temp2,temp2,res[0][3]);
   mpz_sub(det3,temp1,temp2);

/* Compute det4 */

   mpz_mul(temp1,res[0][3],res[3][4]);
   mpz_mul(temp2,res[0][2],res[2][4]);
   mpz_mul(temp3,res[0][1],res[1][4]);
   mpz_add(temp1,temp1,temp2); mpz_add(temp1,temp3,temp1);
   coef = 2;
   mpz_mul_si(temp1,temp1,coef);
   mpz_mul(temp1,temp1,res[0][1]);
   mpz_mul(temp2,res[1][3],res[1][3]);
   mpz_mul(temp3,res[1][2],res[1][2]);
   mpz_add(temp2,temp3,temp2);
   mpz_mul(temp2,temp2,res[0][1]);
   mpz_mul(temp3,res[2][3],diff);
   mpz_sub(temp2,temp3,temp2);
   coef = 4;
   mpz_mul_si(temp2,temp2,coef);
   mpz_add(det4,temp1,temp2);

/* Compute numerator of the radius of the smallest circumsphere of the edge */

   mpz_mul(temp1,det0,det4);
   mpz_mul(temp2,det3,det3);
   mpz_sub(temp2,temp2,temp1);
   mpz_mul(temp1,det2,det2);
   mpz_add(temp2,temp2,temp1);
   mpz_mul(temp1,det1,det1);
   mpz_add(num,temp1,temp2);

/* Compute denominator of the radius of the smallest circumsphere of the edge */

   mpz_mul(den,det0,det0);

/* check if radius is lower than ALPHA         */

   mpz_mul(temp1,den,alp);
   mpz_sub(temp2,num,temp1);

   if(mpz_sgn(temp2) < 0) (*testr)=1;

}
#ifdef __cplusplus
}
#endif

/*------------------------------------------------------------------------*/
/* triangle_radius_gmp:
   This program checks if the radius of the circumsphere of a facet of the
   regular triangulation is smaller than alpha


       Input:

       For the three points a,b,c that form the triangles, the program
       needs as input the following determinants:

       S(i,j)   = Minor(a,b,c,i,j,0)= det | a(i)  a(j)  1 |
                                          | b(i)  b(j)  1 |
                                          | c(i)  c(j)  1 |

       for i in [1,3] and j in [i+1,4]

       and:

       T(i,j) = Minor(a,b,c,i,j,4)=det | a(i) a(j) a(4) |
                                       | b(i) b(j) b(4) |
                                       | c(i) c(j) c(4) |

       and

       Dabc  = Minor(a,b,c,1,2,3)

   Output:

   testr   : flag set to 1 if ALPHA is larger than rho, the radius
        of the circumsphere of the triangle

*/

void F77Name2(triangle_radius_gmp)(double *a, double *b, double *c, double *ra, 
      double *rb, double *rc,
      int *testr, double *scale, double *alpha, int *memory)

{
   int i,j,coef, ivalue;
   double value;

   value = (*alpha)*(*scale); ivalue = (int) floor(value); 
   mpz_set_si(alp,ivalue);

   if(*memory !=1) set_triangle(a, b, c, *ra, *rb, *rc, *scale);

   (*testr) = 0;

   mpz_set_si(temp1,0);
   for (i=1; i<3; i++)
   {
      for (j=i+1; j<4; j++)
      {
         mpz_mul(temp2,S[i][j],S[i][j]);
         mpz_add(temp1,temp1,temp2);
      }
   }

/* Compute det0 */

   coef = 4;
   mpz_mul_si(det0,temp1,coef);

/* Compute det1 */

   mpz_mul(temp1,Dabc,S[2][3]);
   coef = -2;
   mpz_mul_si(temp1,temp1,coef);
   mpz_mul(temp2,S[1][2],S[2][4]);
   mpz_add(temp1,temp2,temp1);
   mpz_mul(temp2,S[1][3],S[3][4]);
   mpz_add(temp1,temp2,temp1);
   mpz_mul_si(det1,temp1,coef);

/* Compute det2 */

   coef = 2;
   mpz_mul(temp1,Dabc,S[1][3]);
   mpz_mul_si(temp1,temp1,coef);
   mpz_mul(temp2,S[2][3],S[3][4]);
   mpz_add(temp1,temp2,temp1);
   mpz_mul(temp2,S[1][2],S[1][4]);
   mpz_sub(temp1,temp2,temp1);
   mpz_mul_si(det2,temp1,coef);

/* Compute det3 */

   mpz_mul(temp1,Dabc,S[1][2]);
   mpz_mul_si(temp1,temp1,coef);
   mpz_mul(temp2,S[1][3],S[1][4]);
   mpz_add(temp1,temp2,temp1);
   mpz_mul(temp2,S[2][3],S[2][4]);
   mpz_add(temp1,temp2,temp1);
   mpz_mul_si(det3,temp1,coef);

/* Compute det4 */

   mpz_mul(temp1,Dabc,Dabc);
   coef = -2;
   mpz_mul_si(temp1,temp1,coef);

   for (i=1; i<3; i++)
   {
      for (j=i+1; j<4; j++)
      {
         mpz_mul(temp2,S[i][j],T[i][j]);
         mpz_add(temp1,temp1,temp2);
      }
   }
   coef = -4;
   mpz_mul_si(det4,temp1,coef);

/* Now compute numerator of the radius of the circumsphere of the triangle */

   mpz_mul(temp1,det0,det4);
   mpz_mul(temp2,det3,det3);
   mpz_sub(temp2,temp2,temp1);
   mpz_mul(temp1,det2,det2);
   mpz_add(temp2,temp2,temp1);
   mpz_mul(temp1,det1,det1);
   mpz_add(num,temp1,temp2);

/* Now compute denominator of the radius of the circumsphere of the triangle */

   mpz_mul(den,det0,det0);

/* Check if radius is lower than ALPHA */

   mpz_mul(temp1,den,alp);
   mpz_sub(temp2,num,temp1);

   if(mpz_sgn(temp2) < 0) (*testr) = 1;

}
#ifdef __cplusplus
}
#endif

/*------------------------------------------------------------------------*/
/* triangle_attach_gmp:
   This program checks if a facet is attached to a vertex

   Input:

   For the three points a,b,c that form the triangles, the program
   needs as input the following determinants:

         S(i,j) = Minor(a,b,c,i,j,0)= det | a(i)  a(j)  1 |
                                          | b(i)  b(j)  1 |
                                          | c(i)  c(j)  1 |
       for all i in [1,3], j in [i+1,4]

       T(i,j) = M(a,b,c,i,j,4)    = Det | a(i) a(j) a(4) |
                                        | b(i) b(j) b(4) |
                                        | c(i) c(j) c(4) |

       for all i in [1,2] and all j in [i+1,3]

       Dabc = Det | a(1) a(2) a(3) |
                  | b(1) b(2) b(3) |
                  | c(1) c(2) c(3) |

       and the coordinates of the fourth vertex d


   testa   : flag set to 1 if the fourth point d is inside the
        circumsphere of (a,b,c)

*/

#ifdef __cplusplus
extern "C" {
#endif
void F77Name2(triangle_attach_gmp)(double *a, double *b, double *c, double *d,
        double *ra, double *rb, double *rc, double *rd, double *scale, 
   int *testa, int *memory)

{
   int i,coef;

   if(*memory !=1) set_triangle(a, b, c, *ra, *rb, *rc, *scale);

   for (i=0; i<3; i++)
   {
      real_to_gmp(d,*scale,i,d_mp[i+1]);
   }
   scalar_to_gmp(*rd,*scale,rd_mp);
   build_weight(d_mp[1],d_mp[2],d_mp[3],rd_mp,d_mp[4]);

/*
       We need to compute:

       det1 = Minor(a,b,c,d,2,3,4,0)
       det2 = Minor(a,b,c,d,1,3,4,0)
       det3 = Minor(a,b,c,d,1,2,4,0)
       det4 = Minor(a,b,c,d,1,2,3,0)
*/
   mpz_mul(temp1,d_mp[2],S[3][4]); mpz_mul(temp2,d_mp[3],S[2][4]); mpz_sub(temp1,temp2,temp1);
   mpz_mul(temp2,d_mp[4],S[2][3]); mpz_sub(temp2,T[2][3],temp2); mpz_add(det1,temp2,temp1);
   mpz_mul(temp1,d_mp[1],S[3][4]); mpz_mul(temp2,d_mp[3],S[1][4]); mpz_sub(temp1,temp2,temp1);
   mpz_mul(temp2,d_mp[4],S[1][3]); mpz_sub(temp2,T[1][3],temp2); mpz_add(det2,temp2,temp1);
   mpz_mul(temp1,d_mp[1],S[2][4]); mpz_mul(temp2,d_mp[2],S[1][4]); mpz_sub(temp1,temp2,temp1);
   mpz_mul(temp2,d_mp[4],S[1][2]); mpz_sub(temp2,T[1][2],temp2); mpz_add(det3,temp2,temp1);
   mpz_mul(temp1,d_mp[1],S[2][3]); mpz_mul(temp2,d_mp[2],S[1][3]); mpz_sub(temp1,temp2,temp1);
   mpz_mul(temp2,d_mp[3],S[1][2]); mpz_sub(temp2,Dabc,temp2); mpz_add(det4,temp2,temp1);

/* check if triangle is attached                               */

   (*testa) = 0;

   mpz_set_si(temp1,1);
/*   for (i=1; i<4; i++)
   {
      mpz_mul(temp2,S[i],S[i]);
      mpz_add(temp1,temp1,temp2);
   }
*/
   mpz_mul(temp2,det4,Dabc);
   coef = -2;
   mpz_mul_si(temp2,temp2,coef);
   mpz_mul(temp3,det3,S[1][2]);
   mpz_add(temp2,temp3,temp2);
   mpz_mul(temp3,det2,S[1][3]);
   mpz_add(temp2,temp3,temp2);
   mpz_mul(temp3,det1,S[2][3]);
   mpz_add(temp2,temp3,temp2);
   mpz_mul(dtest,temp1,temp2);

   if(mpz_sgn(dtest) > 0) (*testa)=1;

/*    Clear local GMP variables */

}
#ifdef __cplusplus
}
#endif

/*------------------------------------------------------------------------*/
/* set_alf_gmp   Version 1 11/24/2000   Patrice Koehl
 
   This procedure initialises all gmp variables that can be used for 
   computing the dual complex
*/

#ifdef __cplusplus
extern "C" {
#endif
void F77Name2(set_alf_gmp)()
{
/*   Initialise local GMP variables */

   int i,j;

   mpz_init(ra2); mpz_init(rb2); mpz_init(dist2);
   mpz_init(dtest);
   mpz_init (r_11); mpz_init (r_22); mpz_init (r_33);
   mpz_init (diff);
   mpz_init (det0);

   for (i= 0; i < 4; i++)
   {
      for (j=0; j < 5; j++)
      {
         mpz_init(res[i][j]);
         mpz_init(res2_c[i][j]);
         mpz_init(Mab[i][j]);
         mpz_init(Mac[i][j]);
         mpz_init(Mbc[i][j]);
         mpz_init(S[i][j]);
         mpz_init(T[i][j]);
      }
   }
   mpz_init (r_14); mpz_init (r_313); mpz_init (r_212);
   mpz_init (det1); mpz_init (det2); mpz_init (det3); mpz_init (det4); 
   mpz_init (wa); mpz_init(wb); mpz_init(wc); mpz_init(wd);

   for (i = 0; i < 5; i++) 
   {
      mpz_init(a_mp[i]);
      mpz_init(b_mp[i]);
      mpz_init(c_mp[i]);
      mpz_init(d_mp[i]);
      mpz_init(Dab[i]);
   }
   mpz_init(ra_mp);mpz_init(rb_mp);
   mpz_init(rc_mp);mpz_init(rd_mp);

   mpz_init (num); mpz_init (den);

   for (i=0; i < 4; i++)
   {
      mpz_init(Sab[i]);
      mpz_init(Sac[i]);
      mpz_init(Sad[i]);
      mpz_init(Sbc[i]);
      mpz_init(Sbd[i]);
      mpz_init(Scd[i]);
      mpz_init(Tab[i]);
      mpz_init(Sa[i]);
      mpz_init(Sb[i]);
      mpz_init(Sc[i]);
      mpz_init(Sd[i]);
      mpz_init(Sam1[i]);
      mpz_init(Sbm1[i]);
      mpz_init(Scm1[i]);
      mpz_init(Sdm1[i]);
      mpz_init(Tc[i]);
      mpz_init(Deter[i]);
   }

   mpz_init(Dabc); mpz_init(Dabd); mpz_init(Dacd); mpz_init(Dbcd);
   mpz_init(alp);mpz_init(Dabcd);

}
#ifdef __cplusplus
}
#endif

/*------------------------------------------------------------------------*/
/* clear_alf_gmp   Version 1 11/24/2000   Patrice Koehl
 
   This procedure clears all gmp variables that were used for 
   computing the dual complex
*/

#ifdef __cplusplus
extern "C" {
#endif
void F77Name2(clear_alf_gmp)()
{

   int i,j;

   mpz_clear(ra2); mpz_clear(rb2); mpz_clear(dist2);
   mpz_clear(dtest);
   mpz_clear (r_11); mpz_clear (r_22); mpz_clear (r_33);
   mpz_clear (diff);
   mpz_clear (det0);

   for (i= 0; i < 4; i++)
   {
      for (j=0; j < 5; j++)
      {
         mpz_clear(res[i][j]);
         mpz_clear(res2_c[i][j]);
         mpz_clear(Mab[i][j]);
         mpz_clear(Mac[i][j]);
         mpz_clear(Mbc[i][j]);
         mpz_clear(S[i][j]);
      }
   }
   for (i= 0; i < 3; i++)
   {
      for (j=0; j < 4; j++)
      {
         mpz_clear(T[i][j]);
      }
   }
   mpz_clear (r_14); mpz_clear (r_313); mpz_clear (r_212);
   mpz_clear (det1); mpz_clear (det2); mpz_clear (det3); mpz_clear (det4); 
   mpz_clear (wa); mpz_clear(wb); mpz_clear(wc); mpz_clear(wd);

   for (i = 0; i < 5; i++) 
   {
      mpz_clear(a_mp[i]);
      mpz_clear(b_mp[i]);
      mpz_clear(c_mp[i]);
      mpz_clear(d_mp[i]);
      mpz_clear(Dab[i]);
   }
   mpz_clear(ra_mp);mpz_clear(rb_mp);
   mpz_clear(rc_mp);mpz_clear(rd_mp);

   mpz_clear (num); mpz_clear (den);

   for (i=0; i < 4; i++)
   {
      mpz_clear(Sab[i]);
      mpz_clear(Sac[i]);
      mpz_clear(Sad[i]);
      mpz_clear(Sbc[i]);
      mpz_clear(Sbd[i]);
      mpz_clear(Scd[i]);
      mpz_clear(Tab[i]);
      mpz_clear(Sa[i]);
      mpz_clear(Sb[i]);
      mpz_clear(Sc[i]);
      mpz_clear(Sd[i]);
      mpz_clear(Sam1[i]);
      mpz_clear(Sbm1[i]);
      mpz_clear(Scm1[i]);
      mpz_clear(Sdm1[i]);
      mpz_clear(Tc[i]);
      mpz_clear(Deter[i]);
   }

   mpz_clear(Dabc);mpz_clear(Dabd); mpz_clear(Dacd); mpz_clear(Dbcd);
   mpz_clear(Dabcd);
   mpz_clear(alp);

}
#ifdef __cplusplus
}
#endif

/*   Sos_minor_gmp.c   Version 1 1/16/2002   Patrice Koehl             */
/*                             */
/*  This is the C version of sos_minor.f, which performs all operations   */
/*  with multi precision arithmetics, using the package GMP               */
/*                             */
/*------------------------------------------------------------------------*/
/*                                                                        */
/* Copyright (C) 2002 Patrice Koehl                                       */
/*                                                                        */
/* This library is free software; you can redistribute it and/or          */
/* modify it under the terms of the GNU Lesser General Public             */
/* License as published by the Free Software Foundation; either           */
/* version 2.1 of the License, or (at your option) any later version.     */
/*                                                                        */
/* This library is distributed in the hope that it will be useful,        */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of         */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU       */
/* Lesser General Public License for more details.                        */
/*                                                                        */
/* You should have received a copy of the GNU Lesser General Public       */
/* License along with this library; if not, write to the Free Software    */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA*/
/*                                                                        */
/*------------------------------------------------------------------------*/
/*                                                                        */
/* Includes :                          */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "uniongmp.h"

#define F77Name1(x) x##_   /* Need to add one underscore to Fortran program */

#if defined intel
#define F77Name2(x) x##_   /* Need to add one underscore to Fortran program for
                                Intel compiler */
#else
#define F77Name2(x) x##__   /* Need to add two underscores to Fortran program for GNU
                                f77 compiler */
#endif

#define round(x) ((x)>=0?(int)((x)+0.5):(int)((x)-0.5))

/*------------------------------------------------------------------------*/

void deter2_gmp(mpz_t deter, mpz_t a, mpz_t b);

void deter3_gmp(mpz_t deter, mpz_t a11, mpz_t a12, mpz_t a21, mpz_t a22,
      mpz_t a31, mpz_t a32);

void deter4_gmp(mpz_t deter, mpz_t a11, mpz_t a12, mpz_t a13, mpz_t a21,
      mpz_t a22, mpz_t a23, mpz_t a31, mpz_t a32,
      mpz_t a33, mpz_t a41, mpz_t a42, mpz_t a43);

void deter5_gmp(mpz_t deter, mpz_t a11, mpz_t a12, mpz_t a13, mpz_t a14,
      mpz_t a21, mpz_t a22, mpz_t a23, mpz_t a24,
      mpz_t a31, mpz_t a32, mpz_t a33, mpz_t a34,
      mpz_t a41, mpz_t a42, mpz_t a43, mpz_t a44,
      mpz_t a51, mpz_t a52, mpz_t a53, mpz_t a54);

void real_to_gmp(double *coord, double scale, int idx, mpz_t val);

void build_weight(mpz_t ax, mpz_t ay, mpz_t az, mpz_t r, mpz_t w);

#ifdef __cplusplus
extern "C"
#endif
void F77Name2(init_sos_gmp)();

#ifdef __cplusplus
extern "C"
#endif
void F77Name2(clear_sos_gmp)();

#ifdef __cplusplus
extern "C"
#endif
void F77Name2(sos_minor2_gmp)(double *coord, double *scale, int *a, int *b, int *i1, int *res);
#ifdef __cplusplus
extern "C"
#endif
void F77Name2(sos_minor3_gmp)(double *coord, double *scale, int *a, int *b, int *c, int *i1, int *i2, int *res);

#ifdef __cplusplus
extern "C"
#endif
void F77Name2(sos_minor4_gmp)(double *coord, double *scale, int *a, int *b, int *c, int *d, int *res);
#ifdef __cplusplus
extern "C"
#endif
void F77Name2(sos_minor5_gmp)(double *coord, double *radius, double *scale, int *a, int *b, int *c, int *d, int *e, int *res);

#ifdef __cplusplus
extern "C"
#endif
void F77Name2(minor4_gmp)(double *coord, double *scale, int *a, int *b, int *c, int *d, int *res);

/*                             */
/*------------------------------------------------------------------------*/
/* init_sos_gmp:                       */
/* This subroutine initializes the GMP variables and arrays
                                                                          */
#ifdef __cplusplus
extern "C" {
#endif
void F77Name2(init_sos_gmp)()
{

   mpz_init(a11_mp);mpz_init(a12_mp); mpz_init(a13_mp); mpz_init(a14_mp);
   mpz_init(a21_mp);mpz_init(a22_mp); mpz_init(a23_mp); mpz_init(a24_mp);
   mpz_init(a31_mp);mpz_init(a32_mp); mpz_init(a33_mp); mpz_init(a34_mp);
   mpz_init(a41_mp);mpz_init(a42_mp); mpz_init(a43_mp); mpz_init(a44_mp);
   mpz_init(a51_mp);mpz_init(a52_mp); mpz_init(a53_mp); mpz_init(a54_mp);

   mpz_init(r1_mp);mpz_init(r2_mp); mpz_init(r3_mp); mpz_init(r4_mp); mpz_init(r5_mp);

   mpz_init(temp1); mpz_init(temp2); mpz_init(temp3); mpz_init(temp4);
   mpz_init(val1);mpz_init(val2); mpz_init(val3);

   mpz_init(c11); mpz_init(c12); mpz_init(c13); mpz_init(c14);
   mpz_init(c21); mpz_init(c22); mpz_init(c23); mpz_init(c24);
   mpz_init(c31); mpz_init(c32); mpz_init(c33); mpz_init(c34);
   mpz_init(c41); mpz_init(c42); mpz_init(c43); mpz_init(c44);
   mpz_init(d1); mpz_init(d2); mpz_init(d3);
   mpz_init(e1); mpz_init(e2); mpz_init(e3);
   mpz_init(f1); mpz_init(f2); mpz_init(f3);
   mpz_init(g1); mpz_init(g2); mpz_init(g3);

}

/*                             */
/*------------------------------------------------------------------------*/
/* real_to_gmp:                       */
/* This subroutine converts one coordinate of one point into mpz_t type
*/

void real_to_gmp(double *coord, double scale, int idx, mpz_t val)
{
   double x;
   int ivalue;

   mpz_set_d(temp3, scale);

   ivalue= (int) coord[idx];
   mpz_set_si(temp1,ivalue);
   mpz_mul(temp1,temp1,temp3);
   x = (coord[idx]-ivalue)*(scale);
   ivalue = (int) round(x);
   mpz_set_si(temp2,ivalue);
   mpz_add(val,temp1,temp2);
}

/*                             */
/*------------------------------------------------------------------------*/
/* clear_sos_gmp:                       */
/* This subroutine clears all pending mpz_t arrays
                                                                          */

#ifdef __cplusplus
extern "C" {
#endif
void F77Name2(clear_sos_gmp)()
{

   mpz_clear(a11_mp);mpz_clear(a12_mp); mpz_clear(a13_mp); mpz_clear(a14_mp);
   mpz_clear(a21_mp);mpz_clear(a22_mp); mpz_clear(a23_mp); mpz_clear(a24_mp);
   mpz_clear(a31_mp);mpz_clear(a32_mp); mpz_clear(a33_mp); mpz_clear(a34_mp);
   mpz_clear(a41_mp);mpz_clear(a42_mp); mpz_clear(a43_mp); mpz_clear(a44_mp);
   mpz_clear(a51_mp);mpz_clear(a52_mp); mpz_clear(a53_mp); mpz_clear(a54_mp);
   mpz_clear(r1_mp);mpz_clear(r2_mp); mpz_clear(r3_mp); mpz_clear(r4_mp); mpz_clear(r5_mp);
   mpz_clear(temp1); mpz_clear(temp2); mpz_clear(temp3); mpz_clear(temp4);
   mpz_clear(val1);mpz_clear(val2);mpz_clear(val3);

   mpz_clear(c11); mpz_clear(c12); mpz_clear(c13); mpz_clear(c14);
   mpz_clear(c21); mpz_clear(c22); mpz_clear(c23); mpz_clear(c24);
   mpz_clear(c31); mpz_clear(c32); mpz_clear(c33); mpz_clear(c34);
   mpz_clear(c41); mpz_clear(c42); mpz_clear(c43); mpz_clear(c44);
   mpz_clear(d1); mpz_clear(d2); mpz_clear(d3);
   mpz_clear(e1); mpz_clear(e2); mpz_clear(e3);
   mpz_clear(f1); mpz_clear(f2); mpz_clear(f3);
   mpz_clear(g1); mpz_clear(g2); mpz_clear(g3);
}

/*                             */
/*------------------------------------------------------------------------*/
/* deter2_gmp:                          */
/*                             */
/*   This subroutine evaluates the determinant:
      D = | b11 1 |
          | b21 1 |

   Input:
      b11, b21
   Output:
      deter
                             */

void deter2_gmp(mpz_t deter, mpz_t b11, mpz_t b21)

{
   mpz_sub(deter,b11,b21);
}

/*                             */
/*------------------------------------------------------------------------*/
/* deter3_gmp:                          */
/*                             */
/*   This subroutine evaluates the determinant:
      D = | b11 b12 1 |
          | b21 b22 1 |
          | b31 b32 1 |

   Input:
      b11, b12, b21, b22, b31, b32
   Output:
      deter3_gmp
                             */

void deter3_gmp(mpz_t deter, mpz_t b11, mpz_t b12, mpz_t b21, 
      mpz_t b22, mpz_t b31, mpz_t b32)

{

   mpz_sub(temp1,b21,b11);
   mpz_sub(temp2,b22,b12);
   mpz_sub(temp3,b31,b11);
   mpz_sub(temp4,b32,b12);

   mpz_mul(val1,temp1,temp4);
   mpz_mul(val2,temp2,temp3);

   mpz_sub(deter,val1,val2);

}

/*                             */
/*------------------------------------------------------------------------*/
/* deter4_gmp:                          */
/*                             */
/*   This subroutine evaluates the determinant:
      D = | b11 b12 b13 1 |
          | b21 b22 b23 1 |
          | b31 b32 b33 1 |
          | b41 b42 b43 1 |

   Input:
      b11, b12, b13, b21, b22, b23, b31, b32, b33
      b41, b42, b43
   Output:
      deter4_gmp
                             */

void deter4_gmp(mpz_t deter, mpz_t b11, mpz_t b12, mpz_t b13, mpz_t b21, 
      mpz_t b22, mpz_t b23, mpz_t b31, mpz_t b32, mpz_t b33, 
      mpz_t b41, mpz_t b42, mpz_t b43)

{
   mpz_sub(c11,b21,b11);mpz_sub(c12,b22,b12);mpz_sub(c13,b23,b13);
   mpz_sub(c21,b31,b11);mpz_sub(c22,b32,b12);mpz_sub(c23,b33,b13);
   mpz_sub(c31,b41,b11);mpz_sub(c32,b42,b12);mpz_sub(c33,b43,b13);

   mpz_mul(temp1,c22,c33);mpz_mul(temp2,c32,c23);mpz_sub(val1,temp1,temp2);
   mpz_mul(temp1,c12,c33);mpz_mul(temp2,c32,c13);mpz_sub(val2,temp1,temp2);
   mpz_mul(temp1,c12,c23);mpz_mul(temp2,c22,c13);mpz_sub(val3,temp1,temp2);

   mpz_mul(temp1,c21,val2);mpz_mul(temp2,c11,val1);mpz_mul(temp3,c31,val3);

   mpz_add(val1,temp2,temp3);
   mpz_sub(deter,temp1,val1);

}

/*                             */
/*------------------------------------------------------------------------*/
/* deter5_gmp:                          */
/*                             */
/*   This subroutine evaluates the determinant:
      D = | b11 b12 b13 b14 1 |
          | b21 b22 b23 b24 1 |
          | b31 b32 b33 b34 1 |
          | b41 b42 b43 b44 1 |
          | b51 b52 b53 b54 1 |

   Input:
      b11, b12, b13, b14, b21, b22, b23, b24, b31, b32, b33, b34
      b41, b42, b43, b44, b51, b52, b53, b54
   Output:
      deter5_gmp
                             */

void deter5_gmp(mpz_t deter, mpz_t b11, mpz_t b12, mpz_t b13, mpz_t b14, 
       mpz_t b21, mpz_t b22, mpz_t b23, mpz_t b24,
       mpz_t b31, mpz_t b32, mpz_t b33, mpz_t b34,
       mpz_t b41, mpz_t b42, mpz_t b43, mpz_t b44,
       mpz_t b51, mpz_t b52, mpz_t b53, mpz_t b54)
{

   mpz_sub(c11,b21,b11); mpz_sub(c12,b22,b12); mpz_sub(c13,b23,b13);
   mpz_sub(c14,b24,b14);
   mpz_sub(c21,b31,b11); mpz_sub(c22,b32,b12); mpz_sub(c23,b33,b13);
   mpz_sub(c24,b34,b14);
   mpz_sub(c31,b41,b11); mpz_sub(c32,b42,b12); mpz_sub(c33,b43,b13);
   mpz_sub(c34,b44,b14);
   mpz_sub(c41,b51,b11); mpz_sub(c42,b52,b12); mpz_sub(c43,b53,b13);
   mpz_sub(c44,b54,b14);

   mpz_mul(temp1,c32,c43); mpz_mul(temp2,c42,c33); mpz_sub(d1,temp1,temp2);
   mpz_mul(temp1,c32,c44); mpz_mul(temp2,c42,c34); mpz_sub(d2,temp1,temp2);
   mpz_mul(temp1,c33,c44); mpz_mul(temp2,c43,c34); mpz_sub(d3,temp1,temp2);

   mpz_mul(temp1,c12,c23); mpz_mul(temp2,c22,c13); mpz_sub(e1,temp1,temp2);
   mpz_mul(temp1,c12,c24); mpz_mul(temp2,c22,c14); mpz_sub(e2,temp1,temp2);
   mpz_mul(temp1,c13,c24); mpz_mul(temp2,c23,c14); mpz_sub(e3,temp1,temp2);

   mpz_mul(temp1,c11,c24); mpz_mul(temp2,c21,c14); mpz_sub(f1,temp1,temp2);
   mpz_mul(temp1,c11,c23); mpz_mul(temp2,c21,c13); mpz_sub(f2,temp1,temp2);
   mpz_mul(temp1,c11,c22); mpz_mul(temp2,c21,c12); mpz_sub(f3,temp1,temp2);

   mpz_mul(temp1,c31,c44); mpz_mul(temp2,c41,c34); mpz_sub(g1,temp1,temp2);
   mpz_mul(temp1,c31,c43); mpz_mul(temp2,c41,c33); mpz_sub(g2,temp1,temp2);
   mpz_mul(temp1,c31,c42); mpz_mul(temp2,c41,c32); mpz_sub(g3,temp1,temp2);
 
   mpz_mul(temp1,e3,g3); mpz_mul(temp2,e2,g2); mpz_sub(temp3,temp1,temp2);
   mpz_mul(temp1,e1,g1); mpz_add(temp3,temp3,temp1);
   mpz_mul(temp1,d3,f3); mpz_add(temp3,temp3,temp1);
   mpz_mul(temp1,d2,f2); mpz_sub(temp3,temp3,temp1);
   mpz_mul(temp1,d1,f1); mpz_add(deter,temp3,temp1);

}

/*                             */
/*------------------------------------------------------------------------*/
/* sos_minor2_gmp:                       */
/*                             
   This subroutine tests the sign of the determinant
      D = | a11 1 |
          | a21 1 |
   If the determinant is found to be 0, then the SoS procedure is used:
   a development of the determinant with respect to a perturbation EPS
   applied to the coordinates in the determinant is computed, and
   the sign of the first non zero term defines the sign of the 
   determinant.
   In the case of a 2x2 determinant, the first term in the expansion
   is the coefficient 1 ...               */
#ifdef __cplusplus
extern "C" {
#endif
void F77Name2(sos_minor2_gmp)(double *coord, double *scale, int *a, int *b, int *ia, int *res)
{
   int icomp,idx_a,idx_b;

/* Initialise local GMP variables */

/* Get coordinates */

   idx_a = (*a)*3 + (*ia) -4;
   idx_b = (*b)*3 + (*ia) -4;
   real_to_gmp(coord,*scale,idx_a,a11_mp);
   real_to_gmp(coord,*scale,idx_b,a21_mp);

/* Compute determinant */

   deter2_gmp(temp1,a11_mp,a21_mp);

   icomp = mpz_sgn(temp1);

   if (icomp != 0) {
      *res = icomp;
   }
   else {
      *res = 1;
   }

/* Clear local GMP variables */

}
#ifdef __cplusplus
}
#endif

/*                             */
/*------------------------------------------------------------------------*/
/* sos_minor3_gmp:                       */
/*                            
   This subroutine tests the sign of the determinant
      D = | a11 a12 1 |
          | a21 a22 1 |
          | a31 a32 1 |
   If the determinant is found to be 0, then the SoS procedure is used:
   a development of the determinant with respect to a perturbation EPS
   applied to the coordinates in the determinant is computed, and
   the sign of the first non zero term defines the sign of the 
   determinant.
   In the case of a 3x3 determinant, the maximum number of terms to be
   checked is 4 ...               */
#ifdef __cplusplus
extern "C" {
#endif
void F77Name2(sos_minor3_gmp)(double *coord, double *scale, int *a, int *b, int *c, int *i1, int *i2, int *res)
{
   int icomp;
   int idx_a1,idx_a2,idx_b1,idx_b2,idx_c1,idx_c2;

/* Initialise local GMP variables */


/* Transfer coordinates to GMP */

   idx_a1 = (*a)*3 + (*i1) -4; idx_a2 = (*a)*3 + (*i2) -4;
   idx_b1 = (*b)*3 + (*i1) -4; idx_b2 = (*b)*3 + (*i2) -4;
   idx_c1 = (*c)*3 + (*i1) -4; idx_c2 = (*c)*3 + (*i2) -4;

   real_to_gmp(coord,*scale,idx_a1,a11_mp);
   real_to_gmp(coord,*scale,idx_a2,a12_mp);
   real_to_gmp(coord,*scale,idx_b1,a21_mp);
   real_to_gmp(coord,*scale,idx_b2,a22_mp);
   real_to_gmp(coord,*scale,idx_c1,a31_mp);
   real_to_gmp(coord,*scale,idx_c2,a32_mp);

/* Compute determinant */

   deter3_gmp(temp1,a11_mp,a12_mp,a21_mp,a22_mp,a31_mp,a32_mp);

   icomp = mpz_sgn(temp1);

/* if major determinant is non 0, return its sign */

   if (icomp != 0) {
      *res = icomp;
      return;
   }

/* Look now at each term in the expansion of the determinant with
   respect to EPS                                               
   The initial determinant is:
      Minor3(i,j,k,1,2,0)            */

/* Term 1: - Minor2(j,k,1,0) */

   deter2_gmp(temp1,a21_mp,a31_mp);
   icomp = mpz_sgn(temp1);

   if (icomp != 0) {
      *res = - icomp;
      return;
   }

/* Term 2: Minor2(j,k,2,0)  */

   deter2_gmp(temp1,a22_mp,a32_mp);
   icomp = mpz_sgn(temp1);

   if (icomp != 0) {
      *res = icomp;
      return;
   }

/* Term 3: Minor2(i,k,1,0)  */

   deter2_gmp(temp1,a11_mp,a31_mp);
   icomp = mpz_sgn(temp1);

   if (icomp != 0) {
      *res = icomp;
      return;
   }

/* Term 4: 1 */

   *res = 1;

}
#ifdef __cplusplus
}
#endif

/*                             */
/*------------------------------------------------------------------------*/
/* sos_minor4_:                          */
/*                             
   This subroutine tests the sign of the determinant
      D = | a11 a12 a13 1 |
          | a21 a22 a23 1 |
          | a31 a32 a33 1 |
          | a41 a42 a43 1 |
   If the determinant is found to be 0, then the SoS procedure is used:
   a development of the determinant with *respect to a perturbation EPS
   applied to the coordinates in the determinant is computed, and
   the sign of the first non zero term defines the sign of the 
   determinant.
   In the case of a 4x4 determinant, the maximum number of terms to be
   checked is 14 ...               */
#ifdef __cplusplus
extern "C" {
#endif
void F77Name2(sos_minor4_gmp)(double *coord, double *scale, int *a, int *b, int *c, int *d, int *res)

{
   int icomp;
   int idx_a1,idx_a2,idx_a3,idx_b1,idx_b2,idx_b3;
   int idx_c1,idx_c2,idx_c3,idx_d1,idx_d2,idx_d3;

/* Transfer coordinates to gmp */

   idx_a1 = (*a)*3 -3; idx_a2 = idx_a1+1; idx_a3 = idx_a2 +1;
   idx_b1 = (*b)*3 -3; idx_b2 = idx_b1+1; idx_b3 = idx_b2 +1;
   idx_c1 = (*c)*3 -3; idx_c2 = idx_c1+1; idx_c3 = idx_c2 +1;
   idx_d1 = (*d)*3 -3; idx_d2 = idx_d1+1; idx_d3 = idx_d2 +1;

   real_to_gmp(coord,*scale,idx_a1,a11_mp);
   real_to_gmp(coord,*scale,idx_a2,a12_mp);
   real_to_gmp(coord,*scale,idx_a3,a13_mp);
   real_to_gmp(coord,*scale,idx_b1,a21_mp);
   real_to_gmp(coord,*scale,idx_b2,a22_mp);
   real_to_gmp(coord,*scale,idx_b3,a23_mp);
   real_to_gmp(coord,*scale,idx_c1,a31_mp);
   real_to_gmp(coord,*scale,idx_c2,a32_mp);
   real_to_gmp(coord,*scale,idx_c3,a33_mp);
   real_to_gmp(coord,*scale,idx_d1,a41_mp);
   real_to_gmp(coord,*scale,idx_d2,a42_mp);
   real_to_gmp(coord,*scale,idx_d3,a43_mp);

/*   printf("a11_mp = %s\n",mpz_get_str(NULL,base,a11_mp)); 
   printf("a12_mp = %s\n",mpz_get_str(NULL,base,a12_mp)); 
   printf("a13_mp = %s\n",mpz_get_str(NULL,base,a13_mp)); 
   printf("a21_mp = %s\n",mpz_get_str(NULL,base,a21_mp)); 
   printf("a22_mp = %s\n",mpz_get_str(NULL,base,a22_mp)); 
   printf("a23_mp = %s\n",mpz_get_str(NULL,base,a23_mp)); 
   printf("a31_mp = %s\n",mpz_get_str(NULL,base,a31_mp)); 
   printf("a32_mp = %s\n",mpz_get_str(NULL,base,a32_mp)); 
   printf("a33_mp = %s\n",mpz_get_str(NULL,base,a33_mp)); 
   printf("a41_mp = %s\n",mpz_get_str(NULL,base,a41_mp)); 
   printf("a42_mp = %s\n",mpz_get_str(NULL,base,a42_mp)); 
   printf("a43_mp = %s\n",mpz_get_str(NULL,base,a43_mp));  */

/* Compute determinant */

   deter4_gmp(temp1,a11_mp,a12_mp,a13_mp,a21_mp,a22_mp,a23_mp,
            a31_mp,a32_mp,a33_mp,a41_mp,a42_mp,a43_mp);

//   printf("Deter4 = %s\n",mpz_get_str(NULL,base,temp1)); 

   icomp = mpz_sgn(temp1);

/* if major determinant is non 0, return its sign 
   (don't forget to clear GMP variables !!)   */

   if (icomp != 0) {
      *res = icomp;
      return;
   }

/* Look now at each term in the expansion of the determinant with
   *respect to EPS                                                
   The initial determinant is:
      Minor4(i,j,k,l,1,2,3,0)            */

/* Term 1:   Minor3(j,k,l,1,2,0)  */

   deter3_gmp(temp1,a21_mp,a22_mp,a31_mp,a32_mp,a41_mp,a42_mp);
   icomp = mpz_sgn(temp1);
   if (icomp != 0) {
      *res = icomp;
      return;
   }

/* Term 2:   -Minor3(j,k,l,1,3,0) */

   deter3_gmp(temp1,a21_mp,a23_mp,a31_mp,a33_mp,a41_mp,a43_mp);
   icomp = mpz_sgn(temp1);
   if (icomp != 0) {
      *res = - icomp;
      return;
   }

/* Term 3:   Minor3(j,k,l,2,3,0) */

   deter3_gmp(temp1,a22_mp,a23_mp,a32_mp,a33_mp,a42_mp,a43_mp);
   icomp = mpz_sgn(temp1);
   if (icomp != 0) {
      *res = icomp;
      return;
   }

/* Term 4:   - Minor3(i,k,l,1,2,0) */

   deter3_gmp(temp1,a11_mp,a12_mp,a31_mp,a32_mp,a41_mp,a42_mp);
   icomp = mpz_sgn(temp1);
   if (icomp != 0) {
      *res = - icomp;
      return;
   }

/* Term 5:   Minor2(k,l,1,0) */

   deter2_gmp(temp1,a31_mp,a41_mp);   
   icomp = mpz_sgn(temp1);
   if (icomp != 0) {
      *res = icomp;
      return;
   }

/* Term 6:   -Minor2(k,l,2,0) */

   deter2_gmp(temp1,a32_mp,a42_mp);   
   icomp = mpz_sgn(temp1);
   if (icomp != 0) {
      *res = - icomp;
      return;
   }

/* Term 7:   Minor3(i,k,l,1,3,0) */

   deter3_gmp(temp1,a11_mp,a13_mp,a31_mp,a33_mp,a41_mp,a43_mp);
   icomp = mpz_sgn(temp1);
   if (icomp != 0) {
      *res =  icomp;
      return;
   }

/* Term 8:   Minor2(k,l,3,0) */

   deter2_gmp(temp1,a33_mp,a43_mp);
   icomp = mpz_sgn(temp1);
   if (icomp != 0) {
      *res =  icomp;
      return;
   }

/* Term 9:   - Minor3(i,k,l,2,3,0) */

   deter3_gmp(temp1,a12_mp,a13_mp,a32_mp,a33_mp,a42_mp,a43_mp);
   icomp = mpz_sgn(temp1);
   if (icomp != 0) {
      *res =  - icomp;
      return;
   }

/* Term 10:   Minor3(i,j,l,1,2,0) */

   deter3_gmp(temp1,a11_mp,a12_mp,a21_mp,a22_mp,a41_mp,a42_mp);
   icomp = mpz_sgn(temp1);
   if (icomp != 0) {
      *res =  icomp;
      return;
   }

/* Term 11:    - Minor2(j,l,1,0) */

   deter2_gmp(temp1,a21_mp,a41_mp);
   icomp = mpz_sgn(temp1);
   if (icomp != 0) {
      *res =  - icomp;
      return;
   }

/* Term 12:   Minor2(j,l,2,0) */

   deter2_gmp(temp1,a22_mp,a42_mp);
   icomp = mpz_sgn(temp1);
   if (icomp != 0) {
      *res =  icomp;
      return;
   }

/* Term 13:   Minor2(i,l,1,0) */

   deter2_gmp(temp1,a11_mp,a41_mp);
   icomp = mpz_sgn(temp1);
   if (icomp != 0) {
      *res =  icomp;
      return;
   }

/* Term 14:   1 */

   *res = 1;

}
#ifdef __cplusplus
}
#endif

/*                             */
/*------------------------------------------------------------------------*/
/* sos_minor5_gmp:                       */
/*                             
   This subroutine tests the sign of the determinant
      D = | a11 a12 a13 a14 1 |
          | a21 a22 a23 a24 1 |
          | a31 a32 a33 a34 1 |
          | a41 a42 a43 a44 1 |
          | a51 a52 a53 a54 1 |
   If the determinant is found to be 0, then the SoS procedure is used:
   a development of the determinant with *respect to a perturbation EPS
   applied to the coordinates in the determinant is computed, and
   the sign of the first non zero term defines the sign of the 
   determinant.
   In the case of a 5x5 determinant, the maximum number of terms to be
   checked is 49 ...               */
#ifdef __cplusplus
extern "C" {
#endif
void F77Name2(sos_minor5_gmp)(double *coord, double *radius, double *scale, int *a, int *b, int *c, int *d, int *e, int *res) 

{
   int icomp;

   int idx_a1,idx_a2,idx_a3,idx_a4,idx_b1,idx_b2,idx_b3,idx_b4;
   int idx_c1,idx_c2,idx_c3,idx_c4,idx_d1,idx_d2,idx_d3,idx_d4;
   int idx_e1,idx_e2,idx_e3,idx_e4;

/*   Initialise local GMP variables */

   idx_a1 = (*a)*3 -3; idx_a2 = idx_a1+1; idx_a3 = idx_a2 +1;
   idx_b1 = (*b)*3 -3; idx_b2 = idx_b1+1; idx_b3 = idx_b2 +1;
   idx_c1 = (*c)*3 -3; idx_c2 = idx_c1+1; idx_c3 = idx_c2 +1;
   idx_d1 = (*d)*3 -3; idx_d2 = idx_d1+1; idx_d3 = idx_d2 +1;
   idx_e1 = (*e)*3 -3; idx_e2 = idx_e1+1; idx_e3 = idx_e2 +1;

   idx_a4 = (*a) -1; idx_b4 = (*b)-1; idx_c4 = (*c)-1; 
   idx_d4 = (*d) -1; idx_e4 = (*e)-1;

   real_to_gmp(coord,*scale,idx_a1,a11_mp);
   real_to_gmp(coord,*scale,idx_a2,a12_mp);
   real_to_gmp(coord,*scale,idx_a3,a13_mp);
   real_to_gmp(coord,*scale,idx_b1,a21_mp);
   real_to_gmp(coord,*scale,idx_b2,a22_mp);
   real_to_gmp(coord,*scale,idx_b3,a23_mp);
   real_to_gmp(coord,*scale,idx_c1,a31_mp);
   real_to_gmp(coord,*scale,idx_c2,a32_mp);
   real_to_gmp(coord,*scale,idx_c3,a33_mp);
   real_to_gmp(coord,*scale,idx_d1,a41_mp);
   real_to_gmp(coord,*scale,idx_d2,a42_mp);
   real_to_gmp(coord,*scale,idx_d3,a43_mp);
   real_to_gmp(coord,*scale,idx_e1,a51_mp);
   real_to_gmp(coord,*scale,idx_e2,a52_mp);
   real_to_gmp(coord,*scale,idx_e3,a53_mp);

   real_to_gmp(radius,*scale,idx_a4,r1_mp);
   real_to_gmp(radius,*scale,idx_b4,r2_mp);
   real_to_gmp(radius,*scale,idx_c4,r3_mp);
   real_to_gmp(radius,*scale,idx_d4,r4_mp);
   real_to_gmp(radius,*scale,idx_e4,r5_mp);

   build_weight(a11_mp,a12_mp,a13_mp,r1_mp,a14_mp);
   build_weight(a21_mp,a22_mp,a23_mp,r2_mp,a24_mp);
   build_weight(a31_mp,a32_mp,a33_mp,r3_mp,a34_mp);
   build_weight(a41_mp,a42_mp,a43_mp,r4_mp,a44_mp);
   build_weight(a51_mp,a52_mp,a53_mp,r5_mp,a54_mp);

/* Compute determinant */

   deter5_gmp(temp1,a11_mp,a12_mp,a13_mp,a14_mp,a21_mp,a22_mp,
     a23_mp,a24_mp,a31_mp,a32_mp,a33_mp,a34_mp,a41_mp,a42_mp,
      a43_mp,a44_mp,a51_mp,a52_mp,a53_mp,a54_mp);

   icomp = mpz_sgn(temp1);

/* if major determinant is non 0, return its sign */

   if (icomp != 0) {
      *res = icomp;
      return;
   }

/* Look now at each term in the expansion of the determinant with
   *respect to EPS                                                
   The initial determinant is:
   Minor5(i,j,k,l,m,1,2,3,4,0)         */

/* Term 1:    -Minor4(j,k,l,m,1,2,3,0) */

   deter4_gmp(temp1,a21_mp,a22_mp,a23_mp,a31_mp,a32_mp,a33_mp,
      a41_mp,a42_mp,a43_mp,a51_mp,a52_mp,a53_mp);
   icomp = mpz_sgn(temp1);
   if (icomp != 0) {
      *res = - icomp;
      return;
   }

/* Term 2:   Minor4(j,k,l,m,1,2,4,0) */
   
   deter4_gmp(temp1,a21_mp,a22_mp,a24_mp,a31_mp,a32_mp,a34_mp,
      a41_mp,a42_mp,a44_mp,a51_mp,a52_mp,a54_mp);
   icomp = mpz_sgn(temp1);
   if (icomp != 0) {
      *res =  icomp;
      return;
   }

/* Term 3:   - Minor4(j,k,l,m,1,3,4,0) */
   
   deter4_gmp(temp1,a21_mp,a23_mp,a24_mp,a31_mp,a33_mp,a34_mp,
      a41_mp,a43_mp,a44_mp,a51_mp,a53_mp,a54_mp);
   icomp = mpz_sgn(temp1);
   if (icomp != 0) {
      *res =  - icomp;
      return;
   }

/* Term 4:   Minor4(j,k,l,m,2,3,4,0) */
   
   deter4_gmp(temp1,a22_mp,a23_mp,a24_mp,a32_mp,a33_mp,a34_mp,
      a42_mp,a43_mp,a44_mp,a52_mp,a53_mp,a54_mp);
   icomp = mpz_sgn(temp1);
   if (icomp != 0) {
      *res =  icomp;
      return;
   }

/* Term 5:   Minor4(i,k,l,m,1,2,3,0) */
   
   deter4_gmp(temp1,a11_mp,a12_mp,a13_mp,a31_mp,a32_mp,a33_mp,
      a41_mp,a42_mp,a43_mp,a51_mp,a52_mp,a53_mp);
   icomp = mpz_sgn(temp1);
   if (icomp != 0) {
      *res =  icomp;
      return;
   }

/* Term 6:   Minor3(k,l,m,1,2,0) */
   
   deter3_gmp(temp1,a31_mp,a32_mp,a41_mp,a42_mp,a51_mp,a52_mp);
   icomp = mpz_sgn(temp1);
   if (icomp != 0) {
      *res =  icomp;
      return;
   }

/* Term 7:   -Minor3(k,l,m,1,3,0) */
   
   deter3_gmp(temp1,a31_mp,a33_mp,a41_mp,a43_mp,a51_mp,a53_mp);
   icomp = mpz_sgn(temp1);
   if (icomp != 0) {
      *res =  -icomp;
      return;
   }

/* Term 8:   Minor3(k,l,m,2,3,0) */
   
   deter3_gmp(temp1,a32_mp,a33_mp,a42_mp,a43_mp,a52_mp,a53_mp);
   icomp = mpz_sgn(temp1);
   if (icomp != 0) {
      *res =  icomp;
      return;
   }

/* Term 9:   -Minor4(i,k,l,m,1,2,4,0) */
   
   deter4_gmp(temp1,a11_mp,a12_mp,a14_mp,a31_mp,a32_mp,a34_mp,
      a41_mp,a42_mp,a44_mp,a51_mp,a52_mp,a54_mp);
   icomp = mpz_sgn(temp1);
   if (icomp != 0) {
      *res =  -icomp;
      return;
   }

/* Term 10:   Minor3(k,l,m,1,4,0) */
   
   deter3_gmp(temp1,a31_mp,a34_mp,a41_mp,a44_mp,a51_mp,a54_mp);
   icomp = mpz_sgn(temp1);
   if (icomp != 0) {
      *res =  icomp;
      return;
   }

/* Term 11:   -Minor3(k,l,m,2,4,0) */
   
   deter3_gmp(temp1,a32_mp,a34_mp,a42_mp,a44_mp,a52_mp,a54_mp);
   icomp = mpz_sgn(temp1);
   if (icomp != 0) {
      *res =  -icomp;
      return;
   }

/* Term 12:   Minor4(i,k,l,m,1,3,4,0) */
   
   deter4_gmp(temp1,a11_mp,a13_mp,a14_mp,a31_mp,a33_mp,a34_mp,
      a41_mp,a43_mp,a44_mp,a51_mp,a53_mp,a54_mp);
   icomp = mpz_sgn(temp1);
   if (icomp != 0) {
      *res =  icomp;
      return;
   }

/* Term 13:   Minor3(k,l,m,3,4,0) */
   
   deter3_gmp(temp1,a33_mp,a34_mp,a43_mp,a44_mp,a53_mp,a54_mp);
   icomp = mpz_sgn(temp1);
   if (icomp != 0) {
      *res =  icomp;
      return;
   }

/* Term 14:   -Minor4(i,k,l,m,2,3,4,0) */
   
   deter4_gmp(temp1,a12_mp,a13_mp,a14_mp,a32_mp,a33_mp,a34_mp,
      a42_mp,a43_mp,a44_mp,a52_mp,a53_mp,a54_mp);
   icomp = mpz_sgn(temp1);
   if (icomp != 0) {
      *res =  -icomp;
      return;
   }

/* Term 15:   -Minor4(i,j,l,m,1,2,3,0) */
   
   deter4_gmp(temp1,a11_mp,a12_mp,a13_mp,a21_mp,a22_mp,a23_mp,
      a41_mp,a42_mp,a43_mp,a51_mp,a52_mp,a53_mp);
   icomp = mpz_sgn(temp1);
   if (icomp != 0) {
      *res =  -icomp;
      return;
   }

/* Term 16:   -Minor3(j,l,m,1,2,0) */
   
   deter3_gmp(temp1,a21_mp,a22_mp,a41_mp,a42_mp,a51_mp,a52_mp);
   icomp = mpz_sgn(temp1);
   if (icomp != 0) {
      *res =  -icomp;
      return;
   }

/* Term 17:   Minor3(j,l,m,1,3,0) */
   
   deter3_gmp(temp1,a21_mp,a23_mp,a41_mp,a43_mp,a51_mp,a53_mp);
   icomp = mpz_sgn(temp1);
   if (icomp != 0) {
      *res =  icomp;
      return;
   }

/* Term 18:   -Minor3(j,l,m,2,3,0) */
   
   deter3_gmp(temp1,a22_mp,a23_mp,a42_mp,a43_mp,a52_mp,a53_mp);
   icomp = mpz_sgn(temp1);
   if (icomp != 0) {
      *res =  -icomp;
      return;
   }

/* Term 19:   Minor3(i,l,m,1,2,0) */
   
   deter3_gmp(temp1,a11_mp,a12_mp,a41_mp,a42_mp,a51_mp,a52_mp);
   icomp = mpz_sgn(temp1);
   if (icomp != 0) {
      *res =  icomp;
      return;
   }

/* Term 20:   -Minor2(l,m,1,0) */

   deter2_gmp(temp1,a41_mp,a51_mp);
   icomp = mpz_sgn(temp1);
   if (icomp != 0) {
      *res =  -icomp;
      return;
   }

/* Term 21:   Minor2(l,m,2,0) */

   deter2_gmp(temp1,a42_mp,a52_mp);
   icomp = mpz_sgn(temp1);
   if (icomp != 0) {
      *res =  icomp;
      return;
   }

/* Term 22:   -Minor3(i,l,m,1,3,0) */
   
   deter3_gmp(temp1,a11_mp,a13_mp,a41_mp,a43_mp,a51_mp,a53_mp);
   icomp = mpz_sgn(temp1);
   if (icomp != 0) {
      *res =  -icomp;
      return;
   }

/* Term 23:   -Minor2(l,m,3,0) */

   deter2_gmp(temp1,a43_mp,a53_mp);
   icomp = mpz_sgn(temp1);
   if (icomp != 0) {
      *res =  -icomp;
      return;
   }

/* Term 24:   Minor3(i,l,m,2,3,0) */
   
   deter3_gmp(temp1,a12_mp,a13_mp,a42_mp,a43_mp,a52_mp,a53_mp);
   icomp = mpz_sgn(temp1);
   if (icomp != 0) {
      *res =  icomp;
      return;
   }

/* Term 25:   Minor4(i,j,l,m,1,2,4,0) */
   
   deter4_gmp(temp1,a11_mp,a12_mp,a14_mp,a21_mp,a22_mp,a24_mp,
      a41_mp,a42_mp,a44_mp,a51_mp,a52_mp,a54_mp);
   icomp = mpz_sgn(temp1);
   if (icomp != 0) {
      *res =  icomp;
      return;
   }

/* Term 26:   -Minor3(j,l,m,1,4,0) */
   
   deter3_gmp(temp1,a21_mp,a24_mp,a41_mp,a44_mp,a51_mp,a54_mp);
   icomp = mpz_sgn(temp1);
   if (icomp != 0) {
      *res =  -icomp;
      return;
   }

/* Term 27:   Minor3(j,l,m,2,4,0) */
   
   deter3_gmp(temp1,a22_mp,a24_mp,a42_mp,a44_mp,a52_mp,a54_mp);
   icomp = mpz_sgn(temp1);
   if (icomp != 0) {
      *res =  icomp;
      return;
   }

/* Term 28:   Minor3(i,l,m,1,4,0) */
   
   deter3_gmp(temp1,a11_mp,a14_mp,a41_mp,a44_mp,a51_mp,a54_mp);
   icomp = mpz_sgn(temp1);
   if (icomp != 0) {
      *res =  icomp;
      return;
   }

/* Term 29:   Minor2(l,m,4,0) */

   deter2_gmp(temp1,a44_mp,a54_mp);
   icomp = mpz_sgn(temp1);
   if (icomp != 0) {
      *res =  icomp;
      return;
   }

/* Term 30:   -Minor3(i,l,m,2,4,0) */
   
   deter3_gmp(temp1,a12_mp,a14_mp,a42_mp,a44_mp,a52_mp,a54_mp);
   icomp = mpz_sgn(temp1);
   if (icomp != 0) {
      *res =  -icomp;
      return;
   }

/* Term 31:   -Minor4(i,j,l,m,1,3,4,0) */
   
   deter4_gmp(temp1,a11_mp,a13_mp,a14_mp,a21_mp,a23_mp,a24_mp,
      a41_mp,a43_mp,a44_mp,a51_mp,a53_mp,a54_mp);
   icomp = mpz_sgn(temp1);
   if (icomp != 0) {
      *res =  -icomp;
      return;
   }

/* Term 32:   -Minor3(j,l,m,3,4,0) */
   
   deter3_gmp(temp1,a23_mp,a24_mp,a43_mp,a44_mp,a53_mp,a54_mp);
   icomp = mpz_sgn(temp1);
   if (icomp != 0) {
      *res =  -icomp;
      return;
   }

/* Term 33:   Minor3(i,l,m,3,4,0) */
   
   deter3_gmp(temp1,a13_mp,a14_mp,a43_mp,a44_mp,a53_mp,a54_mp);
   icomp = mpz_sgn(temp1);
   if (icomp != 0) {
      *res =  icomp;
      return;
   }

/* Term 34:   Minor4(i,j,l,m,2,3,4,0) */
   
   deter4_gmp(temp1,a12_mp,a13_mp,a14_mp,a22_mp,a23_mp,a24_mp,
      a42_mp,a43_mp,a44_mp,a52_mp,a53_mp,a54_mp);
   icomp = mpz_sgn(temp1);
   if (icomp != 0) {
      *res =  icomp;
      return;
   }

/* Term 35:   Minor4(i,j,k,m,1,2,3,0) */
   
   deter4_gmp(temp1,a11_mp,a12_mp,a13_mp,a21_mp,a22_mp,a23_mp,
      a31_mp,a32_mp,a33_mp,a51_mp,a52_mp,a53_mp);
   icomp = mpz_sgn(temp1);
   if (icomp != 0) {
      *res =  icomp;
      return;
   }

/* Term 36:   Minor3(j,k,m,1,2,0) */
   
   deter3_gmp(temp1,a21_mp,a22_mp,a31_mp,a32_mp,a51_mp,a52_mp);
   icomp = mpz_sgn(temp1);
   if (icomp != 0) {
      *res =  icomp;
      return;
   }

/* Term 37:   -Minor3(j,k,m,1,3,0) */
   
   deter3_gmp(temp1,a21_mp,a23_mp,a31_mp,a33_mp,a51_mp,a53_mp);
   icomp = mpz_sgn(temp1);
   if (icomp != 0) {
      *res =  -icomp;
      return;
   }

/* Term 38:   Minor3(j,k,m,2,3,0) */
   
   deter3_gmp(temp1,a22_mp,a23_mp,a32_mp,a33_mp,a52_mp,a53_mp);
   icomp = mpz_sgn(temp1);
   if (icomp != 0) {
      *res =  icomp;
      return;
   }

/* Term 39:   -Minor3(i,k,m,1,2,0) */
   
   deter3_gmp(temp1,a11_mp,a12_mp,a31_mp,a32_mp,a51_mp,a52_mp);
   icomp = mpz_sgn(temp1);
   if (icomp != 0) {
      *res =  -icomp;
      return;
   }

/* Term 40:   Minor2(k,m,1,0) */

   deter2_gmp(temp1,a31_mp,a51_mp);
   icomp = mpz_sgn(temp1);
   if (icomp != 0) {
      *res =  icomp;
      return;
   }

/* Term 41:   -Minor2(k,m,2,0) */

   deter2_gmp(temp1,a32_mp,a52_mp);
   icomp = mpz_sgn(temp1);
   if (icomp != 0) {
      *res =  -icomp;
      return;
   }

/* Term 42:   Minor3(i,k,m,1,3,0) */
   
   deter3_gmp(temp1,a11_mp,a13_mp,a31_mp,a33_mp,a51_mp,a53_mp);
   icomp = mpz_sgn(temp1);
   if (icomp != 0) {
      *res =  icomp;
      return;
   }

/* Term 43:   Minor2(k,m,3,0) */

   deter2_gmp(temp1,a33_mp,a53_mp);
   icomp = mpz_sgn(temp1);
   if (icomp != 0) {
      *res =  icomp;
      return;
   }

/* Term 44:   -Minor3(i,k,m,2,3,0) */
   
   deter3_gmp(temp1,a12_mp,a13_mp,a32_mp,a33_mp,a52_mp,a53_mp);
   icomp = mpz_sgn(temp1);
   if (icomp != 0) {
      *res =  -icomp;
      return;
   }

/* Term 45:   Minor3(i,j,m,1,2,0) */
   
   deter3_gmp(temp1,a11_mp,a12_mp,a21_mp,a22_mp,a51_mp,a52_mp);
   icomp = mpz_sgn(temp1);
   if (icomp != 0) {
      *res =  icomp;
      return;
   }

/* Term 46:   -Minor2(j,m,1,0) */

   deter2_gmp(temp1,a21_mp,a51_mp);
   icomp = mpz_sgn(temp1);
   if (icomp != 0) {
      *res =  -icomp;
      return;
   }

/* Term 47:   Minor2(j,m,2,0) */

   deter2_gmp(temp1,a22_mp,a52_mp);
   icomp = mpz_sgn(temp1);
   if (icomp != 0) {
      *res =  icomp;
      return;
   }

/* Term 48:   Minor2(i,m,1,0) */

   deter2_gmp(temp1,a11_mp,a51_mp);
   icomp = mpz_sgn(temp1);
   if (icomp != 0) {
      *res =  icomp;
      return;
   }

/* Term 49:   1 */

   *res = 1;

}
#ifdef __cplusplus
}
#endif

/*------------------------------------------------------------------------*/
/* minor4_:                          */
/*                             
   This subroutine tests the sign of the determinant
      D = | a11 a12 a13 1 |
          | a21 a22 a23 1 |
          | a31 a32 a33 1 |
          | a41 a42 a43 1 |
   and return 1 if positive, -1 if negative, 0 otherwise   */
#ifdef __cplusplus
extern "C" {
#endif
void F77Name2(minor4_gmp)(double *coord, double *scale, int *a, int *b, int *c, int *d, int *res)

{
   int icomp;
   int idx_a1,idx_a2,idx_a3,idx_b1,idx_b2,idx_b3;
   int idx_c1,idx_c2,idx_c3,idx_d1,idx_d2,idx_d3;

/* Initialise local gmp variables */

/* Transfer coordinates to gmp */

   idx_a1 = (*a)*3 -3; idx_a2 = idx_a1+1; idx_a3 = idx_a2 +1;
   idx_b1 = (*b)*3 -3; idx_b2 = idx_b1+1; idx_b3 = idx_b2 +1;
   idx_c1 = (*c)*3 -3; idx_c2 = idx_c1+1; idx_c3 = idx_c2 +1;
   idx_d1 = (*d)*3 -3; idx_d2 = idx_d1+1; idx_d3 = idx_d2 +1;

   real_to_gmp(coord,*scale,idx_a1,a11_mp);
   real_to_gmp(coord,*scale,idx_a2,a12_mp);
   real_to_gmp(coord,*scale,idx_a3,a13_mp);
   real_to_gmp(coord,*scale,idx_b1,a21_mp);
   real_to_gmp(coord,*scale,idx_b2,a22_mp);
   real_to_gmp(coord,*scale,idx_b3,a23_mp);
   real_to_gmp(coord,*scale,idx_c1,a31_mp);
   real_to_gmp(coord,*scale,idx_c2,a32_mp);
   real_to_gmp(coord,*scale,idx_c3,a33_mp);
   real_to_gmp(coord,*scale,idx_d1,a41_mp);
   real_to_gmp(coord,*scale,idx_d2,a42_mp);
   real_to_gmp(coord,*scale,idx_d3,a43_mp);

/* Compute determinant */

   deter4_gmp(temp1,a11_mp,a12_mp,a13_mp,a21_mp,a22_mp,a23_mp,
            a31_mp,a32_mp,a33_mp,a41_mp,a42_mp,a43_mp);

   icomp = mpz_sgn(temp1);

   *res = icomp;

}
#ifdef __cplusplus
}
#endif
