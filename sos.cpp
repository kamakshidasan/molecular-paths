/***************************************************************************
 *   Copyright (C) 2010 by raghavendra,,,                                  *
 *   raghavendra@incognito                                                 *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include "sos.h"

Sos::Sos()
{
}

Sos::~Sos()
{
}

/*!
    \fn Sos::Deter2(mpz_t deter,mpz_t b11,mpz_t b21)
 */
void Sos::Deter2(mpz_t deter,mpz_t b11,mpz_t b21)
{
        mpz_sub(deter,b11,b21);
}

/*!
    \fn Sos::Deter3(mpz_t deter, mpz_t b11, mpz_t b12, mpz_t b21, mpz_t b22, mpz_t b31, mpz_t b32)
 */
void Sos::Deter3(mpz_t deter, mpz_t b11, mpz_t b12, mpz_t b21, mpz_t b22, mpz_t b31, mpz_t b32)
{
        mpz_t tmp1,tmp2,tmp3,val1,val2,tmp4;

        mpz_init(tmp1); mpz_init(tmp2); mpz_init(tmp3);
        mpz_init(val1); mpz_init(val2); mpz_init(tmp4);

        mpz_sub(tmp1,b21,b11);
        mpz_sub(tmp2,b22,b12);
        mpz_sub(tmp3,b31,b11);
        mpz_sub(tmp4,b32,b12);

        mpz_mul(val1,tmp1,tmp4);
        mpz_mul(val2,tmp2,tmp3);

        mpz_sub(deter,val1,val2);

        mpz_clear(tmp1); mpz_clear(tmp2); mpz_clear(tmp3);
        mpz_clear(val1); mpz_clear(val2); mpz_clear(tmp4);
}

/*!
    \fn Sos::Deter4(mpz_t deter, mpz_t b11, mpz_t b12, mpz_t b13,mpz_t b21, mpz_t b22, mpz_t b23,mpz_t b31, mpz_t b32, mpz_t b33,mpz_t b41, mpz_t b42, mpz_t b43)
 */
void Sos::Deter4(mpz_t deter, mpz_t b11, mpz_t b12, mpz_t b13,mpz_t b21, mpz_t b22, mpz_t b23,mpz_t b31, mpz_t b32, mpz_t b33,mpz_t b41, mpz_t b42, mpz_t b43)
{
        mpz_t c11,c12,c13,c21,c22,c23,c31,c32,c33;
        mpz_t val1,val2,val3,tmp1,tmp2,tmp3;

        mpz_init(c11); mpz_init(c12); mpz_init(c13); mpz_init(c21);
        mpz_init(c22); mpz_init(c23); mpz_init(c31); mpz_init(c32);
        mpz_init(c33);
        mpz_init(tmp1); mpz_init(tmp2); mpz_init(tmp3);
        mpz_init(val1); mpz_init(val2); mpz_init(val3);

        mpz_sub(c11,b21,b11);mpz_sub(c12,b22,b12);mpz_sub(c13,b23,b13);
        mpz_sub(c21,b31,b11);mpz_sub(c22,b32,b12);mpz_sub(c23,b33,b13);
        mpz_sub(c31,b41,b11);mpz_sub(c32,b42,b12);mpz_sub(c33,b43,b13);

        mpz_mul(tmp1,c22,c33);mpz_mul(tmp2,c32,c23);mpz_sub(val1,tmp1,tmp2);
        mpz_mul(tmp1,c12,c33);mpz_mul(tmp2,c32,c13);mpz_sub(val2,tmp1,tmp2);
        mpz_mul(tmp1,c12,c23);mpz_mul(tmp2,c22,c13);mpz_sub(val3,tmp1,tmp2);

        mpz_mul(tmp1,c21,val2);mpz_mul(tmp2,c11,val1);mpz_mul(tmp3,c31,val3);

        mpz_add(val1,tmp2,tmp3);
        mpz_sub(deter,tmp1,val1);

        mpz_clear(tmp1); mpz_clear(tmp2); mpz_clear(tmp3);
        mpz_clear(val1); mpz_clear(val2); mpz_clear(val3);
        mpz_clear(c11); mpz_clear(c12); mpz_clear(c13); mpz_clear(c21);
        mpz_clear(c22); mpz_clear(c23); mpz_clear(c31); mpz_clear(c32);
        mpz_clear(c33);
}

/*!
    \fn Sos::Deter5(mpz_t deter, mpz_t b11, mpz_t b12, mpz_t b13, mpz_t b14,mpz_t b21, mpz_t b22, mpz_t b23, mpz_t b24,mpz_t b31, mpz_t b32, mpz_t b33, mpz_t b34,mpz_t b41, mpz_t b42, mpz_t b43, mpz_t b44,mpz_t b51, mpz_t b52, mpz_t b53, mpz_t b54)
 */
void Sos::Deter5(mpz_t deter, mpz_t b11, mpz_t b12, mpz_t b13, mpz_t b14,mpz_t b21, mpz_t b22, mpz_t b23, mpz_t b24,mpz_t b31, mpz_t b32, mpz_t b33, mpz_t b34,mpz_t b41, mpz_t b42, mpz_t b43, mpz_t b44,mpz_t b51, mpz_t b52, mpz_t b53, mpz_t b54)
{
        mpz_t c11,c12,c13,c14,c21,c22,c23,c24,c31,c32,c33,c34;
        mpz_t c41,c42,c43,c44;
        mpz_t d1,d2,d3,e1,e2,e3,f1,f2,f3,g1,g2,g3;
        mpz_t tmp1,tmp2,tmp3;

        /*FILE * fp = fopen("gmp","w");

        gmp_fprintf(fp,"b11 = %Zd\n",b11);
        gmp_fprintf(fp,"b12 = %Zd\n",b12);
        gmp_fprintf(fp,"b13 = %Zd\n",b13);
        gmp_fprintf(fp,"b14 = %Zd\n",b14);
        gmp_fprintf(fp,"b21 = %Zd\n",b21);
        gmp_fprintf(fp,"b22 = %Zd\n",b22);
        gmp_fprintf(fp,"b23 = %Zd\n",b23);
        gmp_fprintf(fp,"b24 = %Zd\n",b24);
        gmp_fprintf(fp,"b31 = %Zd\n",b31);
        gmp_fprintf(fp,"b32 = %Zd\n",b32);
        gmp_fprintf(fp,"b33 = %Zd\n",b33);
        gmp_fprintf(fp,"b34 = %Zd\n",b34);
        gmp_fprintf(fp,"b41 = %Zd\n",b41);
        gmp_fprintf(fp,"b42 = %Zd\n",b42);
        gmp_fprintf(fp,"b43 = %Zd\n",b43);
        gmp_fprintf(fp,"b44 = %Zd\n",b44);
        gmp_fprintf(fp,"b51 = %Zd\n",b51);
        gmp_fprintf(fp,"b52 = %Zd\n",b52);
        gmp_fprintf(fp,"b53 = %Zd\n",b53);
        gmp_fprintf(fp,"b54 = %Zd\n",b54);

        fclose(fp);*/

        mpz_init(c11); mpz_init(c12); mpz_init(c13); mpz_init(c14);
        mpz_init(c21); mpz_init(c22); mpz_init(c23); mpz_init(c24);
        mpz_init(c31); mpz_init(c32); mpz_init(c33); mpz_init(c34);
        mpz_init(c41); mpz_init(c42); mpz_init(c43); mpz_init(c44);
        mpz_init(d1); mpz_init(d2); mpz_init(d3);
        mpz_init(e1); mpz_init(e2); mpz_init(e3);
        mpz_init(f1); mpz_init(f2); mpz_init(f3);
        mpz_init(g1); mpz_init(g2); mpz_init(g3);
        mpz_init(tmp1); mpz_init(tmp2); mpz_init(tmp3);

        mpz_sub(c11,b21,b11); mpz_sub(c12,b22,b12); mpz_sub(c13,b23,b13);
        mpz_sub(c14,b24,b14);
        mpz_sub(c21,b31,b11); mpz_sub(c22,b32,b12); mpz_sub(c23,b33,b13);
        mpz_sub(c24,b34,b14);
        mpz_sub(c31,b41,b11); mpz_sub(c32,b42,b12); mpz_sub(c33,b43,b13);
        mpz_sub(c34,b44,b14);
        mpz_sub(c41,b51,b11); mpz_sub(c42,b52,b12); mpz_sub(c43,b53,b13);
        mpz_sub(c44,b54,b14);

        mpz_mul(tmp1,c32,c43); mpz_mul(tmp2,c42,c33); mpz_sub(d1,tmp1,tmp2);
        mpz_mul(tmp1,c32,c44); mpz_mul(tmp2,c42,c34); mpz_sub(d2,tmp1,tmp2);
        mpz_mul(tmp1,c33,c44); mpz_mul(tmp2,c43,c34); mpz_sub(d3,tmp1,tmp2);

        mpz_mul(tmp1,c12,c23); mpz_mul(tmp2,c22,c13); mpz_sub(e1,tmp1,tmp2);
        mpz_mul(tmp1,c12,c24); mpz_mul(tmp2,c22,c14); mpz_sub(e2,tmp1,tmp2);
        mpz_mul(tmp1,c13,c24); mpz_mul(tmp2,c23,c14); mpz_sub(e3,tmp1,tmp2);

        mpz_mul(tmp1,c11,c24); mpz_mul(tmp2,c21,c14); mpz_sub(f1,tmp1,tmp2);
        mpz_mul(tmp1,c11,c23); mpz_mul(tmp2,c21,c13); mpz_sub(f2,tmp1,tmp2);
        mpz_mul(tmp1,c11,c22); mpz_mul(tmp2,c21,c12); mpz_sub(f3,tmp1,tmp2);

        mpz_mul(tmp1,c31,c44); mpz_mul(tmp2,c41,c34); mpz_sub(g1,tmp1,tmp2);
        mpz_mul(tmp1,c31,c43); mpz_mul(tmp2,c41,c33); mpz_sub(g2,tmp1,tmp2);
        mpz_mul(tmp1,c31,c42); mpz_mul(tmp2,c41,c32); mpz_sub(g3,tmp1,tmp2);

        mpz_mul(tmp1,e3,g3); mpz_mul(tmp2,e2,g2); mpz_sub(tmp3,tmp1,tmp2);
        mpz_mul(tmp1,e1,g1); mpz_add(tmp3,tmp3,tmp1);
        mpz_mul(tmp1,d3,f3); mpz_add(tmp3,tmp3,tmp1);
        mpz_mul(tmp1,d2,f2); mpz_sub(tmp3,tmp3,tmp1);
        mpz_mul(tmp1,d1,f1); mpz_add(deter,tmp3,tmp1);

        mpz_clear(c11); mpz_clear(c12); mpz_clear(c13); mpz_clear(c14);
        mpz_clear(c21); mpz_clear(c22); mpz_clear(c23); mpz_clear(c24);
        mpz_clear(c31); mpz_clear(c32); mpz_clear(c33); mpz_clear(c34);
        mpz_clear(c41); mpz_clear(c42); mpz_clear(c43); mpz_clear(c44);
        mpz_clear(d1); mpz_clear(d2); mpz_clear(d3);
        mpz_clear(e1); mpz_clear(e2); mpz_clear(e3);
        mpz_clear(f1); mpz_clear(f2); mpz_clear(f3);
        mpz_clear(g1); mpz_clear(g2); mpz_clear(g3);
        mpz_clear(tmp1); mpz_clear(tmp2); mpz_clear(tmp3);
}

/*!
    \fn Sos::Minor2(std::vector<Vertex> & vertexList, int *a, int *b, int *ia, int *res)
 */
void Sos::Minor2(std::vector<Vertex> & vertexList, int *a, int *b, int ia, int *res)
{
        int icomp;

        /* Initialise local GMP variables */

        mpz_t temp1;
        mpz_init(temp1);

        Deter2(temp1,vertexList[*a].V[ia],vertexList[*b].V[ia]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = icomp;
        }
        else
        {
                *res = 1;
        }
        mpz_clear(temp1);
}

/*!
    \fn Sos::Minor3(std::vector<Vertex> & vertexList,  int *a,  int *b,  int *c,  int *i1,  int *i2,  int *res)
 */
void Sos::Minor3(std::vector<Vertex> & vertexList,  int *a,  int *b,  int *c,  int *i1,  int *i2,  int *res)
{
        int icomp;

        mpz_t temp1;
        mpz_init(temp1);

        //Compute determinant:Minor3(i,j,k,1,2,0)
        Deter3(temp1,   vertexList[*a].V[*i1], vertexList[*a].V[*i2],
                        vertexList[*b].V[*i1], vertexList[*b].V[*i2],
                        vertexList[*c].V[*i1], vertexList[*c].V[*i2]);

        // if major determinant is non 0, return its sign
        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = icomp;
                mpz_clear(temp1);
                return;
        }

        //Look now at each term in the expansion of the determinant with respect to EPS

        //Term 1: - Minor2(j,k,1,0)
        Deter2(temp1, vertexList[*b].V[*i1], vertexList[*c].V[*i1]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = -icomp;
                mpz_clear(temp1);
                return;
        }

        //Term 2: Minor2(j,k,2,0)
        Deter2(temp1, vertexList[*b].V[*i2], vertexList[*c].V[*i2]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = icomp;
                mpz_clear(temp1);
                return;
        }

        //Term 3: Minor2(i,k,1,0)
        Deter2(temp1, vertexList[*a].V[*i1], vertexList[*c].V[*i1]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = icomp;
                mpz_clear(temp1);
                return;
        }

        //Term 4: 1
        *res = 1;
        mpz_clear(temp1);
        return;
}

/*!
    \fn Sos::Minor4(std::vector<Vertex> & vertexList,  int *a,  int *b,  int *c,  int *d,  int *res)
 */
void Sos::Minor4(std::vector<Vertex> & vertexList,  int *a,  int *b,  int *c,  int *d,  int *res)
{
        int icomp;

        mpz_t temp1;
        mpz_init(temp1);

        //Compute determinant:Minor4(i,j,k,l,1,2,3,0)
        Deter4(temp1,   vertexList[*a].V[1], vertexList[*a].V[2], vertexList[*a].V[3],
                        vertexList[*b].V[1], vertexList[*b].V[2], vertexList[*b].V[3],
                        vertexList[*c].V[1], vertexList[*c].V[2], vertexList[*c].V[3],
                        vertexList[*d].V[1], vertexList[*d].V[2], vertexList[*d].V[3]);

        // if major determinant is non 0, return its sign
        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = icomp;
                mpz_clear(temp1);
                return;
        }

        //Look now at each term in the expansion of the determinant with *respect to EPS

        //Term 1:Minor3(j,k,l,1,2,0)
        Deter3(temp1,   vertexList[*b].V[1], vertexList[*b].V[2],
                        vertexList[*c].V[1], vertexList[*c].V[2],
                        vertexList[*d].V[1], vertexList[*d].V[2]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = icomp;
                mpz_clear(temp1);
                return;
        }

        //Term 2:-Minor3(j,k,l,1,3,0)
        Deter3(temp1,   vertexList[*b].V[1], vertexList[*b].V[3],
                        vertexList[*c].V[1], vertexList[*c].V[3],
                        vertexList[*d].V[1], vertexList[*d].V[3]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = -icomp;
                mpz_clear(temp1);
                return;
        }

        //Term 3:Minor3(j,k,l,2,3,0)
        Deter3(temp1,   vertexList[*b].V[2], vertexList[*b].V[3],
                        vertexList[*c].V[2], vertexList[*c].V[3],
                        vertexList[*d].V[2], vertexList[*d].V[3]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = icomp;
                mpz_clear(temp1);
                return;
        }

        //Term 4:- Minor3(i,k,l,1,2,0)
        Deter3(temp1,   vertexList[*a].V[1],vertexList[*a].V[2],
                        vertexList[*c].V[1],vertexList[*c].V[2],
                        vertexList[*d].V[1],vertexList[*d].V[2]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = -icomp;
                mpz_clear(temp1);
                return;
        }

        //Term 5:Minor2(k,l,1,0)
        Deter2(temp1,vertexList[*c].V[1],vertexList[*d].V[1]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = icomp;
                mpz_clear(temp1);
                return;
        }

        //Term 6:-Minor2(k,l,2,0)
        Deter2(temp1,vertexList[*c].V[2],vertexList[*d].V[2]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = -icomp;
                mpz_clear(temp1);
                return;
        }

        //Term 7:Minor3(i,k,l,1,3,0)
        Deter3(temp1,   vertexList[*a].V[1],vertexList[*a].V[3],
                        vertexList[*c].V[1],vertexList[*c].V[3],
                        vertexList[*d].V[1],vertexList[*d].V[3]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = icomp;
                mpz_clear(temp1);
                return;
        }

        //Term 8:Minor2(k,l,3,0)
        Deter2(temp1,vertexList[*c].V[3],vertexList[*d].V[3]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = icomp;
                mpz_clear(temp1);
                return;
        }

        //Term 9:- Minor3(i,k,l,2,3,0)
        Deter3(temp1,   vertexList[*a].V[2],vertexList[*a].V[3],
                        vertexList[*c].V[2],vertexList[*c].V[3],
                        vertexList[*d].V[2],vertexList[*d].V[3]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = -icomp;
                mpz_clear(temp1);
                return;
        }

        //Term 10:Minor3(i,j,l,1,2,0)
        Deter3(temp1,   vertexList[*a].V[1],vertexList[*a].V[2],
                        vertexList[*b].V[1],vertexList[*b].V[2],
                        vertexList[*d].V[1],vertexList[*d].V[2]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = icomp;
                mpz_clear(temp1);
                return;
        }

        //Term 11: - Minor2(j,l,1,0)
        Deter2(temp1,vertexList[*b].V[1],vertexList[*d].V[1]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = -icomp;
                mpz_clear(temp1);
                return;
        }

        //Term 12:Minor2(j,l,2,0)
        Deter2(temp1, vertexList[*b].V[2],vertexList[*d].V[2]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = icomp;
                mpz_clear(temp1);
                return;
        }

        //Term 13:Minor2(i,l,1,0)
        Deter2(temp1,vertexList[*a].V[1],vertexList[*d].V[1]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = icomp;
                mpz_clear(temp1);
                return;
        }

        //Term 14:1
        *res = 1;
        mpz_clear(temp1);
        return;
}

/*!
    \fn Sos::Minor5(std::vector<Vertex> & vertexList, int *a, int *b, int *c, int *d, int *e, int *res)
 */
void Sos::Minor5(std::vector<Vertex> & vertexList, int *a, int *b, int *c, int *d, int *e, int *res)
{
        int icomp = 0;
        mpz_t temp1;mpz_init(temp1);

        //Compute determinant:Minor5(i,j,k,l,m,1,2,3,4,0)
        Deter5(temp1,   vertexList[*a].V[1],vertexList[*a].V[2],vertexList[*a].V[3],vertexList[*a].V[4],
                        vertexList[*b].V[1],vertexList[*b].V[2],vertexList[*b].V[3],vertexList[*b].V[4],
                        vertexList[*c].V[1],vertexList[*c].V[2],vertexList[*c].V[3],vertexList[*c].V[4],
                        vertexList[*d].V[1],vertexList[*d].V[2],vertexList[*d].V[3],vertexList[*d].V[4],
                        vertexList[*e].V[1],vertexList[*e].V[2],vertexList[*e].V[3],vertexList[*e].V[4]);

        // if major determinant is non 0, return its sign
        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = icomp;
                mpz_clear(temp1);
                return;
        }

        //Look now at each term in the expansion of the determinant with *respect to EPS

        //Term 1: 	-Minor4(j,k,l,m,1,2,3,0)
        Deter4(temp1,   vertexList[*b].V[1],vertexList[*b].V[2],vertexList[*b].V[3],
                        vertexList[*c].V[1],vertexList[*c].V[2],vertexList[*c].V[3],
                        vertexList[*d].V[1],vertexList[*d].V[2],vertexList[*d].V[3],
                        vertexList[*e].V[1],vertexList[*e].V[2],vertexList[*e].V[3]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = -icomp;
                mpz_clear(temp1);
                return;
        }

        //Term 2:	Minor4(j,k,l,m,1,2,4,0)
        Deter4(temp1,   vertexList[*b].V[1],vertexList[*b].V[2],vertexList[*b].V[4],
                        vertexList[*c].V[1],vertexList[*c].V[2],vertexList[*c].V[4],
                        vertexList[*d].V[1],vertexList[*d].V[2],vertexList[*d].V[4],
                        vertexList[*e].V[1],vertexList[*e].V[2],vertexList[*e].V[4]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = icomp;
                mpz_clear(temp1);
                return;
        }

        //Term 3:	- Minor4(j,k,l,m,1,3,4,0)
        Deter4(temp1,   vertexList[*b].V[1],vertexList[*b].V[3],vertexList[*b].V[4],
                        vertexList[*c].V[1],vertexList[*c].V[3],vertexList[*c].V[4],
                        vertexList[*d].V[1],vertexList[*d].V[3],vertexList[*d].V[4],
                        vertexList[*e].V[1],vertexList[*e].V[3],vertexList[*e].V[4]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = -icomp;
                mpz_clear(temp1);
                return;
        }

        //Term 4:	Minor4(j,k,l,m,2,3,4,0)
        Deter4(temp1,   vertexList[*b].V[2],vertexList[*b].V[3],vertexList[*b].V[4],
                        vertexList[*c].V[2],vertexList[*c].V[3],vertexList[*c].V[4],
                        vertexList[*d].V[2],vertexList[*d].V[3],vertexList[*d].V[4],
                        vertexList[*e].V[2],vertexList[*e].V[3],vertexList[*e].V[4]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = icomp;
                mpz_clear(temp1);
                return;
        }

        //Term 5:	Minor4(i,k,l,m,1,2,3,0)
        Deter4(temp1,   vertexList[*a].V[1],vertexList[*a].V[2],vertexList[*a].V[3],
                        vertexList[*c].V[1],vertexList[*c].V[2],vertexList[*c].V[3],
                        vertexList[*d].V[1],vertexList[*d].V[2],vertexList[*d].V[3],
                        vertexList[*e].V[1],vertexList[*e].V[2],vertexList[*e].V[3]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = icomp;
                mpz_clear(temp1);
                return;
        }

            //Term 6:	Minor3(k,l,m,1,2,0)
        Deter3(temp1,   vertexList[*c].V[1],vertexList[*c].V[2],
                        vertexList[*d].V[1],vertexList[*d].V[2],
                        vertexList[*e].V[1],vertexList[*e].V[2]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = icomp;
                mpz_clear(temp1);
                return;
        }

        //Term 7:	-Minor3(k,l,m,1,3,0)
        Deter3(temp1,   vertexList[*c].V[1],vertexList[*c].V[3],
                        vertexList[*d].V[1],vertexList[*d].V[3],
                        vertexList[*e].V[1],vertexList[*e].V[3]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = -icomp;
                mpz_clear(temp1);
                return;
        }

        //Term 8:	Minor3(k,l,m,2,3,0)
        Deter3(temp1,   vertexList[*c].V[2],vertexList[*c].V[3],
                        vertexList[*d].V[2],vertexList[*d].V[3],
                        vertexList[*e].V[2],vertexList[*e].V[3]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = icomp;
                mpz_clear(temp1);
                return;
        }

        //Term 9:	-Minor4(i,k,l,m,1,2,4,0)
        Deter4(temp1,   vertexList[*a].V[1],vertexList[*a].V[2],vertexList[*a].V[4],
                        vertexList[*c].V[1],vertexList[*c].V[2],vertexList[*c].V[4],
                        vertexList[*d].V[1],vertexList[*d].V[2],vertexList[*d].V[4],
                        vertexList[*e].V[1],vertexList[*e].V[2],vertexList[*e].V[4]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = -icomp;
                mpz_clear(temp1);
                return;
        }

        //Term 10:	Minor3(k,l,m,1,4,0)
        Deter3(temp1,   vertexList[*c].V[1],vertexList[*c].V[4],
                        vertexList[*d].V[1],vertexList[*d].V[4],
                        vertexList[*e].V[1],vertexList[*e].V[4]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = icomp;
                mpz_clear(temp1);
                return;
        }

        //Term 11:	-Minor3(k,l,m,2,4,0)
        Deter3(temp1,   vertexList[*c].V[2],vertexList[*c].V[4],
                        vertexList[*d].V[2],vertexList[*d].V[4],
                        vertexList[*e].V[2],vertexList[*e].V[4]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = -icomp;
                mpz_clear(temp1);
                return;
        }

        //Term 12:	Minor4(i,k,l,m,1,3,4,0)
        Deter4(temp1,   vertexList[*a].V[1],vertexList[*a].V[3],vertexList[*a].V[4],
                        vertexList[*c].V[1],vertexList[*c].V[3],vertexList[*c].V[4],
                        vertexList[*d].V[1],vertexList[*d].V[3],vertexList[*d].V[4],
                        vertexList[*e].V[1],vertexList[*e].V[3],vertexList[*e].V[4]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = icomp;
                mpz_clear(temp1);
                return;
        }

        //Term 13:	Minor3(k,l,m,3,4,0)
        Deter3(temp1,   vertexList[*c].V[3],vertexList[*c].V[4],
                        vertexList[*d].V[3],vertexList[*d].V[4],
                        vertexList[*e].V[3],vertexList[*e].V[4]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = icomp;
                mpz_clear(temp1);
                return;
        }

        //Term 14:	-Minor4(i,k,l,m,2,3,4,0)
        Deter4(temp1,   vertexList[*a].V[2],vertexList[*a].V[3],vertexList[*a].V[4],
                        vertexList[*c].V[2],vertexList[*c].V[3],vertexList[*c].V[4],
                        vertexList[*d].V[2],vertexList[*d].V[3],vertexList[*d].V[4],
                        vertexList[*e].V[2],vertexList[*e].V[3],vertexList[*e].V[4]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = -icomp;
                mpz_clear(temp1);
                return;
        }

        //Term 15:	-Minor4(i,j,l,m,1,2,3,0)
        Deter4(temp1,   vertexList[*a].V[1],vertexList[*a].V[2],vertexList[*a].V[3],
                        vertexList[*b].V[1],vertexList[*b].V[2],vertexList[*b].V[3],
                        vertexList[*d].V[1],vertexList[*d].V[2],vertexList[*d].V[3],
                        vertexList[*e].V[1],vertexList[*e].V[2],vertexList[*e].V[3]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = -icomp;
                mpz_clear(temp1);
                return;
        }

        //Term 16:	-Minor3(j,l,m,1,2,0)
        Deter3(temp1,   vertexList[*b].V[1],vertexList[*b].V[2],
                        vertexList[*d].V[1],vertexList[*d].V[2],
                        vertexList[*e].V[1],vertexList[*e].V[2]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = -icomp;
                mpz_clear(temp1);
                return;
        }

        //Term 17:	Minor3(j,l,m,1,3,0)
        Deter3(temp1,   vertexList[*b].V[1],vertexList[*b].V[3],
                        vertexList[*d].V[1],vertexList[*d].V[3],
                        vertexList[*e].V[1],vertexList[*e].V[3]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = icomp;
                mpz_clear(temp1);
                return;
        }

        //Term 18:	-Minor3(j,l,m,2,3,0)
        Deter3(temp1,   vertexList[*b].V[2],vertexList[*b].V[3],
                        vertexList[*d].V[2],vertexList[*d].V[3],
                        vertexList[*e].V[2],vertexList[*e].V[3]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = -icomp;
                mpz_clear(temp1);
                return;
        }

        //Term 19:	Minor3(i,l,m,1,2,0)
        Deter3(temp1,   vertexList[*a].V[1],vertexList[*a].V[2],
                        vertexList[*d].V[1],vertexList[*d].V[2],
                        vertexList[*e].V[1],vertexList[*e].V[2]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = icomp;
                mpz_clear(temp1);
                return;
        }

        //Term 20:	-Minor2(l,m,1,0)
        Deter2(temp1, vertexList[*d].V[1],vertexList[*e].V[1]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = -icomp;
                mpz_clear(temp1);
                return;
        }

        //Term 21:	Minor2(l,m,2,0)
        Deter2(temp1, vertexList[*d].V[2],vertexList[*e].V[2]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = icomp;
                mpz_clear(temp1);
                return;
        }

        //Term 22:	-Minor3(i,l,m,1,3,0)
        Deter3(temp1,   vertexList[*a].V[1],vertexList[*a].V[3],
                        vertexList[*d].V[1],vertexList[*d].V[3],
                        vertexList[*e].V[1],vertexList[*e].V[3]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = -icomp;
                mpz_clear(temp1);
                return;
        }

        //Term 23:	-Minor2(l,m,3,0)
        Deter2(temp1, vertexList[*d].V[3],vertexList[*e].V[3]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = -icomp;
                mpz_clear(temp1);
                return;
        }

        //Term 24:	Minor3(i,l,m,2,3,0)
        Deter3(temp1,   vertexList[*a].V[2],vertexList[*a].V[3],
                        vertexList[*d].V[2],vertexList[*d].V[3],
                        vertexList[*e].V[2],vertexList[*e].V[3]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = icomp;
                mpz_clear(temp1);
                return;
        }

        //Term 25:	Minor4(i,j,l,m,1,2,4,0)
        Deter4(temp1,   vertexList[*a].V[1],vertexList[*a].V[2],vertexList[*a].V[4],
                        vertexList[*b].V[1],vertexList[*b].V[2],vertexList[*b].V[4],
                        vertexList[*d].V[1],vertexList[*d].V[2],vertexList[*d].V[4],
                        vertexList[*e].V[1],vertexList[*e].V[2],vertexList[*e].V[4]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = icomp;
                mpz_clear(temp1);
                return;
        }

        //Term 26:	-Minor3(j,l,m,1,4,0)
        Deter3(temp1,   vertexList[*b].V[1],vertexList[*b].V[4],
                        vertexList[*d].V[1],vertexList[*d].V[4],
                        vertexList[*e].V[1],vertexList[*e].V[4]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = -icomp;
                mpz_clear(temp1);
                return;
        }

        //Term 27:	Minor3(j,l,m,2,4,0)
        Deter3(temp1,   vertexList[*b].V[2],vertexList[*b].V[4],
                        vertexList[*d].V[2],vertexList[*d].V[4],
                        vertexList[*e].V[2],vertexList[*e].V[4]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = icomp;
                mpz_clear(temp1);
                return;
        }

        //Term 28:	Minor3(i,l,m,1,4,0)
        Deter3(temp1,   vertexList[*a].V[1],vertexList[*a].V[4],
                        vertexList[*d].V[1],vertexList[*d].V[4],
                        vertexList[*e].V[1],vertexList[*e].V[4]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = icomp;
                mpz_clear(temp1);
                return;
        }

        //Term 29:	Minor2(l,m,4,0)
        Deter2(temp1, vertexList[*d].V[4],vertexList[*e].V[4]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = icomp;
                mpz_clear(temp1);
                return;
        }

        //Term 30:	-Minor3(i,l,m,2,4,0)
        Deter3(temp1,   vertexList[*a].V[2],vertexList[*a].V[4],
                        vertexList[*d].V[2],vertexList[*d].V[4],
                        vertexList[*e].V[2],vertexList[*e].V[4]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = -icomp;
                mpz_clear(temp1);
                return;
        }

        //Term 31:	-Minor4(i,j,l,m,1,3,4,0)
        Deter4(temp1,   vertexList[*a].V[1],vertexList[*a].V[3],vertexList[*a].V[4],
                        vertexList[*b].V[1],vertexList[*b].V[3],vertexList[*b].V[4],
                        vertexList[*d].V[1],vertexList[*d].V[3],vertexList[*d].V[4],
                        vertexList[*e].V[1],vertexList[*e].V[3],vertexList[*e].V[4]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = -icomp;
                mpz_clear(temp1);
                return;
        }

        //Term 32:	-Minor3(j,l,m,3,4,0)
        Deter3(temp1,   vertexList[*b].V[3],vertexList[*b].V[4],
                        vertexList[*d].V[3],vertexList[*d].V[4],
                        vertexList[*e].V[3],vertexList[*e].V[4]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = -icomp;
                mpz_clear(temp1);
                return;
        }

        //Term 33:	Minor3(i,l,m,3,4,0)
        Deter3(temp1,   vertexList[*a].V[3],vertexList[*a].V[4],
                        vertexList[*d].V[3],vertexList[*d].V[4],
                        vertexList[*e].V[3],vertexList[*e].V[4]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = icomp;
                mpz_clear(temp1);
                return;
        }

        //Term 34:	Minor4(i,j,l,m,2,3,4,0)
        Deter4(temp1,   vertexList[*a].V[2],vertexList[*a].V[3],vertexList[*a].V[4],
                        vertexList[*b].V[2],vertexList[*b].V[3],vertexList[*b].V[4],
                        vertexList[*d].V[2],vertexList[*d].V[3],vertexList[*d].V[4],
                        vertexList[*e].V[2],vertexList[*e].V[3],vertexList[*e].V[4]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = icomp;
                mpz_clear(temp1);
                return;
        }

        //Term 35:	Minor4(i,j,k,m,1,2,3,0)
        Deter4(temp1,   vertexList[*a].V[1],vertexList[*a].V[2],vertexList[*a].V[3],
                        vertexList[*b].V[1],vertexList[*b].V[2],vertexList[*b].V[3],
                        vertexList[*c].V[1],vertexList[*c].V[2],vertexList[*c].V[3],
                        vertexList[*e].V[1],vertexList[*e].V[2],vertexList[*e].V[3]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = icomp;
                mpz_clear(temp1);
                return;
        }

        //Term 36:	Minor3(j,k,m,1,2,0)
        Deter3(temp1,   vertexList[*b].V[1],vertexList[*b].V[2],
                        vertexList[*c].V[1],vertexList[*c].V[2],
                        vertexList[*e].V[1],vertexList[*e].V[2]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = icomp;
                mpz_clear(temp1);
                return;
        }

        //Term 37:	-Minor3(j,k,m,1,3,0)
        Deter3(temp1,   vertexList[*b].V[1],vertexList[*b].V[3],
                        vertexList[*c].V[1],vertexList[*c].V[3],
                        vertexList[*e].V[1],vertexList[*e].V[3]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = -icomp;
                mpz_clear(temp1);
                return;
        }

        //Term 38:	Minor3(j,k,m,2,3,0)
        Deter3(temp1,   vertexList[*b].V[2],vertexList[*b].V[3],
                        vertexList[*c].V[2],vertexList[*c].V[3],
                        vertexList[*e].V[2],vertexList[*e].V[3]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = icomp;
                mpz_clear(temp1);
                return;
        }

        //Term 39:	-Minor3(i,k,m,1,2,0)
        Deter3(temp1,   vertexList[*a].V[1],vertexList[*a].V[2],
                        vertexList[*c].V[1],vertexList[*c].V[2],
                        vertexList[*e].V[1],vertexList[*e].V[2]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = -icomp;
                mpz_clear(temp1);
                return;
        }

        //Term 40:	Minor2(k,m,1,0)
        Deter2(temp1, vertexList[*c].V[1],vertexList[*e].V[1]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = icomp;
                mpz_clear(temp1);
                return;
        }

        //Term 41:	-Minor2(k,m,2,0)
        Deter2(temp1, vertexList[*c].V[2],vertexList[*e].V[2]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = -icomp;
                mpz_clear(temp1);
                return;
        }

        //Term 42:	Minor3(i,k,m,1,3,0)
        Deter3(temp1,   vertexList[*a].V[1],vertexList[*a].V[3],
                        vertexList[*c].V[1],vertexList[*c].V[3],
                        vertexList[*e].V[1],vertexList[*e].V[3]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = icomp;
                mpz_clear(temp1);
                return;
        }

        //Term 43:	Minor2(k,m,3,0)
        Deter2(temp1, vertexList[*c].V[3],vertexList[*e].V[3]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = icomp;
                mpz_clear(temp1);
                return;
        }

        //Term 44:	-Minor3(i,k,m,2,3,0)
        Deter3(temp1,   vertexList[*a].V[2],vertexList[*a].V[3],
                        vertexList[*c].V[2],vertexList[*c].V[3],
                        vertexList[*e].V[2],vertexList[*e].V[3]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = -icomp;
                mpz_clear(temp1);
                return;
        }

        //Term 45:	Minor3(i,j,m,1,2,0)
        Deter3(temp1,   vertexList[*a].V[1],vertexList[*a].V[2],
                        vertexList[*b].V[1],vertexList[*b].V[2],
                        vertexList[*e].V[1],vertexList[*e].V[2]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = icomp;
                mpz_clear(temp1);
                return;
        }

        //Term 46:	-Minor2(j,m,1,0)
        Deter2(temp1, vertexList[*b].V[1],vertexList[*e].V[1]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = -icomp;
                mpz_clear(temp1);
                return;
        }

        //Term 47:	Minor2(j,m,2,0)
        Deter2(temp1, vertexList[*b].V[2],vertexList[*e].V[2]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = icomp;
                mpz_clear(temp1);
                return;
        }

        //Term 48:	Minor2(i,m,1,0)
        Deter2(temp1, vertexList[*a].V[1],vertexList[*e].V[1]);

        icomp = mpz_sgn(temp1);
        if (icomp != 0)
        {
                *res = icomp;
                mpz_clear(temp1);
                return;
        }

        //Term 49:	1
        *res = 1;
        mpz_clear(temp1);
        return;
}

/*!
    \fn Sos::onlyMinor4(std::vector<Vertex> & vertexList, int *a, int *b, int *c, int *d, int *res)
 */
void Sos::onlyMinor4(std::vector<Vertex> & vertexList, int *a, int *b, int *c, int *d, int *res)
{
        int icomp;

        mpz_t temp1;mpz_init(temp1);

        //Compute determinant
        Deter4(temp1,   vertexList[*a].V[1], vertexList[*a].V[2], vertexList[*a].V[3],
                        vertexList[*b].V[1], vertexList[*b].V[2], vertexList[*b].V[3],
                        vertexList[*c].V[1], vertexList[*c].V[2], vertexList[*c].V[3],
                        vertexList[*d].V[1], vertexList[*d].V[2], vertexList[*d].V[3]);

        icomp = mpz_sgn(temp1);
        *res = icomp;
        mpz_clear(temp1);
        return;
}
