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

#ifndef SIZE_H
#define SIZE_H

#include <cassert>
#include <gmp.h>
#include <vertex.h>
#include <edge.h>
#include <triangle.h>
#include <triangle.h>
#include <deluanaycomplex.h>

class Size
{
    public:
        Size();
        ~Size();

        int AlfRatioCompare(mpz_t* a, mpz_t* b,mpz_t* c, mpz_t* d);
        void VertSize(std::vector<Vertex> & vertexList, int i, mpz_t *p, mpz_t *q);
        int CheckVertex(std::vector<Vertex> & vertexList, int a, int b);
        void CheckEdge(mpz_t *p_mp,mpz_t *q_mp,mpz_t *r_mp,bool *isAttached);
        void EdgeSize(std::vector<Vertex> & vertexList, int a, int b, mpz_t *p, mpz_t *q);
        void TrigSize(std::vector<Vertex> & vertexList, int a, int b, int c, mpz_t *p, mpz_t *q);
        void CheckTriangle(mpz_t *p_mp,mpz_t *q_mp,mpz_t *r_mp,mpz_t *s_mp, bool *isAttached);
        void TetraSize(DeluanayComplex *delcx, std::vector<Vertex> & vertexList, int a, int b, int c, int d, mpz_t *p, mpz_t *q, int index);

    private:
        int ISort3(int *a, int *b, int *c);
        int ISort4(int *a, int *b, int *c, int *d);
};

#endif // SIZE_H
