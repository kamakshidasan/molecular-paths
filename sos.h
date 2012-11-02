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

#ifndef SOS_H
#define SOS_H

#include <gmp.h>
#include <vector>
#include <vertex.h>

class Sos
{
    public:
        Sos();
        ~Sos();

        void Deter2(mpz_t deter,mpz_t b11,mpz_t b21);
        void Deter3(mpz_t deter, mpz_t b11, mpz_t b12, mpz_t b21, mpz_t b22, mpz_t b31, mpz_t b32);
        void Deter4(mpz_t deter, mpz_t b11, mpz_t b12, mpz_t b13,mpz_t b21, mpz_t b22, mpz_t b23,mpz_t b31, mpz_t b32, mpz_t b33,mpz_t b41, mpz_t b42, mpz_t b43);
        void Deter5(mpz_t deter, mpz_t b11, mpz_t b12, mpz_t b13, mpz_t b14,mpz_t b21, mpz_t b22, mpz_t b23, mpz_t b24,mpz_t b31, mpz_t b32, mpz_t b33, mpz_t b34,mpz_t b41, mpz_t b42, mpz_t b43, mpz_t b44,mpz_t b51, mpz_t b52, mpz_t b53, mpz_t b54);
        void Minor2(std::vector<Vertex> & vertexList, int *a, int *b, int ia, int *res);
        void Minor3(std::vector<Vertex> & vertexList,  int *a,  int *b,  int *c,  int *i1,  int *i2,  int *res);
        void Minor4(std::vector<Vertex> & vertexList,  int *a,  int *b,  int *c,  int *d,  int *res);
        void Minor5(std::vector<Vertex> & vertexList, int *a, int *b, int *c, int *d, int *e, int *res);
        void onlyMinor4(std::vector<Vertex> & vertexList, int *a, int *b, int *c, int *d, int *res);

};

#endif // SOS_H
