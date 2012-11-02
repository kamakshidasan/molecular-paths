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

#ifndef FLOW_H
#define FLOW_H

#include <vector>
#include <gmp.h>
#include <vertex.h>

class Flow
{
    private:
        void CheckFacet(mpz_t S[], mpz_t Det1, mpz_t Det2, mpz_t Det3, mpz_t Deter, mpz_t De3, int *testa);

    public:
        Flow();
        ~Flow();

        void CalculateFlow(std::vector<Vertex> & vertexList, int i, int j, int k, int l, int testa[]);
};

#endif // FLOW_H
