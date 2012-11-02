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

#ifndef TLISTNODE_H
#define TLISTNODE_H

#include <gmp.h>

class TlistNode
{
    public:
        unsigned int ftype;
        unsigned int rtype;
        int r;
        int ix;
        int si;
        mpz_t a;
        mpz_t b;

        TlistNode(mpz_t *p,mpz_t *q,int rank, int index, int selin, unsigned int ft, unsigned int rt);
        ~TlistNode();

        bool operator < (const TlistNode& rhs ) const;
        TlistNode operator = (TlistNode tl);
};

#endif // TLISTNODE_H
