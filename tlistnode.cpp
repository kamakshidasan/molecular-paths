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

#include "tlistnode.h"
TlistNode::TlistNode(mpz_t *p,mpz_t *q,int rank, int index, int selin, unsigned int ft, unsigned int rt)
{
        mpz_init(this->a);mpz_init(this->b);
        mpz_set(this->a,*p);
        mpz_set(this->b,*q);
        this->r = rank;
        this->ix = index;
        this->si = selin;
        this->ftype = ft;
        this->rtype = rt;
}

TlistNode::~TlistNode()
{
        //mpz_clear(a);
        //mpz_clear(b);
}

/*!
    \fn TlistNode::operator < (const TlistNode& rhs )
 */
bool TlistNode::operator < (const TlistNode& rhs ) const
{
        mpz_t p;mpz_init(p);
        mpz_t q;mpz_init(q);
        mpz_t r;mpz_init(r);
        int temp;

        mpz_mul(p,this->a,rhs.b);
        mpz_mul(q,this->b,rhs.a);
        mpz_sub(r,p,q);

        temp = mpz_sgn(r);

        if(temp<0)
        {
                mpz_clear(p);
                mpz_clear(q);
                mpz_clear(r);
                return true;
        }
        else
        {

                mpz_clear(p);
                mpz_clear(q);
                mpz_clear(r);
                return false;
        }
}

 TlistNode TlistNode::operator = (TlistNode tl)
{
        ftype = tl.ftype;
        rtype = tl.rtype;
        r = tl.r;
        ix = tl.ix;
        si = tl.si;
        mpz_init(a);mpz_init(b);
        mpz_set(a,tl.a);
        mpz_set(b,tl.b);

        return *this;
}
