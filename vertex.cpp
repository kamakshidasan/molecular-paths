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
#include "vertex.h"
#include <cmath>

Vertex::Vertex(double x, double y, double z, double radius,int index,double scale)
{
        mpz_t temp1,temp2;

        this->Index = index;
        this->Coordinates[1] = this->NormCoordinates[1] = x;
        this->Coordinates[2] = this->NormCoordinates[2] = y;
        this->Coordinates[3] = this->NormCoordinates[3] = z;
        this->Radius = this->BackRadius = radius;
        this->Weight = this->BackWeight = pow(x,2) + pow(y,2) + pow(z,2) - pow(radius,2);
        for(int i=1;i<4;i++)
        {
                mpz_init(this->V[i]);
                mpz_set_d(this->V[i],this->Coordinates[i]*scale);
        }
        mpz_init(temp1);mpz_init(temp2);
        mpz_init(this->V[4]);

        mpz_set_d(temp1,this->Radius*(scale)); mpz_mul(temp1,temp1,temp1);
        mpz_mul(temp2,this->V[3],this->V[3]), mpz_sub(temp1,temp2,temp1);
        mpz_mul(temp2,this->V[2],this->V[2]), mpz_add(temp1,temp2,temp1);
        mpz_mul(temp2,this->V[1],this->V[1]), mpz_add(this->V[4],temp2,temp1);

        mpz_clear(temp1);mpz_clear(temp2);
        this->Redinfo = 0;
        this->ranValue = 0;
        this->Hull = -1;
        this->AlphaStatus = -1;
        this->Coef = 1.0;
        this->Rho = 0;
        this->Mu1 = 0;
        this->Mu2 = 0;
        this->Repeats = -1;
        this->ufKey = -1;
        this->valid = true;
        selected = 0;
}

Vertex::~Vertex()
{
}

/*!
    \fn Vertex::ranSortFunc(const Vertex& lhs,const Vertex& rhs)
 */
bool Vertex::ranSortFunc(const Vertex& lhs,const Vertex& rhs)
{
        return(lhs.ranValue < rhs.ranValue);
}


/*!
    \fn Vertex::operator < (const Vertex& rhs ) const
 */
bool Vertex::operator < (const Vertex& rhs ) const
{
    return(Coordinates[0] < rhs.Coordinates[0]);
}
