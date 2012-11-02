/***************************************************************************
 *   Copyright (C) 2010 by raghavendra, , ,                                   *
 *   raghavendra@incognito                                                 *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License,  or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,        *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not,  write to the                         *
 *   Free Software Foundation,  Inc.,                                        *
 *   59 Temple Place - Suite 330,  Boston,  MA  02111-1307,  USA.             *
 ***************************************************************************/
#include "vector3.h"

Vector3::Vector3(double x, double y, double z)
{
        this->X = x;
        this->Y = y;
        this->Z = z;
}

Vector3::Vector3()
{
        this->X = 0.0;
        this->Y = 0.0;
        this->Z = 0.0;
}

Vector3::~Vector3()
{
}

/*!
    \fn Vector3::DotProduct(Vector3* v1, Vector3 *v2)
 */
void Vector3::DotProduct(Vector3* v1, Vector3 *v2, double *res)
{
        *res = v1->X * v2->X + v1->Y * v2->Y + v1->Z * v2->Z;
}

/*!
    \fn Vector3::CrossProduct(Vector3 *result, Vector3 *u, Vector3 *v)
 */
void Vector3::CrossProduct(Vector3 *result, Vector3 *u, Vector3 *v)
{
        result->X = u->Y * v->Z - u->Z * v->Y;
        result->Y = u->Z * v->X - u->X * v->Z;
        result->Z = u->X * v->Y - u->Y * v->X;
}

/*!
    \fn Vector3::DiffVector(Vector3 *result, Vector3 *u, Vector3 *v)
 */
void Vector3::DiffVector(Vector3 *result, Vector3 *u, Vector3 *v)
{
        result->X = u->X - v->X;
        result->Y = u->Y - v->Y;
        result->Z = u->Z - v->Z;
}
/*!
    \fn Vector3::Normalize(Vector3 *result, Vector3 *v)
 */
void Vector3::Normalize(Vector3 *result, Vector3 *v)
{
        double mod;
        DotProduct(v, v, &mod);
        mod = sqrt(mod);
        result->X = v->X/mod;
        result->Y = v->Y/mod;
        result->Z = v->Z/mod;
}

/*!
    \fn Vector3::Normalize()
 */
void Vector3::Normalize()
{
        double mod;
        DotProduct(this, this, &mod);
        mod = sqrt(mod);
        this->X/=mod;
        this->Y/=mod;
        this->Z/=mod;
}

void Vector3::Sum (Vector3 *result, Vector3 *u, Vector3 *v)
{
    result->X = u->X + v->X;
    result->Y = u->Y + v->Y;
    result->Z = u->Z + v->Z;
}

void Vector3::Scale(Vector3* result, Vector3 *in, double factor)
{
    result->X = in->X * factor;
    result->Y = in->Y * factor;
    result->Z = in->Z * factor;
}

void Vector3::Distance(double *dist, Vector3 * u, Vector3 *v)
{
    double sqrdist = 0.0;
    Vector3 D;
    DiffVector (&D, u, v);
    DotProduct (&D, &D, &sqrdist);
    if(sqrdist < 0.0) sqrdist = 0.0;

    *dist = sqrt(sqrdist);
}
