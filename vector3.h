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
#ifndef VECTOR3_H
#define VECTOR3_H

#include <cmath>

class Vector3
{

    public:
        double X;
        double Y;
        double Z;

        Vector3(double x, double y, double z);
        Vector3();
        ~Vector3();

        float operator[](int i) const
        {
            switch(i)
            {
            case 0:return X;
            case 1:return Y;
            case 2:return Z;
            default:return -1;
            }
        };

        static void DotProduct(Vector3* v1, Vector3 *v2, double *res);
        static void CrossProduct(Vector3 *result, Vector3 *u, Vector3 *v) ;
        static void DiffVector(Vector3 *result, Vector3 *u, Vector3 *v);
        static void Normalize(Vector3 *result, Vector3 *v);
        static void Sum (Vector3 *result, Vector3 *u, Vector3 *v);
        static void Scale(Vector3* result, Vector3 *in, double factor);
        static void Distance(double *dist, Vector3 * u, Vector3 *v);
        void Normalize();
        void Negate(){
            this->X = -this->X;
            this->Y = -this->Y;
            this->Z = -this->Z;
        }
};

#endif // VECTOR3_H
