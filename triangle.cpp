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

#include "triangle.h"

Triangle::Triangle()
{
        this->ReverseLink1 = 0;
        this->ReverseLink2 = 0;
        this->nLink = 0;
        this->isValid = false;

        for (int i = 0; i < 4; i++)
        {
                this->Corners[i] = -1;
                this->TrigLink[i] = -1;
        }
        this->Hull = -1;
        this->AlphaStatus = -1;
        this->IsAttached = false;
        this->Entry = -1;
        this->Rho = 0;
        this->Mu1 = 0;
        this->Mu2 = 0;
        this->Size = -1.0;
        this->ufKey = -1;
        this->PocIndex = -1;
        this->IsMouth = false;
        this->Repeats = -1;
        this->MainRepeats = -1;
        this->Persistence = -1;
        this->AlphaPersistence = 0.0;
        this->Normal = new Vector3(0.0,0.0,0.0);
}

Triangle::Triangle(unsigned int n[3])
{
        this->ReverseLink1 = 0;
        this->ReverseLink2 = 0;
        this->nLink = 0;
        this->isValid = false;

        this->Corners[1] = n[0];this->Corners[2] = n[1];this->Corners[3] = n[2];
        for (int i = 0; i < 4; i++)
        {
                this->TrigLink[i] = -1;
        }
        this->Hull = -1;
        this->AlphaStatus = -1;
        this->IsAttached = false;
        this->Entry = -1;
        this->Rho = 0;
        this->Mu1 = 0;
        this->Mu2 = 0;
        this->Size = -1.0;
        this->ufKey = -1;
        this->PocIndex = -1;
        this->IsMouth = false;
        this->Repeats = -1;
        this->MainRepeats = -1;
        this->Persistence = -1;
        this->AlphaPersistence = 0.0;
        this->Normal = new Vector3(0.0,0.0,0.0);
}

Triangle::~Triangle()
{
}

Vector3 Triangle::getCentroid(std::vector<Vertex> & vertList){
    Vector3 ret;
    Vector3 v1 = vertList[Corners[1]].getCoordVector();
    Vector3 v2 = vertList[Corners[2]].getCoordVector();
    Vector3 v3 = vertList[Corners[3]].getCoordVector();
    ret = v1;
    Vector3::Sum(&ret, &ret, &v2);
    Vector3::Sum(&ret, &ret, &v3);
    Vector3::Scale(&ret, &ret, 1.0/3.0);
    return ret;
}

