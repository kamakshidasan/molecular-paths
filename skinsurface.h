/***************************************************************************
 *   Copyright (C) 2009 by raghavendra,,,                                  *
 *   raghavendra@incognito                                                 *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General License as published by         *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General License for more details.                                 *
 *                                                                         *
 *   You should have received a copy of the GNU General License            *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef SKINSURFACE_HPP
#define SKINSURFACE_HPP

#include <vector>
#include <map>
#include <assert.h>
#include <GL/glew.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include <vertex.h>
#include <triangle.h>
#include "scalarfield.h"

class VertexNormals
{
    public:
        double normals[3];
        VertexNormals(double n[3])
        {
                memcpy(normals,n,sizeof(normals));
        }
};

class SkinSurface
{
    private:
        double scenter[3],sscale;
        int AlphaOrPocket;
        std::vector<Vertex> vertList;
        std::vector<Triangle> triList;
        std::vector< std::vector<int> > verTrimap;
        std::vector<VertexNormals> vertexNormals;

    public:
        SkinSurface(double center[],double scale,int alpha_poc);
        void InitMaterial(int i);
        void Read(const char* vfname,int count);
        void ReadWrite(const char* vfname, const char* outFileName, ScalarField* field, int count);
        void Process();
        void Draw(bool smoothShading,bool skinWireFrame);
        void Draw(bool smoothShading,bool skinWireFrame,int pocIndex);
        void DrawSolid();
        void DrawWithField(ScalarField* field);
};

#endif
