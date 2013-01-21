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

#include <stdio.h>
#include <algorithm>
#include <skinsurface.h>
#include <lightmaterial.h>

static void Cross(double u[3], double v[3], double n[3])
{
	n[0] = u[1]*v[2] - u[2]*v[1];
	n[1] = u[2]*v[0] - u[0]*v[2];
	n[2] = u[0]*v[1] - u[1]*v[0];
}

static void Normalize(double* v)
{
	double l;

	l = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
	v[0] /= l;
	v[1] /= l;
	v[2] /= l;
}

/******************************************************************************************************/
SkinSurface::SkinSurface(double center[],double scale,int alpha_poc)
{
	for(int i = 0; i < 3; i++)
	{
		scenter[i] = center[i];
	}
	sscale = scale;
	vertList.clear();
	triList.clear();
	verTrimap.clear();
        vertexNormals.clear();
        AlphaOrPocket = alpha_poc;
};

/*!
    \fn SkinSurface::InitMaterial(int i)
 */
void SkinSurface::InitMaterial(int i)
{
        if (i == 0)
        {
                glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, LightMaterial::MatAmb[i]);
                glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, LightMaterial::MatDiff[i]);
                glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, LightMaterial::MatSpec[i]);
                glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, LightMaterial::MatShin[i]);
                glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, LightMaterial::MatEmission);
        }
        else
        {
                //int m = i % 9;
                int m = i % 15;
                glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, LightMaterial::MatAmb[m]);
                glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, LightMaterial::MatDiff[m]);
                glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, LightMaterial::MatSpec[m]);
                glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, LightMaterial::MatShin[m]);
                glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, LightMaterial::MatEmission);
        }
}

void SkinSurface::Read(const char* vfname,int count)
{
	FILE* fp = fopen(vfname, "r");
	assert(fp);
	char waste[4];
	int n1 = 0,n2 = 0, i = 0,j;
	unsigned int t[3] = {0,0,0};
	double v[3] = {0,0,0};
	double b = 0;
	int vsize;
	int tsize;
	std::vector<int> temp;
	temp.clear();

	if(count == 0)
	{
		triList.clear();
		vertList.clear();
		verTrimap.clear();
		Triangle tri(t);
		triList.push_back(tri);
                vertexNormals.clear();
		vsize = 0;
		tsize = 0;
	}
	else
	{
		vsize = vertList.size();
		tsize = triList.size()-1;
	}

	fscanf(fp, "%s",waste);
	fscanf(fp, "%d %d %d",&n1,&n2,&j);

	for(i = 0; i < n1; i++)
	{
		fscanf(fp, "%lf %lf %lf", v, v+1, v+2);
		Vertex vert(v[0],v[1],v[2],b,-1,0.0);
		vertList.push_back(vert);
	}

        // translate around center then scale
        /*if(AlphaOrPocket == 0)
        {
            for (i = 0; i < vertList.size(); i++)
            {
                    vertList[i].NormCoordinates[1] -= scenter[0];
                    vertList[i].NormCoordinates[2] -= scenter[1];
                    vertList[i].NormCoordinates[3] -= scenter[2];
                    vertList[i].NormCoordinates[1] *= sscale;
                    vertList[i].NormCoordinates[2] *= sscale;
                    vertList[i].NormCoordinates[3] *= sscale;
            }
        }*/

	for(i = 0; i < n1; i++)
	{
		verTrimap.push_back(temp);
	}

	for(i = 1; i <= n2; i++)
	{
		fscanf(fp, "%d %d %d %d",&j, t, t+1, t+2);
		t[0] += vsize;
		t[1] += vsize;
		t[2] += vsize;
		verTrimap[t[0]].push_back(tsize+i);
		verTrimap[t[1]].push_back(tsize+i);
		verTrimap[t[2]].push_back(tsize+i);
		Triangle tri(t);
                tri.PocIndex = count;
		triList.push_back(tri);
	}

	fclose(fp);
}

void SkinSurface::DrawSolid()
{
    int a,b,c;
    uint i = 0;

    for ( i = 1; i < triList.size(); i++)
    {
        a = triList[i].Corners[1];
        b = triList[i].Corners[2];
        c = triList[i].Corners[3];

        glBegin(GL_TRIANGLES);

            glNormal3dv(vertexNormals[a].normals);
            glVertex3d(vertList[a].NormCoordinates[1],vertList[a].NormCoordinates[2],vertList[a].NormCoordinates[3]);

            glNormal3dv(vertexNormals[b].normals);
            glVertex3d(vertList[b].NormCoordinates[1],vertList[b].NormCoordinates[2],vertList[b].NormCoordinates[3]);

            glNormal3dv(vertexNormals[c].normals);
            glVertex3d(vertList[c].NormCoordinates[1],vertList[c].NormCoordinates[2],vertList[c].NormCoordinates[3]);

        glEnd();
    }
}

#include "scalarfield.h"
void SkinSurface::DrawWithField(ScalarField* field)
{
    int a,b,c;
    uint i = 0;

    Color low(1,0,0);
    Color high(0,0,1);

    for ( i = 1; i < triList.size(); i++)
    {
        a = triList[i].Corners[1];
        b = triList[i].Corners[2];
        c = triList[i].Corners[3];

        glBegin(GL_TRIANGLES);

            float pt[3];
            pt[0] = (float) vertList[a].NormCoordinates[1];
            pt[1] = (float) vertList[a].NormCoordinates[2];
            pt[2] = (float) vertList[a].NormCoordinates[3];
            float val = field->getValueBiLinear(pt);
            Color col = field->getColor2(low, high, val);
//            col = (val<0)? low : high;
            glColor3f(col.r, col.g, col.b);
            glNormal3dv(vertexNormals[a].normals);
            glVertex3d(vertList[a].NormCoordinates[1],vertList[a].NormCoordinates[2],vertList[a].NormCoordinates[3]);

            pt[0] = (float) vertList[b].NormCoordinates[1];
            pt[1] = (float) vertList[b].NormCoordinates[2];
            pt[2] = (float) vertList[b].NormCoordinates[3];
            val = field->getValueBiLinear(pt);
            col = field->getColor2(low, high, val);
//            col = (val<0)? low : high;
            glColor3f(col.r, col.g, col.b);
            glNormal3dv(vertexNormals[b].normals);
            glVertex3d(vertList[b].NormCoordinates[1],vertList[b].NormCoordinates[2],vertList[b].NormCoordinates[3]);

            pt[0] = (float) vertList[c].NormCoordinates[1];
            pt[1] = (float) vertList[c].NormCoordinates[2];
            pt[2] = (float) vertList[c].NormCoordinates[3];
            val = field->getValueBiLinear(pt);
            col = field->getColor2(low, high, val);
//            col = (val<0)? low : high;
            glColor3f(col.r, col.g, col.b);
            glNormal3dv(vertexNormals[c].normals);
            glVertex3d(vertList[c].NormCoordinates[1],vertList[c].NormCoordinates[2],vertList[c].NormCoordinates[3]);

        glEnd();
    }
}

void SkinSurface::Draw(bool smoothShading,bool skinWireFrame)
{
	int a,b,c;
	uint i = 0;

	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, LightMaterial::MatAmb[1]);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, LightMaterial::MatDiff[1]);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, LightMaterial::MatSpec[1]);
    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, LightMaterial::MatShin[1]);
	glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, LightMaterial::MatEmission);

        if(skinWireFrame)
        {
                glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
                glLineWidth(0.5);
        }
        else
        {
                glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
        }

	for ( i = 1; i < triList.size(); i++)
	{
		a = triList[i].Corners[1];
		b = triList[i].Corners[2];
		c = triList[i].Corners[3];

                if(AlphaOrPocket == 0)
                {
                    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, LightMaterial::MatAmb[3]);
                    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, LightMaterial::MatDiff[3]);
                    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, LightMaterial::MatSpec[3]);
                    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, LightMaterial::MatShin[3]);
                    glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, LightMaterial::MatEmission);
                }
                else
                {
                    InitMaterial(triList[i].PocIndex);
                }

		glBegin(GL_TRIANGLES);

                if(!smoothShading)
		{
			glNormal3d(triList[i].Normal->X,triList[i].Normal->Y,triList[i].Normal->Z);
			glVertex3d(vertList[a].NormCoordinates[1],vertList[a].NormCoordinates[2],vertList[a].NormCoordinates[3]);
			glVertex3d(vertList[b].NormCoordinates[1],vertList[b].NormCoordinates[2],vertList[b].NormCoordinates[3]);
			glVertex3d(vertList[c].NormCoordinates[1],vertList[c].NormCoordinates[2],vertList[c].NormCoordinates[3]);
		}
		else
		{
			glNormal3dv(vertexNormals[a].normals);
			glVertex3d(vertList[a].NormCoordinates[1],vertList[a].NormCoordinates[2],vertList[a].NormCoordinates[3]);

			glNormal3dv(vertexNormals[b].normals);
			glVertex3d(vertList[b].NormCoordinates[1],vertList[b].NormCoordinates[2],vertList[b].NormCoordinates[3]);

			glNormal3dv(vertexNormals[c].normals);
			glVertex3d(vertList[c].NormCoordinates[1],vertList[c].NormCoordinates[2],vertList[c].NormCoordinates[3]);
		}

		glEnd();
    }
}

void SkinSurface::Draw(bool smoothShading, bool skinWireFrame, int pocIndex)
{
    int a,b,c;
    uint i = 0;

    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, LightMaterial::MatAmb[1]);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, LightMaterial::MatDiff[1]);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, LightMaterial::MatSpec[1]);
    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, LightMaterial::MatShin[1]);
    glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, LightMaterial::MatEmission);

    if(skinWireFrame)
    {
            glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
    }
    else
    {
            glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
    }

    for ( i = 1; i < triList.size(); i++)
    {
            a = triList[i].Corners[1];
            b = triList[i].Corners[2];
            c = triList[i].Corners[3];

            if(pocIndex!=triList[i].PocIndex)
            {
                continue;
            }

            InitMaterial(triList[i].PocIndex);

            glBegin(GL_TRIANGLES);

            if(!smoothShading)
            {
                    glNormal3d(triList[i].Normal->X,triList[i].Normal->Y,triList[i].Normal->Z);
                    glVertex3d(vertList[a].NormCoordinates[1],vertList[a].NormCoordinates[2],vertList[a].NormCoordinates[3]);
                    glVertex3d(vertList[b].NormCoordinates[1],vertList[b].NormCoordinates[2],vertList[b].NormCoordinates[3]);
                    glVertex3d(vertList[c].NormCoordinates[1],vertList[c].NormCoordinates[2],vertList[c].NormCoordinates[3]);
            }
            else
            {
                    glNormal3dv(vertexNormals[a].normals);
                    glVertex3d(vertList[a].NormCoordinates[1],vertList[a].NormCoordinates[2],vertList[a].NormCoordinates[3]);

                    glNormal3dv(vertexNormals[b].normals);
                    glVertex3d(vertList[b].NormCoordinates[1],vertList[b].NormCoordinates[2],vertList[b].NormCoordinates[3]);

                    glNormal3dv(vertexNormals[c].normals);
                    glVertex3d(vertList[c].NormCoordinates[1],vertList[c].NormCoordinates[2],vertList[c].NormCoordinates[3]);
            }

            glEnd();
    }

}

void SkinSurface::Process()
{
	uint i,j;

	double v1[3],v2[3],v[3];

	//face normals
	for (i = 0; i< triList.size(); i++)
	{
		//a,b,c
		// v1 = b - a
		v1[0] = vertList[triList[i].Corners[2]].Coordinates[1] - vertList[triList[i].Corners[1]].Coordinates[1];
		v1[1] = vertList[triList[i].Corners[2]].Coordinates[2] - vertList[triList[i].Corners[1]].Coordinates[2];
		v1[2] = vertList[triList[i].Corners[2]].Coordinates[3] - vertList[triList[i].Corners[1]].Coordinates[3];
		//v2 = c - a
		v2[0] = vertList[triList[i].Corners[3]].Coordinates[1] - vertList[triList[i].Corners[1]].Coordinates[1];
		v2[1] = vertList[triList[i].Corners[3]].Coordinates[2] - vertList[triList[i].Corners[1]].Coordinates[2];
		v2[2] = vertList[triList[i].Corners[3]].Coordinates[3] - vertList[triList[i].Corners[1]].Coordinates[3];

		// v = v1 X v2
		Cross(v1, v2, v);
		Normalize(v);

		triList[i].Normal->X = v[0];
		triList[i].Normal->Y = v[1];
		triList[i].Normal->Z = v[2];
	}

	//vertex normals
	for(i = 0;i<verTrimap.size();i++)
	{
		v[0] = v[1] = v[2] = 0.0;
		//add all facet normals of which the vertex belongs to
		for(j =0;j< verTrimap[i].size();j++)
		{
			v[0] += triList[verTrimap[i][j]].Normal->X;
			v[1] += triList[verTrimap[i][j]].Normal->Y;
			v[2] += triList[verTrimap[i][j]].Normal->Z;
		}
		//take average
		v[0] /= verTrimap[i].size();
		v[1] /= verTrimap[i].size();
		v[2] /= verTrimap[i].size();
		vertexNormals.push_back(VertexNormals(v));
	}

        /*FILE *fp = fopen("vnorms","w");
	for(i = 0;i<vertexNormals.size();i++)
	{
		fprintf(fp,"%f %f %f\n",vertexNormals[i].normals[0],vertexNormals[i].normals[1],vertexNormals[i].normals[2]);
	}
        fclose(fp);*/
}
