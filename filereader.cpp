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

#include "filereader.h"

static double myMax(double a, double b)
{
        if (b > a) return b;
        return a;
}

/*static double myAbs(double f)
{
        if (f < 0) return -f;
        return f;
}*/

FileReader::FileReader()
{
}

FileReader::~FileReader()
{
}

/*!
    \fn FileReader::ReadVertices(const char* filename,std::vector<Vertex> & vertexList,,double center[],double *Scale)
 */
void FileReader::ReadVertices(const char* filename,std::vector<Vertex> & vertexList,double center[],
                              double *Scale, double min[], double max[])
{
        unsigned int n = 0;
        unsigned int i = 0;
        char *val;
        double r = 0,x = 0,y = 0,z = 0;
        double scale = 100000000.0;
        char line[100];

        FILE* fp = fopen(filename, "r");
        assert(fp);

        val=fgets(line, 100, fp);
        sscanf(line,"%d",&n);
        val=fgets(line, 100, fp);
        Vertex vert(x,y,z,r,0,scale);
        vertexList.push_back(vert);

        FILE *fp1 = fopen("radius.txt","w");

        for(i = 1; i <= n; i++)
        {
                val=fgets(line, 100, fp);
                sscanf(line,"%lf %lf %lf %lf",&x, &y, &z, &r);
                fprintf(fp1,"%d %lf\t",i,r);
                r += 1.4;
                fprintf(fp1,"%lf\n",r);
                Vertex vert(x,y,z,r,i,scale);
                vertexList.push_back(vert);
        }
        fclose(fp1);
        fclose(fp);

        scale = 0;
        double maxx = vertexList[1].Coordinates[1], maxy = vertexList[1].Coordinates[2], maxz = vertexList[1].Coordinates[3];
        double minx = vertexList[1].Coordinates[1], miny = vertexList[1].Coordinates[2], minz = vertexList[1].Coordinates[3];
        double cx, cy, cz, w, h, e;

        for (i = 1; i < vertexList.size(); i++)
        {
                if (vertexList[i].Coordinates[1] > maxx)
                {
                        maxx = vertexList[i].Coordinates[1];
                }
                if (vertexList[i].Coordinates[1] < minx)
                {
                        minx = vertexList[i].Coordinates[1];
                }
                if (vertexList[i].Coordinates[2] > maxy)
                {
                        maxy = vertexList[i].Coordinates[2];
                }
                if (vertexList[i].Coordinates[2] < miny)
                {
                        miny = vertexList[i].Coordinates[2];
                }
                if (vertexList[i].Coordinates[3] > maxz)
                {
                        maxz = vertexList[i].Coordinates[3];
                }
                if (vertexList[i].Coordinates[3] < minz)
                {
                        minz = vertexList[i].Coordinates[3];
                }
        }

        // calculate model width, height, and depth
        w = (maxx) - (minx);
        h = (maxy) - (miny);
        e = (maxz) - (minz);

        // calculate center of the model
        cx = (maxx + minx) / 2.0;
        cy = (maxy + miny) / 2.0;
        cz = (maxz + minz) / 2.0;

        // calculate unitizing scale factor
        //scale = 2.0 / myMax(myMax(w, h), e);
        scale = myMax(myMax(w, h), e);

        center[0] = cx;center[1] = cy;center[2] = cz;
        *Scale = scale;

        min[0] = minx;
        min[1] = minx;
        min[2] = minx;
        max[0] = maxx;
        max[1] = maxy;
        max[2] = maxz;

        // translate around center then scale
        /*for (i = 1; i < vertexList.size(); i++)
        {
                vertexList[i].NormCoordinates[1] -= cx;
                vertexList[i].NormCoordinates[2] -= cy;
                vertexList[i].NormCoordinates[3] -= cz;
                vertexList[i].NormCoordinates[1] *= scale;
                vertexList[i].NormCoordinates[2] *= scale;
                vertexList[i].NormCoordinates[3] *= scale;
        }*/
}

