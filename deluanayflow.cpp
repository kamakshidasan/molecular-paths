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

#include "deluanayflow.h"
#include <iostream>
using namespace std;

DeluanayFlow::DeluanayFlow()
{
        flow = new Flow();
}

DeluanayFlow::~DeluanayFlow()
{
}

/*!
    \fn DeluanayFlow::CalculateDF(DeluanayComplex *delcx, std::vector<Vertex> & vertexList)
 */
void DeluanayFlow::CalculateDF(DeluanayComplex *delcx, std::vector<Vertex> & vertexList)
{
        int i, j, k, l, m, n, p;
        uint idx;
        int iflow[5];
        /*
        cout << "!Calculate DF! " << delcx->DeluanayTet.size() << endl;

        for (idx = 1; idx < delcx->DeluanayTet.size(); idx++)
        {
            Tetrahedron tet = delcx->DeluanayTet[idx];
            for (m = 1; m <= 4; m++)
            {
                cout << tet.Corners[m] << ' ';
            }
            cout << "  ---   ";
            for (m = 1; m <= 4; m++)
            {
                cout << tet.Neighbours[m] << '('<< tet.Nindex[m] << ") ";
            }
            cout << endl;
        }
        */

        for (idx = 1; idx < delcx->DeluanayTet.size(); idx++)
        {
                if (delcx->DeluanayTet[idx].Status == 0) continue;

                i = delcx->DeluanayTet[idx].Corners[1];
                j = delcx->DeluanayTet[idx].Corners[2];
                k = delcx->DeluanayTet[idx].Corners[3];
                l = delcx->DeluanayTet[idx].Corners[4];

                for (m = 1; m <= 4; m++)
                {
                        iflow[m] = 0;
                        n = delcx->DeluanayTet[idx].Neighbours[m];
                        p = delcx->DeluanayTet[idx].Nindex[m];
                        if ((n != 0) && (n < idx))
                        {
                                if (delcx->DeluanayTet[n].TetFlow[p] == 1)
                                        iflow[m] = -1;
                        }
                }

                flow->CalculateFlow(vertexList, i, j, k, l, iflow);

                for (m = 1; m <= 4; m++)
                {
                        if (iflow[m] == 1)
                        {
                                delcx->DeluanayTet[idx].TetFlow[m] = 1;
                        }
                        else
                        {
                                delcx->DeluanayTet[idx].TetFlow[m] = 0;
                        }
                }
        }

        /*FILE *fp = fopen("flowfile","w");
        for (uint idx = 1; idx < delcx->DeluanayTet.size(); idx++)
        {
                fprintf(fp,"%d %d %d %d %d %d\n",idx,delcx->DeluanayTet[idx].Status,delcx->DeluanayTet[idx].TetFlow[1], delcx->DeluanayTet[idx].TetFlow[2],delcx->DeluanayTet[idx].TetFlow[3],delcx->DeluanayTet[idx].TetFlow[4]);
        }
        fclose(fp);*/
}

