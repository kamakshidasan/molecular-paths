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

#include "depth.h"

static int myMax(int a, int b)
{
        if (b > a) return b;
        return a;
}

Depth::Depth()
{
}

Depth::~Depth()
{
}

/*!
    \fn Depth::CalculateDepth(DeluanayComplex *delcx, std::vector <int> &sortedTet)
 */
void Depth::CalculateDepth(DeluanayComplex *delcx, std::vector <int> &sortedTet)
{
        int tetindex, k, dmax;
        int TetSize = delcx->DeluanayTet.size() - delcx->redundantCount;
        int INF = delcx->DeluanayTet.size();

        for (int i = TetSize - 1; i > 0; i--)
        {
                tetindex = sortedTet[i];
                delcx->DeluanayTet[tetindex].Depth = tetindex;

                dmax = delcx->DeluanayTet[tetindex].Depth;

                for (int j = 1; j <= 4; j++)
                {
                        k = delcx->DeluanayTet[tetindex].Neighbours[j];

                        if ((k == 0) && (delcx->DeluanayTet[tetindex].TetFlow[j] == 1))
                        {
                                dmax = INF;
                                break;
                        }

                        if (delcx->DeluanayTet[tetindex].TetFlow[j] == 1)
                        {
                                dmax = myMax(dmax, delcx->DeluanayTet[k].Depth);
                        }
                }

                delcx->DeluanayTet[tetindex].Depth = dmax;
        }

        /*FILE *fp = fopen("depthfile","w");
        for (uint idx = 1; idx < delcx->DeluanayTet.size(); idx++)
        {
                fprintf(fp,"%d %d %d %d %d %d %d\n",idx,delcx->DeluanayTet[idx].Status,delcx->DeluanayTet[idx].TetFlow[1], delcx->DeluanayTet[idx].TetFlow[2],delcx->DeluanayTet[idx].TetFlow[3],delcx->DeluanayTet[idx].TetFlow[4],delcx->DeluanayTet[idx].Depth);
        }
        fclose(fp);*/
}
