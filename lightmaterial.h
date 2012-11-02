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

#ifndef LIGHTMATERIAL_H
#define LIGHTMATERIAL_H

class LightMaterial
{
    public:
        static float DimLightModelAmb[4];
        static float MediumLightModelAmb[4];
        static float BrightLightModelAmb[4];

        static float DimTopLightAmb[4];
        static float DimTopLightDif[4];
        static float DimTopLightPos[4];
        static float DimTopLightDir[4];

        static float DimBottomLightAmb[4];
        static float DimBottomLightDif[4];
        static float DimBottomLightPos[4];
        static float DimBottomLightDir[4];

        static float MediumTopLightAmb[4];
        static float MediumTopLightDif[4];
        static float MediumTopLightPos[4];
        static float MediumTopLightDir[4];

        static float MediumBottomLightAmb[4];
        static float MediumBottomLightDif[4];
        static float MediumBottomLightPos[4];
        static float MediumBottomLightDir[4];

        static float BrightTopLightAmb[4];
        static float BrightTopLightDif[4];
        static float BrightTopLightSpec[4];
        static float BrightTopLightPos[4];
        static float BrightTopLightDir[4];

        static float BrightBottomLightAmb[4];
        static float BrightBottomLightDif[4];
        static float BrightBottomLightSpec[4];
        static float BrightBottomLightPos[4];
        static float BrightBottomLightDir[4];

        /*static float MatAmb[19][4];
        static float MatDiff[19][4];
        static float MatSpec[19][4];
        static float MatShin[19];*/

        static float MatAmb[31][4];
        static float MatDiff[31][4];
        static float MatSpec[31][4];
        static float MatShin[31];

        static float MatAmb2[4];
        static float MatDiff2[4];
        static float MatSpec2[4];
        static float MatShin2[1];

        static float MatEmission[4];
};

#endif // LIGHTMATERIAL_H
