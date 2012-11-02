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

#ifndef DISJOINTSET_H
#define DISJOINTSET_H

#include <vector>
#include <node.h>
#include <stdio.h>
#include <map>
#include <cassert>

class DisJointSet
{
    private:
        int elementCount;
        int setCount;
        std::vector<Node> nodes;

    public:
        DisJointSet();
        ~DisJointSet();
        int ElementCount();
        int SetCount();

        int FindSet(int element);
        void Union(int s1, int s2);
        void Add(int element);

        int FindSet(int element,FILE *fp);
        void Union(int s1, int s2,FILE *fp);
        void Add(int element,FILE *fp);

        void Clear();

        int GetElementAt(int i);
        std::vector<std::vector<int> > Consolidate();
};

#endif // DISJOINTSET_H
