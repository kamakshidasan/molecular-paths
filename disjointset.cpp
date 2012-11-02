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

#include "disjointset.h"

DisJointSet::DisJointSet()
{
        elementCount = setCount = 0;
}

DisJointSet::~DisJointSet()
{
}

/*!
    \fn DisJointSet::ElementCount()
 */
int DisJointSet::ElementCount()
{
    return elementCount;
}

/*!
    \fn DisJointSet::SetCount()
 */
int DisJointSet::SetCount()
{
        return setCount;
}

/*!
    \fn DisJointSet::FindSet(int element)
 */
int DisJointSet::FindSet(int element)
{
        assert (element < elementCount);

        Node *current = &nodes[element];
        while (current->Parent != 0)
        {
                current = current->Parent;
        }
        Node *root = current;

        current = &nodes[element];
        while (current != root)
        {
                Node *next = current->Parent;
                current->Parent = root;
                current = next;
        }
        return root->Index;
}

/*!
    \fn DisJointSet::FindSet(int element,FILE *fp)
 */
int DisJointSet::FindSet(int element,FILE *fp)
{
        assert (element < elementCount);

        fprintf(fp,"element = %d\n",element);

        Node *current = &nodes[element];
        fprintf(fp,"current element = %d ",current->ElementIndex);

        while (current->Parent != 0)
        {
                current = current->Parent;
                fprintf(fp,"current element = %d ",current->ElementIndex);
        }
        Node *root = current;
        fprintf(fp,"\n");

        current = &nodes[element];
        while (current != root)
        {
                Node *next = current->Parent;
                current->Parent = root;
                current = next;
        }
        fprintf(fp,"root element = %d root index = %d\n",root->ElementIndex,root->Index);
        return root->Index;
}

/*!
    \fn DisJointSet::Union(int s1, int s2)
 */
void DisJointSet::Union(int s1, int s2)
{
        assert((s1 < elementCount) && (s2 < elementCount));

        if (s1 == s2) return;

        Node *set1 = &nodes[s1];
        Node *set2 = &nodes[s2];

        if (set1->Rank > set2->Rank)
        {
                set2->Parent = set1;
        }
        else if (set2->Rank > set1->Rank)
        {
                set1->Parent = set2;
        }
        else
        {
                set2->Parent = set1;
                set1->Rank++;
        }
        setCount--;
        assert (setCount != 0);
}

/*!
    \fn DisJointSet::Union(int s1, int s2,FILE *fp)
 */
void DisJointSet::Union(int s1, int s2,FILE *fp)
{
        assert((s1 < elementCount) && (s2 < elementCount));
        fprintf(fp,"s1 = %d s2 = %d\n",s1,s2);

        if (s1 == s2) return;

        Node *set1 = &nodes[s1];
        Node *set2 = &nodes[s2];

        fprintf(fp,"set1 root = %d set2 root = %d\n",set1->ElementIndex,set2->ElementIndex);
        fprintf(fp,"set1 rank = %d set2 rank = %d\n",set1->Rank,set2->Rank);
        if (set1->Rank > set2->Rank)
        {
                set2->Parent = set1;
        }
        else if (set2->Rank > set1->Rank)
        {
                set1->Parent = set2;
        }
        else
        {
                set2->Parent = set1;
                set1->Rank++;
        }
        setCount--;
        fprintf(fp,"setCount = %d\n",setCount);
        assert (setCount != 0);
}

/*!
    \fn DisJointSet::Add(int element)
 */
void DisJointSet::Add(int element)
{
        Node newelement;
        newelement.ElementIndex = element;
        newelement.Index = elementCount;
        newelement.Parent = 0;//null
        newelement.Rank = 0;

        nodes.push_back(newelement);

        elementCount++;
        setCount++;
}

/*!
\fn DisJointSet::Add(int element,FILE *fp)
 */
void DisJointSet::Add(int element,FILE *fp)
{
        Node newelement;
        newelement.ElementIndex = element;
        newelement.Index = elementCount;
        newelement.Parent = 0;//null
        newelement.Rank = 0;

        fprintf(fp,"adding %d at %d\n",element,elementCount);
        nodes.push_back(newelement);

        elementCount++;
        setCount++;
        fprintf(fp,"elementcount = %d setcount = %d\n",elementCount,setCount);
}

/*!
    \fn DisJointSet::Clear()
 */
void DisJointSet::Clear()
{
        if(nodes.size() > 0)
        {
                nodes.clear();
        }
        elementCount = setCount = 0;
}

/*!
    \fn DisJointSet::GetElementAt(int i)
 */
int DisJointSet::GetElementAt(int i)
{
        return nodes[i].ElementIndex;
}
/*!
    \fn DisJointSet::Consolidate()
 */
std::vector<std::vector<int> > DisJointSet::Consolidate()
{
        int k = 0,dest;
        std::vector<std::vector<int> > ConsolidatedSet;
        std::map<int, int> Map;

        for (unsigned int i = 0; i < nodes.size(); i++)
        {
                if (nodes[i].Parent == 0)
                {
                        Map.insert(std::pair<int,int>(nodes[i].Index,k));
                        std::vector <int> temp;
                        ConsolidatedSet.push_back(temp);
                        ConsolidatedSet[k].push_back(nodes[i].ElementIndex);
                        k++;
                }
        }
        for (unsigned int i = 0; i < nodes.size(); i++)
        {
                if (nodes[i].Parent != 0)
                {
                        dest = Map[nodes[FindSet(i)].Index];
                        ConsolidatedSet[dest].push_back(nodes[i].ElementIndex);
                }
        }
        return(ConsolidatedSet);
}
