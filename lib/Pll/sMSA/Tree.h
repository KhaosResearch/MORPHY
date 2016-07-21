/****************************************************************************
 * This file is part of the program sMSA (statistical Multiple Sequence     *
 * Alignment) for joint sampling alignments and mutation parameter values   *
 * of n-given DNA-sequences.                                                *
 *                                                                          *
 *    Copyright (C) 2012  Dimitrios Lyras                                   *
 *    e-mail: dimlyras@gmail.com                                            *
 *                                                                          *
 *  This program is free software; you can redistribute it and/or modify    *
 *  it under the terms of the GNU General Public License as published by    *
 *  the Free Software Foundation; either version 2 of the License, or       *
 *  (at your option) any later version.                                     *
 *                                                                          *
 *  This program is distributed in the hope that it will be useful,         *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *  GNU General Public License for more details.                            *
 *                                                                          *
 *  You should have received a copy of the GNU General Public License       *
 *  along with this program; if not, write to the                           *
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330,             *
 *  Boston, MA  02111-1307  USA                                             *
 ****************************************************************************/

/***************************************************************/
/**                        File Scopes:                        */
/**            Class the represent the 4 leaved tree           */
/***************************************************************/

#ifndef TREE_H
#define TREE_H

/** Header file Inclusions **/
#include <iomanip>
#include "Sequences.h"

using namespace TypeDefinitions;
using namespace Mappings;

class Tree {

    public:
        /** Custom Structs for Nodes Representation: **/
        typedef struct {
            double    Length;  //The branch length between this and its adjacent node
            DoubleMat FPTable; //The Felsenstein Pruning table for the current InnerNode
        }InnerNode;

        struct Node {
            InnerNode N1;
            InnerNode N2;
            InnerNode N3;
        };

        /** Constructors, Deconstructors and Operators Definitions **/
        Tree();                                   //Default Constructor
        Tree(string filename, Sequences* Seqs);   //Constructor that creates the created Tree based on the input Tree file and the Given Sequences
        virtual ~Tree();                          //Virtual Decontructor

        /** Setter And Getter Methods for Private Members **/
        vector<Node>&       Get_TreeNodes()       { return TreeNodes;     }
        const vector<Node>& Get_TreeNodes() const { return TreeNodes;     }
        DoubleVecRef        Get_UniqueLengths()   { return UniqueLengths; }

        /** Extra Functions: **/
        double    GetLengthOfNode(int NodeIdx, int InnerNodeIdx);
        DoubleMat GetFPTableOfNode(int NodeIdx, int InnerNodeIdx);

    protected:

    private:
        /** Class Members **/
        vector<Node>    TreeNodes;      //The Nodes of the tree
        DoubleVec       UniqueLengths;  //The Unique Lengths of the tree

};

/** Definition of the toString Operation **/
ostream& operator<<(ostream& ostr, const Tree& tree);

#endif // TREE_H
