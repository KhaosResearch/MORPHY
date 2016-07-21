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

/**********************************************************************************************/
/**                                       File Scopes:                                        */
/**  Aim: This Class provides with set of useful interactive menus to get user preferences    */
/**********************************************************************************************/

#ifndef MENUS_H
#define MENUS_H

#include <iostream>
#include <stdio.h>

#include "TripleInitializer.h"
#include "PairInitializer.h"
#include "Utils.h"
#include "TypeDefinitions.h"
#include "EvaluationObj.h"

using namespace std;
using namespace Utils;
using namespace TypeDefinitions;

class Menus {
    public:
        Menus();
        virtual ~Menus();

        /** Extra Functions **/
        void ShowMainMenu();                                                    //Function to display the available options to the standard output stream
        void InitTripleInitializer(int SubTree, bool LoadStartingHomologies);   //Function to initialize properly the TripleInit object
        void InitPairInitializer(int Branch, bool LoadStartingHomologies);      //Function to initialize properly the PairInit object
        void ExecuteTriplewiseResampling(int rounds);                           //Function to Perform the left/right subtree triplewise re-alignment
        void ExecutePairwiseResampling(int rounds);                             //Function to Perform pairwise re-alignment
        void ExecuteResampling();                                               //Function to start the re-alignment
        void GetNewHomologyPair(int Branch, IntVecRef RegionStartIdxs, string RegionStartTihlLbl, bool HasExtraStartState, bool HasExtraEndState); //Function to Estimate the new Homology Structure based on the re-sampled evolution history
        void GetNewHomologyTriple(int SubTree, IntVecRef RegionStartIdxs, string RegionStartTihlLbl, bool HasExtraStartState, bool HasExtraEndState); //Function to Estimate the new Homology Structure based on the re-sampled evolution history
        void GetNewAlignment();                                      //Function to Estimate the new Homology Structure based on the newly estimated Homology Structure
        void ExpandNs();                                             //Function to Expand the Ns of the Newly Estimated Evolution History
        void GetAlignmentFromFile(StringVecRef algn, string filename); //Function to read the given file and store the alignment (ignoring the sequence headers) in the given vector of strings
        void LoadNewHomology();                                      //Function to Load a starting Homology from a file
        void SaveNewHomology();                                      //Function to Save the newly Estismated Homology to a file


        void UpdateNewAlignment();                                   //Function to save the Newly estimated Alignment back to the file
        void PrintNewEvolHistory();                                  //Function to Print the Newly sampled Evolution History to the Standard Output Stream
        void PrintNewHomology();                                     //Function to Print the Newly estimated Homology Structure to the Standard Output Stream
        void PrintNewAlignment();                                    //Function to Print the Newly estimated Alignment to the Standard Output Stream
        void PrintIntervals(IntervalVecRef itvs);                    //Function to Print the Created Intervals
        void PrintAlignmentRegion(IntVecRef StartIdxs, IntVecRef EndIdxs, string StartTihlLbl, string StartSoanLbl, bool HasExtraStartState, bool HasExtraEndState, IntVecRef Dimensions, DoubleVecRef Lengths);     //Function to Print the specification for a given alignment region

        int GetIdxOfAllHomologousColumn(IntVecRef SubTreeIdxs, IntMatRef OverallHomologies, int StartPos, int Direction);

    protected:

    private:
        /** Class Members **/
        TripleInitializer*   TripleInit;              //An instance that will initialize all the program prerequisites for the TripleWise Alignment
        PairInitializer*     PairInit;                //An instance that will initialize all the program prerequisites for the PairWise Alignment
        StringVec            NewEvolHistory;          //A Vector containing the final (after re-aligning all sites) re-sampled evolution history
        IntMat               NewHomology;             //A Matrix representing the new Homology structure of the newly re-sampled evolution history
        StringVec            NewAlignment;            //A Vector representing the new Alignment of the input sequences after the a re-sampling step has been completed
};

#endif // MENUS_H
