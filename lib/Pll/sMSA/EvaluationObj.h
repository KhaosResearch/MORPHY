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

/*********************************************************************************************************/
/**                                         Class Scopes:                                                */
/**        Aim: Class responsible to assess an alignment compared to a reference alignment using         */
/**                various assessment metrics (e.g. MCC, TC score, SP score etc)                         */
/*********************************************************************************************************/

#ifndef EVALUATIONOBJ_H
#define EVALUATIONOBJ_H

/** Header file Inclusions **/
#include <iostream>
#include <iomanip>
#include <cmath>
#include "omp.h"

#include "TypeDefinitions.h"
#include "Utils.h"

using namespace std;
using namespace TypeDefinitions;
using namespace Utils;

class EvaluationObj {

    public:
        /** Constructors, Deconstructors and Operators Definitions **/
        EvaluationObj();                                              //Default Constructor
        virtual ~EvaluationObj();                                     //Virtual Decontructor

         /** Custom Structs for Scoring Representation: **/
        struct MCCScoreSummary {
            double    MCCScoreValue;    //The pairwise MCC score of the examined alignments
            double    SensitivityScore; //The sensitivity Score of the alignment (i.e. TP/(TP+FN)
            double    SelectivityScore; //The selectivity Score of the alignment (i.e. TP/(TP+FP)
        };

        struct SPSimScoreSummary {
            double    SPScoreValue;     //The pairwise SP score of the alignment
            int       SPSimScore;       //The representing the sum of similar pairs found within the alignment
            int       TotalNumOfPairs;  //The total number of pairs found in the alignment
        };

        struct SPDistScoreSummary {
            double    SPScoreValue;     //The pairwise SP score of the alignment
            int       SPDistScore;      //The representing the sum of different pairs found within the alignment
            int       TotalNumOfPairs;  //The total number of pairs found in the alignment
        };

        struct TCScoreSummary {
            double    TCScoreValue;     //The columnwise TC score of the examined alignments
            int       SharedColumns;    //The number of columns exactly shared by both alignments
            int       RefColumns;       //The total number of alignment columns in the reference alignment
        };


        /** Extra Functions: **/
        void     ConvertAlignmentToMatrix(StringVecRef alignment, IntMatRef mat); //Function to convert a Given alignment into a matrix of indices representing the alignment
        string   GetMCCScore(IntMatRef predMat, IntMatRef refMat);                //Function to return the MCC score after comparing the predicted and the reference alignment at a pairwise level
        string   GetSPSimScore(StringVecRef align);                               //Function to return the similarity-based SP (Sum of Pairs) score at a pairwise level. Shows how close the aligned sequences are to each other
        string   GetSPDistScore(StringVecRef align);                              //Function to return the distance-based SP (Sum of Pairs) score at a pairwise level. Shows how distant the aligned sequences are to each other
        string   GetTCScore(StringVecRef pred, StringVecRef ref);                 //Function to return the TC score after comparing the predicted and the reference alignment at a columnwise level

    protected:

    private:
};

/** Definition of the toString Operation **/
ostream& operator<<(ostream& ostr, const EvaluationObj& evalObj);

#endif // EVALUATIONOBJ_H
