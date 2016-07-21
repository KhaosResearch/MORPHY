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
/** Aim: Class to represent transition probability between the states while alignining sequences         */
/*********************************************************************************************************/


#ifndef PROBABILITYOBJECT_H
#define PROBABILITYOBJECT_H

/** Header file Inclusions **/
#include <iomanip>

#include "Utils.h"
#include "TypeDefinitions.h"

using namespace std;
using namespace Utils;
using namespace TypeDefinitions;

/** Typedef Definitions **/


class ProbabilityObject {
    public:
        /** Constructors, Deconstructors and Operators Definitions **/
        ProbabilityObject();                                 //Default (empty Constructor)
        virtual ~ProbabilityObject();                        //Virtual Decontructor
        bool operator==  (const ProbabilityObject&);         //Definition of The equality   (==) operator
        bool operator!=  (const ProbabilityObject&);         //Definition of The inequality (!=) operator
        bool operator<   (const ProbabilityObject&) const;   //Overloading of the comparison operators
        bool operator<=  (const ProbabilityObject&) const;   //Overloading of the comparison operators
        bool operator>   (const ProbabilityObject&) const;   //Overloading of the comparison operators
        bool operator>=  (const ProbabilityObject&) const;   //Overloading of the comparison operators

         /** Setter And Getter Methods for Private Members **/
         void Set_OverallProbability       (double val                             ) { OverallProbability       = val; }
         void Set_Contributors             (vector<ProbabilityObject*> val         ) { Contributors             = val; }
         void Set_PosIn_1D                 (LongInt val                            ) { PosIn_1D                 = val; }
         void Set_PosIn_ND                 (vector<int> val                        ) { PosIn_ND                 = val; }
         void Set_PosIn_ND                 (LongInt val, const vector<int>& Dimensions) { GetIndexesForPos(val, Dimensions, PosIn_ND); }
         void Set_EndState                 (bool        val                        ) { EndState                 = val; }
         void Set_SoanLbl                  (string      val                        ) { SoanLbl                  = val; }
         void Set_ContributorTihls         (vector<string> val                     ) { ContributorTihls         = val; }
         void Set_ContributorProbabilities (vector<double> val                     ) { ContributorProbabilities = val; }

         double                            Get_OverallProbability()       const { return OverallProbability;       }
         const vector<ProbabilityObject*>& Get_Contributors()             const { return Contributors;             }
         LongInt                           Get_PosIn_1D()                 const { return PosIn_1D;                 }
         const vector<int>&                Get_PosIn_ND()                 const { return PosIn_ND;                 }
         bool                              Is_EndState()                  const { return EndState;                 }
         string                            Get_SoanLbl()                  const { return SoanLbl;                  }
         const vector<string>&             Get_ContributorTihls()         const { return ContributorTihls;         }
         const vector<double>&             Get_ContributorProbabilities() const { return ContributorProbabilities; }

         /** Extra Functions: **/
         void               AddContributor(ProbabilityObject* contr);              //Function that adds a new contributor to the collection (vector) of contributos
         ProbabilityObject* GetContributorAtPos(int pos) const;                    //Function that returns the contributor whose position equals to the given -as argument- position
         void               AddContributorTihl(string tihl);                       //Function that adds a new string corresponding to a contributor
         string             GetContributorTihlAtPos(int pos) const;                //Function that returns the tihl of the contributor whose position equals to the given -as argument- position
         void               AddContributorTihlExt(string tihl);                    //Function that adds a new string corresponding to a contributor
         string             GetContributorTihlExtAtPos(int pos) const;             //Function that returns the tihl of the contributor whose position equals to the given -as argument- position
         void               AddContributorProbability(double val);                 //Function that adds a new probability value corresponding to a contributor
         double             GetContributorProbabilityAtPos(int pos) const;         //Function that returns the contributed probability of the contributor whose position equals to the given -as argument- position
         void               IncreaseOverallProbabilityBy(double val);              //Function to increase the OverallProbability by the given value
         string             GetBriefDescription() const;                           //Function to return a brief description (i.e. without accoubting for its contributors) of current ProbabilityObject
         void               AddContributorWingFold(bool val);                      //Function that adds true or false according to whether the given contributor is wingfolded or not
         bool               GetContributorWingFoldAtPos(int pos) const;            //Function that returns whether the contributor atb the given position is wing folded or not
         void               AddContributorDirectTransProb(double val);             //Function that adds the direct transition probability corresponding to a contributor
         double             GetContributorDirectTransProbAtPos(int pos) const;     //Function that returns the direct transition probability corresponding to the given contributor
         void               AddContributorInDirectTransProb(double val);           //Function that adds the indirect transition probability corresponding to a contributor
         double             GetContributorInDirectTransProbAtPos(int pos) const;   //Function that returns the indirect transition probability corresponding to the given contributor
         void               AddContributorSelfTransProb(double val);               //Function that adds the self transition probability (to the non emitting state) corresponding to a contributor
         double             GetContributorSelfTransProbAtPos(int pos) const;       //Function that returns the self transition probability (to the non emitting state) corresponding to the given contributor

    protected:

    private:
        /** Class Members **/
        LongInt                    PosIn_1D;                     //The index of the current ProbabilityObject transformed in the 1D plane.
        vector<int>                PosIn_ND;                     //The indices of the current ProbabilityObject transformed in the ND plane.
        bool                       EndState;                     //Checks if current state is an end state (i.e. all indexes of current state equal to the lengths of the sequences)
        string                     SoanLbl;                      //The current label of this State (Soan)
        double                     OverallProbability;           //The overall probability formed as the sum of the probabilities of its contributors
        vector<ProbabilityObject*> Contributors;                 //A vector containing all individual contributors leading to formning this overall probability
        vector<string>             ContributorTihls;             //A vector with all the corresponding tihls for each Contributor
        vector<string>             ContributorTihlsExt;          //A vector with an extensive description of the corresponding tihls (i.e. expressing the wing folding history) for each Contributor
        vector<double>             ContributorProbabilities;     //A vector with all the corresponding probabilities for each Contributor
        vector<bool>               ContributorWingFolds;         //A vector with boolean values representing if the corresponding contributor is wing folded or not
        vector<double>             ContibutorDirectTransition;   //A vector which for each contibutor stores the probability for direct transition to current state (for wing folded cases)
        vector<double>             ContibutorInDirectTransition; //A vector which for each contibutor stores the probability for indirect transition to current state (for wing folded cases)
        vector<double>             ContibutorSelfTransition;     //A vector which for each contibutor stores the probability for self-transitioning to the non-emitting state (for wing folded cases)
};

/** Definition of the toString Operation **/
ostream& operator<<(ostream& ostr, const ProbabilityObject& prob);

#endif // PROBABILITYOBJECT_H
