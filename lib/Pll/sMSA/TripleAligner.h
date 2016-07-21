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
/**    Aim: Class responsible to perform the statistical Multiple Sequence Alignment for a given Tree    */
/*********************************************************************************************************/

#ifndef TRIPLEALIGNER_H
#define TRIPLEALIGNER_H

/** Header file Inclusions **/
#include <float.h>
#include <queue>
#include <map>
#include <iomanip>
#include <cmath>

#include "Parameters.h"
#include "Triplet.h"
#include "ProbabilityObject.h"
#include "BoostRandomGenerator.h"
#include "F84Model.h"
#include "TN93Model.h"
#include "GTRModel.h"
#include "TypeDefinitions.h"
#include "Utils.h"

using namespace std;
using namespace Mappings;
using namespace Utils;
using namespace TypeDefinitions;

/** Typedef Definitions **/
typedef vector< vector<Triplet* > >             TripletMat;
typedef map<LongInt, ProbabilityObject*>        ProbMap;

class TripleAligner {
    public:
        /** Constructors, Deconstructors and Operators Definitions **/
        TripleAligner();                              //Default Constructor
        virtual ~TripleAligner();                     //Virtual Decontructor
        TripleAligner(const TripleAligner& other);          //Copy Constructor to Clone objects
        bool operator==  (const TripleAligner&);      //Definition of The equality   (==) operator
        bool operator!=  (const TripleAligner&);      //Definition of The inequality (!=) operator

        /** Setter And Getter Methods for Private Members **/
        Parameters*             Get_Parameters()       const { return Params;           }
        FIDModel*               Get_FidModel()         const { return Fid;              }
        ProbMap                 Get_Probabilities()    const { return Probabilities;    }
        IntCVecRef              Get_Dimensions()       const { return Dimensions;       }
        DoubleCVecRef           Get_Lengths()          const { return Lengths;          }
        DoubleCHypRef           Get_FPTables()         const { return FPTables;           }
        const F84Model*         Get_F84()              const { return f84;              }
        const TN93Model*        Get_TN93()             const { return tn93;             }
        const GTRModel*         Get_GTR()              const { return gtr;              }
        IntCMatRef              Get_HomologiesIdxs()   const { return HomologiesIdxs;   }

        void Set_Parameters(Parameters* val)          { Params           = val; }
        void Set_FidModel(FIDModel* val)              { Fid              = val; }
        void Set_Probabilities(ProbMap val)           { Probabilities    = val; }
        void Set_Dimensions(IntVecRef val)            { Dimensions       = val; }
        void Set_Lengths(DoubleVecRef val)            { Lengths          = val; }
        void Set_FPTables(DoubleHypRef val)           { FPTables         = val; }
        void Set_F84(F84Model* val)                   { f84              = val; }
        void Set_TN93(TN93Model* val)                 { tn93             = val; }
        void Set_GTR(GTRModel* val)                   { gtr              = val; }
        void Set_HomologiesIdxs(IntCMatRef val)       { HomologiesIdxs   = val; }
        void Set_EndPosIn_ND(IntVecRef val)           { EndPosIn_ND      = val; }
        void Set_StartPosIn_ND(IntVecRef val)         { StartPosIn_ND    = val; }
        void Set_SaveToLogFile(bool val, string fName){ SaveToLogFile    = val; LogFileName = fName; }
        void Set_RandomGenerator(BoostRandomGenerator* val){ MyRand      = val; }

        /** Extra Functions: **/
        //Starter Functions

        //Initialization Functions
        void  InitParameters(Parameters* val);                               //Function to Initialize all Program's Parameters
        void  InitSubModel(F84Model* val1, TN93Model* val2, GTRModel* val3); //Function to Initialize the Substitution Model
        void  InitFidModel();                                                //Function to Initialize the FID Model
        void  InitHomologyIndices(IntCMatRef Homologs);                      //Function to Initialize the homologies structure from the Given HomologyTable and  Indices to be used
        void  InitLengths(DoubleCVecRef val);                                //Function to Initialize the corresponding branch lengths
        void  InitDimensions(IntCVecRef val);                                //Function to Initialize the Dimensions (equal to the sequence sizes)
        void  InitFPTables(DoubleCHypRef val);                               //Function to Initialize the Felsenstein Pruning Tables of the Tree according to the given values
        void  InitTriplets();                                                //Function to Initialize all Possible triplets
        void  InitProbabilities(string startSoanLbl, string startTihlLbl, IntVecRef startIdxs);   //Function to intialize all ProbabilityObjects with default values
        void  InitSubMatrices(DoubleCMatRef val);                            //Function to Initialize the Substitution Matrices for the given branch lengths


        //Assisting Functions for In-Between calculations
        void     GetAllReachableSoans(vector<string> &soans, const vector<string> &tihls);            //Function to Create All Reachable Soans, based on all possible tihls and an initial (all H-labeled) soan
        void     GetAllPossibleTriplets(const vector<string>& soans, const vector<string>& tihls);    //Function to Create All Possible, based on all possible tihls and all reachable soans
        bool     CanTihlBeAppliedToSoan(string soan, string tihl);                                    //Function that examines whether a given tihl may be applied to a given soan or not
        string   ApplyTihlToSoan(string soan, string tihl);                                           //Function that Applies a given Tihl to a given Soan
        string   ApplyTihlToSoan(string soan, string tihl, bool PairWise);                            //Function that Applies a given Tihl to a given Soan (for pairwise alignment)
        Triplet* Get_ToNonEmmitingTriplet(int r_idx) const;                                           //Function to return the triplet corresponding to the given index, that leads to the non emmitting state
        Triplet* Get_FromNonEmmitingTriplet(string sAfter, IntCVecRef idxs) const;                    //Function to return the triplet that leaves from the non emmitting state and transits to a soan labeled according to the given argument
        int      IndexOfSoan(string soan);                                                            //Function that returns the correspondinx index (from 0 - soans.size()-1) that the given soan corresponds to in the vector with all possible soans (e.g. HHH is at index 0, HHB is at index 1,... ebB is at index 40....)
        bool     IsValidTransition(IntVecRef NextIdxs, string tihl);                                  //Function that checks wheteher a Transition is Valid according to the given Homologies
        double   GetSubstitutionProbability(string tihlLbl, IntCVecRef Idxs, int stx) const;          //Function to get the Corresponding Substitution Probability for the given transition to the next state
        Triplet* GetEndStateTriplet(string StartSoanLbl);                                             //Function to Get the Triplet starting with the given SoanLbl and transitioning to the End State (i.e. HHH..H) via the corresponding all Homologous Tihl

        //Probabilities Estimation Functions
        void   EstimateProbabilities(int choice, bool HasExtraEndState);                //Function to calculate all the probabilities describing all evolutionary histories for the alignement of the given sequences. Argument defines which MakeTransition Function should be called
        void   MakeTransition(IntVecRef CurIdxs, LongInt CurPosIn_1D, int stx, bool HasExtraEndState); //Given Homologies: Overloaded Function to make the transition from current state to the next state
        void   TransitToNonEmittingState(); //Function to Transit to the Non Emitting State for the cases where the End State is (H)EEE (i.e. state 26)

        //Functions to Sample Alignments
        void   GetEvolHistories(StringMatRef evolHist, DoubleVecRef probs); //Function to sample evolutionary histories based on the calculated Probabilities
        void   GetEvolHistory(StringVecRef hist, double& prob, bool FixNs); //Function to sample only one evolutionary history based on the calculated Probabilities
        void   GetEvolHistory(StringVecRef hist, double& prob);             //Function to sample only one evolutionary history based on the calculated Probabilities
        void   GetEvolHistory(StringVecRef hist, double& prob, StringVecRef soans, StringVecRef tihls, IntMatRef idxs); //Function to sample only one evolutionary history based on the calculated Probabilities and also "return" the sequences od soans as well as their corresponding indices

        //Functions to output results
        void SaveResultsToFile(string filename); //Function to save the estimated probabilities to file
        void PrintParameters();                  //Function to Output the Program Parameters to the standard output stream
        void PrintSubstitutionMatrices();        //Function to Output the Substitution Matrices of interest to the standard output stream
        void PrintLengths();                     //Function to Print the subtree lengths
        void PrintFPTables();                    //Function to Print the Felsenstein Pruning Tables for the Current Fragment to the standard output stream
        void PrintDimensions();                  //Function to Print the Dimensions for the current fragment to the standard output stream
        void PrintTriplets();                    //Function to Print the Estimated Triplets
        void PrintModels();                      //Function to Print the Substitution and Fid Models
        void PrintBoundaries();                  //Function to Print the Alignment Boundaries
        void PrintStatus();                      //Function to call the above printing functions in order to output the curretn aligner instance setup
        void SaveFPTables();                     //Fucntion to save the employed FPTables to an output file

    protected:

    private:
        /** Class Members **/
        Parameters*    Params;               //The Aligner Parameters
        FIDModel*      Fid;                  //The Fid Model
        TripletMat     Triplets;             //All the available triplets organized per each distinct reachable soan (i.e. 41 vectors of triplets, each vector representing the triplets per each reachable soan)
        ProbMap        Probabilities;        //All calculated Probabilities (size equal to the product of the (sequences lengths + 1) multiplied by the number of reachable soans -1 (not to account for the non emitting state since wing folding will be considered)
        IntVec         Dimensions;           //The Lengths (i.e. number of nucleotides) per sequence and the number of reachable Soans
        DoubleVec      Lengths;              //The corresponding Branch lengths for the sequences to be aligned
        int            idx_R;                //The row index in the Triplets vector of the non emitting triplet
        int            idx_C;                //The column index in the Triplets vector of the non emitting triplet
        F84Model*      f84;                  //The F84 Model
        TN93Model*     tn93;                 //The TN93 Model
        GTRModel*      gtr;                  //The GTR Model
        DoubleHyp      FPTables;             //The hypermatrix (3D) representing the felsenstein probabilities per each Sequence
        DoubleMat      SubMats;              //A Matrix representing the substitution matrices for each branch length
        IntMat         HomologiesIdxs;       //A Matrix representing the homology structures between the sequences
        IntMat         Anchors;              //A Matrix representing the homology structures between the sequences converted in the form of Anchor points with regards to the root
        LongInt        StartPosIn_1D;        //The starting position (index) of the Markov Chain in 1D
        IntVec         StartPosIn_ND;        //The starting positions (indices) of the Markov Chain in ND
        string         StartTihlLbl;         //The starting tihl label
        LongInt        EndPosIn_1D;          //The ending position (index) of the Markov Chain in 1D
        IntVec         EndPosIn_ND;          //The ending positions (indices) of the Markov Chain in ND
        bool           SaveToLogFile;        //A boolean variable to define whether a Log of the visited states will be kept or not
        string         LogFileName;          //The name of the file that will be used to log the visited states
        BoostRandomGenerator* MyRand;        //The Random Generator Object (so as not to create a new one per each run)

};

/** Definition of the toString Operation **/
ostream& operator<<(ostream& ostr, const TripleAligner& aligner);

#endif // TRIPLEALIGNER_H
