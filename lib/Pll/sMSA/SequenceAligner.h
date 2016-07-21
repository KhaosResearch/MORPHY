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

#ifndef SEQUENCEALIGNER_H
#define SEQUENCEALIGNER_H

/** Header file Inclusions **/
#include <float.h>
#include <queue>
#include <map>

#include "Parameters.h"
#include "Triplet.h"
#include "ProbabilityObject.h"
#include "BoostRandomGenerator.h"
#include "F84Model.h"
#include "TN93Model.h"
#include "GTRModel.h"
#include "TypeDefinitions.h"

using namespace std;
using namespace Mappings;
using namespace Utils;
using namespace TypeDefinitions;

/** Typedef Definitions **/
typedef vector< vector<Triplet* > >         TripletMat;
typedef vector <vector< vector<double > > > DoubleMat;
typedef vector<ProbabilityObject*>          ProbMat;
typedef vector< vector<string > >&          StringMatRef;
typedef vector< vector<int > >              IntMat;
typedef vector< vector<int > >&             IntMatRef;
typedef vector<int>                         IntVec;
typedef vector<int>&                        IntVecRef;

class SequenceAligner {
    public:
        /** Constructors, Deconstructors and Operators Definitions **/
        SequenceAligner();                              //Default Constructor
        virtual ~SequenceAligner();                     //Virtual Decontructor
        SequenceAligner(const SequenceAligner& other);  //Copy Constructor to Clone objects
        bool operator==  (const SequenceAligner&);      //Definition of The equality   (==) operator
        bool operator!=  (const SequenceAligner&);      //Definition of The inequality (!=) operator

        /** Setter And Getter Methods for Private Members **/
        Parameters*             Get_Parameters()       const { return Params;           }
        FIDModel*               Get_FidModel()         const { return Fid;              }
        ProbMat                 Get_Probabilities()    const { return Probabilities;    }
        const vector<int>&      Get_Dimensions()       const { return Dimensions;       }
        const vector<string>&   Get_Sequences()        const { return Sequences;        }
        const vector<double>&   Get_Lengths()          const { return Lengths;          }
        const DoubleMat&        Get_SeqMat()           const { return SeqMat;           }

        void Set_Parameters(Parameters* val)          { Params           = val; }
        void Set_FidModel(FIDModel* val)              { Fid              = val; }
        void Set_Probabilities(ProbMat val)           { Probabilities    = val; }
        void Set_Dimensions(vector<int> val)          { Dimensions       = val; }
        void Set_Sequences(vector<string> val)        { Sequences        = val; }
        void Set_Lengths(vector<double> val)          { Lengths          = val; }
        void Set_SeqMat(DoubleMat val)                { SeqMat           = val; }

        /** Extra Functions: **/
        //Starter Functions
        void          AlignPlainSequences();                            //Function to align plain sequences
        void          AlignProbSequences();                             //Function to align sequences provided in the form of probability matrices
        void          SampleOnGivenAlignment();                         //Function to estimate all probabilities for a given (pre-defined) alignment
        void          SampleOnGivenHomologies();                        //Function to estimate all probabilities for a given (pre-defined) set of Homologies between the sequences

        //Initialization Functions
        void          InitParameters(string filename);                  //Function to Initialize all Program's Parameters
        void          InitFidModel();                                   //Function to Initialize the FID Model
        void          InitSequences();                                  //Function to Initialize the tree (i.e. Sequences to be aligned as well as equivalent branch lengths) (plain sequences form)
        void          InitSequences(string filename);                   //Function to Initialize the tree (i.e. Sequences to be aligned as well as equivalent branch lengths) (plain sequences form - read from file)
        void          InitSequences(int x);                             //Function to Initialize the tree (i.e. Sequences to be aligned as well as equivalent branch lengths) (felsenstein form)
        void          InitTriplets();                                   //Function to Initialize all Possible triplets
        void          InitSubModel();                                   //Function to Initialize the Substitution Model
        void          InitSubMatrices();                                //Function to Initialize the Substitution Matrices for the given branch lengths
        void          InitProbabilities();                              //Function to intialize all ProbabilityObjects with default values
        void          InitAlignedSequences();                           //Function to Get a given alignment from the user
        void          InitAlignedSequences(string filename);            //Function to Get a given alignment from a file
        void          InitHomologyIndices();                            //Function to Get a given homology from the user
        void          InitHomologyIndices(string filename);             //Function to Get a given homology from a file

        //Assisting Functions for In-Between calculations
        string        ApplyTihlToSoan(string soan, string tihl);        //Function that Applies a given Tihl to a given Soan
        bool          CanTihlBeAppliedToSoan(string soan, string tihl); //Function that examines whether a given tihl may be applied to a given soan or not
        void          GetAllReachableSoans(vector<string> &soans, const vector<string> &tihls);  //Function to Create All Reachable Soans, based on all possible tihls and an initial (all H-labeled) soan
        void          GetAllPossibleTriplets(const vector<string>& soans, const vector<string>& tihls);    //Function to Create All Possible, based on all possible tihls and all reachable soans
        int           IndexOfSoan(string soan);                         //Function that returns the correspondinx index (from 0 - soans.size()-1) that the given soan corresponds to in the vector with all possible soans (e.g. HHH is at index 0, HHB is at index 1,... ebB is at index 40....)
        Triplet*      Get_ToNonEmmitingTriplet(int r_idx) const;        //Function to return the triplet corresponding to the given index, that leads to the non emmitting state
        Triplet*      Get_FromNonEmmitingTriplet(string sAfter, const vector<int>& idxs) const;  //Function to return the triplet that leaves from the non emmitting state and transits to a soan labeled according to the given argument
        double        GetSubstitutionProbability(string tihlLbl, const vector<int>& Idxs) const; //Function to get the Corresponding Substitution Probability for the given transition to the next state
        double        GetSubstitutionProbability(string tihlLbl, const vector<int>& Idxs, int x) const; //Overloaded Function to get the Corresponding Substitution Probability for the given transition to the next state
        void          GetIndicesFromAlignment(vector<string>& alignment, IntMatRef idxs);  //Function to get the corresponding representation of the given alignment in the form of sequence indices
        void          GetIndicesFromEvolHistory(vector<string>& evolHist, IntMatRef idxs); //Function to get the corresponding representation of the given evolutionary history in the form of sequence indices
        bool          IsValidTransition(IntVecRef NextIdxs, string tihl);                  //Function that checks wheteher a Transition is Valid according to the given Homologies

        //Probabilities Estimation Functions
        void          EstimateProbabilities(int choice);                //Function to calculate all the probabilities describing all evolutionary histories for the alignement of the given sequences. Argument defines which MakeTransition Function should be called
        void          MakeTransition(vector<int>& CurIdxs, int CurPosIn_1D); //Plain Sequences: Function to make the transition from current state to the next state
        void          MakeTransition(vector<int>& CurIdxs, int CurPosIn_1D, int x);      //Felsenstein Matrices: Overloaded Function to make the transition from current state to the next state
        void          MakeTransition(vector<int>& CurIdxs, int CurPosIn_1D, float x);    //Given Alignment: Overloaded Function to make the transition from current state to the next state
        void          MakeTransition(vector<int>& CurIdxs, int CurPosIn_1D, bool reset); //Given Homologies: Overloaded Function to make the transition from current state to the next state
        void          TransitToEndState();                              //Function to make the final transition to the End State (i.e. the all H(omologous) edges labeled state)

        //Functions to Sample Alignments
        void          GetEvolHistories(StringMatRef evolHist, vector<double>&  probs);//Function to sample evolutionary histories based on the calculated Probabilities
        void          GetAlignmentForEvolHistory(vector<string >& evolHist, vector<string >& align); //Function that returns the corresponding alignment for the given evolutionary history
        void          AssessAlignment();                                //Function to assess the probability of an alignment given from the keyboard based on the computed Probability Objects
        double        AssessAlignment(vector<string>& alignment);       //Function to assess the probability of a given alignment based on the computed Probability Objects

        //Functions to output results
        void          SaveResultsToFile();                              //Function to save derived results to file
        void          ShowAllSampledAlignments(StringMatRef evolHistories, vector<double>& probs);//Function to print to screeen all Sampled Alignments alongside with their corresponding evolutionary histories
        void          PrintSequenceProbabilities();                     //Function to output the initial probabilities of the Sequences
        void          PrintSubstitutionMatrices();                      //Function to output to the screen the Values of the Substitution matrices according to the substitution model and the corresponding branch lengths
        string        GetVectorInPrintableForm(string description, IntVecRef vec);  //Function to return a string with the vector's values in a printable for (e.g. <1 2 3>)

    protected:

    private:
        /** Class Members **/
        Parameters*    Params;               //The Aligner Parameters
        FIDModel*      Fid;                  //The Fid Model
        TripletMat     Triplets;             //All the available triplets organized per each distinct reachable soan (i.e. 41 vectors of triplets, each vector representing the triplets per each reachable soan)
        ProbMat        Probabilities;        //All calculated Probabilities (size equal to the product of the (sequences lengths + 1) multiplied by the number of reachable soans -1 (not to account for the non emitting state since wing folding will be considered)
        vector<int>    Dimensions;           //The Lengths (i.e. number of nucleotides) per sequence and the number of reachable Soans
        vector<string> Sequences;            //The Given Sequences to be aligned
        vector<double> Lengths;              //The corresponding Branch lengths for the sequences to be aligned
        int            idx_R;                //The row index in the Triplets vector of the non emitting triplet
        int            idx_C;                //The column index in the Triplets vector of the non emitting triplet
        F84Model*      f84;                  //The F84 Model
        TN93Model*     tn93;                 //The TN93 Model
        GTRModel*      gtr;                  //The GTR Model
        DoubleMat      SeqMat;               //The vector representing the felsenstein probabilities per each Sequence
        vector< vector<double> > SubMats;    //A Vector of Vectors representing the substitution matrices for each branch length
        vector<string> AlignedSequences;     //A Vector for storing the already aligned sequences
        vector<int>    AlignedSequencesIdxs; //A Vector representing an alignment in the form of Indices
        IntMat         HomologiesIdxs;       //A Vector representing the homology structures between the sequences
        IntMat         Anchors;              //A Vector representing the homology structures between the sequences converted in the form of Anchor points with regards to the root
};

/** Definition of the toString Operation **/
ostream& operator<<(ostream& ostr, const SequenceAligner& aligner);

#endif // SEQUENCEALIGNER_H
