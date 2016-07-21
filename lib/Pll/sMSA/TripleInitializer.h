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

/***********************************************************************************/
/**                                 File Scopes:                                   */
/**   Aim: This Class provides with objects that will encapsulate all the required */
/**       by the program memebers (e.g. Parameters, SubMats, Homologies etc)       */
/***********************************************************************************/

#ifndef TRIPLEINITIALIZER_H
#define TRIPLEINITIALIZER_H

/** Header file Inclusions **/
#include <iostream>
#include <stdio.h>
#include "Utils.h"
#include "TypeDefinitions.h"
#include "Parameters.h"
#include "Sequences.h"
#include "Tree.h"
#include "F84Model.h"
#include "TN93Model.h"
#include "GTRModel.h"
#include "TripleAligner.h"

class TripleInitializer {
    public:
        /** Constructors, Deconstructors and Operators Definitions **/
        TripleInitializer();
        virtual ~TripleInitializer();

        /** Setter And Getter Methods for Private Members **/
        Sequences*            Get_Sequences()                          { return Seqs;                  }
        Parameters*           Get_Parameters()                         { return Params;                }
        Tree*                 Get_Tree()                               { return FTree;                 }
        F84Model*             Get_F84Model()                           { return f84;                   }
        TN93Model*            Get_TN93Model()                          { return tn93;                  }
        GTRModel*             Get_GTRModel()                           { return gtr;                   }
        DoubleMatRef          Get_SubstitutionMatrices()               { return SubMats;               }
        DoubleMatRef          Get_OverallSubstitutionMatrices()        { return OverallSubMats;        }
        DoubleVecRef          Get_OverallSubstitutionMatricesLengths() { return OverallSubMatsLengths; }
        StringVecRef          Get_Soans3Leaved()                       { return Soans3Leaved;          }
        StringVecRef          Get_Soans4Leaved()                       { return Soans4Leaved;          }
        BoostRandomGenerator* Get_RandomGenerator()                    { return RandGen;               }
        int                   Get_NumOfSteps()                         { return NumOfSteps;            }
        IntVecRef             Get_SubTreeIdxs()                        { return SubTreeIdxs;           }
        int                   Get_Subtree()                            { return SubTree;               }
        DoubleVecRef          Get_Lengths()                            { return Lengths;               }
        DoubleHypRef          Get_FPTables()                           { return FPTables;              }
        DoubleMatRef          Get_SubMatrices()                        { return SubMats;               }
        IntVecRef             Get_Dimensions()                         { return Dimensions;            }
        IntMatRef             Get_RegionHomologies()                   { return RegionHomologies;      }
        IntMatRef             Get_AlignmentHomologies()                { return AlignmentHomologies;   }
        IntMatRef             Get_SampledHomologies()                  { return SampledHomologies;     }
        IntMatRef             Get_Homologies()                         { return SampledHomologies;     }
        StringVecRef          Get_EvolutionHistory()                   { return EvolHistory;           }
        StringVecRef          Get_Soans()                              { return Soans;                 }
        StringVecRef          Get_Tihls()                              { return Tihls;                 }
        IntMatRef             Get_SoansIdxs()                          { return SoansIdxs;             }

        /** Extra Functions: **/
        //Initialization Functions
        void InitSequences();                                //Function to Read Sequences and Alignment from File and Initiliaze the Corresponding Sequences Object. The subtree argument is needed in order to delete all the "All-1" Columns from the OverallHomologies Matrix according to the specified subtree
        void InitParameters();                               //Function to Read the Desired parameters (e.g. ã, ë rates, Prior Nucleotide frequencies etc) and Initialize the Parameters instance
        void InitTree();                                     //Function to Initialize the Felsentein Tree Structure for the Given Sequences
        void InitOverallSubstitutionMatrices();              //Function to Initialize the substitution matrices for each branch length
        void InitTreeFPTables();                             //Function to Initialize the Felsenstein Pruning Tables on the Inner Nodes of the Tree
        void InitSoans(int treeSize);                        //Function to Initialize and store to the given vector of string all possible Tihls for the given tree size
        void InitRandomGenerator();                          //Function to Initialize the (Global) Random generator Object to generate random numbers
        void InitNumOfSteps();                               //Function to Initialize how many times the left-right subtree re-alignment procedure will take place
        void ResampleHistoryGivenHomologies();               //Function to Resample a new Homology Stucture for the given Subtree based on the initial homologies (of the initial alignment)
        void InitSubtree(int Subtree);                       //Function to Initialize the Subtree
        void InitSubTreeIdxs();                              //Function to Initialize the SubtreeIdxs
        void InitLengths();                                  //Function to Initialize the Lengths for the current Subtree
        void InitFPTables();                                 //Function to Initialize the Felsenstein Pruning tables for the current Subtree
        void InitAlignmentHomologies(IntMatRef Homologies);  //Function to Initialize the Alignment Homologies for the current Subtree
        void InitRegionHomologies();                         //Function to Initialize the Region Homologies for the current Subtree
        void InitStartAndEndIndices();                       //Function to Initialize the Start and End Soan Indices for the current Subtree
        void InitSubstitutionMatrices();                     //Function to Initialize the Substitution Matrices needed for the current Subtree
        void InitDimensions();                               //Function to Initialize the Dimensions needed for the current Subtree alignment region
        void MakeProperInitializationsAll(int subtree, IntMatRef Homologies); //Function to Call the above Initialization Functions in the Proper Order so as to ensure that all inter-dependencies between successive calls are met. The argument defines which subtree will be initialized

        //Soans, Tihls and SoansIdxs Getter Functions
        int       GetFirstIdxofSoanWhereRootIs(int rootIdx);        //Function that parses the Soan Indices Matrix to find the FIRST soan whose index corresponding to the root has a value equal to the one provided as argument. When such a soan is found, the corresponding index is returned or -1 otherwise
        int       GetLastIdxofSoanWhereRootIs(int rootIdx);         //Function that parses the Soan Indices Matrix to find the FIRST soan whose index corresponding to the root has a value equal to the one provided as argument. When such a soan is found, the corresponding index is returned or -1 otherwise
        string    GetTihlAtPosition(int pos, bool includeRoot);     //Function that returns the Tihl lying at the given position within the Tihls vector. The boolean argument represents whether the tihl corresponding to the root will also be included at the returned result or not
        string    GetSoanAtPosition(int pos, bool includeRoot);     //Function that returns the Soan lying at the given position within the Soans vector. The boolean argument represents whether the tihl corresponding to the root will also be included at the returned result or not
        IntVecRef GetSoanIdxsAtPosition(int pos);                   //Function that returns the Soan Indices at the given position within the SoansIdxs matrix. The boolean argument represents whether the tihl corresponding to the root will also be included at the returned result or not
        int       GetRootMaxValue();                                //Function that returns the maximum value of the soan index corresponding to the root of the current subtree


        //Functions to output results
        void PrintSequences();                   //Function to Output the Formed Sequence Objects to the standard output stream
        void PrintParameters();                  //Function to Output the Program Parameters to the standard output stream
        void PrintOverallSubstitutionMatrices(); //Function to Output the Estimated Substitution Matrices to the standard output stream
        void PrintSubstitutionMatrices();        //Function to Output the Substitution Matrices of interest to the standard output stream
        void PrintTree();                        //Function to Output the Estimated Felsenstein Tree Structure of the Sequences to the standard output stream
        void SaveTreeToFile();                   //Function to Save the Estimated Felsenstein Tree Structure of the Sequences to an Output File
        void PrintSubtree();                     //Function to Print which subtree (left or right) is currently being examined
        void PrintSubtreeIdxs();                 //Function to Print the SubTreeidxs of the current subtree
        void PrintLengths();                     //Function to Print the subtree lengths
        void PrintFPTables();                    //Function to Print the Felsenstein Pruning Tables for the Current Fragment to the standard output stream
        void PrintDimensions();                  //Function to Print the Dimensions for the current fragment to the standard output stream
        void PrintAlignmentHomologies();         //Function to Print the AlignmentHomologies to the Standard Output Stream
        void PrintRegionHomologies();            //Function to Print the RegionHomologies to the Standard Output Stream
        void PrintSampledHomologies();           //Function to Print the SampledHomologies to the Standard Output Stream
        void PrintHomologies();                  //Function to Print the Homologies (i.e. SampledHomologies) to the Standard Output Stream
        void SaveAlignmentHomologiesToFile();    //Function to Print the AlignmentHomologies to the Standard Output Stream
        void SaveRegionHomologiesToFile();       //Function to Print the RegionHomologies to an output file
        void SaveSampledHomologiesToFile();      //Function to Print the SampledHomologies to an output file
        void SaveHomologiesToFile();             //Function to Print the Homologies (i.e. SampledHomologies) to an output file
        void PrintEvolHistory();                 //Function to Print the Sampled Evolution History to the Standard Output Stream
        void PrintEvolHistoryReport();           //Function to Print the Sampled Soans, Soan Indices and Tihls to the Standard Output Stream


    protected:

    private:
        /** Private Members: **/
        Sequences*            Seqs;                    //The Input Sequences and Alignment
        Parameters*           Params;                  //The program Parameters
        Tree*                 FTree;                   //The Felsenstein Tree structure of the Sequences
        F84Model*             f84;                     //The F84 Model
        TN93Model*            tn93;                    //The TN93 Model
        GTRModel*             gtr;                     //The GTR Model
        DoubleMat             OverallSubMats;          //A Vector of Vectors representing the substitution matrices for each branch length
        DoubleVec             OverallSubMatsLengths;   //A Vector representing the lengths for each substitution matrix
        StringVec             Soans3Leaved;            //A vector containing all possible soans for a 3 leaved tree
        StringVec             Soans4Leaved;            //A vector containing all possible soans for a 4 leaved tree
        BoostRandomGenerator* RandGen;                 //The (global) RandomGenerator Object to be passed to the Aligner Objects
        int                   NumOfSteps;              //The number of times that the left-right subtree re-alignment procedure will be executed
        IntMat                SampledHomologies;       //A Matrix to store the complete matrix of the homologies between the given sequences (and their inner nodes) as defined by the "Resampling considering the Given Homologies Step)
        IntMat                RegionHomologies;        //A Matrix to store either the complete matrix of the homologies for this region (i.e. igonring all places where the root is -1)
        IntMat                AlignmentHomologies;     //Like Homologies, but straight as calculated from the Given Sequences Alignment (i.e. without skipping the "All-1" Columns)

        IntVec                SubTreeIdxs;             //The Subtree indices depending on whether the left or thr=e right subtree is being re-aligned
        int                   SubTree;                 //An integer value representing the subtree to be re-aligned. 0 stands for the right subtree whereas 1 stands for the left subtree
        DoubleMat             SubMats;                 //A Vector of Vectors representing the substitution matrices for each branch length
        DoubleVec             Lengths;                 //The Brach lengths for the edges of the current alignment region
        DoubleHyp             FPTables;                //The Felsenstein Pruning Tables for this alignment region
        IntVec                Dimensions;              //The Dimensions for this alignment region
        StringVec             EvolHistory;             //The Deriving Evolution History after resampling considering the given homologies
        StringVec             Soans;                   //The Deriving sequence of soans for the sampled Evolution History
        StringVec             Tihls;                   //The Deriving sequence of tihls for the sampled Evolution History
        IntMat                SoansIdxs;               //The Deriving seqeunce of Soans expressed by their corresponding indices
        IntVec                StartIdxs;               //The Start Indices for the current subtree
        IntVec                EndIdxs;                 //The End Indices for the current subtree


};

#endif // TRIPLEINITIALIZER_H
