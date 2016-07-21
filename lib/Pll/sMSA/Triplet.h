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
/**                                 Class Scopes:                                  */
/**   Aim: Class defining the properties needed to represent a Triplet object      */
/**          (i.e. an object representing a soan->tihl->soan transition)           */
/***********************************************************************************/

#ifndef TRIPLET_H
#define TRIPLET_H

/** Header file Inclusions **/
#include <iostream>
#include <vector>
#include <string>

#include "FIDModel.h"
#include "Mappings.h"
#include "Utils.h"

using namespace std;
using namespace Mappings;
using namespace Utils;

class Triplet {
    public:
        /** Constructors, Deconstructors and Operators Definitions **/
        Triplet();                                                           //Default (empty Constructor)
        virtual ~Triplet();                                                  //Virtual Decontructor
        Triplet(const Triplet& other);                                       //Copy Constructor to Clone objects
        bool operator== (const Triplet&);                                    //Definition of The equality    (==) operator
        bool operator!= (const Triplet&);                                    //Definition of The inequality  (!=) operator


        /** Setter And Getter Methods for Private Members **/
        void Set_SoanBefore  (string val)  { SoanBefore  = val; }
        void Set_Tihl        (string val)  { Tihl        = val; }
        void Set_SoanAfter   (string val)  { SoanAfter   = val; }
        void Set_Length      (double    val)  { Length      = val; }
        void Set_Probability (double    val)  { Probability = val; }
        void Set_Probability (FIDModel* fid, const vector<double>& lengths); //Function to Calculate the Transition Probabiltiy of this Triplet based on the provided fid model and vector of branch lengths
        void Set_Probability (FIDModel* fid, const vector<double>& lengths, bool pairwise); //Function to Calculate the Transition Probabiltiy of this Triplet based on the provided fid model and vector of branch lengths for pairwise alignment
        void Set_Valid       (bool      val)  { Valid    = val; }

        string Get_SoanBefore()  const { return SoanBefore;  }
        string Get_Tihl()        const { return Tihl;        }
        string Get_SoanAfter()   const { return SoanAfter;   }
        double Get_Length()      const { return Length;      }
        double Get_Probability() const { return Probability; }
        bool   Is_Valid()        const { return Valid;       }

        /** Extra Functions: **/
        const vector<int> GetTripletIdxChanges()  const;      //Function that returns a vector having +1 for every H, B or N event, or 0 for every E, or - event
        int GetTripletNumberOfAlignmentColumns()  const;      //Function that returns the number of columns that are required in order to convert this evolution history into an alignment (e.g. B and N require one column each, whereas H and E can coexist in one column)
    protected:

    private:
        /** Class Members: **/
        string SoanBefore;   //To represent the state before applying the tihl
        string Tihl;         //The tihl to be applied to current soan (SoanBefore)
        string SoanAfter;    //The resulting soan after the tihl has been applied to cuurent state
        double Length;       //The branch length - used to estimate the transition probability
        double Probability;  //The transition probability of this event to happen (int time = length)
        bool   Valid;        //A boolean variable signifying if this transition is valid or not
};

/** Definition of the toString Operation **/
ostream& operator<<(ostream& ostr, const Triplet& trip);

#endif // TRIPLET_H
