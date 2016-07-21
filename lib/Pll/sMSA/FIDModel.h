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

/**************************************************************************************/
/**                                   Class Scopes:                                   */
/**          Aim: This class aims at defining the properties needed to represent      */
/**                                   the FID Model                                   */
/**************************************************************************************/

#ifndef FIDMODEL_H
#define FIDMODEL_H

/** Header file Inclusions **/
#include <vector>
#include <stdlib.h>
#include <cmath>

#include "Mappings.h"

using namespace std;
using namespace Mappings;

/** Typedef Definitions **/
typedef unsigned int uint;
typedef vector< vector<double> > RealMatrix;

class FIDModel {
    public:
        /** Constructors, Deconstructors and Operators Definitions **/
        FIDModel();                                 //Default (empty Constructor)
        FIDModel(uint GammaVal, double LambdaVal);  //Constructor that builds a FID object based on the given Lambda and Gamma parameters
        virtual ~FIDModel();                        //Virtual Decontructor
        FIDModel(const FIDModel& other);            //Copy Constructor to Clone objects
        bool operator== (const FIDModel&);          //Definition of The equality (==) operator
        bool operator!= (const FIDModel&);          //Definition of The inequality (!=) operator


        /** Setter And Getter Methods for Private Members **/
        void Set_Gamma (uint val) { Gamma = val;  }
        void Set_Lambda(uint val) { Lambda = val; }

        uint   Get_Gamma()  const { return Gamma;  }
        double Get_Lambda() const { return Lambda; }

        /** Extra Functions: **/
        void EvaluateTransitionsMatrix(double length);                                    //Function that evaluates all the transition probabilities for a given branch length
        double GetProbability(EventType1 prevEvent, EventType1 nextEvent, double length); //Function that returns the transition probaility of being at state prevEvent and transitioning to state nextEvent at time length
        double GetProbability(EventType2 prevEvent, EventType2 nextEvent, double length); //Function that returns the transition probaility of being at state prevEvent and transitioning to state nextEvent at time length
        const RealMatrix& Get_TransitionsMatrix();                                        //Function that returns the already estimated transition probabilities in a matrix form
        const RealMatrix& Get_TransitionsMatrix(double length);                           //Function that estimates the transition probabilities for a given length and returns the already estimated transition probabilities in a matrix form
    protected:

    private:
        /** Class Members: **/
        uint   Gamma;                       //The expected fragment length (>=1)
        double Lambda;                      //The indel rate per site      (>=0)
        double ProbabilitiesMatrix[3][3];   //A 3x3 Probabilities Matrix to store the transition probabilities between the different states
};

/** Definition of the toString Operation **/
ostream& operator<<(ostream& ostr, const FIDModel& model);

#endif // FIDMODEL_H
