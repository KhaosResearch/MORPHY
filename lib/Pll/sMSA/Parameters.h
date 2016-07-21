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
/**  Aim: This Class provides with set of useful functions that help initialize    */
/**  all input parameters either from the command prompt or from an input file     */
/***********************************************************************************/

#ifndef PARAMETERS_H
#define PARAMETERS_H

/** Header file Inclusions **/
#include <string>
#include <iostream>
#include <fstream>
#include <stdlib.h>

#include "Mappings.h"

using namespace std;
using namespace Mappings;

class Parameters {
    public:
        /** Constructors, Deconstructors and Operators Definitions **/
        Parameters();                                   //Default (empty Constructor)
        virtual ~Parameters();                          //Virtual Decontructor
        Parameters(const Parameters& other);            //Copy Constructor to Clone objects
        bool operator== (const Parameters&);            //Definition of The equality (==) operator
        bool operator!= (const Parameters&);            //Definition of The inequality (!=) operator

        /** Setter And Getter Methods for Private Members **/
        void Set_BaseFreq       (double val[], size_t size) { for (unsigned int i=0;i<size;i++) { BaseFreq[i] = val[i]; } }
        void Set_Lambda         (double val)                { Lambda          = val; }
        void Set_Gamma          (double val)                { Gamma           = val; }
        void Set_SampleSize     (int val   )                { SampleSize      = val; }
        void Set_AlignFragSize  (int val   )                { AlignFragSize   = val; }
        void Set_AlignmentOutput(string val)                { AlignmentOutput = val; }
        void Set_SubModel       (string val)                { SubModel        = val; }
        void Set_Kappa          (double val)                { k               = val; }
        void Set_AG_Rate        (double val)                { k1              = val; }
        void Set_CT_Rate        (double val)                { k2              = val; }
        void Set_X1_Param       (double val)                { x1              = val; }
        void Set_X2_Param       (double val)                { x2              = val; }
        void Set_X3_Param       (double val)                { x3              = val; }
        void Set_X4_Param       (double val)                { x4              = val; }
        void Set_X5_Param       (double val)                { x5              = val; }
        void Set_X6_Param       (double val)                { x6              = val; }

        const double* Get_BaseFreq()        const { return BaseFreq;        }
        double        Get_Lambda()          const { return Lambda;          }
        double        Get_Gamma()           const { return Gamma;           }
        int           Get_SampleSize()      const { return SampleSize;      }
        int           Get_AlignFragSize()   const { return AlignFragSize;   }
        string        Get_AlignmentOutput() const { return AlignmentOutput; }
        string        Get_SubModel()        const { return SubModel;        }
        double        Get_Kappa()           const { return k;               }
        double        Get_AG_Rate()         const { return k1;              }
        double        Get_CT_Rate()         const { return k2;              }
        double        Get_X1_Param()        const { return x1;              }
        double        Get_X2_Param()        const { return x2;              }
        double        Get_X3_Param()        const { return x3;              }
        double        Get_X4_Param()        const { return x4;              }
        double        Get_X5_Param()        const { return x5;              }
        double        Get_X6_Param()        const { return x6;              }

        /** Extra Functions: **/
        void ReadDataFromFile(string filename);  //Function to Read Parameters from an input file
        void ReadDataFromCmd();                  //Function to Read Parameters from the CMD
        bool HaveBeenInitialized();              //Function to Check wheter the Parameters have been initiliazed

    protected:

    private:
        /** Class Members **/
        double BaseFreq[4];         //The prior Base Frequencies
        double Lambda;              //The indel rate per site (should be >=0)
        double Gamma;               //The expected fragment length (should be >=1)
        int    SampleSize;          //The Desired Sample size when sampling sequences
        int    AlignFragSize;       //The Alignment Fragment size when sampling sequences
        string AlignmentOutput;     //The Desired Alignment type that should be used for the output
        string SubModel;            //The Desired Substitution Model to be employed
        double k;                   //the kappa parameter that distinguishes between the rate of transitions and transversions
        double k1;                  //The A<->G rate
        double k2;                  //The C<->T rate
        double x1;                  //The x1 substitution parameter
        double x2;                  //The x2 substitution parameter
        double x3;                  //The x3 substitution parameter
        double x4;                  //The x4 substitution parameter
        double x5;                  //The x5 substitution parameter
        double x6;                  //The x6 substitution parameter
};

/** Definition of the toString Operation **/
ostream& operator<<(ostream& ostr, const Parameters& p);

#endif // PARAMETERS_H
