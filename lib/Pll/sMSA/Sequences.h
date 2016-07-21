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

/************************************************************************/
/**                         File Scopes:                                */
/**    Class the represent the input Sequences and Alignments           */
/************************************************************************/

#ifndef SEQUENCES_H
#define SEQUENCES_H

/** Header file Inclusions **/
#include <iomanip>

#include "TypeDefinitions.h"
#include "Mappings.h"
#include "Utils.h"

using namespace TypeDefinitions;
using namespace Mappings;
using namespace Utils;

class Sequences {

    public:
        /** Constructors, Deconstructors and Operators Definitions **/
        Sequences();                                   //Default Constructor
        Sequences(int size);                           //Constructor to init all vector with default values according to the given size
        virtual ~Sequences();                          //Virtual Decontructor

        /** Setter And Getter Methods for Private Members **/
        void Set_BareSeqs (StringVecRef val) { BareSequences  = val; }
        void Set_Alignment(StringVecRef val) { AlignSequences = val; }
        void Set_SeqsIds  (StringVecRef val) { SequenceIds    = val; }

        StringCVecRef Get_BareSeqs () const { return BareSequences;  }
        StringCVecRef Get_AlignSeqs() const { return AlignSequences; }
        StringCVecRef Get_SeqsIds  () const { return SequenceIds;    }
        IntCMatRef    Get_IndicesFromAlignment() const { return IndicesFromAlignment; }

        /** Extra Functions: **/
        void AddBareSequence   (string val) { BareSequences.push_back(val);  }
        void AddAlignedSequence(string val) { AlignSequences.push_back(val); }
        void AddSequenceId     (string val) { SequenceIds.push_back(val);    }

        void SetBareSequenceAt   (string val, int pos) { BareSequences[pos]  = (val);  }
        void SetAlignedSequenceAt(string val, int pos) { AlignSequences[pos] = (val);  }
        void SetSequenceIdAt     (string val, int pos) { SequenceIds[pos]    = (val);  }

        string GetBareSequenceAt   (int pos) const { return BareSequences[pos];  }
        string GetAlignedSequenceAt(int pos) const { return AlignSequences[pos]; }
        string GetSequenceIdAt     (int pos) const { return SequenceIds[pos];    }

        void   SetIndicesFromAlignment();   //Function that converts the Alignemnt into a matrix of integer values representing the correspondences between the sequences' sites

        //string GetBareSequenceWithId   (string SeqId) const { int pos = GetPosOfElement(SequenceIds, SeqId);  return BareSequences[pos];  }
        //string GetAlignedSequenceWithId(string SeqId) const { int pos = GetPosOfElement(SequenceIds, SeqId);  return AlignSequences[pos]; }
        //string GetIdOfBareSequence     (string Seq  ) const { int pos = GetPosOfElement(BareSequences, Seq);  return SequenceIds[pos];    }
        //string GetIdOfAlignedSequence  (string Seq  ) const { int pos = GetPosOfElement(AlignSequences, Seq); return SequenceIds[pos];    }


    protected:

    private:
         /** Class Members **/
         StringVec  BareSequences;          //A Vector of strings representing the Unaligned sequences at the tips of the tree
         StringVec  AlignSequences;         //A Vector of strings representing the Aligned sequences at the tips of the tree
         StringVec  SequenceIds;            //A Vector of strings representing the ids sequences at the tips of the tree
         IntMat     IndicesFromAlignment;   //A Matrix of indices representing the Homologies between the sequences based on their alignment

};

/** Definition of the toString Operation **/
ostream& operator<<(ostream& ostr, const Sequences& seqs);

#endif // SEQUENCES_H
