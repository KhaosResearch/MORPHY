#include "Triplet.h"

//Default (empty Constructor)
Triplet::Triplet()  { }

//Virtual Decontructor
Triplet::~Triplet() { }

//Copy Constructor to Clone objects
Triplet::Triplet(const Triplet& other){
    SoanBefore  = other.Get_SoanBefore();
    Tihl        = other.Get_Tihl();
    SoanAfter   = other.Get_SoanAfter();
    Length      = other.Get_Length();
    Probability = other.Get_Probability();
}

//Definition of The equality (==) operator
bool Triplet::operator==(const Triplet &other) {
    if ( this->Get_SoanBefore().compare(other.Get_SoanBefore())==0  &&
         this->Get_Tihl().compare(other.Get_Tihl())==0              &&
         this->Get_SoanAfter().compare(other.Get_SoanAfter())==0    &&
         this->Get_Length() == other.Get_Length()                   &&
         this->Get_Probability() == other.Get_Probability()
        )
         return true;
    else
        return false;
}

//Definition of The inequality (!=) operator
bool Triplet::operator!=(const Triplet &other) {
    return !( *this == other);
}

//Function to Calculate the Transition Probabiltiy of this Triplet based on the provided
//fid model and vector of branch lengths
void Triplet::Set_Probability(FIDModel* fid, const vector<double>& lengths) {
    double trans_prob = 1;
    int Pos_Of_B_Node = FirstIndexOf(Tihl, "B");

    for (unsigned int k=0;k<Tihl.size();k++){
        if (k<Pos_Of_B_Node && Pos_Of_B_Node!=-1) {
            if      (SoanBefore[k]=='H'){ trans_prob *= (1 - fid->GetProbability(BB, _B, lengths[k])); }
            else if (SoanBefore[k]=='B'){ trans_prob *= (1 - fid->GetProbability(_B, _B, lengths[k])); }
        }

        if      (Tihl[k]=='H'){
            if      (SoanBefore[k]=='H') { trans_prob *= fid->GetProbability(BB, BB, lengths[k]);  }
            else if (SoanBefore[k]=='B') { trans_prob *= fid->GetProbability(_B, BB, lengths[k]);  }
            else if (SoanBefore[k]=='h') { trans_prob *= ( fid->GetProbability(BB, BB, lengths[k]) / (1-fid->GetProbability(BB, _B, lengths[k])) );  }
            else if (SoanBefore[k]=='b') { trans_prob *= ( fid->GetProbability(_B, BB, lengths[k]) / (1-fid->GetProbability(_B, _B, lengths[k])) );  }
            else if (SoanBefore[k]=='e') { trans_prob *= ( fid->GetProbability(B_, BB, lengths[k]) / (1-fid->GetProbability(B_, _B, lengths[k])) );  }
        }else if (Tihl[k]=='B'){
            if      (SoanBefore[k]=='H') { trans_prob *= fid->GetProbability(BB, _B, lengths[k]);  }
            else if (SoanBefore[k]=='B') { trans_prob *= fid->GetProbability(_B, _B, lengths[k]);  }
            else                         { trans_prob *= 0;                                                      }
        }else if (Tihl[k]=='N'){
            if      (SoanBefore[k]=='H') { trans_prob *= ( fid->GetProbability(BB, B_, lengths[k]) * fid->GetProbability(B_, _B, lengths[k])  );  }
            else if (SoanBefore[k]=='B') { trans_prob *= ( fid->GetProbability(_B, B_, lengths[k]) * fid->GetProbability(B_, _B, lengths[k])  );  }
            else if (SoanBefore[k]=='h') { trans_prob *= ( fid->GetProbability(BB, B_, lengths[k]) * fid->GetProbability(B_, _B, lengths[k])  ) / (1-fid->GetProbability(BB, _B, lengths[k]));  }
            else if (SoanBefore[k]=='b') { trans_prob *= ( fid->GetProbability(_B, B_, lengths[k]) * fid->GetProbability(B_, _B, lengths[k])  ) / (1-fid->GetProbability(_B, _B, lengths[k]));  }
            else if (SoanBefore[k]=='e') { trans_prob *= ( fid->GetProbability(B_, B_, lengths[k]) * fid->GetProbability(B_, _B, lengths[k])  ) / (1-fid->GetProbability(B_, _B, lengths[k]));  }
        }else if (Tihl[k]=='E'){
            if      (SoanBefore[k]=='H') { trans_prob *= ( fid->GetProbability(BB, B_, lengths[k]) * (1-fid->GetProbability(B_, _B, lengths[k]))  );  }
            else if (SoanBefore[k]=='B') { trans_prob *= ( fid->GetProbability(_B, B_, lengths[k]) * (1-fid->GetProbability(B_, _B, lengths[k]))  );  }
            else if (SoanBefore[k]=='h') { trans_prob *= ( fid->GetProbability(BB, B_, lengths[k]) * (1-fid->GetProbability(B_, _B, lengths[k]))  ) / (1-fid->GetProbability(BB, _B, lengths[k]));  }
            else if (SoanBefore[k]=='b') { trans_prob *= ( fid->GetProbability(_B, B_, lengths[k]) * (1-fid->GetProbability(B_, _B, lengths[k]))  ) / (1-fid->GetProbability(_B, _B, lengths[k]));  }
            else if (SoanBefore[k]=='e') { trans_prob *= ( fid->GetProbability(B_, B_, lengths[k]));  }
        }
    }

    Probability = trans_prob;
}

//Function to Calculate the Transition Probabiltiy of this Triplet based on the provided fid model and vector of branch lengths for pairwise alignment
void Triplet::Set_Probability (FIDModel* fid, const vector<double>& lengths, bool pairwise) {
    double trans_prob = 1;
    unsigned int k;

    if (Tihl.size()==1) {
        k = 0;
    }else{
        k = 1;
    }
    for (k;k<Tihl.size();k++){
        if      (Tihl[k]=='H'){
            if      (SoanBefore[k]=='H') { trans_prob *= fid->GetProbability(BB, BB, lengths[k]);  }
            else if (SoanBefore[k]=='B') { trans_prob *= fid->GetProbability(_B, BB, lengths[k]);  }
            else if (SoanBefore[k]=='e') { trans_prob *= ( fid->GetProbability(B_, BB, lengths[k]) / (1-fid->GetProbability(B_, _B, lengths[k])) );  }
        }else if (Tihl[k]=='B'){
            if      (SoanBefore[k]=='H') { trans_prob *= fid->GetProbability(BB, _B, lengths[k]);  }
            else if (SoanBefore[k]=='B') { trans_prob *= fid->GetProbability(_B, _B, lengths[k]);  }
            else                         { trans_prob *= 0;                                                      }
        }else if (Tihl[k]=='N'){
            if      (SoanBefore[k]=='H') { trans_prob *= ( fid->GetProbability(BB, B_, lengths[k]) * fid->GetProbability(B_, _B, lengths[k])  );  }
            else if (SoanBefore[k]=='B') { trans_prob *= ( fid->GetProbability(_B, B_, lengths[k]) * fid->GetProbability(B_, _B, lengths[k])  );  }
            else if (SoanBefore[k]=='e') { trans_prob *= ( ( fid->GetProbability(B_, B_, lengths[k]) * fid->GetProbability(B_, _B, lengths[k])  ) / (1-fid->GetProbability(B_, _B, lengths[k])) );  }
        }else if (Tihl[k]=='E'){
            if      (SoanBefore[k]=='H') { trans_prob *= ( fid->GetProbability(BB, B_, lengths[k]) * (1-fid->GetProbability(B_, _B, lengths[k]))  );  }
            else if (SoanBefore[k]=='B') { trans_prob *= ( fid->GetProbability(_B, B_, lengths[k]) * (1-fid->GetProbability(B_, _B, lengths[k]))  );  }
            else if (SoanBefore[k]=='e') { trans_prob *= ( fid->GetProbability(B_, B_, lengths[k]) );  }
        }
    }

    Probability = trans_prob;
}

//Function that returns a vector of size equal to the labeling.substr(1) of the SoanAfter soan, having +1 if an
//H or B has come up, or 0 if an h, b, or e has come up
const vector<int> Triplet::GetTripletIdxChanges() const {
    vector<int> idxs(Tihl.size()+1, 0);

    for (unsigned int i=0;i<Tihl.size();i++){
        if (Tihl[i]=='H' || Tihl[i]=='B' || Tihl[i]=='N') idxs[i] = 1;
    }

    return idxs;
}

//Function that returns the number of columns that are required in order to convert this evolution history into
//an alignment (e.g. B and N require one column each, whereas H and E can coexist in one column)
int Triplet::GetTripletNumberOfAlignmentColumns() const {
    int colCount = 1; //minimum 1 column is required
    for (unsigned int i=0;i<Tihl.size();i++){
        if (Tihl[i]=='B' || Tihl[i]=='N') colCount++;
    }
    return colCount;
}

//Definition of the toString Operation
ostream& operator<<(ostream& ostr, const Triplet& trip) {
    ostr << trip.Get_SoanBefore() << " -> ";
    ostr << trip.Get_Tihl()       << " -> ";
    ostr << trip.Get_SoanAfter();

    return ostr;
}
