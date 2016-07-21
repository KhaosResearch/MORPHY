#include "Sequences.h"

//Default Constructor
Sequences::Sequences()  { }

//Constructor to init all vector with default values according to the given size
Sequences::Sequences(int size)  {
    BareSequences.resize(size, "");
    AlignSequences.resize(size, "");
    SequenceIds.resize(size, "");
}

//Virtual Decontructor
Sequences::~Sequences() {
    vector<string>().swap(BareSequences);
    vector<string>().swap(AlignSequences);
    vector<string>().swap(SequenceIds);
}

//Function that converts the Alignment into a matrix of integer values representing the correspondences between the sequences' sites
void Sequences::SetIndicesFromAlignment() {
    IndicesFromAlignment.resize(6, IntVec(AlignSequences[0].size(),-1));
    IntVec counter(6, 1);

    for (int j=0;j<AlignSequences[0].size();j++){
        for (int i=0;i<AlignSequences.size();i++){
            if (AlignSequences[i][j]!='-'){
                IndicesFromAlignment[i][j] = counter[i]++;
            }
        }

        if ( (AlignSequences[0][j]!='-' && AlignSequences[1][j]!='-')  ||
             (AlignSequences[0][j]!='-' && (AlignSequences[2][j]!='-' || AlignSequences[3][j]!='-')) ||
             (AlignSequences[1][j]!='-' && (AlignSequences[2][j]!='-' || AlignSequences[3][j]!='-'))  ){
            IndicesFromAlignment[4][j] = counter[4]++;
        }

        if ( (AlignSequences[2][j]!='-' && AlignSequences[3][j]!='-')  ||
             (AlignSequences[2][j]!='-' && (AlignSequences[0][j]!='-' || AlignSequences[1][j]!='-')) ||
             (AlignSequences[3][j]!='-' && (AlignSequences[0][j]!='-' || AlignSequences[1][j]!='-'))  ){
            IndicesFromAlignment[5][j] = counter[5]++;
        }
    }
}

//Definition of the toString Operation
ostream& operator<<(ostream& ostr, const Sequences& seqs) {
    ostr<<"Sequence Headers:\n-----------------\n";
    for (int i=0;i<seqs.Get_SeqsIds().size();i++){
        ostr << seqs.GetSequenceIdAt(i) <<endl;
    }

    ostr<<endl;
    ostr<<"Unaligned Sequences:\n--------------------\n";
    for (int i=0;i<seqs.Get_BareSeqs().size();i++){
        ostr << seqs.GetBareSequenceAt(i) <<endl;
    }

    ostr<<endl;

    ostr<<"Alignment:\n----------\n";
    for (int i=0;i<seqs.Get_AlignSeqs().size();i++){
        ostr << seqs.GetAlignedSequenceAt(i) <<endl;
    }

    ostr<<endl;

    return ostr;
}
