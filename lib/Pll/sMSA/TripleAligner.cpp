#include "TripleAligner.h"

//Default Constructor
TripleAligner::TripleAligner()  {
    SaveToLogFile = false;
}

//Virtual Decontructor
TripleAligner::~TripleAligner() {
    /*
    delete Params;
    delete Fid;
    for (int i=0;i<Triplets.size();i++){
        for (int j = 0;j<Triplets[i].size();j++) {
            delete Triplets[i][j];
        }
    }

    Triplets.clear();
    */

    for( map<LongInt, ProbabilityObject* >::iterator ii=Probabilities.begin(); ii!=Probabilities.end(); ++ii) {
        delete (*ii).second;
    }
    Dimensions.clear();
    //Lengths.clear();
    //delete f84;
    //delete tn93;
    //delete gtr;
    FPTables.clear();
    //SubMats.clear();
    HomologiesIdxs.clear();
    Anchors.clear();
}

/** Starter Functions **/ //{

//}

/** Initialization Functions **/ //{
//Function to Initialize all Program's Parameters
void TripleAligner::InitParameters(Parameters* val) {
    Params = val;
}

//Function to Initialize the Substitution Model
void TripleAligner::InitSubModel(F84Model* val1, TN93Model* val2, GTRModel* val3) {
    f84  = val1;
    tn93 = val2;
    gtr  = val3;
}

//Function to Initialize the FID Model
void TripleAligner::InitFidModel() {
    Fid =  new FIDModel(Params->Get_Gamma(), Params->Get_Lambda());
}

//Function to Initialize the homologies structure from the Given HomologyTable and  Indices to be used
void TripleAligner::InitHomologyIndices(IntCMatRef Homologies) {
    HomologiesIdxs = Homologies;

    int rootLength = FPTables[0][0].size();
    Anchors.resize(HomologiesIdxs[0].size(), vector<int>(HomologiesIdxs.size(), -1));
    for (int i=0; i<HomologiesIdxs.size(); i++) {
        for (int j=0; j<HomologiesIdxs[0].size(); j++) {
            Anchors[j][i] = HomologiesIdxs[i][j];
        }
    }

    /*For Debugging*/
    /*cout<<"\n\nHomologies:\n-----------\n";
    for (int i=0;i<HomologiesIdxs.size();i++){
        for (int j=0;j<HomologiesIdxs[i].size();j++){
            cout<<right<<setw(4)<<HomologiesIdxs[i][j];
        }
        cout<<endl;
    }

    cout<<"\n\nAnchors:\n--------\n";
    for (int i=0;i<Anchors.size();i++){
        for (int j=0;j<Anchors[i].size();j++){
            cout<<right<<setw(4)<<Anchors[i][j];
        }
        cout<<endl;
    }*/
}

//Function to Initialize the corresponding branch lengths
void TripleAligner::InitLengths(DoubleCVecRef val) {
    Lengths = val;
}

//Function to Initialize the Dimensions (equal to the sequence sizes)
void TripleAligner::InitDimensions(IntCVecRef val) {
    Dimensions = val;
}

//Function to Initialize the Felsenstein Pruning Tables of the Tree according to the given values
void TripleAligner::InitFPTables(DoubleCHypRef val) {
    FPTables = val;
}

//Function to Initialize all Triplets
void TripleAligner::InitTriplets() {
    /** Create all Possible Tihls **/
    vector<string> tihlLbls;
    vector<string> Tihls;
    if (Lengths.size()>1) {
        tihlLbls.push_back("H");
        tihlLbls.push_back("B");
        tihlLbls.push_back("N");
        tihlLbls.push_back("E");
        tihlLbls.push_back("-");
        GetAllCombinations("", Lengths.size(), tihlLbls.size(), tihlLbls, "TihlTree", true, Tihls);
    }else {
        Tihls.push_back("H");
        Tihls.push_back("B");
        Tihls.push_back("N");
        Tihls.push_back("E");
    }

    /** Create all Reachable Soans **/
    vector<string> Soans;
    GetAllReachableSoans(Soans, Tihls);

    /*To Store Reachable Soans in Files*/
    /*string filename;
    if (Lengths.size()==2) {
        filename = ReachableSoans2LeavedFile;
    }else if (Lengths.size()==3) {
        filename = ReachableSoans3LeavedFile;
    }else if (Lengths.size()==4){
        filename = ReachableSoans4LeavedFile;
    }
    WriteToFile(filename, "", false, true);
    stringstream ss; ss.str("");
    for (int i=0;i<Soans.size();i++) {
        ss<<Soans[i]<<endl;
    }
    WriteToFile(filename, ss.str(), false, false);*/

    /** Create all Possible Triplets organized per each reachable soan **/
    GetAllPossibleTriplets(Soans, Tihls);

    //For Debugging
    /*int counter = 0;
    double sum = 0;
    for (int i=0;i<Triplets.size();i++){
        for (int j=0;j<Triplets[i].size();j++){
            cout<< *(Triplets[i][j])<<endl;
            counter++;
            sum+=Triplets[i][j]->Get_Probability();
        }
    }
    cout<<"Total: " <<counter <<" triplets."<<endl;
    cout<<"Sum: " <<sum<<endl;
    Pause();*/
}

//Function to intialize all ProbabilityObject with default values
void TripleAligner::InitProbabilities(string startSoanLbl, string startTihlLbl, IntVecRef startIdxs) {
    EndPosIn_1D = 1;
    Dimensions.push_back(Triplets.size());

    for (unsigned int i=0; i<Dimensions.size(); i++) {
        EndPosIn_1D *= (Dimensions[i]);
    }

    ProbabilityObject* prob = new ProbabilityObject();
    prob->Set_OverallProbability(log(1.0));
    prob->Set_SoanLbl(startSoanLbl);
    StartPosIn_ND = startIdxs;
    StartTihlLbl  = startTihlLbl;
    prob->Set_PosIn_ND(StartPosIn_ND);
    StartPosIn_1D = GetPosForIndexes(startIdxs, Dimensions);
    prob->Set_PosIn_1D(StartPosIn_1D);
    prob->Set_EndState(false);
    prob->AddContributor(new ProbabilityObject());
    prob->AddContributorProbability(log(0.0));
    prob->AddContributorTihl(startTihlLbl);
    prob->AddContributorTihlExt(startTihlLbl);
    prob->AddContributorWingFold(false);
    Probabilities[StartPosIn_1D] = prob;

    EndPosIn_1D = GetPosForIndexes(EndPosIn_ND, Dimensions);
}

//Function to Initialize the Substitution Matrices for the given branch lengths
void TripleAligner::InitSubMatrices(DoubleCMatRef val) {
    SubMats = val;
}

//}

/** Assisting Functions for In-Between calculations **/ //{
//Function to Create All Reachable Soans, based on all possible tihls and an initial (all H-labeled) soan
void TripleAligner::GetAllReachableSoans(vector<string>& soans, const vector<string>& tihls) {
    stringstream ss;
    ss.str("");
    for (unsigned int i=0; i<Dimensions.size(); i++) {
        ss<<"H";
    }

    queue<string> soansQueue;
    soansQueue.push(ss.str());
    while (!soansQueue.empty()) {
        string topSoan = soansQueue.front();
        soansQueue.pop();
        //Apply every tihl to current soan
        for (unsigned int i=0; i < tihls.size(); i++) {
            if (!CanTihlBeAppliedToSoan(topSoan, tihls[i])) continue;
            string resultingSoan;
            if (Lengths.size()<=2) {
                resultingSoan = ApplyTihlToSoan(topSoan, tihls[i], true);
            } else {
                resultingSoan = ApplyTihlToSoan(topSoan, tihls[i]);
            }

            //cout<<topSoan<<"\t"<<tihls[i]<<"\t"<<resultingSoan<<endl;
            if ( std::find(soans.begin(), soans.end(), resultingSoan) != soans.end()) {
                continue;
            }

            soans.push_back(resultingSoan);
            soansQueue.push(resultingSoan);
        }
    }
}

//Function to Create All Possible, based on all possible tihls and all reachable soans
void TripleAligner::GetAllPossibleTriplets(const vector<string>& soans, const vector<string>& tihls) {
    idx_R = -1;
    idx_C = -1;
    int pos;
    for (unsigned int i=0; i<soans.size(); i++) {
        vector<Triplet *> tripvec;
        for (unsigned int j=0; j<tihls.size(); j++) {
            if (!CanTihlBeAppliedToSoan(soans[i], tihls[j])) continue;
            Triplet* trip = new Triplet();
            trip->Set_SoanBefore(soans[i]);
            trip->Set_Tihl(tihls[j]);
            if (Lengths.size()<=2) {
                trip->Set_Probability(Fid, Lengths, true);
                trip->Set_SoanAfter(ApplyTihlToSoan(soans[i], tihls[j], true));
            } else {
                trip->Set_Probability(Fid, Lengths);
                trip->Set_SoanAfter(ApplyTihlToSoan(soans[i], tihls[j]));
            }

            //if (!CanTihlBeAppliedToSoan(soans[i], tihls[j])) {
            //    trip->Set_Valid(false);
            // } else {
            trip->Set_Valid(true);
            // }

            tripvec.push_back(trip);

            if (StringCount(tihls[j], "E")==tihls[j].size()) {
                pos = tripvec.size()-1;
            }
        }
        Triplets.push_back(tripvec);

        if (idx_R !=-1) continue;
        if (StringCount(soans[i], "e")==soans[i].size()) {
            idx_R  = i;
            idx_C  = pos;
        }
    }
}

//Function that examines whether a given tihl may be applied to a given soan or not
bool TripleAligner::CanTihlBeAppliedToSoan(string soan, string tihl) {
    for (unsigned int i=0; i<soan.size(); i++) {
        if ( (soan.at(i)=='e' || soan.at(i)=='h' || soan.at(i)=='b') && (tihl.at(i)=='B'))
            return false;
    }
    return true;
}

//Function that Applies a given Tihl to a given Soan
string TripleAligner::ApplyTihlToSoan(string soan, string tihl) {
    stringstream ss;
    ss.str("");
    for (unsigned int i=0; i<tihl.size(); i++) {
        switch (soan.at(i)) {
        case 'H':
            switch(tihl.at(i)) {
            case 'H':
                ss<<"H";
                break;
            case 'B':
                ss<<"B";
                break;
            case 'N':
                ss<<"B";
                break;
            case 'E':
                ss<<"e";
                break;
            case '-':
                if (i<FirstIndexOf(tihl, "B")) {
                    ss<<"h";
                } else                           {
                    ss<<"H";
                }
            }
            break;

        case 'B':
            switch(tihl.at(i)) {
            case 'H':
                ss<<"H";
                break;
            case 'B':
                ss<<"B";
                break;
            case 'N':
                ss<<"B";
                break;
            case 'E':
                ss<<"e";
                break;
            case '-':
                if (i<FirstIndexOf(tihl, "B")) {
                    ss<<"b";
                } else                           {
                    ss<<"B";
                }
            }
            break;

        case 'h':
            switch(tihl.at(i)) {
            case 'H':
                ss<<"H";
                break;
            case 'N':
                ss<<"B";
                break;
            case 'E':
                ss<<"e";
                break;
            case '-':
                ss<<"h";
                break;
            }
            break;

        case 'b':
            switch(tihl.at(i)) {
            case 'H':
                ss<<"H";
                break;
            case 'N':
                ss<<"B";
                break;
            case 'E':
                ss<<"e";
                break;
            case '-':
                ss<<"b";
                break;
            }
            break;

        case 'e':
            switch(tihl.at(i)) {
            case 'H':
                ss<<"H";
                break;
            case 'N':
                ss<<"B";
                break;
            case 'E':
                ss<<"e";
                break;
            case '-':
                ss<<"e";
                break;
            }
            break;

        }
    }
    return ss.str();
}

//Function that Applies a given Tihl to a given Soan (for pairwise alignment)
string TripleAligner::ApplyTihlToSoan(string soan, string tihl, bool PairWise) {
    stringstream ss;
    ss.str("");
    for (unsigned int i=0; i<tihl.size(); i++) {
        switch (soan.at(i)) {
        case 'H':
            switch(tihl.at(i)) {
            case 'H':
                ss<<"H";
                break;
            case 'B':
                ss<<"B";
                break;
            case 'N':
                ss<<"B";
                break;
            case 'E':
                ss<<"e";
                break;
            case '-':
                ss<<"-";
            }
            break;

        case 'B':
            switch(tihl.at(i)) {
            case 'H':
                ss<<"H";
                break;
            case 'B':
                ss<<"B";
                break;
            case 'N':
                ss<<"B";
                break;
            case 'E':
                ss<<"e";
                break;
            case '-':
                ss<<"-";
            }
            break;

        case 'h':
            switch(tihl.at(i)) {
            case 'H':
                ss<<"H";
                break;
            case 'N':
                ss<<"B";
                break;
            case 'E':
                ss<<"e";
                break;
            case '-':
                ss<<"-";
                break;
            }
            break;

        case 'b':
            switch(tihl.at(i)) {
            case 'H':
                ss<<"H";
                break;
            case 'N':
                ss<<"B";
                break;
            case 'E':
                ss<<"e";
                break;
            case '-':
                ss<<"-";
                break;
            }
            break;

        case 'e':
            switch(tihl.at(i)) {
            case 'H':
                ss<<"H";
                break;
            case 'N':
                ss<<"B";
                break;
            case 'E':
                ss<<"e";
                break;
            case '-':
                ss<<"-";
                break;
            }
            break;

        case '-':
            switch(tihl.at(i)) {
            case 'H':
                ss<<"H";
                break;
            case 'N':
                ss<<"B";
                break;
            case 'E':
                ss<<"e";
                break;
            case '-':
                ss<<"-";
                break;
            }
            break;

        }
    }
    return ss.str();
}

//Function to return the triplet corresponding to the given index, that leads to the non emmitting state
Triplet* TripleAligner::Get_ToNonEmmitingTriplet(int r_idx) const {
    for (int i=0; i<Triplets[r_idx].size(); i++) {
        if ( StringCount(Triplets[r_idx][i]->Get_SoanAfter(), "e")==Triplets[r_idx][i]->Get_SoanAfter().size() ) return Triplets[r_idx][i];
    }
    return NULL;
}

//Function to return the triplet that leaves from the non emmitting state and transits to a soan labeled according to the given argument
Triplet* TripleAligner::Get_FromNonEmmitingTriplet(string sAfter, IntCVecRef idxs) const {
    for (int i=0; i<Triplets[idx_R].size(); i++) {
        if ( (Triplets[idx_R][i]->Get_SoanAfter().compare(sAfter)==0 ) &&
                (Triplets[idx_R][i]->GetTripletIdxChanges() == idxs )
           ) return Triplets[idx_R][i];
    }
    return NULL;
}

//Function that returns the corresponding index (from 0 - soans.size()-1) that the given soan
//corresponds to in the vector with all possible soans (e.g. HHH is at index 0, HHB is at index 1,...
//ebB is at index 40....)
int TripleAligner::IndexOfSoan(string soan) {
    for (unsigned int i=0; i<Triplets.size(); i++) {
        if ( Triplets[i][0]->Get_SoanBefore().compare(soan) == 0 )return i;
    }
    return -1;
}

//Function that checks wheteher a Transition is Valid according to the given Homologies
bool TripleAligner::IsValidTransition(IntVecRef NextIdxs, string tihl) {

    for (int j=0; j<NextIdxs.size()-1; j++) {
        int pos = -1;
        for (int i=0; i<Anchors.size(); i++) {
            if (Anchors[i][j]==NextIdxs[j]) {
                pos = i;
                break;
            }
        }

        if (pos!=-1) { //NextIdxs[j] was found in the pos-th Anchor point
            for (int k=0; k<Anchors[pos].size(); k++) {
                if (j==k || Anchors[pos][k]< NextIdxs[k] || (Anchors[pos][k]==NextIdxs[k] && (tihl.at(k)=='B' || tihl.at(k)=='-' || tihl.at(k)=='E') ) ) {
                    continue;
                }

                if (Anchors[pos][k]!= NextIdxs[k] && Anchors[pos][k]!=-1) {
                    return false;
                } else if (Anchors[pos][k]!= -1 && tihl.at(k)!='H') {
                    return false;
                }
            }
        } else {        //NextIdxs[j] corresponds to -1 in the Anchors vector
            if (tihl.at(j)=='H') {
                return false;
            }
        }
    }

    return true;
}

//Felsenstein Tables: Overloaded Function to get the Corresponding Substitution Probability for the given transition to the next state
double TripleAligner::GetSubstitutionProbability(string tihlLbl, const vector<int>& Idxs, int stx) const {
    //Argument (stx) is used in order to define the starting position for the H case. More specifically,
    //it can be either set to 0 (when every H contained in the tihl will be normally examined - disregarding the given Homologies)
    // or to 1 for the case when the given Homologies are considered and thus the first H (e.g. HHEH) should be negelected from
    //the Substitution Probability estimation

    double subProb = 1;
    for (int i=stx; i<tihlLbl.size(); i++) {
        switch(tihlLbl.at(i)) {
        case '-':
        case 'E':
            break;
        case 'H': {
            vector<int> H_Idxs;
            GetIndicesOf(H_Idxs, tihlLbl, "H");
            double temp = 0;
            for (int j=0; j<4; j++) { //For every nucleotide A,C,G,T
                double tempProd = Params->Get_BaseFreq()[j];
                for (int k=stx; k<H_Idxs.size(); k++) { //Parse every H position
                    double tempSum = 0;
                    for (int l=0; l<4; l++) { //For every nucleotide A,C,G,T
                        tempSum += SubMats[H_Idxs[k]][GetAlphabetPair(GetAlphabet(j), GetAlphabet(l))] * FPTables[H_Idxs[k]][l][Idxs[H_Idxs[k]]];
                    }
                    tempProd *= tempSum;
                }
                temp += tempProd;
            }

            subProb *= temp;

            for (int k=0; k<H_Idxs.size(); k++) { //here always start from 0 (and not from stx) since even the initial H (when considering the given Holomologies) should be converted in -
                tihlLbl.at(H_Idxs[k]) = '-';
            }

            break;
        }
        case 'B':
        case 'N': {
            double tempSum = 0;

            for (int j=0; j<4; j++) { //For every nucleotide A,C,G,T
                tempSum += Params->Get_BaseFreq()[j] * FPTables[i][j][Idxs[i]];
            }
            subProb *= tempSum;
            break;
        }
        }
    }
    return subProb;
}

//Function to Get the Triplet starting with the given SoanLbl and transitioning to the End State (i.e. HHH..H) via the corresponding all Homologous Tihl
Triplet* TripleAligner::GetEndStateTriplet(string StartSoanLbl) {
    for (int i=0; i<Triplets.size(); i++) {
        if (Triplets[i][0]->Get_SoanBefore().compare(StartSoanLbl)!=0) continue;
        return Triplets[i][0];
    }
}

//}

/** Probabilities Estimation Functions **/ //{
//Function to calculate all the probabilities describing all evolutionary histories for the alignement of the given sequences. Argument defines which MakeTransition Function should be called
void TripleAligner::EstimateProbabilities(int choice, bool HasExtraEndState) {
    if (SaveToLogFile) {
        WriteToFile(LogFileName, "", false, true);
    }

    for( map<LongInt, ProbabilityObject* >::iterator ii=Probabilities.begin(); ii!=Probabilities.end(); ++ii) {
        LongInt CurPosIn_1D = (*ii).second->Get_PosIn_1D();
        vector<int> CurIdxs = (*ii).second->Get_PosIn_ND();

        if (CurPosIn_1D==EndPosIn_1D) continue;

        /** 2. Transit to the next state **/
        switch (choice) {
        case 0: //Probabilities Ignoring Given Homologies - Forward Algorithm (For triplewise alignment)
            MakeTransition(CurIdxs, CurPosIn_1D, 0, HasExtraEndState);
            break;

        case 1: //Probabilities Considering Given Homologies - Forward Algorithm (For triplewise alignment)
            MakeTransition(CurIdxs, CurPosIn_1D, 1, HasExtraEndState);
            break;
        }
    }

    //If the Current end state is the non-emitting one then make a transition to the non-emitting state (wing folding fix)
    if (EndPosIn_ND[EndPosIn_ND.size()-1]==26){
        TransitToNonEmittingState();
    }
}

//Triplewise Alignment: Function to make the transition from current state to the next state
void TripleAligner::MakeTransition(IntVecRef CurIdxs, LongInt CurPosIn_1D, int stx, bool HasExtraEndState) {
    /** Log Visiting States */
    stringstream logStream;
    if (SaveToLogFile) {
        logStream.str("");
        logStream<<endl<<"Moving from State: [";
        for (long i=0; i<CurIdxs.size(); i++) {
            if (i<CurIdxs.size()-1) {
                logStream<<left<<setw(3)<<CurIdxs[i]<<",";
            } else {
                logStream<<left<<setw(3)<<CurIdxs[i];
            }
        }
        logStream<<"] - "<<setw(4)<<CurPosIn_1D << " - " << Probabilities[CurPosIn_1D]->Get_SoanLbl() <<" : "<<endl;
    }

    /** 2. Start Applying each possible tihl to this state **/
    Triplet* nonEmittingTriplet = Triplets[idx_R][idx_C];

    for (unsigned int i=0; i<Triplets[CurIdxs[CurIdxs.size()-1]].size(); i++) {
        Triplet* currentTriplet  = Triplets[CurIdxs[CurIdxs.size()-1]][i];

        if (!currentTriplet->Is_Valid()) continue;

        //if (CurPosIn_1D==737) { cout<<"Current Tihl: "<<currentTriplet->Get_Tihl()<<endl;  }

        //When Considering the given Homologies (i.e. stx==1) Skip tihls that do not start with 'H' or '-' -> i.e. for the root node to either have a B or not
        if (stx==1) {
            if (!StartsWith(currentTriplet->Get_Tihl(), "H") && !StartsWith(currentTriplet->Get_Tihl(), "-")) {
                //if (CurPosIn_1D==737) { cout<<"Skipped - Not Starting with H nor -"<<endl;  }
                continue;
            }
        }

        //Skip non-emitting state - non-emitting tihl since it is considered at the emitting tihls afterwards
        if (StringCount(currentTriplet->Get_SoanAfter(), "e")==currentTriplet->Get_SoanAfter().size()) {
            //if (CurPosIn_1D==737) { cout<<"Skipped - Is non emitting state"<<endl;  }
            continue;
        }

        Triplet* toNonEmmitingTriplet   = NULL;
        Triplet* fromNonEmmitingTriplet = NULL;

        if (!Contains(currentTriplet->Get_Tihl(), "B")) {
            toNonEmmitingTriplet   = Get_ToNonEmmitingTriplet(CurIdxs[CurIdxs.size()-1]);
            fromNonEmmitingTriplet = Get_FromNonEmmitingTriplet(currentTriplet->Get_SoanAfter(), currentTriplet->GetTripletIdxChanges());
        }

        //2a. Create the corresponding indices (both in 1D and in ND) for the next state
        //(i.e. the state to go, after applying current Tihl to current state - i.e. states[0])
        IntVec NextIdxs(CurIdxs);
        NextIdxs += currentTriplet->GetTripletIdxChanges();
        NextIdxs[NextIdxs.size()-1] = IndexOfSoan(currentTriplet->Get_SoanAfter());
        LongInt NextPosIn_1D = GetPosForIndexes(NextIdxs, Dimensions);

        //2b. Check if this transition leads to creating sequence(s) larger than the ones provided. If so, then skip this tihl
        if (NextPosIn_1D > EndPosIn_1D ) {
            //if (CurPosIn_1D==737) { cout<<"Skipped - NextPosIn_1D > EndPosIn_1D ("<<NextPosIn_1D<<">"<<EndPosIn_1D<<")"<<endl;  }
            continue;
        }

        bool LeadsToLargerSequence = false;

        for (unsigned int j=0; j<NextIdxs.size()-1; j++) {
            if (NextIdxs[j]>=(Dimensions[j])) {
                LeadsToLargerSequence = true;    //i.e. this tihl leads to the emission of more sites for the i-th sequence than those that the i-th sequence has
                //if (CurPosIn_1D==737) { cout<<"Skipped - Leads to Larger Sequence"<<endl; cout<<Dimensions<<endl; }
                break;
            }
        }
        if (LeadsToLargerSequence) continue;

        /*
        if (CurPosIn_1D==2748472) {
            cout<<"NextIds: [";
            for (int k=0;k<NextIdxs.size();k++){
                cout<<setw(2)<<NextIdxs[k];
                if (k<NextIdxs.size()-1) cout<<",";
            }
            cout<<"]"<<endl;
        }*/

        //2c. Check if the NextIdxs are in accordance with the given H(omologies) <=> Anchor Points
        //check is performed only when stx==1, i.e. given homolgies are considered, otherwise this check is ignored (i.e. when diregarding the given homologies)
        if (stx==1) {
            if (!IsValidTransition(NextIdxs, currentTriplet->Get_Tihl())) {
                //if (CurPosIn_1D==737) { cout<<"Skipped - Invalid according to anchor points"<<endl;  }
                continue;
            }
        }

        //2d. Add a new contributor to the ProbabilityObject of the Next State (defined by NextIdxs)
        //First Create an Empty ProbabilityObject at the NextPos if no such obejct exists
        if (Probabilities.count(NextPosIn_1D)==0) {
            ProbabilityObject* prob = new ProbabilityObject();
            prob->Set_OverallProbability(log(0.0));
            prob->Set_SoanLbl("");
            prob->Set_PosIn_1D(NextPosIn_1D);
            prob->Set_PosIn_ND(NextPosIn_1D, Dimensions);
            prob->Set_EndState(false);
            Probabilities[NextPosIn_1D] = prob;
        }

        Probabilities[NextPosIn_1D]->AddContributor(Probabilities[CurPosIn_1D]);
        if (fromNonEmmitingTriplet!=NULL) {
            stringstream st;
            st.str("");
            st << "(" << currentTriplet->Get_Tihl() << ") AND (" << toNonEmmitingTriplet->Get_Tihl() <<" -> " << fromNonEmmitingTriplet->Get_Tihl() << ")";
            Probabilities[NextPosIn_1D]->AddContributorTihlExt(st.str());
        } else {
            Probabilities[NextPosIn_1D]->AddContributorTihlExt(currentTriplet->Get_Tihl());
        }
        Probabilities[NextPosIn_1D]->AddContributorTihl(currentTriplet->Get_Tihl());

        double contributedProbability = 0;

        if (fromNonEmmitingTriplet==NULL) {
            double subProb;
            if (HasExtraEndState) {
                if (NextPosIn_1D==EndPosIn_1D) {
                    subProb = 1;
                } else {
                    subProb = ( GetSubstitutionProbability(currentTriplet->Get_Tihl(), NextIdxs, stx) );
                }
            } else {
                subProb = ( GetSubstitutionProbability(currentTriplet->Get_Tihl(), NextIdxs, stx) );
            }
            contributedProbability = Probabilities[CurPosIn_1D]->Get_OverallProbability() + log(currentTriplet->Get_Probability() ) + log(subProb);

            Probabilities[NextPosIn_1D]->AddContributorWingFold(false);
            Probabilities[NextPosIn_1D]->AddContributorDirectTransProb(contributedProbability);
            Probabilities[NextPosIn_1D]->AddContributorInDirectTransProb(0.0f);
            Probabilities[NextPosIn_1D]->AddContributorSelfTransProb(0.0f);
        } else {
            double subProb;
            if (HasExtraEndState) {
                if (NextPosIn_1D==EndPosIn_1D) {
                    subProb = 1;
                } else {
                    subProb = ( GetSubstitutionProbability(currentTriplet->Get_Tihl(), NextIdxs, stx) );
                }
            } else {
                subProb = ( GetSubstitutionProbability(currentTriplet->Get_Tihl(), NextIdxs, stx) );
            }

            double direct   = ( Probabilities[CurPosIn_1D]->Get_OverallProbability() ) + log(currentTriplet->Get_Probability() ) + log (subProb );
            double indirect = ( Probabilities[CurPosIn_1D]->Get_OverallProbability() ) + log(toNonEmmitingTriplet->Get_Probability() * (fromNonEmmitingTriplet->Get_Probability()/(1-nonEmittingTriplet->Get_Probability())) ) + log(subProb);
            contributedProbability = AddLog(direct, indirect);

            Probabilities[NextPosIn_1D]->AddContributorWingFold(true);
            Probabilities[NextPosIn_1D]->AddContributorDirectTransProb(direct);
            Probabilities[NextPosIn_1D]->AddContributorInDirectTransProb(indirect);
            Probabilities[NextPosIn_1D]->AddContributorSelfTransProb(log(nonEmittingTriplet->Get_Probability()));
        }

        Probabilities[NextPosIn_1D]->AddContributorProbability(contributedProbability);
        Probabilities[NextPosIn_1D]->Set_PosIn_1D(NextPosIn_1D);
        Probabilities[NextPosIn_1D]->Set_PosIn_ND(NextIdxs);
        Probabilities[NextPosIn_1D]->IncreaseOverallProbabilityBy(contributedProbability);
        Probabilities[NextPosIn_1D]->Set_SoanLbl(currentTriplet->Get_SoanAfter());

        if (SaveToLogFile) {
            logStream<<"\t     to State: [";
            for (int i=0; i<NextIdxs.size(); i++) {
                if (i<NextIdxs.size()-1) {
                    logStream<<left<<setw(3)<<NextIdxs[i]<<",";
                } else {
                    logStream<<left<<setw(3)<<NextIdxs[i];
                }
            }
            logStream<<"] - "<<setw(4)<<NextPosIn_1D << " - " << Probabilities[NextPosIn_1D]->Get_SoanLbl()<<" via Tihl: " << Probabilities[NextPosIn_1D]->GetContributorTihlExtAtPos(Probabilities[NextPosIn_1D]->Get_Contributors().size()-1)<<"  Prob:"<<contributedProbability<<endl;
        }

        //2e. Examine if the Next State is a (pre)End State. Although the actual end state will be the one with all edges H(omomlogous)
        //if the Next State defined until now has indexing equal to the lengths of the input sequences then it is a
        //pre-final state and does not have to be visited again apart from transitioning to the all H-state
        bool isEndState = true;
        for (unsigned int j=0; j<NextIdxs.size()-1; j++) {
            if (NextIdxs[j]<(Dimensions[j]-1)) {
                isEndState = false;
                break;
            }
        }

        if (isEndState)    {
            Probabilities[NextPosIn_1D]->Set_EndState(true);
            continue;
        } else               {
            Probabilities[NextPosIn_1D]->Set_EndState(false);
        }
        //if (CurPosIn_1D==737) { cout<<"Examined Thoroughly"<<endl;  }
    }
    if (SaveToLogFile) {
        WriteToFile(LogFileName, logStream.str(), true, false);
    }
}

//Function to Transit to the Non Emitting State for the cases where the End State is (H)EEE (i.e. state 26) since this is never achieved during the normal forward computation step because of the wing-folding
void TripleAligner::TransitToNonEmittingState() {
    ProbabilityObject* prob = new ProbabilityObject();
    prob->Set_OverallProbability(log(0.0));
    prob->Set_PosIn_1D(EndPosIn_1D);
    prob->Set_PosIn_ND(EndPosIn_ND);
    prob->Set_EndState(true);
    Probabilities[EndPosIn_1D] = prob;

    IntVec MinPosIn_ND(EndPosIn_ND);
    MinPosIn_ND[MinPosIn_ND.size()-1] = 0;
    LongInt MinPosIn_1D = GetPosForIndexes(MinPosIn_ND, Dimensions);
    int counter = 0;

    //Now form the end state's contributors
    for( map<LongInt, ProbabilityObject* >::reverse_iterator ii=Probabilities.rbegin(); ii!=Probabilities.rend(); ++ii) {
        if (counter ==0) { counter++; continue; }
        ProbabilityObject* cur_prob = (*ii).second;

        if (!cur_prob->Is_EndState()             ) continue;
        if ( cur_prob->Get_PosIn_1D() < MinPosIn_1D) break;

        IntVec CurIdxs = cur_prob->Get_PosIn_ND();
        Triplet* toNonEmmitingTriplet = Get_ToNonEmmitingTriplet(CurIdxs[CurIdxs.size()-1]);
        double contributedProbability = 1 * cur_prob->Get_OverallProbability();
        Probabilities[EndPosIn_1D]->AddContributor(cur_prob);
        Probabilities[EndPosIn_1D]->AddContributorTihl(toNonEmmitingTriplet->Get_Tihl());
        Probabilities[EndPosIn_1D]->AddContributorTihlExt(toNonEmmitingTriplet->Get_Tihl());
        Probabilities[EndPosIn_1D]->AddContributorProbability(contributedProbability);
        Probabilities[EndPosIn_1D]->IncreaseOverallProbabilityBy(contributedProbability);
        Probabilities[EndPosIn_1D]->Set_SoanLbl(toNonEmmitingTriplet->Get_SoanAfter());
        Probabilities[EndPosIn_1D]->AddContributorWingFold(false);
        Probabilities[EndPosIn_1D]->AddContributorDirectTransProb(contributedProbability);
        Probabilities[EndPosIn_1D]->AddContributorInDirectTransProb(0.0f);
        Probabilities[EndPosIn_1D]->AddContributorSelfTransProb(0.0f);
    }
}

//}

/** Functions to Sample Alignments **/ //{
//Function to sample evolutionary histories based on the calculated Probabilities
void TripleAligner::GetEvolHistories(StringMatRef evolHist, DoubleVecRef probs) {
    //For each one of the possible evolution histories....
    for (int i=0; i<Params->Get_SampleSize(); i++) {
        StringVec tempStrV;
        double tempProb;
        GetEvolHistory(tempStrV, tempProb);
        evolHist.push_back(tempStrV);
        probs.push_back(tempProb);
    }
}

//Function to sample only one evolutionary history based on the calculated Probabilities
//This Function expands the Ns and also accounts for the root sites

/* To Go back only at most probable previous state (to MAX)*/
/*
void TripleAligner::GetEvolHistory(StringVecRef hist, double& prob, bool FixNs) {
    prob = log(0);

    //Start from the End state...
    ProbabilityObject* state = Probabilities[EndPosIn_1D];

    //Init the string vector representing the evolution history of each sequence
    hist.resize(Dimensions.size(), "");

    //Move Backwards until the start state has been reached
    //while(state != Probabilities[0]) {
    while(state != NULL) {
        string tihl;
        int j;

        //if it is the "default" extra start state skip the contributors finding process
        if (state->Get_PosIn_1D()!=StartPosIn_1D) {
            //Get a random number and select a contributor according to this random value...
            double max = state->GetContributorProbabilityAtPos(0);
            for (int k=0;k<state->Get_Contributors().size();k++){
                if (state->GetContributorProbabilityAtPos(k)>max){
                    max = state->GetContributorProbabilityAtPos(k);
                }
            }

            for (j=0; j<state->Get_Contributors().size(); j++) {
                if (state->GetContributorProbabilityAtPos(j)==max){
                    break;
                }
            }


            //Get the contributor's tihl and form corresponding evolution history and alignment for each of the sequences
            tihl = state->GetContributorTihlAtPos(j);

            if (Contains(tihl, "B")) {
                tihl.insert(0, "-");
            }else{
                tihl.insert(0, "H");
            }
        }else{
            tihl = StartTihlLbl;
            if (Contains(tihl, "B")) {
                tihl.insert(0, "-");
            }else{
                tihl.insert(0, "H");
            }
        }


        //Continue back to the contributor (after having considered the wing folded case)
        //Convert Ns to Distinct Columns
        vector<string> tihlVec;
        if (Contains(tihl, "N")) {
            int pos = 0;
            while( (pos = LastIndexOf(tihl, "N"))>=0) {
                string s = string(tihl.size(), '-');
                s.at(pos) = 'N';
                tihl.at(pos) = '-';
                tihlVec.push_back(s);
            }
        }

        tihlVec.push_back(tihl);

        for (int k=0; k<tihlVec.size(); k++) {
            for (int p=0; p<tihlVec[k].size(); p++) {
                stringstream st;
                st.str("");
                st << tihlVec[k].at(p) << hist[p];
                hist[p] = st.str();
            }
        }
        if (state->Get_PosIn_1D()!=StartPosIn_1D) { prob = AddLog(prob, state->GetContributorDirectTransProbAtPos(j)); }
        if (state == Probabilities[StartPosIn_1D]) {
            state = NULL;
        }else{
            state = state->GetContributorAtPos(j);
        }
    }
}
*/

/* To Sample previous Step (randomly sampled)*/
void TripleAligner::GetEvolHistory(StringVecRef hist, double& prob, bool FixNs) {
    prob = log(0);

    //Start from the End state...
    ProbabilityObject* state = Probabilities[EndPosIn_1D];

    //Init the string vector representing the evolution history of each sequence
    hist.resize(Dimensions.size(), "");

    //Move Backwards until the start state has been reached
    //while(state != Probabilities[0]) {
    while(state != NULL) {
        string tihl;
        int j;

        //if it is the "default" extra start state skip the contributors finding process
        if (state->Get_PosIn_1D()!=StartPosIn_1D) {
            //Get a random number and select a contributor according to this random value...
            double randProb = (*MyRand).GetRandomDouble(0.0f, 1.0f);
            for (j=0; randProb>0; j++) {
                randProb -= exp(state->GetContributorProbabilityAtPos(j) - state->Get_OverallProbability());
            }
            j--;

            //Get the contributor's tihl and form corresponding evolution history and alignment for each of the sequences
            tihl = state->GetContributorTihlAtPos(j);

            if (Contains(tihl, "B")) {
                tihl.insert(0, "-");
            } else {
                tihl.insert(0, "H");
            }

            //Sample non-emitting states
            if (state->GetContributorWingFoldAtPos(j)) {
                randProb = (*MyRand).GetRandomDouble(0.0f, 1.0f);
                randProb -= exp(state->GetContributorDirectTransProbAtPos(j) - state->GetContributorProbabilityAtPos(j));
                if ( randProb > 0 ) {
                    prob = AddLog(prob, state->GetContributorInDirectTransProbAtPos(j));
                    for (int j=0; j<tihl.size(); j++) {
                        stringstream st;
                        st.str("");
                        if(j==0) {
                            st << "H" << hist[j];
                        } else {
                            st << "E" << hist[j];
                        }
                        hist[j] = st.str();
                    }

                    while( (*MyRand).GetRandomDouble(0.0f, 1.0f) < exp(state->GetContributorSelfTransProbAtPos(j))) {
                        prob = AddLog(prob, state->GetContributorSelfTransProbAtPos(j));
                        for (int j=0; j<tihl.size(); j++) {
                            stringstream st;
                            st.str("");

                            if(j==0) {
                                st << "H" << hist[j];
                            } else {
                                st << "E" << hist[j];
                            }
                            hist[j] = st.str();
                        }
                    }
                }
            }
        } else {
            tihl = StartTihlLbl;
            if (Contains(tihl, "B")) {
                tihl.insert(0, "-");
            } else {
                tihl.insert(0, "H");
            }
        }


        //Continue back to the contributor (after having considered the wing folded case)
        //Convert Ns to Distinct Columns
        vector<string> tihlVec;
        if (Contains(tihl, "N")) {
            int pos = 0;
            while( (pos = LastIndexOf(tihl, "N"))>=0) {
                string s = string(tihl.size(), '-');
                s.at(pos) = 'N';
                tihl.at(pos) = '-';
                tihlVec.push_back(s);
            }
        }

        tihlVec.push_back(tihl);

        for (int k=0; k<tihlVec.size(); k++) {
            for (int p=0; p<tihlVec[k].size(); p++) {
                stringstream st;
                st.str("");
                st << tihlVec[k].at(p) << hist[p];
                hist[p] = st.str();
            }
        }
        if (state->Get_PosIn_1D()!=StartPosIn_1D) {
            prob = AddLog(prob, state->GetContributorDirectTransProbAtPos(j));
        }
        if (state == Probabilities[StartPosIn_1D]) {
            state = NULL;
        } else {
            state = state->GetContributorAtPos(j);
        }
    }
}

//Function to sample only one evolutionary history based on the calculated Probabilities
void TripleAligner::GetEvolHistory(StringVecRef hist, double& prob) {
    prob = log(0);

    //Start from the End state...
    ProbabilityObject* state = Probabilities[EndPosIn_1D];

    //Init the string vector representing the evolution history of each sequence
    hist.resize(Dimensions.size(), "");

    //Move Backwards until the start state has been reached
    //while(state != Probabilities[0]) {
    while(state != NULL) {
        string tihl;
        int j;

        //if it is the "default" extra start state skip the contributors finding process
        if (state->Get_PosIn_1D()!=StartPosIn_1D) {
            //Get a random number and select a contributor according to this random value...
            double randProb = (*MyRand).GetRandomDouble(0.0f, 1.0f);
            for (j=0; randProb>0; j++) {
                randProb -= exp(state->GetContributorProbabilityAtPos(j) - state->Get_OverallProbability());
            }
            j--;

            //Get the contributor's tihl and form corresponding evolution history and alignment for each of the sequences
            tihl = state->GetContributorTihlAtPos(j);

            if (Contains(tihl, "B")) {
                tihl.insert(0, "-");
            } else {
                tihl.insert(0, "H");
            }

            //Sample non-emitting states
            if (state->GetContributorWingFoldAtPos(j)) {
                randProb = (*MyRand).GetRandomDouble(0.0f, 1.0f);
                randProb -= exp(state->GetContributorDirectTransProbAtPos(j) - state->GetContributorProbabilityAtPos(j));
                if ( randProb > 0 ) {
                    prob = AddLog(prob, state->GetContributorInDirectTransProbAtPos(j));
                    for (int j=0; j<tihl.size(); j++) {
                        stringstream st;
                        st.str("");
                        if(j==0) {
                            st << "H" << hist[j];
                        } else {
                            st << "E" << hist[j];
                        }
                        hist[j] = st.str();
                    }

                    while( (*MyRand).GetRandomDouble(0.0f, 1.0f) < exp(state->GetContributorSelfTransProbAtPos(j))) {
                        prob = AddLog(prob, state->GetContributorSelfTransProbAtPos(j));
                        for (int j=0; j<tihl.size(); j++) {
                            stringstream st;
                            st.str("");

                            if(j==0) {
                                st << "H" << hist[j];
                            } else {
                                st << "E" << hist[j];
                            }
                            hist[j] = st.str();
                        }
                    }
                }
            }
        } else {
            tihl = StartTihlLbl;
            if (Contains(tihl, "B")) {
                tihl.insert(0, "-");
            } else {
                tihl.insert(0, "H");
            }
        }

        //Continue back to the contributor (after having considered the wing folded case)
        for (int p=0; p<tihl.size(); p++) {
            stringstream st;
            st.str("");
            st << tihl.at(p) << hist[p];
            hist[p] = st.str();
        }

        if (state->Get_PosIn_1D()!=StartPosIn_1D) {
            prob = AddLog(prob, state->GetContributorDirectTransProbAtPos(j));
        }

        if (state == Probabilities[StartPosIn_1D]) {
            state = NULL;
        } else {
            state = state->GetContributorAtPos(j);
        }
    }
}

//Function to sample only one evolutionary history based on the calculated Probabilities and
//also "return" the sequences of soans as well as their corresponding indices
void TripleAligner::GetEvolHistory(StringVecRef hist, double& prob, StringVecRef soans, StringVecRef tihls, IntMatRef idxs) {
    for (int i=0; i<idxs.size(); i++) {
        idxs[i].clear();
    }
    idxs.clear();

    prob = log(0);
    //Start from the End state...
    ProbabilityObject* state = Probabilities[EndPosIn_1D];

    //Init the string vector representing the evolution history of each sequence
    hist.resize(Dimensions.size()-1, "");

    //Move Backwards until the start state has been reached
    //while(state != Probabilities[0]) {
    while(state != NULL) {
        string tihl;
        int j;

        //if it is the "default" extra start state skip the contributors finding process
        if (state->Get_PosIn_1D()!=StartPosIn_1D) {
            //Get a random number and select a contributor according to this random value...
            double randProb = (*MyRand).GetRandomDouble(0.0f, 1.0f);
            for (j=0; randProb>0; j++) {
                randProb -= exp(state->GetContributorProbabilityAtPos(j) - state->Get_OverallProbability());
            }
            j--;

            //Get the contributor's tihl and form corresponding evolution history and alignment for each of the sequences
            tihl = state->GetContributorTihlAtPos(j);

            //Sample non-emitting states
            if (state->GetContributorWingFoldAtPos(j)) {
                randProb = (*MyRand).GetRandomDouble(0.0f, 1.0f);
                randProb -= exp(state->GetContributorDirectTransProbAtPos(j) - state->GetContributorProbabilityAtPos(j));
                if ( randProb > 0 ) {
                    prob = AddLog(prob, state->GetContributorInDirectTransProbAtPos(j));

                    for (int j=0; j<tihl.size(); j++) {
                        stringstream st;
                        st.str("");
                        st << "E" << hist[j];
                        hist[j] = st.str();
                    }
                    string nullSoan(tihl.size(), 'e');
                    nullSoan.at(0) = 'H';
                    soans.push_back(nullSoan);

                    string nullTihl(tihl.size(), 'E');
                    nullTihl.at(0) = 'H';
                    tihls.push_back(nullTihl);

                    idxs.push_back(state->Get_PosIn_ND());


                    while( (*MyRand).GetRandomDouble(0.0f, 1.0f) < exp(state->GetContributorSelfTransProbAtPos(j))) {
                        prob = AddLog(prob, state->GetContributorSelfTransProbAtPos(j));
                        for (int j=0; j<tihl.size(); j++) {
                            stringstream st;
                            st.str("");
                            st << "E" << hist[j];
                            hist[j] = st.str();
                        }
                        soans.push_back(nullSoan);
                        tihls.push_back(nullTihl);
                        idxs.push_back(state->Get_PosIn_ND());
                    }
                }
            }
        } else {
            tihl = StartTihlLbl;
        }

        //Continue back to the contributor (after having considered the wing folded case)
        for (int p=0; p<tihl.size(); p++) {
            stringstream st;
            st.str("");
            st << tihl.at(p) << hist[p];
            hist[p] = st.str();
        }

        soans.push_back(state->Get_SoanLbl());
        tihls.push_back(tihl);
        idxs.push_back(state->Get_PosIn_ND());

        if (state->Get_PosIn_1D()!=StartPosIn_1D) {
            prob = AddLog(prob, state->GetContributorProbabilityAtPos(j));
        }

        if (state == Probabilities[StartPosIn_1D]) {
            state = NULL;
        } else {
            state = state->GetContributorAtPos(j);
        }
    }

    reverse(soans.begin(),soans.end());
    reverse(tihls.begin(),tihls.end());
    reverse(idxs.begin(),idxs.end());
}

//}

/** Functions to Output results **/ //{
//Function to save the estimated probabilities to file
void TripleAligner::SaveResultsToFile(string filename) {
    WriteToFile(filename, "", true, true);
    for( map<LongInt, ProbabilityObject* >::iterator ii=Probabilities.begin(); ii!=Probabilities.end(); ++ii) {
        ProbabilityObject* temp = (*ii).second;

        if (temp->Get_OverallProbability()<0) {
            stringstream ss;
            ss.str("");
            ss<<"-"<<*(temp);
            WriteToFile(filename, ss.str(), true, false);
        } else if (temp->Get_OverallProbability()>0) {
            cout<<"Error - probability " << (*ii).first << " has a POSITIVE value!"<<endl;
            stringstream ss;
            ss.str("");
            ss<<*(temp);
            WriteToFile(filename, ss.str(), true, false);
        }
    }
}

//Function to Output the Program Parameters to the standard output stream
void TripleAligner::PrintParameters() {
    if (Params==NULL) {
        cout<<"\n Please specify program parameters first.\n\n";
        return;
    }
    cout<<" \n"<< *(Params)<<endl;
}

//Function to Print the subtree lengths
void TripleAligner::PrintLengths() {
    cout<<"\nCurrent Subtree Lengths:\n------------------------\n<";
    for (int i=0; i<Lengths.size(); i++) {
        cout<<setw(4)<<Lengths[i];
    }
    cout<<">\n\n";
}

//Function to Print the Substitution Matrices corresponding to the current fragment to the standard output stream
void TripleAligner::PrintSubstitutionMatrices() {
    cout<<"\nFragment's Substitution Matrices:\n---------------------------------\n";
    cout.precision(5);
    for (int i=0; i<SubMats.size(); i++) {
        cout<<"For Length: "<<Lengths[i]<<endl;
        cout<<"\t"<<setw(15)<<"A"<<setw(12)<<"C"<<setw(12)<<"G"<<setw(12)<<"T"<<setw(11)<<"SUM:"<<endl<<"\t";
        //vector<double> sumC(4, 0.0);

        for (int j=0; j<4; j++) {
            double sumL = 0;
            cout<<setw(3)<<GetAlphabet_String(j);
            for (int k=0; k<4; k++) {
                double tempVal = SubMats[i][GetAlphabetPair(j,k)];
                sumL += tempVal;
                cout<<setw(12)<<tempVal;
                //sumC[k] += tempVal;
            }
            cout<<"\t"<<setw(4)<<sumL<<endl<<"\t";
        }
        //cout<<"SUM:";
        //for (int k=0;k<4;k++) {
        //    cout<<setw(11)<<sumC[k]; }
        cout<<endl<<endl;
    }
}

//Function to Print the Felsenstein Pruning Tables for the Current Fragment to the standard output stream
void TripleAligner::PrintFPTables() {
    cout<<"\nFragment's FPTables:\n--------------------\n";

    //cout<<"(Sizes: #of FPTables:"<<FPTables.size()<<"  #of rows in FPTable[0]:"<<FPTables[0].size()<<"  #of sites in FPTable[0][0]:"<<FPTables[0][0].size()<<")"<<endl;
    cout.precision(5);
    for (int i=0; i<FPTables.size(); i++) {
        cout<<"FPTable - "<<i<<endl<<"-----------"<<endl;
        for(int j=0; j<FPTables[i].size(); j++) {
            cout<<"\t";
            for (int k=0; k<FPTables[i][j].size(); k++) {
                cout<<left<<setw(12)<<FPTables[i][j][k];
            }
            cout<<endl;
        }
        cout<<endl;
    }
}

//Function to Print the Dimensions for the current fragment to the standard output stream
void TripleAligner::PrintDimensions() {
    cout<<"\nDimensions:\n-----------\n<";
    for (int i=0; i<Dimensions.size(); i++) {
        cout<<setw(5)<<Dimensions[i];
    }
    cout<<">\n\n";
}

//Function to Print the Estimated Triplets
void TripleAligner::PrintTriplets() {
    cout<<"\nTriplets:\n---------\n";
    for (int i=0; i<Triplets.size(); i++) {
        double sum = 0;
        for (int j=0; j<Triplets[i].size(); j++) {
            sum += Triplets[i][j]->Get_Probability();
            cout<< *(Triplets[i][j])<<"\tProb:"<<Triplets[i][j]->Get_Probability()<<endl;
        }
        cout<<"\t\t SUM:"<<sum<<endl<<endl;
    }
}

//Function to Print the Substitution and Fid Models
void TripleAligner::PrintModels() {
    cout<<"\n"<< (*Fid)<<endl;
}

//Function to Print the Alignment Boundaries
void TripleAligner::PrintBoundaries() {
    cout<<"\nAlignment Region:\n-----------------\n";
    cout<<"Start State:  <";
    for (int i=0; i<StartPosIn_ND.size(); i++) {
        cout<<setw(5)<<StartPosIn_ND[i];
        if (i<StartPosIn_ND.size()-1) {
            cout<<",";
        }
    }
    cout<<">  (In 1D:"<<StartPosIn_1D<<")"<<endl;

    cout<<"  End State:  <";
    for (int i=0; i<EndPosIn_ND.size(); i++) {
        cout<<setw(5)<<EndPosIn_ND[i];
        if (i<EndPosIn_ND.size()-1) {
            cout<<",";
        }
    }
    cout<<">  (In 1D:"<<EndPosIn_1D<<")"<<endl;
}

//Function to call the above printing functions in order to output the curretn TripleAligner instance setup
void TripleAligner::PrintStatus() {
    PrintBoundaries();
    PrintDimensions();
    PrintLengths();
    //PrintFPTables();
    SaveFPTables();
    PrintSubstitutionMatrices();
    PrintParameters();
    PrintModels();
    //PrintTriplets();
}

//Fucntion to save the employed FPTables to an output file
void TripleAligner::SaveFPTables() {
    stringstream ss;
    ss.str("");
    ss<<"\nFragment's FPTables:\n--------------------\n";

    //cout<<"(Sizes: #of FPTables:"<<FPTables.size()<<"  #of rows in FPTable[0]:"<<FPTables[0].size()<<"  #of sites in FPTable[0][0]:"<<FPTables[0][0].size()<<")"<<endl;
    ss.precision(5);
    for (int i=0; i<FPTables.size(); i++) {
        ss<<"FPTable - "<<i<<endl<<"-----------"<<endl;
        for(int j=0; j<FPTables[i].size(); j++) {
            ss<<"\t";
            for (int k=0; k<FPTables[i][j].size(); k++) {
                ss<<left<<setw(12)<<FPTables[i][j][k];
            }
            ss<<endl;
        }
        ss<<endl;
    }

    WriteToFile(FPTablesFile, ss.str(), false, true);
}

//}

//Definition of the toString Operation
ostream& operator<<(ostream& ostr, const TripleAligner& aligner) {
    ostr << "Not defined yet." << endl;
    return ostr;
}
