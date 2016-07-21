#include "SequenceAligner.h"

//Default Constructor
SequenceAligner::SequenceAligner()  {
    //AlignPlainSequences();
    //AlignProbSequences();
    //SampleOnGivenAlignment();
    //SampleOnGivenHomologies();
}

//Virtual Decontructor
SequenceAligner::~SequenceAligner() { }

/** Starter Functions **/ //{
//Function to align plain sequences
void SequenceAligner::AlignPlainSequences() {
    InitParameters("Parameters-9.txt");     //Initialize Program's Parameters
    InitSubModel();                         //Initialize the Substitution Model
    InitFidModel();                         //Initialize the FID Model
    InitSequences("Demo-9.txt");            //Initialize the tree (i.e. Sequences to be aligned as well as equivalent branch lengths)
    InitTriplets();                         //Initialize all possible Triplets (organized per each reachable soan)
    InitProbabilities();                    //Initialize all ProbabilityObject with default values
    InitSubMatrices();                      //Initialize the Substitution Matrices for the given branch lengths

    cout<<"\t\tPARAMETERS"<<endl;
    cout<<"\t\t=========="<<endl<<(*Params);

    cout<<endl<<endl<<"\t\tSEQUENCES"<<endl;
    cout<<"\t\t========="<<endl;
    for (int i=0; i<Sequences.size(); i++) {
        cout<<"Seq-"<<(i+1)<<": "<<left<<setw(30)<<Sequences[i]<<"   Length: "<<left<<setw(9)<< Lengths[i] <<"   # of Sites: " << Sequences[i].size()<<right<<endl;
    }
    PrintSubstitutionMatrices();

    //start counting execution time
    unsigned t_start=clock(),t_end;
    cout<<"\nEstimating Transition Probabilities...";

    EstimateProbabilities(0);
    TransitToEndState();

    t_end=clock()-t_start;
    cout<<" ->  Time elapsed: " << (((float)t_end)/CLOCKS_PER_SEC) <<" seconds."<<endl;
    //SaveResultsToFile();

    cout<<"Sampling " << Params->Get_SampleSize() << " evolutionary histories...";
    t_start=clock();

    //Sample Evolutionary Histories...
    vector< vector<string> > evolHist;
    vector<double> evolHistProbs;
    GetEvolHistories(evolHist, evolHistProbs);

    t_end=clock()-t_start;
    cout<<" ->  Time elapsed: " << (((float)t_end)/CLOCKS_PER_SEC) <<" seconds."<<endl;

    //ShowAllSampledAlignments(evolHist, evolHistProbs);

    //Most Probable Evolutionary History
    int evolPos = GetPosOfMaxElement(evolHistProbs);
    vector<string> alignForEvolHist;
    GetAlignmentForEvolHistory(evolHist[evolPos], alignForEvolHist);

    cout<<"\n\n\nMost Probable Evolutionary History:\n";
    for (int i=0; i<alignForEvolHist.size(); i++) {
        cout<<alignForEvolHist[i]<<"\t"<<evolHist[evolPos][i]<<endl;
    }
    cout<<"\tAlign. Log-Likelihood:"<<AssessAlignment(alignForEvolHist)<<"\tHist. Log-Likelihood:"<<evolHistProbs[evolPos]<<endl;

    AssessAlignment();
}

//Function to align sequences provided in the form of probability matrices
void SequenceAligner::AlignProbSequences() {
    InitParameters("Parameters-1.txt");       //Initialize Program's Parameters
    InitSubModel();                         //Initialize the Substitution Model
    InitFidModel();                         //Initialize the FID Model
    InitSequences(5);                       //Initilaize the tree (i.e. Sequences to be aligned as well as equivalent branch lengths)
    InitTriplets();                         //Initilaize all possible Triplets (organized per each reachable soan)
    InitProbabilities();                    //Intialize all ProbabilityObject with default values
    InitSubMatrices();                      //Initialize the Substitution Matrices for the given branch lengths

    cout<<"\t\tPARAMETERS"<<endl;
    cout<<"\t\t=========="<<endl<<(*Params);
    PrintSequenceProbabilities();
    PrintSubstitutionMatrices();

    //start counting execution time
    unsigned t_start=clock(),t_end;
    cout<<"\nEstimating Transition Probabilities...";

    EstimateProbabilities(1);
    TransitToEndState();

    t_end=clock()-t_start;
    cout<<" ->  Time elapsed: " << (((float)t_end)/CLOCKS_PER_SEC) <<" seconds."<<endl;
    //SaveResultsToFile();

    cout<<"Sampling " << Params->Get_SampleSize() << " alignments...";
    t_start=clock();

    //Sample Evolutionary Histories and Alignments...
    vector< vector<string> > evolHist;
    vector<double> probs;
    GetEvolHistories(evolHist, probs);
    t_end=clock()-t_start;
    cout<<" ->  Time elapsed: " << (((float)t_end)/CLOCKS_PER_SEC) <<" seconds."<<endl;

    int pos = GetPosOfMaxElement(probs);

    cout<<"\n\n\nMost Probable Evolutionary History:\n";
    for (int i=0; i<evolHist[pos].size(); i++) {
        cout<<evolHist[pos][i]<<endl;
    }
    cout<<"\tHist. Log-Likelihood:"<<probs[pos]<<endl;
}

//Function to estimate all probabilities for a given (pre-defined) alignment
void SequenceAligner::SampleOnGivenAlignment() {
    InitParameters("Parameters-1.txt");       //Initialize Program's Parameters
    InitSubModel();                           //Initialize the Substitution Model
    InitFidModel();                           //Initialize the FID Model
    InitAlignedSequences("Alignment-1.txt");  //Initialize the tree (i.e. Sequences to be aligned as well as equivalent branch lengths)
    InitSubMatrices();                        //Initialize the Substitution Matrices for the given branch lengths

    cout<<"\t\tPARAMETERS"<<endl;
    cout<<"\t\t=========="<<endl<<(*Params);

    cout<<endl<<endl<<"\t\tGIVEN ALIGNMENT"<<endl;
    cout<<"\t\t=============="<<endl;
    for (int i=0; i<AlignedSequences.size(); i++) {
        cout<<"Seq-"<<(i+1)<<": "<< AlignedSequences[i] <<endl;
    }

    cout<<endl<<endl<<"\t\tSEQUENCES"<<endl;
    cout<<"\t\t========="<<endl;
    for (int i=0; i<Sequences.size(); i++) {
        cout<<"Seq-"<<(i+1)<<": "<<left<<setw(30)<<Sequences[i]<<"   Length: "<<left<<setw(9)<< Lengths[i] <<"   # of Sites: " << Sequences[i].size()<<right<<endl;
    }
    PrintSubstitutionMatrices();

    //start counting execution time
    unsigned t_start=clock(),t_end;
    cout<<"\nEstimating Transition Probabilities...";

    EstimateProbabilities(2);
    TransitToEndState();

    t_end=clock()-t_start;
    cout<<" ->  Time elapsed: " << (((float)t_end)/CLOCKS_PER_SEC) <<" seconds."<<endl;
    //SaveResultsToFile();

    cout<<"\n\n\nGiven Alignment:\n";
    for (int i=0; i<AlignedSequences.size(); i++) {
        cout<<"Seq-"<<(i+1)<<": "<<AlignedSequences[i]<<endl;
    }
    cout<<"\tAlign. Log-Likelihood:"<<AssessAlignment(AlignedSequences)<<endl;

    //Sample Evolutionary Histories...
    cout<<"\n\nSampling " << Params->Get_SampleSize() << " evolutionary histories...";
    t_start=clock();

    vector< vector<string> > evolHist;
    vector<double> evolHistProbs;
    GetEvolHistories(evolHist, evolHistProbs);

    t_end=clock()-t_start;
    cout<<" ->  Time elapsed: " << (((float)t_end)/CLOCKS_PER_SEC) <<" seconds."<<endl;

    //ShowAllSampledAlignments(evolHist, evolHistProbs);

    //Most Probable Evolutionary History
    int evolPos = GetPosOfMaxElement(evolHistProbs);
    vector<string> alignForEvolHist;
    GetAlignmentForEvolHistory(evolHist[evolPos], alignForEvolHist);

    cout<<"\n\n\nMost Probable Evolutionary History:\n";
    for (int i=0; i<alignForEvolHist.size(); i++) {
        cout<<alignForEvolHist[i]<<"\t"<<evolHist[evolPos][i]<<endl;
    }
    cout<<"\tAlign. Log-Likelihood:"<<AssessAlignment(alignForEvolHist)<<"\tHist. Log-Likelihood:"<<evolHistProbs[evolPos]<<endl;

}

//Function to estimate all probabilities for a given (pre-defined) set of Homologies between the sequences
void SequenceAligner::SampleOnGivenHomologies() {
    InitParameters("Parameters-1.txt");       //Initialize Program's Parameters
    InitSubModel();                           //Initialize the Substitution Model
    InitFidModel();                           //Initialize the FID Model
    InitHomologyIndices("Homology-3.txt");    //Initialize the provided homologies between the sequences
    InitTriplets();                           //Initialize all possible Triplets (organized per each reachable soan)
    InitProbabilities();                      //Initialize all ProbabilityObject with default values
    InitSubMatrices();                        //Initialize the Substitution Matrices for the given branch lengths

    cout<<"\t\tPARAMETERS"<<endl;
    cout<<"\t\t=========="<<endl<<(*Params);

    cout<<endl<<endl<<"\t\tGIVEN HOMOLOGIES"<<endl;
    cout<<"\t\t================"<<endl;
    for (int i=0; i<HomologiesIdxs.size(); i++) {
        cout<<"Root to Seq-"<<(i+1)<<": ";
        for (int j=0; j<HomologiesIdxs[i].size(); j++) {
            cout<<right<<setw(3)<<HomologiesIdxs[i][j] <<" ";
        }
        cout<<endl;
    }

    cout<<endl<<endl<<"\t\tANCHOR POINTS"<<endl;
    cout<<"\t\t============="<<endl;
    for (int i=0; i<Anchors.size(); i++) {
        cout<<"Anchor Point-"<<(i+1)<<": ";
        for (int j=0; j<Anchors[i].size(); j++) {
            cout<<right<<setw(3)<<Anchors[i][j] <<" ";
        }
        cout<<endl;
    }

    cout<<endl<<endl<<"\t\tSEQUENCES"<<endl;
    cout<<"\t\t========="<<endl;
    for (int i=0; i<Sequences.size(); i++) {
        if (i==0) {
            cout<<" Root: "<<left<<setw(30)<<Sequences[i]<<"   Length: "<<left<<setw(9)<< Lengths[i] <<"   # of Sites: " << Sequences[i].size()<<right<<endl;
        } else {
            cout<<"Seq-"<<(i)<<": "<<left<<setw(30)<<Sequences[i]<<"   Length: "<<left<<setw(9)<< Lengths[i] <<"   # of Sites: " << Sequences[i].size()<<right<<endl;
        }
    }
    PrintSubstitutionMatrices();

    //start counting execution time
    unsigned t_start=clock(),t_end;
    cout<<"\nEstimating Transition Probabilities...";

    EstimateProbabilities(3);
    TransitToEndState();

    t_end=clock()-t_start;
    cout<<" ->  Time elapsed: " << (((float)t_end)/CLOCKS_PER_SEC) <<" seconds."<<endl;
    SaveResultsToFile();

    //Sample Evolutionary Histories...
    cout<<"Sampling " << Params->Get_SampleSize() << " evolutionary histories...";
    t_start=clock();

    vector< vector<string> > evolHist;
    vector<double> evolHistProbs;
    GetEvolHistories(evolHist, evolHistProbs);

    t_end=clock()-t_start;
    cout<<" ->  Time elapsed: " << (((float)t_end)/CLOCKS_PER_SEC) <<" seconds."<<endl;

    //ShowAllSampledAlignments(evolHist, evolHistProbs);

    //Most Probable Evolutionary History
    int evolPos = GetPosOfMaxElement(evolHistProbs);
    vector<string> alignForEvolHist;
    GetAlignmentForEvolHistory(evolHist[evolPos], alignForEvolHist);

    cout<<"\n\nMost Probable Evolutionary History:\n";
    for (int i=0; i<alignForEvolHist.size(); i++) {
        cout<<alignForEvolHist[i]<<"\t"<<evolHist[evolPos][i]<<endl;
    }
    cout<<"\tAlign. Log-Likelihood:"<<AssessAlignment(alignForEvolHist)<<"\tHist. Log-Likelihood:"<<evolHistProbs[evolPos]<<endl;
}

//}

/** Initialization Functions **/ //{

//Function to Initialize all Program's Parameters
void SequenceAligner::InitParameters(string filename) {
    Params = new Parameters();
    Params->ReadDataFromFile(filename);
}

//Function to Initialize the FID Model
void SequenceAligner::InitFidModel() {
    Fid =  new FIDModel(Params->Get_Gamma(), Params->Get_Lambda());
}

//Function to Initialize the Substitution Model
void SequenceAligner::InitSubModel() {
    if (Params->Get_SubModel().compare("F84")==0) {
        f84 = new F84Model(Params->Get_BaseFreq(), Params->Get_Kappa());
        /* For Debugging */
        /*f84->SetLength(0.1);
        vector<double> colsSum(4, 0.0);
        for (int i=0;i<4;i++){
            double sum = 0;
            for (int j=0;j<4;j++){
                double prob = f84->GetProbability(GetAlphabet(i), GetAlphabet(j));
                sum += prob;
                colsSum[j] += prob;
                cout<<  prob << "\t";
            }
            cout<<sum<<endl;
        }

        for (int i=0;i<colsSum.size();i++){
            cout<<colsSum[i]<<"\t";
        }
        cout<<endl;
        Pause();
        */
    } else if (Params->Get_SubModel().compare("TN93")==0) {
        tn93 = new TN93Model(Params->Get_BaseFreq(), Params->Get_AG_Rate(), Params->Get_CT_Rate());
        /* For Debugging */
        /*tn93->SetLength(0.1);
        vector<double> colsSum(4, 0.0);
        for (int i=0;i<4;i++){
            double sum = 0;
            for (int j=0;j<4;j++){
                double prob = tn93->GetProbability(GetAlphabet(i), GetAlphabet(j));
                sum += prob;
                colsSum[j] += prob;
                cout<<  prob << "\t";
            }
            cout<<sum<<endl;
        }

        for (int i=0;i<colsSum.size();i++){
            cout<<colsSum[i]<<"\t";
        }
        cout<<endl;
        Pause();
        */
    } else if (Params->Get_SubModel().compare("GTR")==0) {
        gtr = new GTRModel(Params->Get_BaseFreq(), Params->Get_X1_Param(), Params->Get_X2_Param(), Params->Get_X3_Param(), Params->Get_X4_Param(), Params->Get_X5_Param(), Params->Get_X6_Param());
        /*For Debugging*/
        /*gtr->SetLength(0.1);
        vector<double> colsSum(4, 0.0);
        for (int i=0;i<4;i++){
            double sum = 0;
            for (int j=0;j<4;j++){
                double prob = gtr->GetProbability(GetAlphabet(i), GetAlphabet(j));
                sum += prob;
                colsSum[j] += prob;
                cout<<  prob << "\t";
            }
            cout<<sum<<endl;
        }

        for (int i=0;i<colsSum.size();i++){
            cout<<colsSum[i]<<"\t";
        }
        cout<<endl;
        Pause();
        */
    }
}

//Function to Initialize the Substitution Matrices for the given branch lengths
void SequenceAligner::InitSubMatrices() {
    vector<double> tempVec(16, 0.0);

    for (int i=0; i<Lengths.size(); i++) {
        double* subM = NULL;
        if (Params->Get_SubModel().compare("F84")==0) {
            subM = f84->GetSubstitutionMatrix(Lengths[i]);
        } else if (Params->Get_SubModel().compare("TN93")==0) {
            subM = tn93->GetSubstitutionMatrix(Lengths[i]);
        } else if (Params->Get_SubModel().compare("GTR")==0) {
            subM = gtr->GetSubstitutionMatrix(Lengths[i]);
        }

        for (int j=0; j<tempVec.size(); j++) {
            tempVec[j] = subM[j];
        }
        SubMats.push_back(tempVec);
    }
}

//Function to Initialize the tree (i.e. Sequences to be aligned as well as equivalent branch lengths) (plain sequences form)
void SequenceAligner::InitSequences() {
    /** Init lengths **/
    int numOfChildren;
    cout<<"Enter the number of leaf nodes of the tree -> ";
    cin>>numOfChildren;
    cout<<endl;

    for (int i=0; i<numOfChildren; i++) {
        double len;
        cout<<"\tEnter the length of Child - " << (i+1) <<": -> ";
        cin>>len;

        Lengths.push_back(len);
    }

    cout<<endl;

    /** Init Sequences **/
    for (int i=0; i<numOfChildren; i++) {
        string seq = "";
        cout<<"\tEnter the Sequence for the " << (i+1) << " Child: -> ";
        cin>>seq;
        Sequences.push_back(seq);
        Dimensions.push_back(seq.size());
    }
}

//Function to Initialize the tree (i.e. Sequences to be aligned as well as equivalent branch lengths) (plain sequences form - read from file)
void SequenceAligner::InitSequences(string filename) {
    /** Init lengths **/
    int numOfChildren;
    ifstream inputstr(filename.c_str());
    if (!inputstr) {
        cerr << "cannot open file " << filename.c_str() << endl;
        exit(1);
    }
    inputstr >> numOfChildren;

    for (int i=0; i<numOfChildren; i++) {
        double len;
        inputstr >> len;

        Lengths.push_back(len);
    }

    /** Init Sequences **/
    for (int i=0; i<numOfChildren; i++) {
        string seq = "";
        inputstr >> seq;
        Sequences.push_back(seq);
        Dimensions.push_back(seq.size());
    }

    inputstr.close();
}

//Function to Initialize the tree (i.e. Sequences to be aligned as well as equivalent branch lengths) (felsenstein form)
void SequenceAligner::InitSequences(int x) {
    /** To Generate a Set of Random Probabilities **/
    //vector<int> seqLens(3,0); seqLens[0] = 11; seqLens[1] = 11; seqLens[2] = 11;
    //Get_N_RandomSequences(3, seqLens); exit(0);

    ifstream inputstr("Sequences.txt");
    if (!inputstr) {
        cerr << "cannot open file " << "Sequences.txt" << endl;
        exit(1);
    }
    int numOfChildren;
    inputstr >> numOfChildren;

    for (int i=0; i<numOfChildren; i++) {
        double len;
        inputstr >> len;
        Lengths.push_back(len);
    }

    for (int i=0; i<numOfChildren; i++) {
        string line, token;
        vector <vector <double> > seqValues;
        for (int j=0; j<4; j++) {
            getline(inputstr, line);
            if (i==0 && j==0) getline(inputstr, line); //reads an extra line because a \n is left from before
            stringstream ss;
            ss.str("");
            ss<<line;
            vector<double> lineValues;
            while( getline(ss, token, '\t') ) {
                lineValues.push_back(atof(token.c_str()));
            }
            seqValues.push_back(lineValues);
        }
        SeqMat.push_back(seqValues);
        Dimensions.push_back(seqValues[0].size());
    }

    inputstr.close();
}

//Function to Initialize all Triplets
void SequenceAligner::InitTriplets() {
    /** Create all Possible Tihls **/
    vector<string> tihlLbls;
    tihlLbls.push_back("H");
    tihlLbls.push_back("B");
    tihlLbls.push_back("N");
    tihlLbls.push_back("E");
    tihlLbls.push_back("-");
    vector<string> Tihls;
    GetAllCombinations("", Lengths.size(), tihlLbls.size(), tihlLbls, "TihlTree", true, Tihls);

    /** Create all Reachable Soans **/
    vector<string> Soans;
    GetAllReachableSoans(Soans, Tihls);

    /** Create all Possible Triplets organized per each reachable soan **/
    GetAllPossibleTriplets(Soans, Tihls);
    /* //For Debugging
    int counter = 0;
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
    */
}

//Function to intialize all ProbabilityObject with default values
void SequenceAligner::InitProbabilities() {
    int product = 1;
    Dimensions.push_back(Triplets.size()-1);   //Because a +1 is later added through the FixPos vector
    vector<int> FixPos(Dimensions.size(), 1);
    Dimensions = Dimensions + FixPos;

    for (unsigned int i=0; i<Dimensions.size(); i++) {
        product *= (Dimensions[i]);
    }

    for (unsigned int i=0; i<product+1; i++) { //an extra state for the end state (i.e. all H-labeled with indices equal to m+1, n+1, k+1 etc (m, n, k,... being the sequences' lengths)
        ProbabilityObject* prob = new ProbabilityObject();
        if (i==0) {
            prob->Set_OverallProbability(1.0);
            prob->Set_SoanLbl("HHH");
        } else      {
            prob->Set_OverallProbability(0.0);
        }
        prob->Set_PosIn_1D(i);
        prob->Set_PosIn_ND(i, Dimensions);
        prob->Set_EndState(false);
        Probabilities.push_back(prob);
    }
}

//Function to Get a given alignment from the user
void SequenceAligner::InitAlignedSequences() {
    /** Init lengths **/
    int numOfChildren;
    cout<<"Enter the number of leaf nodes of the tree -> ";
    cin>>numOfChildren;
    cout<<endl;

    for (int i=0; i<numOfChildren; i++) {
        double len;
        cout<<"\tEnter the length of Child - " << (i+1) <<": -> ";
        cin>>len;

        Lengths.push_back(len);
    }

    cout<<endl;

    /** Init Sequences **/
    for (int i=0; i<numOfChildren; i++) {
        string seq = "";
        cout<<"\tEnter the Sequence for the " << (i+1) << " Child: -> ";
        cin>>seq;
        AlignedSequences.push_back(seq);
        ReplaceAll(seq, "-", "");
        Sequences.push_back(seq);
        Dimensions.push_back(seq.size());
    }

    InitTriplets();        //Initialize all possible Triplets (organized per each reachable soan)
    InitProbabilities();   //Initialize all ProbabilityObject with default values

    IntMat idxs_ND;
    GetIndicesFromAlignment(AlignedSequences, idxs_ND);

    for (int i=0; i<idxs_ND.size(); i++) {
        idxs_ND[i].push_back(0);
        AlignedSequencesIdxs.push_back(GetPosForIndexes(idxs_ND[i], Dimensions));
    }
}

//Function to Get a given alignment from a file
void SequenceAligner::InitAlignedSequences(string filename) {
    /** Init lengths **/
    int numOfChildren;
    ifstream inputstr(filename.c_str());
    if (!inputstr) {
        cerr << "cannot open file " << filename.c_str() << endl;
        exit(1);
    }
    inputstr >> numOfChildren;

    for (int i=0; i<numOfChildren; i++) {
        double len;
        inputstr >> len;

        Lengths.push_back(len);
    }

    /** Init Sequences **/
    for (int i=0; i<numOfChildren; i++) {
        string seq = "";
        inputstr >> seq;
        AlignedSequences.push_back(seq);
        ReplaceAll(seq, "-", "");
        Sequences.push_back(seq);
        Dimensions.push_back(seq.size());
    }

    inputstr.close();

    InitTriplets();        //Initialize all possible Triplets (organized per each reachable soan)
    InitProbabilities();   //Initialize all ProbabilityObject with default values

    IntMat idxs_ND;
    GetIndicesFromAlignment(AlignedSequences, idxs_ND);

    for (int i=0; i<idxs_ND.size(); i++) {
        idxs_ND[i].push_back(0);
        AlignedSequencesIdxs.push_back(GetPosForIndexes(idxs_ND[i], Dimensions));
    }
}

//Function to Get a given homology from the user
void SequenceAligner::InitHomologyIndices() {
    /** Init lengths **/
    int numOfChildren;
    cout<<"Enter the number of leaf nodes of the tree -> ";
    cin>>numOfChildren;
    cout<<endl;

    //the +1 is used in order to also input the root sequence
    for (int i=0; i<(numOfChildren+1); i++) {
        double len;
        if (i==0) {
            cout<<"\tEnter the length of the root: -> ";
        } else {
            cout<<"\tEnter the length of Child - " << (i+1) <<": -> ";
        }
        cin>>len;

        Lengths.push_back(len);
    }

    cout<<endl;

    /** Init Sequences **/
    //the +1 is used in order to also input the root sequence
    for (int i=0; i<(numOfChildren+1); i++) {
        string seq = "";
        if (i==0) {
            cout<<"\tEnter the Sequence for the  root: -> ";
        } else {
            cout<<"\tEnter the Sequence for the " << (i+1) << " Child: -> ";
        }
        cin>>seq;
        Sequences.push_back(seq);
        Dimensions.push_back(seq.size());
    }

    int rootLength = Sequences[0].length();
    HomologiesIdxs.resize(Sequences.size()-1, vector<int>(rootLength, -1));
    Anchors.resize(rootLength, vector<int>(Sequences.size(), -1));
    for (int i=0; i<numOfChildren; i++) {
        cout<<"\tFor sequence-"<<(i+1)<<":"<<endl;
        for (int j=0; j<rootLength; j++) {
            int idx;
            cout<<"\t\tPos-"<<(j+1)<<" is Homologous which pos of seq-"<<(i+1)<<" ? -> ";
            cin>>idx;
            HomologiesIdxs[i][j] = idx;
            Anchors[j][0] = j+1;
            Anchors[j][i+1] = idx;
        }
    }
}

//Function to Get a given homology from a file
void SequenceAligner::InitHomologyIndices(string filename) {
    /** Init lengths **/
    int numOfChildren;
    ifstream inputstr(filename.c_str());
    if (!inputstr) {
        cerr << "cannot open file " << filename.c_str() << endl;
        exit(1);
    }
    inputstr >> numOfChildren;

    //the +1 is used in order to also input the root sequence
    for (int i=0; i<(numOfChildren+1); i++) {
        double len;
        inputstr >> len;

        Lengths.push_back(len);
    }

    /** Init Sequences **/
    //the +1 is used in order to also input the root sequence
    for (int i=0; i<(numOfChildren+1); i++) {
        string seq = "";
        inputstr >> seq;
        Sequences.push_back(seq);
        Dimensions.push_back(seq.size());
    }

    int rootLength = Sequences[0].length();
    HomologiesIdxs.resize(Sequences.size()-1, vector<int>(rootLength, -1));
    Anchors.resize(rootLength, vector<int>(Sequences.size(), -1));
    for (int i=0; i<numOfChildren; i++) {
        for (int j=0; j<rootLength; j++) {
            int idx;
            inputstr>>idx;
            HomologiesIdxs[i][j] = idx;
            Anchors[j][0] = j+1;
            Anchors[j][i+1] = idx;
        }
    }

    inputstr.close();
}

//}

/** Assisting Functions for In-Between calculations **/ //{
//Function to Create All Possible, based on all possible tihls and all reachable soans
void SequenceAligner::GetAllPossibleTriplets(const vector<string>& soans, const vector<string>& tihls) {
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
            trip->Set_SoanAfter(ApplyTihlToSoan(soans[i], tihls[j]));
            trip->Set_Probability(Fid, Lengths);
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

//Function to Create All Reachable Soans, based on all possible tihls and an initial (all H-labeled) soan
void SequenceAligner::GetAllReachableSoans(vector<string>& soans, const vector<string>& tihls) {
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
            string resultingSoan = ApplyTihlToSoan(topSoan, tihls[i]);
            //cout<<topSoan<<"\t"<<tihls[i]<<"\t"<<resultingSoan<<endl;
            if ( std::find(soans.begin(), soans.end(), resultingSoan) != soans.end()) {
                continue;
            }

            soans.push_back(resultingSoan);
            soansQueue.push(resultingSoan);
        }
    }
}

//Function that examines whether a given tihl may be applied to a given soan or not
bool SequenceAligner::CanTihlBeAppliedToSoan(string soan, string tihl) {
    for (unsigned int i=0; i<soan.size(); i++) {
        if ( (soan.at(i)=='e' || soan.at(i)=='h' || soan.at(i)=='b') && (tihl.at(i)=='B'))
            return false;
    }
    return true;
}

//Function that Applies a given Tihl to a given Soan
string SequenceAligner::ApplyTihlToSoan(string soan, string tihl) {
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

//Plain Sequences: Function to get the Corresponding Substitution Probability for the given transition to the next state
double SequenceAligner::GetSubstitutionProbability(string tihlLbl, const vector<int>& Idxs) const {
    double subProb = 1;
    for (int i=0; i<tihlLbl.size(); i++) {
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
                for (int k=0; k<H_Idxs.size(); k++) { //Parse every H position
                    tempProd *= SubMats[H_Idxs[k]][GetAlphabetPair(GetAlphabet(j), GetAlphabet(Sequences[H_Idxs[k]].at(Idxs[H_Idxs[k]]-1)))];
                }
                temp += tempProd;
            }

            subProb *= temp;

            for (int k=0; k<H_Idxs.size(); k++) {
                tihlLbl.at(H_Idxs[k]) = '-';
            }

            break;
        }
        case 'B':
        case 'N':
            subProb *= Params->Get_BaseFreq()[GetAlphabet(Sequences[i].at(Idxs[i]-1))];
            break;
        }
    }
    return subProb;
}

//Felsenstein Tables: Overloaded Function to get the Corresponding Substitution Probability for the given transition to the next state
double SequenceAligner::GetSubstitutionProbability(string tihlLbl, const vector<int>& Idxs, int x) const {
    double subProb = 1;
    for (int i=0; i<tihlLbl.size(); i++) {
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
                for (int k=0; k<H_Idxs.size(); k++) { //Parse every H position
                    double tempSum = 0;
                    for (int l=0; l<4; l++) { //For every nucleotide A,C,G,T
                        tempSum += SubMats[H_Idxs[k]][GetAlphabetPair(GetAlphabet(j), GetAlphabet(l))] * SeqMat[H_Idxs[k]][l][Idxs[H_Idxs[k]]-1];
                        //For Debugging
                        //cout<<endl<<"P(l"<<(H_Idxs[k]+1)<<")"<<GetAlphabet_String(j) <<"->"<<GetAlphabet_String(l) <<": "<< SubMats[H_Idxs[k]][GetAlphabetPair(GetAlphabet(j), GetAlphabet(l))]<<endl;
                        //cout<<"      p"<<GetAlphabet_String(l) << Idxs[H_Idxs[k]] <<": "<< SeqMat[H_Idxs[k]][l][Idxs[H_Idxs[k]]-1]<<endl;
                        //Pause();
                    }
                    tempProd *= tempSum;
                }
                temp += tempProd;
            }

            subProb *= temp;

            for (int k=0; k<H_Idxs.size(); k++) {
                tihlLbl.at(H_Idxs[k]) = '-';
            }

            break;
        }
        case 'B':
        case 'N': {
            double tempSum = 0;

            for (int j=0; j<4; j++) { //For every nucleotide A,C,G,T
                tempSum += Params->Get_BaseFreq()[j] * SeqMat[i][j][Idxs[i]-1];
            }
            subProb *= tempSum;
            break;
        }
        }
    }
    return subProb;
}

//Function to return the triplet corresponding to the given index, that leads to the non emmitting state
Triplet* SequenceAligner::Get_ToNonEmmitingTriplet(int r_idx) const {
    for (int i=0; i<Triplets[r_idx].size(); i++) {
        if ( StringCount(Triplets[r_idx][i]->Get_SoanAfter(), "e")==Triplets[r_idx][i]->Get_SoanAfter().size() ) return Triplets[r_idx][i];
    }
    return NULL;
}

//Function to return the triplet that leaves from the non emmitting state and transits to a soan labeled according to the given argument
Triplet* SequenceAligner::Get_FromNonEmmitingTriplet(string sAfter, const vector<int>& idxs) const {
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
int SequenceAligner::IndexOfSoan(string soan) {
    for (unsigned int i=0; i<Triplets.size(); i++) {
        if ( Triplets[i][0]->Get_SoanBefore().compare(soan) == 0 )return i;
    }
    return -1;
}

//Function to get the corresponding representation of the given alignment in the form of sequence indices
void SequenceAligner::GetIndicesFromAlignment(vector<string>& alignment, IntMatRef idxs) {
    //Convert Alignment to Matrix of Indices
    idxs.clear();
    idxs.resize(alignment[0].size(), vector<int> ( alignment.size(),0));
    for (int i=0; i<alignment.size(); i++) {
        int sitesCounter = 0;
        for (int j=0; j<alignment[i].size(); j++) {
            if (alignment[i].at(j)!='-') {
                sitesCounter++;
            }
            idxs[j][i] = sitesCounter;
        }
    }
}

//Function to get the corresponding representation of the given evolutionary history in the form of sequence indices
void SequenceAligner::GetIndicesFromEvolHistory(vector<string>& evolHist, IntMatRef idxs) {
    //Convert Alignment to Matrix of Indices
    idxs.clear();
    idxs.resize(evolHist[0].size(), vector<int> ( evolHist.size(),0));
    for (int i=0; i<evolHist.size(); i++) {
        int sitesCounter = 0;
        for (int j=0; j<evolHist[i].size(); j++) {
            if (evolHist[i].at(j)!='-' && evolHist[i].at(j)!='E') {
                sitesCounter++;
            }
            idxs[j][i] = sitesCounter;
        }
    }
}

//Function that checks wheteher a Transition is Valid according to the given Homologies
bool SequenceAligner::IsValidTransition(IntVecRef NextIdxs, string tihl) {
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
                if (j==k || Anchors[pos][k]< NextIdxs[k] || (Anchors[pos][k]==NextIdxs[k] && (tihl.at(k)=='B' || tihl.at(k)=='-')) ){
                    continue;
                }

                if (Anchors[pos][k]!= NextIdxs[k] && Anchors[pos][k]!=-1) {
                    return false;
                }else if (Anchors[pos][k]!= -1 && tihl.at(k)!='H') {
                    return false;
                }
            }
        }else {        //NextIdxs[j] corresponds to -1 in the Anchors vector
            if (tihl.at(j)=='H') {
                return false;
            }
        }
    }

    return true;
}

//}

/** Probabilities Estimation Functions **/ //{
//Function to calculate all the probabilities describing all evolutionary histories for the alignement of the given sequences. Argument defines which MakeTransition Function should be called
void SequenceAligner::EstimateProbabilities(int choice) {
    WriteToFile("Log.txt", "", false, true);

    for (unsigned int i=0; i<Probabilities.size()-1-Dimensions[Dimensions.size()-1]; i++) {
        if (Probabilities[i]->Get_OverallProbability()==0) continue;

        int CurPosIn_1D = i;
        vector<int> CurIdxs = Probabilities[i]->Get_PosIn_ND();

        /** 2. Transit to the next state **/
        switch (choice) {
        case 0: //The plain Sequences case
            MakeTransition(CurIdxs, CurPosIn_1D);
            break;

        case 1: //The Felsenstein Tables Case
            MakeTransition(CurIdxs, CurPosIn_1D, 5);
            break;

        case 2: //The Given Alignment Case
            MakeTransition(CurIdxs, CurPosIn_1D, 5.0f);
            break;

        case 3: //The Given Homologies Case
            MakeTransition(CurIdxs, CurPosIn_1D, true);
            break;

        }
    }
}

//Plain Sequences: Function to make the transition from current state to the next state
void SequenceAligner::MakeTransition(vector<int>& CurIdxs, int CurPosIn_1D) {
    /** 2. Start Applying each possible tihl to this state **/
    Triplet* nonEmittingTriplet = Triplets[idx_R][idx_C];

    for (unsigned int i=0; i<Triplets[CurIdxs[CurIdxs.size()-1]].size(); i++) {
        Triplet* currentTriplet         = Triplets[CurIdxs[CurIdxs.size()-1]][i];
        if (StringCount(currentTriplet->Get_SoanAfter(), "e")==currentTriplet->Get_SoanAfter().size()) {
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
        vector<int> NextIdxs(CurIdxs);
        NextIdxs += currentTriplet->GetTripletIdxChanges();
        NextIdxs[NextIdxs.size()-1] = IndexOfSoan(currentTriplet->Get_SoanAfter());
        int NextPosIn_1D = GetPosForIndexes(NextIdxs, Dimensions);

        //2b. Check if this transition leads to creating sequence(s) larger than the ones provided
        //If so, then skip this tihl
        bool LeadsToLargerSequence = false;

        for (unsigned int j=0; j<NextIdxs.size()-1; j++) {
            if (NextIdxs[j]>=(Dimensions[j])) {
                LeadsToLargerSequence = true;    //i.e. this tihl leads to the emission of more sites for the i-th sequence than those that the i-th sequence has
                break;
            }
        }
        if (LeadsToLargerSequence) continue;

        //2c. Add a new contributor to the ProbabilityObject of the Next State (defined by NextIdxs)
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
            contributedProbability = (Probabilities[CurPosIn_1D]->Get_OverallProbability() * currentTriplet->Get_Probability() ) *
                                     ( GetSubstitutionProbability(currentTriplet->Get_Tihl(), NextIdxs) )  ;
            Probabilities[NextPosIn_1D]->AddContributorWingFold(false);
            Probabilities[NextPosIn_1D]->AddContributorDirectTransProb(contributedProbability);
            Probabilities[NextPosIn_1D]->AddContributorSelfTransProb(0.0f);
        } else {
            /*contributedProbability = (
                                      ( (Probabilities[CurPosIn_1D]->Get_OverallProbability() * currentTriplet->Get_Probability() )
                                        +
                                        (Probabilities[CurPosIn_1D]->Get_OverallProbability() * toNonEmmitingTriplet->Get_Probability() * (fromNonEmmitingTriplet->Get_Probability()/(1-nonEmittingTriplet->Get_Probability())) )
                                      ) * ( GetSubstitutionProbability(currentTriplet->Get_Tihl(), NextIdxs) )
                                     );
            */
            double subProb = ( GetSubstitutionProbability(currentTriplet->Get_Tihl(), NextIdxs) );
            double direct = (Probabilities[CurPosIn_1D]->Get_OverallProbability() * currentTriplet->Get_Probability() ) * subProb;
            double indirect = (Probabilities[CurPosIn_1D]->Get_OverallProbability() * toNonEmmitingTriplet->Get_Probability() * (fromNonEmmitingTriplet->Get_Probability()/(1-nonEmittingTriplet->Get_Probability())) ) * subProb;
            contributedProbability = direct + indirect;

            Probabilities[NextPosIn_1D]->AddContributorWingFold(true);
            Probabilities[NextPosIn_1D]->AddContributorDirectTransProb(direct);
            Probabilities[NextPosIn_1D]->AddContributorSelfTransProb(nonEmittingTriplet->Get_Probability());
        }

        Probabilities[NextPosIn_1D]->AddContributorProbability(contributedProbability);
        Probabilities[NextPosIn_1D]->Set_PosIn_1D(NextPosIn_1D);
        Probabilities[NextPosIn_1D]->Set_PosIn_ND(NextIdxs);
        Probabilities[NextPosIn_1D]->IncreaseOverallProbabilityBy(contributedProbability);
        Probabilities[NextPosIn_1D]->Set_SoanLbl(currentTriplet->Get_SoanAfter());

        //2d. Examine if the Next State is a (pre)End State. Although the actual end state will be the one with all edges H(omomlogous)
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
    }
}

//Felsenstein Tables: Overloaded Function to make the transition from current state to the next state
void SequenceAligner::MakeTransition(vector<int>& CurIdxs, int CurPosIn_1D, int x) {
    /** 2. Start Applying each possible tihl to this state **/
    Triplet* nonEmittingTriplet = Triplets[idx_R][idx_C];

    for (unsigned int i=0; i<Triplets[CurIdxs[CurIdxs.size()-1]].size(); i++) {
        Triplet* currentTriplet         = Triplets[CurIdxs[CurIdxs.size()-1]][i];
        if (StringCount(currentTriplet->Get_SoanAfter(), "e")==currentTriplet->Get_SoanAfter().size()) {
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
        vector<int> NextIdxs(CurIdxs);
        NextIdxs += currentTriplet->GetTripletIdxChanges();
        NextIdxs[NextIdxs.size()-1] = IndexOfSoan(currentTriplet->Get_SoanAfter());
        int NextPosIn_1D = GetPosForIndexes(NextIdxs, Dimensions);

        //2b. Check if this transition leads to creating sequence(s) larger than the ones provided
        //If so, then skip this tihl
        bool LeadsToLargerSequence = false;

        for (unsigned int j=0; j<NextIdxs.size()-1; j++) {
            if (NextIdxs[j]>=(Dimensions[j])) {
                LeadsToLargerSequence = true;    //i.e. this tihl leads to the emission of more sites for the i-th sequence than those that the i-th sequence has
                break;
            }
        }
        if (LeadsToLargerSequence) continue;

        //2c. Add a new contributor to the ProbabilityObject of the Next State (defined by NextIdxs)
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
            contributedProbability = (Probabilities[CurPosIn_1D]->Get_OverallProbability() * currentTriplet->Get_Probability() ) *
                                     ( GetSubstitutionProbability(currentTriplet->Get_Tihl(), NextIdxs, x) )  ;
            Probabilities[NextPosIn_1D]->AddContributorWingFold(false);
            Probabilities[NextPosIn_1D]->AddContributorDirectTransProb(contributedProbability);
            Probabilities[NextPosIn_1D]->AddContributorSelfTransProb(0.0f);
        } else {
            /*contributedProbability = (
                                      ( (Probabilities[CurPosIn_1D]->Get_OverallProbability() * currentTriplet->Get_Probability() )
                                        +
                                        (Probabilities[CurPosIn_1D]->Get_OverallProbability() * toNonEmmitingTriplet->Get_Probability() * (fromNonEmmitingTriplet->Get_Probability()/(1-nonEmittingTriplet->Get_Probability())) )
                                      ) * ( GetSubstitutionProbability(currentTriplet->Get_Tihl(), NextIdxs, x) )
                                     );
            */
            double subProb = ( GetSubstitutionProbability(currentTriplet->Get_Tihl(), NextIdxs, x) );
            double direct = (Probabilities[CurPosIn_1D]->Get_OverallProbability() * currentTriplet->Get_Probability() ) * subProb;
            double indirect = (Probabilities[CurPosIn_1D]->Get_OverallProbability() * toNonEmmitingTriplet->Get_Probability() * (fromNonEmmitingTriplet->Get_Probability()/(1-nonEmittingTriplet->Get_Probability())) ) * subProb;
            contributedProbability = direct + indirect;

            Probabilities[NextPosIn_1D]->AddContributorWingFold(true);
            Probabilities[NextPosIn_1D]->AddContributorDirectTransProb(direct);
            Probabilities[NextPosIn_1D]->AddContributorSelfTransProb(nonEmittingTriplet->Get_Probability());
        }

        Probabilities[NextPosIn_1D]->AddContributorProbability(contributedProbability);
        Probabilities[NextPosIn_1D]->Set_PosIn_1D(NextPosIn_1D);
        Probabilities[NextPosIn_1D]->Set_PosIn_ND(NextIdxs);
        Probabilities[NextPosIn_1D]->IncreaseOverallProbabilityBy(contributedProbability);
        Probabilities[NextPosIn_1D]->Set_SoanLbl(currentTriplet->Get_SoanAfter());

        //2d. Examine if the Next State is a (pre)End State. Although the actual end state will be the one with all edges H(omomlogous)
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
    }
}

//Given Alignment: Function to make the transition from current state to the next state
void SequenceAligner::MakeTransition(vector<int>& CurIdxs, int CurPosIn_1D, float x) {
    /** Log Visiting States */
    stringstream log;
    log.str("");
    log<<endl<<"Moving from State: [";
    for (int i=0; i<CurIdxs.size(); i++) {
        if (i<CurIdxs.size()-1) {
            log<<left<<setw(3)<<CurIdxs[i]<<",";
        } else {
            log<<left<<setw(3)<<CurIdxs[i];
        }
    }
    log<<"] - "<<setw(4)<<CurPosIn_1D<<" : "<<endl;

    /** 2. Start Applying each possible tihl to this state **/
    Triplet* nonEmittingTriplet = Triplets[idx_R][idx_C];

    for (unsigned int i=0; i<Triplets[CurIdxs[CurIdxs.size()-1]].size(); i++) {
        Triplet* currentTriplet         = Triplets[CurIdxs[CurIdxs.size()-1]][i];

        //Skip state that leads to alignment that needs more columns that then given ones
        if (currentTriplet->GetTripletNumberOfAlignmentColumns()>1) continue;

        //Skip non-emitting state
        if (StringCount(currentTriplet->Get_SoanAfter(), "e")==currentTriplet->Get_SoanAfter().size()) {
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
        vector<int> NextIdxs(CurIdxs);
        NextIdxs += currentTriplet->GetTripletIdxChanges();
        NextIdxs[NextIdxs.size()-1] = IndexOfSoan(currentTriplet->Get_SoanAfter());
        int NextPosIn_1D = GetPosForIndexes(NextIdxs, Dimensions);


        //2b. Check if the indices of the current transition are within the acceptable (as defined by the given
        //alignment range. If not, then skip this tihl
        bool FoundWithinRange = false;
        for (int j=0; j<AlignedSequencesIdxs.size(); j++) {
            if ( AlignedSequencesIdxs[j]<=NextPosIn_1D &&  NextPosIn_1D<(AlignedSequencesIdxs[j]+Triplets.size()) ) {
                FoundWithinRange = true;
                break;
            }
        }
        if (!FoundWithinRange) continue;


        //2c. Check if this transition leads to creating sequence(s) larger than the ones provided
        //If so, then skip this tihl
        bool LeadsToLargerSequence = false;

        for (unsigned int j=0; j<NextIdxs.size()-1; j++) {
            if (NextIdxs[j]>=(Dimensions[j])) {
                LeadsToLargerSequence = true;    //i.e. this tihl leads to the emission of more sites for the i-th sequence than those that the i-th sequence has
                break;
            }
        }
        if (LeadsToLargerSequence) continue;

        //2d. Add a new contributor to the ProbabilityObject of the Next State (defined by NextIdxs)
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
            contributedProbability = (Probabilities[CurPosIn_1D]->Get_OverallProbability() * currentTriplet->Get_Probability() ) *
                                     ( GetSubstitutionProbability(currentTriplet->Get_Tihl(), NextIdxs) )  ;
            Probabilities[NextPosIn_1D]->AddContributorWingFold(false);
            Probabilities[NextPosIn_1D]->AddContributorDirectTransProb(contributedProbability);
            Probabilities[NextPosIn_1D]->AddContributorSelfTransProb(0.0f);
        } else {
            double subProb = ( GetSubstitutionProbability(currentTriplet->Get_Tihl(), NextIdxs) );
            double direct = (Probabilities[CurPosIn_1D]->Get_OverallProbability() * currentTriplet->Get_Probability() ) * subProb;
            double indirect = (Probabilities[CurPosIn_1D]->Get_OverallProbability() * toNonEmmitingTriplet->Get_Probability() * (fromNonEmmitingTriplet->Get_Probability()/(1-nonEmittingTriplet->Get_Probability())) ) * subProb;
            contributedProbability = direct + indirect;

            Probabilities[NextPosIn_1D]->AddContributorWingFold(true);
            Probabilities[NextPosIn_1D]->AddContributorDirectTransProb(direct);
            Probabilities[NextPosIn_1D]->AddContributorSelfTransProb(nonEmittingTriplet->Get_Probability());
        }

        Probabilities[NextPosIn_1D]->AddContributorProbability(contributedProbability);
        Probabilities[NextPosIn_1D]->Set_PosIn_1D(NextPosIn_1D);
        Probabilities[NextPosIn_1D]->Set_PosIn_ND(NextIdxs);
        Probabilities[NextPosIn_1D]->IncreaseOverallProbabilityBy(contributedProbability);
        Probabilities[NextPosIn_1D]->Set_SoanLbl(currentTriplet->Get_SoanAfter());

        log<<"\tto State: [";
        for (int i=0; i<NextIdxs.size(); i++) {
            if (i<NextIdxs.size()-1) {
                log<<left<<setw(3)<<NextIdxs[i]<<",";
            } else {
                log<<left<<setw(3)<<NextIdxs[i];
            }
        }
        log<<"] - "<<setw(4)<<NextPosIn_1D<<" via Tihl: " << Probabilities[NextPosIn_1D]->GetContributorTihlExtAtPos(Probabilities[NextPosIn_1D]->Get_Contributors().size()-1)<<endl;

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
    }
    WriteToFile("Log.txt", log.str(), true, false);
}

//Given Homologies: Function to make the transition from current state to the next state
void SequenceAligner::MakeTransition(vector<int>& CurIdxs, int CurPosIn_1D, bool x) {
    /** Log Visiting States */
    stringstream log;
    log.str("");
    log<<endl<<"Moving from State: [";
    for (int i=0; i<CurIdxs.size(); i++) {
        if (i<CurIdxs.size()-1) {
            log<<left<<setw(3)<<CurIdxs[i]<<",";
        } else {
            log<<left<<setw(3)<<CurIdxs[i];
        }
    }
    log<<"] - "<<setw(4)<<CurPosIn_1D << " - " << Probabilities[CurPosIn_1D]->Get_SoanLbl() <<" : "<<endl;

    /** 2. Start Applying each possible tihl to this state **/
    Triplet* nonEmittingTriplet = Triplets[idx_R][idx_C];

    for (unsigned int i=0; i<Triplets[CurIdxs[CurIdxs.size()-1]].size(); i++) {
        Triplet* currentTriplet         = Triplets[CurIdxs[CurIdxs.size()-1]][i];

        //Skip tihls that do not start with 'H' or '-' -> i.e. for the root node to either have a B or not
        if (!StartsWith(currentTriplet->Get_Tihl(), "H") && !StartsWith(currentTriplet->Get_Tihl(), "-")) {
            continue;
        }

        //Skip non-emitting state
        if (StringCount(currentTriplet->Get_SoanAfter(), "e")==currentTriplet->Get_SoanAfter().size()) {
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
        int NextPosIn_1D = GetPosForIndexes(NextIdxs, Dimensions);

        //2b. Check if the NextIdxs are in accordance with the given H(omologies) <=> Anchor Points
        if (!IsValidTransition(NextIdxs, currentTriplet->Get_Tihl())) {
            continue;
        }

        //2c. Check if this transition leads to creating sequence(s) larger than the ones provided
        //If so, then skip this tihl
        bool LeadsToLargerSequence = false;

        for (unsigned int j=0; j<NextIdxs.size()-1; j++) {
            if (NextIdxs[j]>=(Dimensions[j])) {
                LeadsToLargerSequence = true;    //i.e. this tihl leads to the emission of more sites for the i-th sequence than those that the i-th sequence has
                break;
            }
        }
        if (LeadsToLargerSequence) continue;

        //2d. Add a new contributor to the ProbabilityObject of the Next State (defined by NextIdxs)
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
            contributedProbability = (Probabilities[CurPosIn_1D]->Get_OverallProbability() * currentTriplet->Get_Probability() ) *
                                     ( GetSubstitutionProbability(currentTriplet->Get_Tihl(), NextIdxs) )  ;
            Probabilities[NextPosIn_1D]->AddContributorWingFold(false);
            Probabilities[NextPosIn_1D]->AddContributorDirectTransProb(contributedProbability);
            Probabilities[NextPosIn_1D]->AddContributorSelfTransProb(0.0f);
        } else {
            double subProb = ( GetSubstitutionProbability(currentTriplet->Get_Tihl(), NextIdxs) );
            double direct = (Probabilities[CurPosIn_1D]->Get_OverallProbability() * currentTriplet->Get_Probability() ) * subProb;
            double indirect = (Probabilities[CurPosIn_1D]->Get_OverallProbability() * toNonEmmitingTriplet->Get_Probability() * (fromNonEmmitingTriplet->Get_Probability()/(1-nonEmittingTriplet->Get_Probability())) ) * subProb;
            contributedProbability = direct + indirect;

            Probabilities[NextPosIn_1D]->AddContributorWingFold(true);
            Probabilities[NextPosIn_1D]->AddContributorDirectTransProb(direct);
            Probabilities[NextPosIn_1D]->AddContributorSelfTransProb(nonEmittingTriplet->Get_Probability());
        }

        Probabilities[NextPosIn_1D]->AddContributorProbability(contributedProbability);
        Probabilities[NextPosIn_1D]->Set_PosIn_1D(NextPosIn_1D);
        Probabilities[NextPosIn_1D]->Set_PosIn_ND(NextIdxs);
        Probabilities[NextPosIn_1D]->IncreaseOverallProbabilityBy(contributedProbability);
        Probabilities[NextPosIn_1D]->Set_SoanLbl(currentTriplet->Get_SoanAfter());

        log<<"\t     to State: [";
        for (int i=0; i<NextIdxs.size(); i++) {
            if (i<NextIdxs.size()-1) {
                log<<left<<setw(3)<<NextIdxs[i]<<",";
            } else {
                log<<left<<setw(3)<<NextIdxs[i];
            }
        }
        log<<"] - "<<setw(4)<<NextPosIn_1D << " - " << Probabilities[NextPosIn_1D]->Get_SoanLbl()<<" via Tihl: " << Probabilities[NextPosIn_1D]->GetContributorTihlExtAtPos(Probabilities[NextPosIn_1D]->Get_Contributors().size()-1)<<endl;

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
    }
    WriteToFile("Log.txt", log.str(), true, false);
}

//Function to make the final transition to the End State (i.e. the all H(omologous) edges labeled state)
void SequenceAligner::TransitToEndState() {
    int counter = 0;
    for (unsigned int i=Probabilities.size()-1-Dimensions[Dimensions.size()-1] ; i<Probabilities.size()-1; i++) {
        if (Probabilities[i]->Get_OverallProbability()==0) continue;

        Triplet* endStateTriplet = Triplets[counter++][0];

        double contributedProbability = 1 * Probabilities[i]->Get_OverallProbability();
        Probabilities[Probabilities.size()-1]->AddContributor(Probabilities[i]);
        Probabilities[Probabilities.size()-1]->AddContributorTihl(endStateTriplet->Get_Tihl());
        Probabilities[Probabilities.size()-1]->AddContributorTihlExt(endStateTriplet->Get_Tihl());
        Probabilities[Probabilities.size()-1]->AddContributorProbability(contributedProbability);
        Probabilities[Probabilities.size()-1]->Set_PosIn_1D(Probabilities.size()-1);
        Probabilities[Probabilities.size()-1]->Set_PosIn_ND(GetIndexesForPos(Probabilities.size()-1, Dimensions));
        Probabilities[Probabilities.size()-1]->IncreaseOverallProbabilityBy(1 * Probabilities[i]->Get_OverallProbability());
        Probabilities[Probabilities.size()-1]->Set_SoanLbl(endStateTriplet->Get_SoanAfter());
        Probabilities[Probabilities.size()-1]->AddContributorWingFold(false);
        Probabilities[Probabilities.size()-1]->AddContributorDirectTransProb(contributedProbability);
        Probabilities[Probabilities.size()-1]->AddContributorSelfTransProb(0.0f);


    }
}

//}

/** Functions to Sample Alignments **/ //{
//Function to sample evolutionary histories based on the calculated Probabilities
void SequenceAligner::GetEvolHistories(StringMatRef evolHist, vector<double>& probs) {
    //Init the random Generator
    BoostRandomGenerator rand;
    rand.SetSeed(-1); //to seed with random values

    //For each one of the possible evolution histories....
    for (int i=0; i<Params->Get_SampleSize(); i++) {
        double prob = 0;

        //Start from the End state...
        ProbabilityObject* state = Probabilities[Probabilities.size()-1];

        //Init the string vector representing the evolution history of each sequence
        vector<string> hist(Dimensions.size()-1, "");

        //Move Backwards until the start state has been reached
        while(state != Probabilities[0]) {
            //Get a random number and select a contributor according to this random value...
            double randProb = rand.GetRandomDouble(0.0f, 1.0f)*state->Get_OverallProbability();
            //double randProb = rand.GetRandomDouble(0, state->Get_OverallProbability());
            int j;
            for (j=0; randProb>=0; j++) {
                randProb -= state->GetContributorProbabilityAtPos(j);
            }
            j--;

            //To skip trying to form an alignment for the transition to the end state...
            if (state == Probabilities[Probabilities.size()-1]) {
                state = state->GetContributorAtPos(j);
                continue;
            }

            //Get the contributor's tihl and form corresponding evolution history and alignment for each of the sequences
            string tihl = state->GetContributorTihlAtPos(j);

            //Sample non-emitting states
            randProb = rand.GetRandomDouble(0.0f, 1.0f)*state->GetContributorProbabilityAtPos(j);

            if (state->GetContributorWingFoldAtPos(j) && randProb > state->GetContributorDirectTransProbAtPos(j)) {
                prob += log(state->GetContributorProbabilityAtPos(j)-state->GetContributorDirectTransProbAtPos(j));
                for (int j=0; j<tihl.size(); j++) {
                    stringstream st;
                    st.str("");
                    st << "E" << hist[j];
                    hist[j] = st.str();
                }

                while(rand.GetRandomDouble(0.0f, 1.0f) < state->GetContributorSelfTransProbAtPos(j)) {
                    prob += log(state->GetContributorSelfTransProbAtPos(j));
                    for (int j=0; j<tihl.size(); j++) {
                        stringstream st;
                        st.str("");
                        st << "E" << hist[j];
                        hist[j] = st.str();
                    }
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
            prob += log(state->GetContributorProbabilityAtPos(j));
            state = state->GetContributorAtPos(j);
        }
        evolHist.push_back(hist);
        probs.push_back(prob);
    }
}

//Function that returns the corresponding alignment for the given evolutionary history
void SequenceAligner::GetAlignmentForEvolHistory(vector<string >& evolHist, vector<string >& align) {
    stringstream st;
    for (int i=0; i<evolHist.size(); i++) {
        int pos = 0;
        st.str("");
        for (int j=0; j<evolHist[i].size(); j++) {
            switch(evolHist[i].at(j)) {
            case 'H':
            case 'N':
            case 'B':
                st << Sequences[i].at(pos++);
                break;
            case 'E':
            case '-':
                st << "-";
                break;
            }
        }
        align.push_back(st.str());
    }
}

//Function to assess the probability of a given alignment based on the computed Probability Objects
void SequenceAligner::AssessAlignment() {
    vector<string> givenAligns(Lengths.size(), "");
    cout<<endl<<endl;

    for (int i=0; i<Lengths.size(); i++) {
        cout<<"Enter alignment for the " << (i+1) <<" sequence... -> ";
        cin>>givenAligns[i];
    }
    cout<<"\tAlign.Log-Likelihood: "<<AssessAlignment(givenAligns)<<endl;
}

//Function to assess the probability of a given alignment based on the computed Probability Objects
double SequenceAligner::AssessAlignment(vector<string>& alignment) {
    //Convert Alignment to Matrix of Indices
    vector< vector<int> > idxs;
    GetIndicesFromAlignment(alignment, idxs);

    vector<double> probs(idxs.size(), 0.0f);
    for (int i=0; i<idxs.size(); i++) {
        for (int j=0; j<Triplets.size(); j++) {
            vector<int> tempIdx = idxs[i];
            tempIdx.push_back(j);
            int pos = GetPosForIndexes(tempIdx, Dimensions);
            probs[i] += Probabilities[pos]->Get_OverallProbability();
        }
    }

    double overallProbability = 0.0f;
    for (int i=0; i<probs.size(); i++) {
        overallProbability += log(probs[i]);
    }

    return overallProbability;
}

//}

/** Functions to Output results **/ //{
//Function to save derived results to file
void SequenceAligner::SaveResultsToFile() {
    WriteToFile("Probs.txt", "", true, true);
    for (unsigned int i=0; i<Probabilities.size(); i++) {
        if (Probabilities[i]->Get_OverallProbability()>0) {
            stringstream ss;
            ss.str("");
            ss<<*(Probabilities[i]);
            WriteToFile("Probs.txt", ss.str(), true, false);
        } else if (Probabilities[i]->Get_OverallProbability()<0) {
            cout<<"Error - probability " << i << " has a NEGATIVE value!"<<endl;
        }
    }
}

//Function to print to screeen all Sampled Alignments alongside with their corresponding evolutionary histories
void SequenceAligner::ShowAllSampledAlignments(StringMatRef evolHistories, vector<double>& probs) {
    cout<<endl;

    for (int i=0; i<evolHistories.size(); i++) {
        vector<string> align;
        GetAlignmentForEvolHistory(evolHistories[i], align);

        //Uncomment to show only histories containing non-emitting states
        //for (int k=0;k<evolHistories[i][0].size();k++){
        //    if (evolHistories[i][0].at(k)=='E' && evolHistories[i][1].at(k)=='E' && evolHistories[i][2].at(k)=='E' ) {

        cout<<"Sample: " << (i+1) <<endl;
        for (int j=0; j<evolHistories[i].size(); j++) {
            cout<<"\t"<<align[j] << "\t" << evolHistories[i][j] <<endl;
        }
        cout<< "\t\tProb:"<<probs[i]<<"\n\n";

        if ( ((i+1) %10) == 0) Pause();
        // } }
    }
}

//Function to output the initial probabilities of the Sequences
void SequenceAligner::PrintSequenceProbabilities() {
    cout<<"\t\tSEQUENCE PROBABILITIES"<<endl;
    cout<<"\t\t======================"<<endl;

    for (int i=0; i<SeqMat.size(); i++) {
        cout<<"Seq-"<<(i+1)<<", of Length :"<< Lengths[i] <<", # of Sites: " << SeqMat[i][0].size() <<endl;
        for (int j=0; j<SeqMat[i].size(); j++) {
            cout<<setw(5)<<GetAlphabet_String(j)<<" ";
            for (int k=0; k<SeqMat[i][j].size(); k++) {
                cout<<setw(5)<<SeqMat[i][j][k]<<"  ";
            }
            cout<<endl;
        }
        cout<<endl;
    }
}

//Function to output to the screen the Values of the Substitution matrices according to the substitution model and the corresponding branch lengths
void SequenceAligner::PrintSubstitutionMatrices() {
    cout<<endl<<endl<<"\t\tSUBSTITUTION MATRICES"<<endl;
    cout<<"\t\t====================="<<endl;

    for (int i=0; i<SubMats.size(); i++) {
        cout<<"For Seq-"<<(i+1)<< ", of Length: "<<Lengths[i]<<endl;
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

//Function to return a string with the vector's values in a printable for (e.g. <1 2 3>)
string SequenceAligner::GetVectorInPrintableForm(string description, IntVecRef vec) {
    stringstream ss;
    ss.str("");

    ss<<description<<": <";
    for (int i=0; i<vec.size(); i++) {
        if (i<vec.size()-1) {
            ss<<setw(3)<<vec[i]<<" ";
        } else {
            ss<<setw(3)<<vec[i];
        }
    }
    ss<<">";

    return ss.str();
}

//}

//Definition of the toString Operation
ostream& operator<<(ostream& ostr, const SequenceAligner& aligner) {
    ostr << "Not defined yet." << endl;
    return ostr;
}
