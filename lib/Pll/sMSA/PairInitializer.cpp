#include "PairInitializer.h"

PairInitializer::PairInitializer() {
    Seqs       = NULL;
    Params     = NULL;
    FTree      = NULL;
    f84        = NULL;
    tn93       = NULL;
    gtr        = NULL;
    NumOfSteps = 5;
}

PairInitializer::~PairInitializer() {
    delete Params;
    delete Seqs;
    delete FTree;
    delete f84;
    delete tn93;
    delete gtr;
    delete RandGen;
}

/** Initialization Functions **/ //{
//Function to Read Sequences and Alignment from File and Initiliaze the Corresponding Sequences Object.
//The subtree argument is needed in order to delete all the "All-1" Columns from the OverallHomologies Matrix according to the specified subtree
void PairInitializer::InitSequences() {
    ClearInputStream();

    Seqs = new Sequences();

    /** Read Bare Sequences and Sequence Ids **/
    string line;
    ifstream myfile (BareSequencesFile.c_str());
    string filename = BareSequencesFile;

    stringstream seqId;
    stringstream seqText;
    bool firstLine = true;

    if (myfile.is_open()) {
        while ( myfile.good() ) {
            getline (myfile,line);

            if (line.compare("")==0 || line.compare("\n")==0 || line.compare(" ")==0) continue; //skip empty lines
            if (StartsWith(line, ">") || StartsWith(line,";")) {
                if (!firstLine) {
                    Seqs->AddBareSequence(seqText.str());
                    Seqs->AddSequenceId(seqId.str());
                    seqId.str("");
                    seqText.str("");
                }
                seqId << line;
                firstLine = false;
            } else {
                seqText << line;
            }

        }
        Seqs->AddBareSequence(seqText.str());
        Seqs->AddSequenceId(seqId.str());
    } else cerr << " Unable to open file: "<<filename<<" for reading.";
    myfile.close();

    /** Read Alignment **/
    string line1;
    ifstream myfile1 (AlignedSequencesFile.c_str());
    string filename1 = AlignedSequencesFile;

    stringstream seqText1;
    bool firstLine1 = true;

    if (myfile1.is_open()) {
        while ( myfile1.good() ) {
            getline (myfile1,line1);

            if (line1.compare("")==0 || line1.compare("\n")==0 || line1.compare(" ")==0) continue; //skip empty lines

            if (StartsWith(line1, ">") || StartsWith(line1,";")) {
                if (!firstLine1) {
                    Seqs->AddAlignedSequence(seqText1.str());
                    seqText1.str("");
                }
                firstLine1 = false;
            } else {
                seqText1 << line1;
            }

        }
        Seqs->AddAlignedSequence(seqText1.str());
    } else cerr << " Unable to open file: "<<filename1<<" for reading.";
    myfile1.close();

    Seqs->SetIndicesFromAlignment();
}

//Function to Read the Desired parameters (e.g. ã, ë rates, Prior Nucleotide frequencies etc) and Initialize the Parameters instance
void PairInitializer::InitParameters() {
    Params = new Parameters();
    Params->ReadDataFromFile(ParametersFile);
}

//Function to Initialize the Felsentein Tree Structure for the Given Sequences
void PairInitializer::InitTree() {
    FTree = new Tree(TreeFile, Seqs);
}

//Function to Initialize the substitution matrices for each branch length
void PairInitializer::InitOverallSubstitutionMatrices() {
    //Init SubModel
    if (Params->Get_SubModel().compare("F84")==0) {
        f84 = new F84Model(Params->Get_BaseFreq(), Params->Get_Kappa());
    } else if (Params->Get_SubModel().compare("TN93")==0) {
        tn93 = new TN93Model(Params->Get_BaseFreq(), Params->Get_AG_Rate(), Params->Get_CT_Rate());
    } else if (Params->Get_SubModel().compare("GTR")==0) {
        gtr = new GTRModel(Params->Get_BaseFreq(), Params->Get_X1_Param(), Params->Get_X2_Param(), Params->Get_X3_Param(), Params->Get_X4_Param(), Params->Get_X5_Param(), Params->Get_X6_Param());
    }

    //Compute Substitution Matrices according to the Selected Substitution Model
    DoubleVec tempVec(16, 0.0);
    for (int i=0; i<FTree->Get_UniqueLengths().size(); i++) {
        double* subM = NULL;
        if (Params->Get_SubModel().compare("F84")==0) {
            subM = f84->GetSubstitutionMatrix(FTree->Get_UniqueLengths()[i]);
        } else if (Params->Get_SubModel().compare("TN93")==0) {
            subM = tn93->GetSubstitutionMatrix(FTree->Get_UniqueLengths()[i]);
        } else if (Params->Get_SubModel().compare("GTR")==0) {
            subM = gtr->GetSubstitutionMatrix(FTree->Get_UniqueLengths()[i]);
        }

        for (int j=0; j<tempVec.size(); j++) {
            tempVec[j] = subM[j];
        }
        OverallSubMats.push_back(tempVec);
        OverallSubMatsLengths.push_back(FTree->Get_UniqueLengths()[i]);
    }
}

//Function to initialize the Felsenstein Pruning Tables on the Inner Nodes of the Tree
//This is not done simulataneously with the tree creation because the Substitution Matrices are needed for the
//Felsenstein Probablilities estimation, and for these, the corresponding branch lengths are required (i.e. the tree has to be fully created first)
void PairInitializer::InitTreeFPTables() {
    //Nodes A, B, C and D
    {
        for (int k=0;k<4;k++){
            FTree->Get_TreeNodes()[k].N1.FPTable.resize(4, DoubleVec(Seqs->Get_BareSeqs()[k].size(),0));
            for (int i=0; i<Seqs->Get_BareSeqs()[k].size(); i++) {
                FTree->Get_TreeNodes()[k].N1.FPTable[GetAlphabet(Seqs->Get_BareSeqs()[k][i])][i]= 1.0;
            }
            FTree->Get_TreeNodes()[k].N2.FPTable.resize(4, DoubleVec(Seqs->Get_BareSeqs()[k].size(),0));
            FTree->Get_TreeNodes()[k].N3.FPTable.resize(4, DoubleVec(Seqs->Get_BareSeqs()[k].size(),0));

        }

    }

    //Node E
    {
        //InnerNode N1 conditioned on Nodes A and B
        double len1 = FTree->Get_TreeNodes()[4].N2.Length;
        int idx1 = GetPosOfElement(OverallSubMatsLengths, len1, 0);

        double len2 = FTree->Get_TreeNodes()[4].N3.Length;
        int idx2 = GetPosOfElement(OverallSubMatsLengths, len2, 0);

        FTree->Get_TreeNodes()[4].N1.FPTable.resize(4, DoubleVec(0,0));
        for (int j=0; j<AlignmentHomologies[4].size(); j++) {
            if (AlignmentHomologies[4][j]==-1) {
                continue;
            }
            for (int i=0; i<4; i++) { //For each possibility, meaning A,C,G or T
                double prob = 1;
                if (AlignmentHomologies[0][j]!=-1) {
                    prob *= OverallSubMats[idx1][GetAlphabetPair(i, Seqs->Get_BareSeqs()[0][AlignmentHomologies[0][j]-1])];

                }

                if (AlignmentHomologies[1][j]!=-1) {
                    prob *= OverallSubMats[idx2][GetAlphabetPair(i, Seqs->Get_BareSeqs()[1][AlignmentHomologies[1][j]-1])];

                }
                FTree->Get_TreeNodes()[4].N1.FPTable[i].push_back(prob);
            }
        }
    }

    //Node 6
    {
        //InnerNode N1 conditioned on Nodes C and D
        double len1 = FTree->Get_TreeNodes()[5].N2.Length;
        int idx1 = GetPosOfElement(OverallSubMatsLengths, len1, 0);

        double len2 = FTree->Get_TreeNodes()[5].N3.Length;
        int idx2 = GetPosOfElement(OverallSubMatsLengths, len2, 0);

        FTree->Get_TreeNodes()[5].N1.FPTable.resize(4, DoubleVec(0,0));
        for (int j=0; j<AlignmentHomologies[5].size(); j++) {
            if (AlignmentHomologies[5][j]==-1) {
                continue;
            }
            for (int i=0; i<4; i++) { //For each possibility, meaning A,C,G or T
                double prob = 1;
                if (AlignmentHomologies[2][j]!=-1) {
                    prob *= OverallSubMats[idx1][GetAlphabetPair(i, Seqs->Get_BareSeqs()[2][AlignmentHomologies[2][j]-1])];

                }

                if (AlignmentHomologies[3][j]!=-1) {
                    prob *= OverallSubMats[idx2][GetAlphabetPair(i, Seqs->Get_BareSeqs()[3][AlignmentHomologies[3][j]-1])];

                }
                FTree->Get_TreeNodes()[5].N1.FPTable[i].push_back(prob);
            }
        }
    }


    //Node 5 - InnerNodes N2 and N3
    {
        //InnerNode N2 conditioned on Nodes B, C and D
        {
            double len1 = FTree->Get_TreeNodes()[4].N3.Length;
            int idx1 = GetPosOfElement(OverallSubMatsLengths, len1, 0);

            double len2 = FTree->Get_TreeNodes()[4].N1.Length;
            int idx2 = GetPosOfElement(OverallSubMatsLengths, len2, 0);

            FTree->Get_TreeNodes()[4].N2.FPTable.resize(4, DoubleVec(0,0));
            for (int j=0; j<AlignmentHomologies[4].size(); j++) {
                if (AlignmentHomologies[4][j]==-1) {
                    continue;
                }

                for (int i=0; i<4; i++) { //For each possibility, meaning A,C,G or T
                    double prob = 1;
                    if (AlignmentHomologies[1][j]!=-1) {
                        prob *= OverallSubMats[idx1][GetAlphabetPair(i, Seqs->Get_BareSeqs()[1][AlignmentHomologies[1][j]-1])];

                    }

                    //if (FTree->Get_TreeNodes()[4]->Get_Homologies()[5][j]!=-1) {
                    if (AlignmentHomologies[5][j]!=-1) {
                        double sum = 0;
                        for (int k=0; k<4; k++) {
                            sum += OverallSubMats[idx2][GetAlphabetPair(i, k)]* (FTree->Get_TreeNodes()[5].N1.FPTable[k][AlignmentHomologies[5][j]-1]);
                        }
                        prob *= sum;
                    }

                    FTree->Get_TreeNodes()[4].N2.FPTable[i].push_back(prob);
                }
            }
        }

        //InnerNode N3 conditioned on Nodes A, C and D
        {
            double len1 = FTree->Get_TreeNodes()[4].N2.Length;
            int idx1 = GetPosOfElement(OverallSubMatsLengths, len1, 0);

            double len2 = FTree->Get_TreeNodes()[4].N1.Length;
            int idx2 = GetPosOfElement(OverallSubMatsLengths, len2, 0);

            FTree->Get_TreeNodes()[4].N3.FPTable.resize(4, DoubleVec(0,0));
            for (int j=0; j<AlignmentHomologies[4].size(); j++) {
                if (AlignmentHomologies[4][j]==-1) {
                    continue;
                }

                for (int i=0; i<4; i++) { //For each possibility, meaning A,C,G or T
                    double prob = 1;
                    if (AlignmentHomologies[0][j]!=-1) {
                        prob *= OverallSubMats[idx1][GetAlphabetPair(i, Seqs->Get_BareSeqs()[0][AlignmentHomologies[0][j]-1])];

                    }

                    if (AlignmentHomologies[5][j]!=-1) {
                        double sum = 0;
                        for (int k=0; k<4; k++) {
                            sum += OverallSubMats[idx2][GetAlphabetPair(i, k)]* (FTree->Get_TreeNodes()[5].N1.FPTable[k][AlignmentHomologies[5][j]-1]);
                        }
                        prob *= sum;
                    }

                    FTree->Get_TreeNodes()[4].N3.FPTable[i].push_back(prob);
                }
            }
        }

    }


    //Node 6 - InnerNodes N2 and N3
    {
        //InnerNode N2 conditioned on Nodes A, B and D
        {
            double len1 = FTree->Get_TreeNodes()[5].N3.Length;
            int idx1 = GetPosOfElement(OverallSubMatsLengths, len1, 0);

            double len2 = FTree->Get_TreeNodes()[5].N1.Length;
            int idx2 = GetPosOfElement(OverallSubMatsLengths, len2, 0);

            FTree->Get_TreeNodes()[5].N2.FPTable.resize(4, DoubleVec(0,0));
            for (int j=0; j<AlignmentHomologies[5].size(); j++) {
                if (AlignmentHomologies[5][j]==-1) {
                    continue;
                }

                for (int i=0; i<4; i++) { //For each possibility, meaning A,C,G or T
                    double prob = 1;
                    if (AlignmentHomologies[3][j]!=-1) {
                        prob *= OverallSubMats[idx1][GetAlphabetPair(i, Seqs->Get_BareSeqs()[3][AlignmentHomologies[3][j]-1])];

                    }

                    if (AlignmentHomologies[4][j]!=-1) {
                        double sum = 0;
                        for (int k=0; k<4; k++) {
                            sum += OverallSubMats[idx2][GetAlphabetPair(i, k)]* (FTree->Get_TreeNodes()[4].N1.FPTable[k][AlignmentHomologies[4][j]-1]);
                        }
                        prob *= sum;
                    }

                    FTree->Get_TreeNodes()[5].N2.FPTable[i].push_back(prob);
                }
            }

        }

        //InnerNode N3 conditioned on Nodes A, B and C
        {
            double len1 = FTree->Get_TreeNodes()[5].N2.Length;
            int idx1 = GetPosOfElement(OverallSubMatsLengths, len1, 0);

            double len2 = FTree->Get_TreeNodes()[5].N1.Length;
            int idx2 = GetPosOfElement(OverallSubMatsLengths, len2, 0);

            FTree->Get_TreeNodes()[5].N3.FPTable.resize(4, DoubleVec(0,0));
            for (int j=0; j<AlignmentHomologies[5].size(); j++) {
                if (AlignmentHomologies[5][j]==-1) {
                    continue;
                }

                for (int i=0; i<4; i++) { //For each possibility, meaning A,C,G or T
                    double prob = 1;
                    if (AlignmentHomologies[2][j]!=-1) {
                        prob *= OverallSubMats[idx1][GetAlphabetPair(i, Seqs->Get_BareSeqs()[2][AlignmentHomologies[2][j]-1])];

                    }

                    if (AlignmentHomologies[4][j]!=-1) {
                        double sum = 0;
                        for (int k=0; k<4; k++) {
                            sum += OverallSubMats[idx2][GetAlphabetPair(i, k)]* (FTree->Get_TreeNodes()[4].N1.FPTable[k][AlignmentHomologies[4][j]-1]);
                        }
                        prob *= sum;
                    }

                    FTree->Get_TreeNodes()[5].N3.FPTable[i].push_back(prob);
                }
            }
        }

    }

}

//Function to Initialize the (Global) Random generator Object to generate random numbers
void PairInitializer::InitRandomGenerator() {
    RandGen = new BoostRandomGenerator();
    //RandGen->SetSeed(5);  //Is Set to 5 for repeatable experiments. Set to -1 for a random seed
    RandGen->SetSeed(-1);  //Is Set to 5 for repeatable experiments. Set to -1 for a random seed
}

//Function to Initialize how many times the left-right subtree re-alignment procedure will take place
void PairInitializer::InitNumOfSteps() {
    NumOfSteps = 5;
}

//Function to Initiliaze the Homology Stuctur for the given Subtree of interest. The initialization is performed by resampling the whole initial alignment using the given homologies
void PairInitializer::ResampleHistoryGivenHomologies() {
    //Resample all the region considering the given Homologies //{
    double prob;

    PairAligner* aligner = new PairAligner();
    aligner->Set_RandomGenerator(RandGen);
    aligner->InitParameters(Params);
    aligner->InitSubModel(f84, tn93, gtr);
    aligner->InitFidModel();
    aligner->InitLengths(Lengths);
    aligner->InitDimensions(Dimensions);
    aligner->InitFPTables(FPTables);
    aligner->InitHomologyIndices(RegionHomologies);
    aligner->InitTriplets();
    aligner->Set_EndPosIn_ND(EndIdxs);
    aligner->InitProbabilities("HH", "HH", StartIdxs);
    aligner->InitSubMatrices(SubMats);
//aligner->Set_SaveToLogFile(true, LogHomologsFile);
    aligner->EstimateProbabilities(1, true);
//aligner->SaveResultsToFile(ProbabilitiesHomologsFile);
    aligner->GetEvolHistory(EvolHistory, prob, Soans, Tihls, SoansIdxs);

    /*For Debugging*/
    //PrintEvolHistory();
    //PrintEvolHistoryReport();

    delete aligner;
    //}
}

//Function to Initialize the Branch
void PairInitializer::InitBranch(int branch) {
    Branch = branch;
}

//Function to Initialize the SubtreeIdxs
void PairInitializer::InitBranchIdxs() {
    switch(Branch){
        case 1:
        {
            //Branch 1 (E-F)
            BranchIdxs.push_back(4); //Node E is the root
            BranchIdxs.push_back(5); //Node F is the 1st Child
            break;
        }

        case 2:
        {
            //Branch 2 (E-A)
            BranchIdxs.push_back(4); //Node E is the root
            BranchIdxs.push_back(0); //Node A is the 1st Child
            break;
        }

        case 3:
        {
            //Branch 3 (E-B)
            BranchIdxs.push_back(4); //Node E is the root
            BranchIdxs.push_back(1); //Node B is the 1st Child
            break;
        }

        case 4:
        {
            //Branch 4 (F-C)
            BranchIdxs.push_back(5); //Node F is the root
            BranchIdxs.push_back(2); //Node C is the 1st Child
            break;
        }

        case 5:
        {
            //Branch 5 (F-D)
            BranchIdxs.push_back(5); //Node F is the root
            BranchIdxs.push_back(3); //Node D is the 1st Child
            break;
        }

        case 6:
        {
            //Branch 6 (F-E)
            BranchIdxs.push_back(5); //Node F is the root
            BranchIdxs.push_back(4); //Node E is the 1st Child
            break;
        }
        default: cout<<"Error branch index!"<<endl; exit(1);
    }
}

//Function to Initialize the Lengths for the current Subtree
void PairInitializer::InitLengths() {
    Lengths.push_back(FTree->GetLengthOfNode(BranchIdxs[0], 1));
    Lengths.push_back(FTree->GetLengthOfNode(BranchIdxs[1], 1));
}

//Function to Initialize the Felsenstein Pruning tables for the current Subtree
void PairInitializer::InitFPTables() {
    if (Branch==1){
        FPTables.push_back(FTree->GetFPTableOfNode(BranchIdxs[0], 1));
        FPTables.push_back(FTree->GetFPTableOfNode(BranchIdxs[1], 1));
    }else if (Branch==2 || Branch==4){
        FPTables.push_back(FTree->GetFPTableOfNode(BranchIdxs[0], 2));
        FPTables.push_back(FTree->GetFPTableOfNode(BranchIdxs[1], 1));
    }else if (Branch==3 || Branch==5){
        FPTables.push_back(FTree->GetFPTableOfNode(BranchIdxs[0], 3));
        FPTables.push_back(FTree->GetFPTableOfNode(BranchIdxs[1], 1));
    }


    //Add an 1 Column for the Extra Start State and another one for the Extra End State
    for (int i=0;i<FPTables.size();i++) {
        for (int j=0;j<FPTables[i].size();j++){
            FPTables[i][j].insert(FPTables[i][j].begin(), 1);
            FPTables[i][j].insert(FPTables[i][j].end(), 1);
        }
    }
}

//Function to Initialize the Alignment Homologies for the current Subtree
void PairInitializer::InitAlignmentHomologies(IntMatRef Homologies) {
    if (Homologies.empty()){
        AlignmentHomologies.clear();
        AlignmentHomologies = Seqs->Get_IndicesFromAlignment();
    }else{
        AlignmentHomologies = Homologies;
    }
}

//Function to Initialize the Region Homologies for the current Subtree
void PairInitializer::InitRegionHomologies() {
    RegionHomologies.clear();
    RegionHomologies.resize(2, IntVec(1,0));
    for (int i=0;i<BranchIdxs.size();i++){
        IntVec temp(1, 0);
        for (int j=0;j<AlignmentHomologies[0].size();j++) {
            if (AlignmentHomologies[BranchIdxs[0]][j]==-1){ continue; }
            if (AlignmentHomologies[BranchIdxs[i]][j]==-1){
                RegionHomologies[i].push_back(-1);
            }else{
                RegionHomologies[i].push_back(AlignmentHomologies[BranchIdxs[i]][j]);
            }
        }
        RegionHomologies[i].push_back(GetFirstNonNegativeValueReverse(AlignmentHomologies[BranchIdxs[i]], AlignmentHomologies[BranchIdxs[i]].size()-1)+1);
    }
}

//Function to Initialize the Start and End Soan Indices for the current Subtree
void PairInitializer::InitStartAndEndIndices() {
    StartIdxs.clear();
    EndIdxs.clear();

    StartIdxs.resize(3, 0);
    for (int i=0;i<RegionHomologies.size();i++) {
        EndIdxs.push_back(RegionHomologies[i][RegionHomologies[i].size()-1]);
    }
    EndIdxs.push_back(0);
}

//Function to Initialize the Substitution Matrices needed for the current Subtree
void PairInitializer::InitSubstitutionMatrices() {
    for (int i=0;i<Lengths.size();i++){
        SubMats.push_back(OverallSubMats[GetPosOfElement(OverallSubMatsLengths, Lengths[i], 0)]);
    }
}

//Function to Initialize the Dimensions needed for the current Subtree alignment region
void PairInitializer::InitDimensions() {
    for (int i=0;i<EndIdxs.size()-1;i++){
        Dimensions.push_back(EndIdxs[i]+1);
    }
}

//Function to Call the above Initialization Function in the Proper Order so as to assure that all inter-dependencies between successive calls are met. The argument defines which branch will be initialized
void PairInitializer::MakeProperInitializations(int branch, IntMatRef Homologies) {
    //A. Make General Initializations in proper order to ensure proper tree construction and parameters intialization //{
        //1. Initialize the Parameters from the Parameters File
        InitParameters();

        //2. Initialize the Rand Number Generator Engine by creating a new instance of the BoostRandomGenerator class
        InitRandomGenerator();

        //3. Initialize the Number of Steps that the left-right subtree resampling will be performed (from the previously created Parameters instance)
        InitNumOfSteps();

        //4. Initialize the Sequences (i.e. Bare Seqeunces, Initial Alignment, Sequence Headers etc) from the Corresponding Input Files
        InitSequences();

        //5. Initialize the AlignmentHomologies based on the previously read initial alignment
        InitAlignmentHomologies(Homologies);

        //6. Initialize (create) the Felsenstein Tree structure based on the red sequences and parameters
        InitTree();

        //7. Now that the tree structure is ready, Initialize ALL the Substitution Matrices for EVERY possible length of the Tree (and based on the defined - by the parameters - substitution model)
        InitOverallSubstitutionMatrices();

        //8.In a same way as before, now that the tree structure is ready, Initialize ALL the Felsenstein Pruning Tables for EVERY Node of the tree
        InitTreeFPTables();

    //}

    //B. Proceed to Make Proper Initializations for the variables and parameters concerning ONLY the subtree of interest (i.e. either the left or the right subtree) //{
        //1. Initialize the Subtree parameter to be aware of which subtree is being resampled
        InitBranch(branch);

        //2. Initialize the corresponding SubtreeIdxs for the current subtree
        InitBranchIdxs();

        //3. Initialize the Lengths corresponding only to the current subtree
        InitLengths();

        //4. Initialize the FPTables corresponding only to the current subtree
        InitFPTables();

        //5. Initialize (calculate) the Corresponding RegionHomologies based on the AlignmentHomologies corresponding only to the current fragment
        InitRegionHomologies();

        //6. Initialize the Start and End Indices for the current subtree
        InitStartAndEndIndices();

        //7. Initialize the Substitution Matrices corresponding only to the current subtree (and equivalently to the current subtree's lengths)
        InitSubstitutionMatrices();

        //8. Initialize the Dimensions (based on the previously estimated EndIdxs) corresponding  to the current subtree
        InitDimensions();

        //9. Now that everything is properly initialized, proceed to resample the whole alignment region (considering the given homologies) in order to come up with an evolution history, the Soans, the SoansIdxs and the Tihls for EVERY positon of the Root node (defined by the current subtree)
        ResampleHistoryGivenHomologies();

    //}

}

//}

/** Assisting Functions for In-Between calculations - Getters **/ //{

//Function that parses the Soan Indices Matrix to find the FIRST soan whose index corresponding to the root has a value equal to the one provided as argument. When such a soan is found, the corresponding index is returned or -1 otherwise
int PairInitializer::GetFirstIdxofSoanWhereRootIs(int rootIdx) {
    for (int i=0;i<SoansIdxs.size();i++) {
        if (SoansIdxs[i][0]==rootIdx) return i;
    }
    return -1;
}

//Function that parses the Soan Indices Matrix to find the FIRST soan whose index corresponding to the root has a value equal to the one provided as argument. When such a soan is found, the corresponding index is returned or -1 otherwise
int PairInitializer::GetLastIdxofSoanWhereRootIs(int rootIdx) {
    for (int i=SoansIdxs.size()-1;i>=0;i--) {
        if (SoansIdxs[i][0]==rootIdx) return i;
    }
    return -1;
}

//Function that returns the Tihl lying at the given position within the Tihls vector. The boolean argument represents whether the tihl corresponding to the root will also be included at the returned result or not
string PairInitializer::GetTihlAtPosition(int pos, bool includeRoot) {
    if (includeRoot) { return Tihls[pos];           }
    else             { return Tihls[pos].substr(1); }
}

//Function that returns the Soan lying at the given position within the Soans vector. The boolean argument represents whether the tihl corresponding to the root will also be included at the returned result or not
string PairInitializer::GetSoanAtPosition(int pos, bool includeRoot) {
    if (includeRoot) { return Soans[pos];           }
    else             { return Soans[pos].substr(1); }
}

//Function that returns the Soan Indices at the given position within the SoansIdxs matrix.
IntVecRef PairInitializer::GetSoanIdxsAtPosition(int pos) {
    //Fix the indices of the soan to be "rescaled" to the corresponding indices of the 3-leaved tree
    return SoansIdxs[pos];
}

//Function that returns the maximum value of the soan index corresponding to the root of the current subtree
int PairInitializer::GetRootMaxValue(){
    return SoansIdxs[SoansIdxs.size()-1][0];
}

//}

/** Functions to Output results **/ //{
//Function to Output the Formed Sequence Objects to the standard output stream
void PairInitializer::PrintSequences() {
    if (Seqs==NULL || Seqs->Get_AlignSeqs().size()==0) {
        cout<<"\n Please specify input sequences first.\n\n";
        return;
    }
    cout << "\n" << *Seqs <<endl;
}

//Function to Output the Program Parameters to the standard output stream
void PairInitializer::PrintParameters() {
    if (Params==NULL) {
        cout<<"\n Please specify program parameters first.\n\n";
        return;
    }
    cout<<" \n"<< *(Params)<<endl;
}

//Function to Output the Estimated Substitution Matrices to the standard output stream
void PairInitializer::PrintOverallSubstitutionMatrices() {
    if (OverallSubMats.size()==0) {
        cout<<"\n Please init sequences and program parameters first (option 2).\n\n";
        return;
    }

    cout<<endl<<endl<<"\t\tOVERALL SUBSTITUTION MATRICES"<<endl;
    cout<<"\t\t============================="<<endl;

    for (int i=0; i<OverallSubMats.size(); i++) {
        cout<<"For Length: "<<OverallSubMatsLengths[i]<<endl;
        cout<<"\t"<<setw(15)<<"A"<<setw(12)<<"C"<<setw(12)<<"G"<<setw(12)<<"T"<<setw(11)<<"SUM:"<<endl<<"\t";
        //vector<double> sumC(4, 0.0);

        for (int j=0; j<4; j++) {
            double sumL = 0;
            cout<<setw(3)<<GetAlphabet_String(j);
            for (int k=0; k<4; k++) {
                double tempVal = OverallSubMats[i][GetAlphabetPair(j,k)];
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

//Function to Output the Estimated Felsenstein Tree Structure of the Sequences to the standard output stream
void PairInitializer::PrintTree() {
    if (FTree==NULL){
        cout<<"\n Please init sequences and program parameters first (option 2).\n\n";
        return;
    }
    cout<< *FTree <<endl;
}

//Function to Save the Estimated Felsenstein Tree Structure of the Sequences to an Output File
void PairInitializer::SaveTreeToFile() {
    if (FTree ==  NULL) {
        cout<<"\n Please init sequences and program parameters first (option 2).\n\n";
        return;
    }
    stringstream ss;
    ss.str("");
    ss<< *FTree <<endl;
    WriteToFile(ExportTreeFile, ss.str(), false, true);
}

//Function to Print which Branch is currently being examined
void PairInitializer::PrintBranch() {
    if (Branch==1){
        cout<<"\n\nBranch 1 (E-A), E being the root.\n";
    }else if (Branch==2){
        cout<<"\n\nBranch 2 (E-B), E being the root.\n";
    }else if (Branch==3){
        cout<<"\n\nBranch 3 (E-F), E being the root.\n";
    }else if (Branch==4){
        cout<<"\n\nBranch 4 (F-C), F being the root.\n";
    }else if (Branch==5){
        cout<<"\n\nBranch 5 (F-D), F being the root.\n";
    }else if (Branch==6){
        cout<<"\n\nBranch 6 (F-E), F being the root.\n";
    }
}

//Function to Print the subtree lengths
void PairInitializer::PrintLengths() {
    cout<<"\nCurrent Subtree Lengths:\n------------------------\n";
    for (int i=0;i<Lengths.size();i++){
        cout<<setw(4)<<Lengths[i];
    }
    cout<<"\n\n";
}

//Function to Print the BranchIdxs of the current realignment region
void PairInitializer::PrintBranchIdxs() {
    cout<<"\nBranch Indices:\n----------------\n\t<";
    for (int i=0;i<BranchIdxs.size();i++){
        cout<<setw(4)<<BranchIdxs[i];
    }
    cout<<">\n\n";
}

//Function to Print the Substitution Matrices corresponding to the current fragment to the standard output stream
void PairInitializer::PrintSubstitutionMatrices() {
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
void PairInitializer::PrintFPTables() {
    cout<<"\nFragment's FPTables:\n--------------------\n";

    //cout<<"(Sizes: #of FPTables:"<<FPTables.size()<<"  #of rows in FPTable[0]:"<<FPTables[0].size()<<"  #of sites in FPTable[0][0]:"<<FPTables[0][0].size()<<")"<<endl;
    cout.precision(5);
    for (int i=0;i<FPTables.size();i++){
        cout<<"FPTable - "<<i<<endl<<"-----------"<<endl;
        for(int j=0;j<FPTables[i].size();j++){
            cout<<"\t";
            for (int k=0;k<FPTables[i][j].size();k++){
                cout<<left<<setw(12)<<FPTables[i][j][k];
            }
            cout<<endl;
        }
        cout<<endl;
    }
}

//Function to Print the Dimensions for the current fragment to the standard output stream
void PairInitializer::PrintDimensions() {
    cout<<"\nFragment Dimensions:\n--------------------\n";
    for (int i=0;i<Dimensions.size();i++){
        cout<<setw(4)<<Dimensions[i];
    }
    cout<<">\n\n";
}

//Function to Print the AlignmentHomologies to the Standard Output Stream
void PairInitializer::PrintAlignmentHomologies() {
    stringstream ss; ss.str("");
    for (int i=0;i<AlignmentHomologies.size();i++){
        if (i==0)  {
            for (int j=0;j<AlignmentHomologies[i].size();j++){
                ss<<setw(5)<<j;
            }
            ss<<endl;
        }
        for (int j=0;j<AlignmentHomologies[i].size();j++){
            ss<<setw(5)<<AlignmentHomologies[i][j];
        }
        ss<<endl;
    }
    ss<<endl;
    cout<<ss.str()<<endl;
}

//Function to Print the RegionHomologies to the Standard Output Stream
void PairInitializer::PrintRegionHomologies() {
    stringstream ss; ss.str("");
    for (int i=0;i<RegionHomologies.size();i++){
        if (i==0)  {
            for (int j=0;j<RegionHomologies[i].size();j++){
                ss<<setw(5)<<j;
            }
            ss<<endl;
        }
        for (int j=0;j<RegionHomologies[i].size();j++){
            ss<<setw(5)<<RegionHomologies[i][j];
        }
        ss<<endl;
    }
    ss<<endl;
    cout<<ss.str()<<endl;
}

//Function to Print the SampledHomologies to the Standard Output Stream
void PairInitializer::PrintSampledHomologies() {
    stringstream ss; ss.str("");
    for (int i=0;i<SampledHomologies.size();i++){
        if (i==0)  {
            for (int j=0;j<SampledHomologies[i].size();j++){
                ss<<setw(5)<<j;
            }
            ss<<endl;
        }
        for (int j=0;j<SampledHomologies[i].size();j++){
            ss<<setw(5)<<SampledHomologies[i][j];
        }
        ss<<endl;
    }
    ss<<endl;
    cout<<ss.str()<<endl;
}

//Function to Print the Homologies (i.e. SampledHomologies) to the Standard Output Stream
void PairInitializer::PrintHomologies() {
    stringstream ss; ss.str("");
    for (int i=0;i<SampledHomologies.size();i++){
        if (i==0)  {
            for (int j=0;j<SampledHomologies[i].size();j++){
                ss<<setw(5)<<j;
            }
            ss<<endl;
        }
        for (int j=0;j<SampledHomologies[i].size();j++){
            ss<<setw(5)<<SampledHomologies[i][j];
        }
        ss<<endl;
    }
    ss<<endl;
    cout<<ss.str()<<endl;
}

//Function to Save the Estimated AlignmentHomologies to an Output File
void PairInitializer::SaveAlignmentHomologiesToFile() {
    stringstream ss; ss.str("");
    for (int i=0;i<AlignmentHomologies.size();i++){
        if (i==0)  {
            for (int j=0;j<AlignmentHomologies[i].size();j++){
                ss<<setw(5)<<j;
            }
            ss<<endl;
        }
        for (int j=0;j<AlignmentHomologies[i].size();j++){
            ss<<setw(5)<<AlignmentHomologies[i][j];
        }
        ss<<endl;
    }
    ss<<endl;
    WriteToFile(AlignmentHomologiesFile, ss.str(), false, true);
}

//Function to Save the Estimated RegionHomologies to an Output File
void PairInitializer::SaveRegionHomologiesToFile() {
    stringstream ss; ss.str("");
    for (int i=0;i<RegionHomologies.size();i++){
        if (i==0)  {
            for (int j=0;j<RegionHomologies[i].size();j++){
                ss<<setw(5)<<j;
            }
            ss<<endl;
        }
        for (int j=0;j<RegionHomologies[i].size();j++){
            ss<<setw(5)<<RegionHomologies[i][j];
        }
        ss<<endl;
    }
    ss<<endl;
    WriteToFile(RegionHomologiesFile, ss.str(), false, true);
}

//Function to Save the Estimated SampledHomologies to an Output File
void PairInitializer::SaveSampledHomologiesToFile() {
    stringstream ss; ss.str("");
    for (int i=0;i<SampledHomologies.size();i++){
        if (i==0)  {
            for (int j=0;j<SampledHomologies[i].size();j++){
                ss<<setw(5)<<j;
            }
            ss<<endl;
        }
        for (int j=0;j<SampledHomologies[i].size();j++){
            ss<<setw(5)<<SampledHomologies[i][j];
        }
        ss<<endl;
    }
    ss<<endl;
    WriteToFile(HomologiesFile, ss.str(), false, true);
}

//Function to Save the Estimated Homologies (i.e. SampledHomologies) to an Output File
void PairInitializer::SaveHomologiesToFile() {
    stringstream ss; ss.str("");
    ss<<setw(7)<<"Tihl"<<setw(7)<<"Soan"<<left<<setw(30)<<"          Soan Indices"<<right<<endl;
    ss<<setw(7)<<"===="<<setw(7)<<"===="<<left<<setw(30)<<"  ==========================="<<right<<endl;
    for (int i=0;i<SoansIdxs.size();i++) {
        ss<<setw(7)<<Tihls[i]<<setw(7)<<Soans[i]<<"  <";
        for (int j=0;j<SoansIdxs[i].size();j++){
            ss<<setw(5)<<SoansIdxs[i][j];
        }
        ss<<">"<<endl;
    }
    ss<<endl;
    WriteToFile(HomologiesFile, ss.str(), false, true);
}

//Function to Print the Sampled Evolution History to the Standard Output Stream
void PairInitializer::PrintEvolHistory() {
    cout<<"Sampled Evolution History:\n--------------------------\n"<<EvolHistory<<endl;
}

//Function to Print the Sampled Soans, Soan Indices and Tihls to the Standard Output Stream
void PairInitializer::PrintEvolHistoryReport() {
    cout<<setw(7)<<"Tihl"<<setw(7)<<"Soan"<<left<<setw(30)<<"          Soan Indices"<<right<<endl;
    cout<<setw(7)<<"===="<<setw(7)<<"===="<<left<<setw(30)<<"  ==========================="<<right<<endl;
    for (int i=0;i<SoansIdxs.size();i++) {
        cout<<setw(7)<<Tihls[i]<<setw(7)<<Soans[i]<<"  <";
        for (int j=0;j<SoansIdxs[i].size();j++){
            cout<<setw(5)<<SoansIdxs[i][j];
        }
        cout<<">"<<endl;
    }
}


//}

