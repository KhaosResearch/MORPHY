#include "Menus.h"

//Class Constructor
Menus::Menus() {
    TripleInit = new TripleInitializer();
    PairInit   = new PairInitializer();
}

//Virtual Deconstructor
Menus::~Menus() {
    delete TripleInit;
    delete PairInit;
}

void Menus::ShowMainMenu() {
    while (true) {
        ClearScreen();

        cout<<" *====================================================================*\n";
        cout<<" *            sMSA - Statistical  Multiple Sequence Aligner           *\n";
        cout<<" *====================================================================*\n";
        cout<<"\n";

        cout<<"     1. Define  Program (Parameters, Sequences and Alignment)\n";
        cout<<"     2. Init    Program (Parameters, Sequences and Alignment)\n";
        cout<<"     3. Preview Program Parameters\n";
        cout<<"     4. Preview Sequences and Alignment\n";
        cout<<"     5. Preview Substitution Matrices\n";
        cout<<"     6. Export  Tree to File\n";
        cout<<"     7. Align   Fragment of the Sequences\n";
        cout<<"     H. Show    Help\n";
        cout<<"     X. EXIT    (Terminate Program)\n\n\n";
        cout<<" Please type your choice -> ";
        char choice;
        cin>>choice;

        switch (toupper(choice)) {
        case '1':
            ExecuteProcess("java -jar \"Jar Executables/sMSA.jar\" Input");
            break;
        case '2': {
            InitTripleInitializer(1, false);
            //InitPairInitializer(1);
            cout<<"\n Parameters and Sequences initialized successfully!\n";
            Pause();
            break;
        }
        case '3':
            TripleInit->PrintParameters();
            Pause();
            break;
        case '4':
            TripleInit->PrintSequences();
            Pause();
            break;
        case '5':
            TripleInit->PrintSubstitutionMatrices();
            Pause();
            break;
        case '6': {
            TripleInit->SaveTreeToFile();
            TripleInit->SaveAlignmentHomologiesToFile();
            TripleInit->SaveHomologiesToFile();
            Pause();
            break;
        }
        case '7': {
            if (TripleInit->Get_Tree() ==  NULL) {
                cout<<"\n Please init sequences and program parameters first (option 2).\n\n";
                Pause();
                break;
            }
            ExecuteResampling();
            Pause();
            break;
        }
        case '?':
        case 'H': {
            cout<<"\n\nHelp not available yet\n\n";
            Pause();
            break;
        }
        case 'Q':
        case 'X':
            exit(0);
            break;
        default:
            cout<<"\n\nUnrecognised Command\n\n";
            Pause();
            break;
        }
    }
}

//Function to initialize properly the TripleInit object
void Menus::InitTripleInitializer(int SubTree, bool LoadStartingHomologies) {
    if (LoadStartingHomologies){
        LoadNewHomology();
    }
    TripleInit->MakeProperInitializationsAll(SubTree, NewHomology);
}

//Function to initialize properly the PairInit object
void Menus::InitPairInitializer(int Branch, bool LoadStartingHomologies) {
    if (LoadStartingHomologies){
        LoadNewHomology();
    }
    PairInit->MakeProperInitializations(Branch, NewHomology);
}

//Function to start the re-alignment
void Menus::ExecuteResampling() {
    //Check if the Status File exists and if it exists then delete it
    if (FileExists(AlignerStatusFile)) {
        DeleteFile(AlignerStatusFile);
    }

    ExecuteTriplewiseResampling(10);
    ExecutePairwiseResampling(5);
    ExecuteTriplewiseResampling(3);
    ExecutePairwiseResampling(3);
    //Reset Homologies
    NewHomology.clear();
    ExecuteTriplewiseResampling(9);
    ExecutePairwiseResampling(5);
    ExecuteTriplewiseResampling(10);
    ExecutePairwiseResampling(5);
    ExecuteTriplewiseResampling(8);

    //Now that the alignment is completed write result to Status File
    stringstream ss; ss.str("");
    ss<<"Alignment Complete - Rate: "<<TripleInit->Get_Parameters()->Get_Lambda()<<endl;

    time_t t = time(0);   // get time now
    struct tm * now = localtime( & t );
    ss<<"Timestamp: " <<  now->tm_mday << '/'
        << (now->tm_mon + 1) << '/'
        << (now->tm_year + 1900) << '\t'
        <<  now->tm_hour << ':'
        <<  now->tm_min  << ':' << now->tm_sec
        << endl;
    ss<<"\n\nSequences:\n"<<(*(TripleInit->Get_Sequences()))<<endl;
    WriteToFile(AlignerStatusFile, ss.str(), false, true);
}

//Function to Perform pairwise re-alignment
void Menus::ExecutePairwiseResampling(int rounds) {
    bool DebugMode              = false;
    bool VerifySeqs             = true;
    bool LogProbs               = false;
    bool StopAtIntervals        = false;
    bool StopAtBranches         = false;
    bool ShowNewEvolHistory     = true;
    bool ShowNewHomology        = false;
    bool ShowNewAlignment       = true;
    bool SaveNewAlignment       = true;
    bool ShowEvaluation         = false;
    bool UseCompare             = false;
    bool SaveStartingHomologies = true;
    bool LoadStartingHomologies = false;

    //Execute the predefined number of Steps to re-align sequences
    for (int k=1; k<=rounds; k++) {
        cout<<"\n\nRunning Pairwise for Nodes: ";
        if (k==1) {
            cout<<"E-F (E is root)"<<endl;
        } else if (k==2) {
            cout<<"E-A (E is root)"<<endl;
        } else if (k==3) {
            cout<<"E-B (E is root)"<<endl;
        } else if (k==4) {
            cout<<"F-C (F is root)"<<endl;
        } else if (k==5) {
            cout<<"F-D (F is root)"<<endl;
        } else if (k==6) {
            cout<<"F-E (F is root)"<<endl;
        }

        //1. Initialize the "Initializer" instance containing all the relevant information for the alignment region //{
        delete PairInit;
        PairInit  = new PairInitializer();
        InitPairInitializer(k, LoadStartingHomologies);
        //}

        //2. Find the Start and End Value of the region to be resampled //{
        int   RegionEnd = PairInit->GetRootMaxValue();
        int RegionStart = max(PairInit->Get_RandomGenerator()->GetRandomInt(0, PairInit->Get_Parameters()->Get_AlignFragSize()-2), 0);
        //RegionStart = 9;
        cout<<"RegionStart:"<<RegionStart<<endl;
        //}

        //3. Initialize some variables used to keep track of the newly estimated alignment //{
        //Also, for each convergence step, clear the newly estimated Homology structure
        NewEvolHistory.clear();
        NewHomology.clear();
        NewAlignment.clear();

        //And create a new EvolHistory and Homology map object
        map<int, StringVec> EvolHistory;
        //}

        //4. Form the corresponding intervals that should be examined (i.e. segment the alignment region in intervals according to the specified fragment size) //{
        //Create the Corresponding Fragments based on the Alignment and the Desired Fragment Size
        int FragmentSize = PairInit->Get_Parameters()->Get_AlignFragSize();
        FragmentSize = 200;
        IntervalVec intervals;

        for (int i=RegionStart; i<RegionEnd; i+=FragmentSize) {
            Interval itv;
            itv.start = i;

            if ( (i+FragmentSize-1)>RegionEnd-3 ) {
                itv.end = RegionEnd;
            } else {
                itv.end = (i+FragmentSize-1);
            }
            intervals.push_back(itv);

            if (itv.end==RegionEnd) break;
        }

        if (intervals.size()==0) {
            cout<<"\n\n Error:No intervals were created!"<<endl;
            return;
        }

        int    RegionStartPosIdx = PairInit->GetFirstIdxofSoanWhereRootIs(intervals[0].start);
        IntVec RegionStartIdxs;
        bool RegionHasExtraStartState, RegionHasExtraEndState;
        if (RegionStart==0) {
            RegionStartIdxs.resize(4, 1);
            RegionHasExtraStartState = true;
        } else {
            RegionStartIdxs  = PairInit->GetSoanIdxsAtPosition(RegionStartPosIdx);
            RegionHasExtraStartState = false;
        }

        if (RegionEnd==PairInit->GetRootMaxValue()) {
            RegionHasExtraEndState = true;
        } else {
            RegionHasExtraEndState = false;
        }
        string RegionStartTihlLbl = PairInit->GetTihlAtPosition(RegionStartPosIdx, true);

        if (DebugMode) {
            PrintIntervals(intervals);
            PairInit->SaveRegionHomologiesToFile();
            PairInit->SaveAlignmentHomologiesToFile();
            PairInit->SaveHomologiesToFile();
            PairInit->SaveTreeToFile();
        }
        //}

        //5. Resample each interval //{
        //omp_set_num_threads( 4 );
        #pragma omp parallel for
        for (int i=0; i<intervals.size(); i++) {
            //Proceed to sample the region of interest ignoring the given homologies
            StringVec IntervalEvolHist;
            double prob;

            //Some initializations concerning the current interval //{
            //a. The start and end positions of the interval
            int StartPos = PairInit->GetFirstIdxofSoanWhereRootIs(intervals[i].start);
            int   EndPos = PairInit->GetLastIdxofSoanWhereRootIs(intervals[i].end);

            //b. The starting tihl label and soan label
            string StartTihlLbl = PairInit->GetTihlAtPosition(StartPos, false);
            string StartSoanLbl = PairInit->GetSoanAtPosition(StartPos, false);

            //c. The starting and ending soan indices
            IntVec StartIdxs  = PairInit->GetSoanIdxsAtPosition(StartPos);
            IntVec   EndIdxs  = PairInit->GetSoanIdxsAtPosition(EndPos);

            //d. The dimensions corresponding to the current interval
            IntVec Dimensions(EndIdxs.begin(), EndIdxs.end()-1) ;
            for (int j=0; j<Dimensions.size(); j++) {
                Dimensions[j]++;
            }

            //e. The lengths of the current interval
            DoubleVec Lengths(PairInit->Get_Lengths().begin(), PairInit->Get_Lengths().end()) ;

            //f. The Felsenstein pruning tables of the current interval
            DoubleHyp FPTables = PairInit->Get_FPTables();

            //g. The substitution matrices corresponding to the lengths of the current interval
            DoubleMat SubMats = PairInit->Get_SubMatrices();

            //h. Definition of whether this interval has extra start/end states
            bool HasExtraStartState = false;
            bool HasExtraEndState   = false;

            if (intervals[i].start==0)         {
                HasExtraStartState = true;
            }
            if (intervals[i].end  ==RegionEnd) {
                HasExtraEndState   = true;
            }


            if (DebugMode) {
                PrintAlignmentRegion(StartIdxs, EndIdxs, StartTihlLbl, StartSoanLbl, HasExtraStartState, HasExtraEndState, Dimensions, Lengths);
            }
            //}

            PairAligner* aligner = new PairAligner();
            if (LogProbs) {
                aligner->Set_SaveToLogFile(true, LogFile);
            }
            aligner->Set_RandomGenerator(PairInit->Get_RandomGenerator());
            aligner->InitParameters(PairInit->Get_Parameters());
            aligner->InitSubModel(PairInit->Get_F84Model(), PairInit->Get_TN93Model(), PairInit->Get_GTRModel());
            aligner->InitFidModel();
            aligner->InitLengths(Lengths);
            aligner->InitDimensions(Dimensions);
            aligner->InitFPTables(FPTables);
            if (DebugMode) {
                aligner->SaveFPTables();
            }
            aligner->InitTriplets();
            aligner->Set_EndPosIn_ND(EndIdxs);
            aligner->InitProbabilities(StartSoanLbl, StartTihlLbl, StartIdxs);
            aligner->InitSubMatrices(SubMats);
            unsigned t_start, t_end;
            if (LogProbs) {
                cout<<"\n\n Resampling new Alignment...";
                t_start=clock();
            }
            //aligner->PrintStatus();
            //Pause();

            aligner->EstimateProbabilities(0, HasExtraEndState);

            if (LogProbs) {
                t_end=clock()-t_start;
                cout<<" -> took: " << (((float)t_end)/CLOCKS_PER_SEC) <<" seconds."<<endl<<endl;
                aligner->SaveResultsToFile(ProbabilitiesFile);
            }
            aligner->GetEvolHistory(IntervalEvolHist, prob, true);

            if (DebugMode) {
                cout<<"\n\nSAMPLED REGION:\n"<<IntervalEvolHist<<endl;
            }

            if (HasExtraStartState) {
                if (DebugMode) {
                    cout<<"Deleted 1st Column -> HasExtraStartState"<<endl;
                }
                for (int j=0; j<IntervalEvolHist.size(); j++) {
                    IntervalEvolHist[j].erase(IntervalEvolHist[j].begin());
                }
            }

            if (HasExtraEndState) {
                if (DebugMode) {
                    cout<<"Deleted last Column -> HasExtraEndState"<<endl;
                }
                for (int j=0; j<IntervalEvolHist.size(); j++) {
                    IntervalEvolHist[j].erase(IntervalEvolHist[j].end()-1);
                }
            }


            if (DebugMode) {
                cout<<"\nKEPT REGION:\n"<<IntervalEvolHist<<endl;
            }

            EvolHistory[i] = IntervalEvolHist;
            delete aligner;

            //For Debugging: Tests to Check whether the deriving alignment leads to more/less sites than it should //{
            if (VerifySeqs) {
                IntVec ExpectedSites;

                for (int l=0; l<EndIdxs.size()-1; l++) {
                    if (HasExtraEndState) {
                        EndIdxs[l]--;
                    }

                    if (HasExtraStartState) {
                        StartIdxs[l]++;
                    }
                    ExpectedSites.push_back(EndIdxs[l] - StartIdxs[l] +1);
                }

                if (StartTihlLbl.compare("H")==0 || StartTihlLbl.compare("E")==0  || StartTihlLbl.compare("N")==0 ) {
                    StartTihlLbl.insert(StartTihlLbl.begin(), 'H');
                }else{
                    StartTihlLbl.insert(StartTihlLbl.begin(), '-');
                }

                for (int l=0; l<StartTihlLbl.size(); l++) {
                    if (StartTihlLbl.at(l)=='E' || StartTihlLbl.at(l)=='-') {
                        ExpectedSites[l]--;
                    }
                }


                for (int l=0; l<IntervalEvolHist.size(); l++) {
                    int Hs = StringCount(IntervalEvolHist[l], "H");
                    int Ns = StringCount(IntervalEvolHist[l], "N");
                    int Bs = StringCount(IntervalEvolHist[l], "B");
                    int sum = Hs + Ns + Bs;

                    if (sum!=ExpectedSites[l]) {
                        cout<<"Inconistency detected for current fragment -Node:"<<l<<endl;
                        cout<<"Expected " <<ExpectedSites[l]<<" sites but resampled "<< sum<<" sites!"<<endl;
                        exit(1);
                    }

                }
            }
            //}

            if (StopAtIntervals) {
                if (i<intervals.size()-1) {
                    cout<<"\n\nProceed with the next interval ("<<(i+1)<<")?";
                    Pause();
                }
            }
        }
        //}

        //6. After all intervals have been properly re-aligned, create the new Evolution History, Homology Structure and Alignment and then Update the new Alignment to the corresponding file //{
        //For Debugging - print each part of the newly sampled EvolHistory:
        /*
        for( map<int, StringVec>::iterator ii=EvolHistory.begin(); ii!=EvolHistory.end(); ++ii) {
                cout<< (*ii).first <<"\t"<< ((*ii).second)[0].size() <<endl;
                cout<< (*ii).second <<endl;
         }
         PrintIntervals(intervals);
        */

        int counter = 0;
        for( map<int, StringVec>::iterator ii=EvolHistory.begin(); ii!=EvolHistory.end(); ++ii) {
            StringVecRef temp = (*ii).second;
            if (counter++==0) {
                NewEvolHistory.resize(temp.size(), "");
            }
            for (int j=0; j<temp.size(); j++) {
                stringstream ss;
                ss.str("");
                ss<<NewEvolHistory[j]<<temp[j];
                NewEvolHistory[j] = ss.str();
            }
        }

        if (ShowNewEvolHistory) {
            cout<<"\n\nNew Evolution History:\n----------------------\n"<<NewEvolHistory<<endl;
        }

        GetNewHomologyPair(k, RegionStartIdxs, RegionStartTihlLbl, RegionHasExtraStartState, RegionHasExtraEndState);

        if (ShowNewHomology) {
            PrintNewHomology();
        }

        if (SaveStartingHomologies) {
            SaveNewHomology();
        }

        GetNewAlignment();
        if (ShowNewAlignment) {
            PrintNewAlignment();
        }

        //For Debugging: Tests to Check whether the deriving alignment leads to more/less sites than it should //{
        if (VerifySeqs) {
            for (int i=0; i<NewAlignment.size(); i++) {
                string tempS = NewAlignment[i];
                ReplaceAll(tempS, "-", "");
                if (tempS.compare(PairInit->Get_Sequences()->Get_BareSeqs()[i])!=0) {
                    cout<<"\n\nPROBLEM: The Newly Estimated Alignment leads to Different Sequences than the Original ones! (Sequence:"<<(i+1)<<")"<<endl;
                    cout<<"Original Sequence: "<<PairInit->Get_Sequences()->Get_BareSeqs()[i]<<endl;
                    cout<<"     New Sequence: "<<tempS<<endl;
                    exit(1);
                }
            }
        }
        //}

        if (SaveNewAlignment) {
            UpdateNewAlignment();
        }

        if (ShowEvaluation) {
            StringVec algn,    refr,    extr;
            IntMat    algnMat, refrMat, extrMat;

            GetAlignmentFromFile(algn, AlignedSequencesFile);
            GetAlignmentFromFile(refr, ReferenceSequencesFile);
            GetAlignmentFromFile(extr, ExtraAlignSequencesFile);

            EvaluationObj eval;
            eval.ConvertAlignmentToMatrix(algn, algnMat);
            eval.ConvertAlignmentToMatrix(refr, refrMat);
            if (UseCompare) { eval.ConvertAlignmentToMatrix(extr, extrMat); }

            cout<<"Predicted Alignment against Reference Alignment:"<<endl;
            cout<<"\t"<<eval.GetMCCScore(algnMat, refrMat);
            cout<<"\t\n"<<eval.GetTCScore( algn,    refr);
            cout<<"\t\nPredicted Alignment:"<<endl;
            cout<<"\t"<<eval.GetSPSimScore(algn);
            cout<<"\t"<<eval.GetSPDistScore(algn);
            cout<<"\t\nReference Alignment:"<<endl;
            cout<<"\t"<<eval.GetSPSimScore(refr);
            cout<<"\t"<<eval.GetSPDistScore(refr);

            if (UseCompare) {
                cout<<"Extra Alignment against Reference Alignment:"<<endl;
                cout<<"\t"<<eval.GetMCCScore(algnMat, refrMat);
                cout<<"\t\n"<<eval.GetTCScore( algn,    refr);
                cout<<"\t\nExtra Alignment:"<<endl;
                cout<<"\t"<<eval.GetSPSimScore(algn);
                cout<<"\t"<<eval.GetSPDistScore(algn);
                cout<<"\t\nReference Alignment:"<<endl;
                cout<<"\t"<<eval.GetSPSimScore(refr);
                cout<<"\t"<<eval.GetSPDistScore(refr);
            }
        }

        //}

        if (StopAtBranches) {
            if (k<=PairInit->Get_NumOfSteps()) {
                cout<<"\n\nProceed with the next Branch ("<<(k+1)<<")?";
                Pause();
            }
        }
    }
}

//Function to Perform the left/right subtree triplewise re-alignment
void Menus::ExecuteTriplewiseResampling(int rounds) {
    bool DebugMode              = false;
    bool VerifySeqs             = true;
    bool LogProbs               = false;
    bool StopAtIntervals        = false;
    bool StopAtSubtrees         = false;
    bool ShowNewEvolHistory     = true;
    bool ShowNewHomology        = false;
    bool ShowNewAlignment       = true;
    bool SaveNewAlignment       = true;
    bool ShowEvaluation         = false;
    bool UseCompare             = false;
    bool SaveStartingHomologies = true;
    bool LoadStartingHomologies = false;


    //Execute the predefined number of Steps to re-align sequences
    for (int k=1; k<=rounds; k++) {
        cout<<"\n\nRunning Step: "<<k;
        if (k%2==1) {
            cout<<" - Left Subtree"<<endl;
        } else {
            cout<<" - Right Subtree"<<endl;
        }

        //1. Initialize the "Initializer" instance containing all the relevant information for the alignment region //{
        delete TripleInit;//comment this line for the Batch Benchmarker program
        TripleInit  = new TripleInitializer();
        InitTripleInitializer((k%2), LoadStartingHomologies);
        //}

        //2. Find the Start and End Value of the region to be resampled //{
        int   RegionEnd = TripleInit->GetRootMaxValue();
        //int RegionStart = TripleInit->Get_RandomGenerator()->GetRandomInt(0, (int)(RegionEnd-3));
        int RegionStart = max(TripleInit->Get_RandomGenerator()->GetRandomInt(0, TripleInit->Get_Parameters()->Get_AlignFragSize()-2), 0);
        //int RegionStart = TripleInit->Get_RandomGenerator()->GetRandomInt(0, 4);
        //RegionStart = 79;
        cout<<"RegionStart:"<<RegionStart<<endl;
        //}

        //3. Initialize some variables used to keep track of the newly estimated alignment //{
        //Also, for each convergence step, clear the newly estimated Homology structure
        NewEvolHistory.clear();
        NewHomology.clear();
        NewAlignment.clear();

        //And create a new EvolHistory and Homology map object
        map<int, StringVec> EvolHistory;

        //}

        //4. Form the corresponding intervals that should be examined (i.e. segment the alignment region in intervals according to the specified fragment size) //{
        //Create the Corresponding Fragments based on the Alignment and the Desired Fragment Size
        int FragmentSize = TripleInit->Get_Parameters()->Get_AlignFragSize();
        FragmentSize = 8;
        IntervalVec intervals;

        for (int i=RegionStart; i<RegionEnd; i+=FragmentSize) {
            Interval itv;
            itv.start = i;

            if ( (i+FragmentSize-1)>RegionEnd-3 ) {
                itv.end = RegionEnd;
            } else {
                itv.end = (i+FragmentSize-1);
            }
            intervals.push_back(itv);

            if (itv.end==RegionEnd) break;
        }

        if (intervals.size()==0) {
            cout<<"\n\n Error:No intervals were created!"<<endl;
            return;
        }

        int    RegionStartPosIdx = TripleInit->GetFirstIdxofSoanWhereRootIs(intervals[0].start);
        IntVec RegionStartIdxs;
        bool RegionHasExtraStartState, RegionHasExtraEndState;
        if (RegionStart==0) {
            RegionStartIdxs.resize(4, 1);
            RegionHasExtraStartState = true;
        } else {
            RegionStartIdxs  = TripleInit->GetSoanIdxsAtPosition(RegionStartPosIdx);
            RegionHasExtraStartState = false;
        }

        if (RegionEnd==TripleInit->GetRootMaxValue()) {
            RegionHasExtraEndState = true;
        } else {
            RegionHasExtraEndState = false;
        }
        string RegionStartTihlLbl = TripleInit->GetTihlAtPosition(RegionStartPosIdx, true);

        if (DebugMode) {
            PrintIntervals(intervals);
            TripleInit->SaveRegionHomologiesToFile();
            TripleInit->SaveAlignmentHomologiesToFile();
            TripleInit->SaveHomologiesToFile();
            TripleInit->SaveTreeToFile();
        }
        //}

        //5. Resample each interval //{
        //omp_set_num_threads( 4 );
        #pragma omp parallel for
        for (int i=0; i<intervals.size(); i++) {
            //Proceed to sample the region of interest ignoring the given homologies
            StringVec IntervalEvolHist;
            double prob;

            //Some initializations concerning the current interval //{
            //a. The start and end positions of the interval
            int StartPos = TripleInit->GetFirstIdxofSoanWhereRootIs(intervals[i].start);
            int   EndPos = TripleInit->GetLastIdxofSoanWhereRootIs(intervals[i].end);

            //b. The starting tihl label and soan label
            string StartTihlLbl = TripleInit->GetTihlAtPosition(StartPos, false);
            string StartSoanLbl = TripleInit->GetSoanAtPosition(StartPos, false);

            //c. The starting and ending soan indices
            IntVec StartIdxs  = TripleInit->GetSoanIdxsAtPosition(StartPos);
            IntVec   EndIdxs  = TripleInit->GetSoanIdxsAtPosition(EndPos);

            //d. The dimensions corresponding to the current interval
            IntVec Dimensions(EndIdxs.begin()+1, EndIdxs.end()-1) ;
            for (int j=0; j<Dimensions.size(); j++) {
                Dimensions[j]++;
            }

            //e. The lengths of the current interval
            DoubleVec Lengths(TripleInit->Get_Lengths().begin()+1, TripleInit->Get_Lengths().end()) ;

            //f. The Felsenstein pruning tables of the current interval
            DoubleHyp FPTables = TripleInit->Get_FPTables();
            FPTables.erase(FPTables.begin());

            //g. The substitution matrices corresponding to the lengths of the current interval
            DoubleMat SubMats = TripleInit->Get_SubMatrices();
            SubMats.erase(SubMats.begin());

            //h. Definition of whether this interval has extra start/end states
            bool HasExtraStartState = false;
            bool HasExtraEndState   = false;

            if (intervals[i].start==0)         {
                HasExtraStartState = true;
            }
            if (intervals[i].end  ==RegionEnd) {
                HasExtraEndState   = true;
            }

            //i. Fix for the Start/End indices (deleting the 1st index corresponding to the root)
            StartIdxs.erase(StartIdxs.begin());
            EndIdxs.erase(EndIdxs.begin());

            if (DebugMode) {
                PrintAlignmentRegion(StartIdxs, EndIdxs, StartTihlLbl, StartSoanLbl, HasExtraStartState, HasExtraEndState, Dimensions, Lengths);
            }
            //}

            TripleAligner* aligner = new TripleAligner();
            if (LogProbs) {
                aligner->Set_SaveToLogFile(true, LogFile);
            }
            aligner->Set_RandomGenerator(TripleInit->Get_RandomGenerator());
            aligner->InitParameters(TripleInit->Get_Parameters());
            aligner->InitSubModel(TripleInit->Get_F84Model(), TripleInit->Get_TN93Model(), TripleInit->Get_GTRModel());
            aligner->InitFidModel();
            aligner->InitLengths(Lengths);
            aligner->InitDimensions(Dimensions);
            aligner->InitFPTables(FPTables);
            if (DebugMode) {
                aligner->SaveFPTables();
            }
            aligner->InitTriplets();
            aligner->Set_EndPosIn_ND(EndIdxs);
            aligner->InitProbabilities(StartSoanLbl, StartTihlLbl, StartIdxs);
            aligner->InitSubMatrices(SubMats);
            unsigned t_start, t_end;
            if (LogProbs) {
                cout<<"\n\n Resampling new Alignment...";
                t_start=clock();
            }

            //aligner->PrintStatus();
            //Pause();

            aligner->EstimateProbabilities(0, HasExtraEndState);

            if (LogProbs) {
                t_end=clock()-t_start;
                cout<<" -> took: " << (((float)t_end)/CLOCKS_PER_SEC) <<" seconds."<<endl<<endl;
                aligner->SaveResultsToFile(ProbabilitiesFile);
            }
            aligner->GetEvolHistory(IntervalEvolHist, prob, true);

            if (DebugMode) {
                cout<<"\n\nSAMPLED REGION:\n"<<IntervalEvolHist<<endl;
            }

            if (HasExtraStartState) {
                if (DebugMode) {
                    cout<<"Deleted 1st Column -> HasExtraStartState"<<endl;
                }
                for (int j=0; j<IntervalEvolHist.size(); j++) {
                    IntervalEvolHist[j].erase(IntervalEvolHist[j].begin());
                }
            }

            if (HasExtraEndState) {
                if (DebugMode) {
                    cout<<"Deleted last Column -> HasExtraEndState"<<endl;
                }
                for (int j=0; j<IntervalEvolHist.size(); j++) {
                    IntervalEvolHist[j].erase(IntervalEvolHist[j].end()-1);
                }
            }


            if (DebugMode) {
                cout<<"\nKEPT REGION:\n"<<IntervalEvolHist<<endl;
            }

            EvolHistory[i] = IntervalEvolHist;
            delete aligner;

            //For Debugging: Tests to Check whether the deriving alignment leads to more/less sites than it should //{
            if (VerifySeqs) {
                IntVec ExpectedSites;
                for (int l=0; l<EndIdxs.size()-1; l++) {
                    if (HasExtraEndState) {
                        EndIdxs[l]--;
                    }

                    if (HasExtraStartState) {
                        StartIdxs[l]++;
                    }
                    ExpectedSites.push_back(EndIdxs[l] - StartIdxs[l] +1);
                }

                for (int l=0; l<StartTihlLbl.size(); l++) {
                    if (StartTihlLbl.at(l)=='E' || StartTihlLbl.at(l)=='-') {
                        ExpectedSites[l]--;
                    }
                }

                for (int l=1; l<IntervalEvolHist.size(); l++) {
                    int Hs = StringCount(IntervalEvolHist[l], "H");
                    int Ns = StringCount(IntervalEvolHist[l], "N");
                    int Bs = StringCount(IntervalEvolHist[l], "B");
                    int sum = Hs + Ns + Bs;

                    if (sum!=ExpectedSites[l-1]) {
                        char node;
                        if (k%2==0) {
                            if      (l==1) {
                                node = 'C';
                            } else if (l==2) {
                                node = 'D';
                            } else if (l==3) {
                                node = 'E';
                            }
                        } else if (k%2==1) {
                            if      (l==1) {
                                node = 'A';
                            } else if (l==2) {
                                node = 'B';
                            } else if (l==3) {
                                node = 'F';
                            }
                        }
                        cout<<"Inconistency detected for current fragment -Node:"<<node<<endl;
                        cout<<"Expected " <<ExpectedSites[l-1]<<" sites but resampled "<< sum<<" sites!"<<endl;
                        exit(1);
                    }

                }
            }
            //}

            if (StopAtIntervals) {
                if (i<intervals.size()-1) {
                    cout<<"\n\nProceed with the next interval ("<<(i+1)<<")?";
                    Pause();
                }
            }
        }
        //}

        //6. After all intervals have been properly re-aligned, create the new Evolution History, Homology Structure and Alignment and then Update the new Alignment to the corresponding file //{
        //For Debugging - print each part of the newly sampled EvolHistory:
        /*
        for( map<int, StringVec>::iterator ii=EvolHistory.begin(); ii!=EvolHistory.end(); ++ii) {
                cout<< (*ii).first <<"\t"<< ((*ii).second)[0].size() <<endl;
                cout<< (*ii).second <<endl;
         }
         PrintIntervals(intervals);
        */

        int counter = 0;
        for( map<int, StringVec>::iterator ii=EvolHistory.begin(); ii!=EvolHistory.end(); ++ii) {
            StringVecRef temp = (*ii).second;
            if (counter++==0) {
                NewEvolHistory.resize(temp.size(), "");
            }
            for (int j=0; j<temp.size(); j++) {
                stringstream ss;
                ss.str("");
                ss<<NewEvolHistory[j]<<temp[j];
                NewEvolHistory[j] = ss.str();
            }
        }

        if (ShowNewEvolHistory) {
            cout<<"\n\nNew Evolution History:\n----------------------\n"<<NewEvolHistory<<endl;
        }
        GetNewHomologyTriple((k)%2, RegionStartIdxs, RegionStartTihlLbl, RegionHasExtraStartState, RegionHasExtraEndState);
        if (ShowNewHomology) {
            PrintNewHomology();
        }

        if (SaveStartingHomologies) {
            SaveNewHomology();
        }

        GetNewAlignment();
        if (ShowNewAlignment) {
            PrintNewAlignment();
        }

        //For Debugging: Tests to Check whether the deriving alignment leads to more/less sites than it should //{
        if (VerifySeqs) {
            for (int i=0; i<NewAlignment.size(); i++) {
                string tempS = NewAlignment[i];
                ReplaceAll(tempS, "-", "");
                if (tempS.compare(TripleInit->Get_Sequences()->Get_BareSeqs()[i])!=0) {
                    cout<<"\n\nPROBLEM: The Newly Estimated Alignment leads to Different Sequences than the Original ones! (Sequence:"<<(i+1)<<")"<<endl;
                    cout<<"Original Sequence: "<<TripleInit->Get_Sequences()->Get_BareSeqs()[i]<<endl;
                    cout<<"     New Sequence: "<<tempS<<endl;
                    exit(1);
                }
            }
        }
        //}

        if (SaveNewAlignment) {
            UpdateNewAlignment();
        }

        if (ShowEvaluation) {
            StringVec algn,    refr,    extr;
            IntMat    algnMat, refrMat, extrMat;

            GetAlignmentFromFile(algn, AlignedSequencesFile);
            GetAlignmentFromFile(refr, ReferenceSequencesFile);
            GetAlignmentFromFile(extr, ExtraAlignSequencesFile);

            EvaluationObj eval;
            eval.ConvertAlignmentToMatrix(algn, algnMat);
            eval.ConvertAlignmentToMatrix(refr, refrMat);
            if (UseCompare) { eval.ConvertAlignmentToMatrix(extr, extrMat); }

            cout<<"Predicted Alignment against Reference Alignment:"<<endl;
            cout<<"\t"<<eval.GetMCCScore(algnMat, refrMat);
            cout<<"\t\n"<<eval.GetTCScore( algn,    refr);
            cout<<"\t\nPredicted Alignment:"<<endl;
            cout<<"\t"<<eval.GetSPSimScore(algn);
            cout<<"\t"<<eval.GetSPDistScore(algn);
            cout<<"\t\nReference Alignment:"<<endl;
            cout<<"\t"<<eval.GetSPSimScore(refr);
            cout<<"\t"<<eval.GetSPDistScore(refr);

            if (UseCompare) {
                cout<<"Extra Alignment against Reference Alignment:"<<endl;
                cout<<"\t"<<eval.GetMCCScore(algnMat, refrMat);
                cout<<"\t\n"<<eval.GetTCScore( algn,    refr);
                cout<<"\t\nExtra Alignment:"<<endl;
                cout<<"\t"<<eval.GetSPSimScore(algn);
                cout<<"\t"<<eval.GetSPDistScore(algn);
                cout<<"\t\nReference Alignment:"<<endl;
                cout<<"\t"<<eval.GetSPSimScore(refr);
                cout<<"\t"<<eval.GetSPDistScore(refr);
            }
        }

        //}

        if (StopAtSubtrees) {
            if (k<=TripleInit->Get_NumOfSteps()) {
                cout<<"\n\nProceed with the next Subtree (";
                if (k%2==0) {
                    cout<<"Left)?";
                } else {
                    cout<<"Right)?";
                }

                Pause();
            }
        }
    }
}

//Function to Expand the Ns of the Newly Estimated Evolution History
void Menus::ExpandNs() {
    stringstream ss;
    for (int j=0; j<NewEvolHistory[0].size(); j++) {
        ss.str("");
        ss<<NewEvolHistory[0][j];
        ss<<NewEvolHistory[1][j];
        ss<<NewEvolHistory[2][j];

        if (StringCount(ss.str(), "N")>0) {
            if (ss.str().compare("NHH")==0) {
                NewEvolHistory[0][j] = '-';
                NewEvolHistory[0].insert(j+1, "N");
                NewEvolHistory[1].insert(j+1, "-");
                NewEvolHistory[2].insert(j+1, "-");
            } else if (ss.str().compare("HNH")==0) {
                NewEvolHistory[1][j] = '-';
                NewEvolHistory[0].insert(j+1, "-");
                NewEvolHistory[1].insert(j+1, "N");
                NewEvolHistory[2].insert(j+1, "-");
            } else if (ss.str().compare("HHN")==0) {
                NewEvolHistory[2][j] = '-';
                NewEvolHistory[0].insert(j+1, "-");
                NewEvolHistory[1].insert(j+1, "-");
                NewEvolHistory[2].insert(j+1, "N");
            } else if (ss.str().compare("NNH")==0) {
                NewEvolHistory[0][j] = '-';
                NewEvolHistory[1][j] = '-';
                NewEvolHistory[0].insert(j+1, "N-");
                NewEvolHistory[1].insert(j+1, "-N");
                NewEvolHistory[2].insert(j+1, "--");
            } else if (ss.str().compare("NHN")==0) {
                NewEvolHistory[0][j] = '-';
                NewEvolHistory[2][j] = '-';
                NewEvolHistory[0].insert(j+1, "N-");
                NewEvolHistory[1].insert(j+1, "--");
                NewEvolHistory[2].insert(j+1, "-N");
            } else if (ss.str().compare("HNN")==0) {
                NewEvolHistory[1][j] = '-';
                NewEvolHistory[2][j] = '-';
                NewEvolHistory[0].insert(j+1, "--");
                NewEvolHistory[1].insert(j+1, "N-");
                NewEvolHistory[2].insert(j+1, "-N");
            } else if (ss.str().compare("NNN")==0) {
                NewEvolHistory[0][j] = '-';
                NewEvolHistory[1][j] = '-';
                NewEvolHistory[2][j] = '-';
                NewEvolHistory[0].insert(j+1, "N--");
                NewEvolHistory[1].insert(j+1, "-N-");
                NewEvolHistory[2].insert(j+1, "--N");
            }
        }
    }
}

//Function to Estimate the new Homology Structure based on the re-sampled evolution history
void Menus::GetNewHomologyTriple(int SubTree, IntVecRef RegionStartIdxs, string RegionStartTihlLbl, bool HasExtraStartState, bool HasExtraEndState) {
    for (int i=0; i<RegionStartTihlLbl.size(); i++) {
        if (RegionStartTihlLbl.at(i)=='E' || RegionStartTihlLbl.at(i)=='-') {
            RegionStartIdxs[i]++;
        }
    }

    switch(SubTree) {
    case 1: {   //Left Subtree
        //1. Convert the sampled Evolution History to the Corresponding Homology indices //{
        NewHomology.resize(6, IntVec(0,0));
        if (!HasExtraStartState) {
            int beforePos = TripleInit->GetLastIdxofSoanWhereRootIs(RegionStartIdxs[0]-1);

            for (int k=1; k<=beforePos; k++) {
                string tihl = TripleInit->GetTihlAtPosition(k, true);
                for (int j=0; j<tihl.size(); j++) {
                    int H_idx;
                    if      (j==0) {
                        H_idx = 4;
                    } else if (j==1) {
                        H_idx = 0;
                    } else if (j==2) {
                        H_idx = 1;
                    } else if (j==3) {
                        H_idx = 5;
                    }

                    if (tihl.at(j)=='E' || tihl.at(j)=='-') {
                        NewHomology[H_idx].push_back(-1);
                    } else {
                        NewHomology[H_idx].push_back(TripleInit->GetSoanIdxsAtPosition(k)[j]);
                    }
                }
            }
        }

        for (int i=0; i<NewEvolHistory.size(); i++) {
            for (int j=0; j<NewEvolHistory[i].size(); j++) {
                int H_idx;
                if      (i==0) {
                    H_idx = 4;
                } else if (i==1) {
                    H_idx = 0;
                } else if (i==2) {
                    H_idx = 1;
                } else if (i==3) {
                    H_idx = 5;
                }

                if (NewEvolHistory[i].at(j)=='E' || NewEvolHistory[i].at(j)=='-') {
                    NewHomology[H_idx].push_back(-1);
                } else {
                    NewHomology[H_idx].push_back(RegionStartIdxs[i]);
                    RegionStartIdxs[i]++;
                }
            }
        }

        if (!HasExtraEndState) {
            int afterPos = TripleInit->GetFirstIdxofSoanWhereRootIs(RegionStartIdxs[0]+1);

            for (int k=afterPos; k<TripleInit->Get_SoansIdxs().size()-1; k++) {
                string tihl = TripleInit->GetTihlAtPosition(k, true);
                for (int j=0; j<tihl.size(); j++) {
                    int H_idx;
                    if      (j==0) {
                        H_idx = 4;
                    } else if (j==1) {
                        H_idx = 0;
                    } else if (j==2) {
                        H_idx = 1;
                    } else if (j==3) {
                        H_idx = 5;
                    }

                    if (tihl.at(j)=='E' || tihl.at(j)=='-') {
                        NewHomology[H_idx].push_back(-1);
                    } else {
                        NewHomology[H_idx].push_back(TripleInit->GetSoanIdxsAtPosition(k)[j]);
                    }
                }
            }
        }
        //}

        //2. Form the "Skipped nodes" indices based on their parental node //{
        NewHomology[2].resize(NewHomology[5].size(), -1);
        NewHomology[3].resize(NewHomology[5].size(), -1);

        for (int j=0; j<NewHomology[5].size(); j++) {
            if (NewHomology[5][j]!=-1) {
                int pos = GetPosOfElement(TripleInit->Get_AlignmentHomologies()[5], NewHomology[5][j], 0);
                int val1 = TripleInit->Get_AlignmentHomologies()[2][pos];
                NewHomology[2][j] = val1;

                int val2 = TripleInit->Get_AlignmentHomologies()[3][pos];
                NewHomology[3][j] = val2;
            }
        }
        //}

        //3. Add the singleton values of the skipped nodes  //{
        //Node C Singletons
        //int CurValue = GetFirstNonNegativeValueReverse(TripleInit->Get_AlignmentHomologies()[2], TripleInit->Get_AlignmentHomologies()[2].size()-1);
        int CurValue = TripleInit->Get_Sequences()->GetBareSequenceAt(2).size();
        for (int j=NewHomology[2].size()-1; j>=0; j--) {
            if (NewHomology[2][j]==-1)       {
                continue;
            }
            if (NewHomology[2][j]==CurValue) {
                CurValue--;
                continue;
            }

            for (int k=CurValue; k>NewHomology[2][j]; k--) {
                NewHomology[0].insert(NewHomology[0].begin()+j+1, -1);
                NewHomology[1].insert(NewHomology[1].begin()+j+1, -1);
                NewHomology[2].insert(NewHomology[2].begin()+j+1,  k);
                NewHomology[3].insert(NewHomology[3].begin()+j+1, -1);
                NewHomology[4].insert(NewHomology[4].begin()+j+1, -1);
                NewHomology[5].insert(NewHomology[5].begin()+j+1, -1);
                CurValue--;
            }
            CurValue--;
        }

        for (int k=CurValue; k>0; k--) {
            NewHomology[0].insert(NewHomology[0].begin(), -1);
            NewHomology[1].insert(NewHomology[1].begin(), -1);
            NewHomology[2].insert(NewHomology[2].begin(),  k);
            NewHomology[3].insert(NewHomology[3].begin(), -1);
            NewHomology[4].insert(NewHomology[4].begin(), -1);
            NewHomology[5].insert(NewHomology[5].begin(), -1);
            CurValue--;
        }

        //Node D Singletons
        //CurValue = GetFirstNonNegativeValueReverse(TripleInit->Get_AlignmentHomologies()[3], TripleInit->Get_AlignmentHomologies()[3].size()-1);
        CurValue = TripleInit->Get_Sequences()->GetBareSequenceAt(3).size();
        for (int j=NewHomology[3].size()-1; j>=0; j--) {
            if (NewHomology[3][j]==-1)       {
                continue;
            }
            if (NewHomology[3][j]==CurValue) {
                CurValue--;
                continue;
            }

            for (int k=CurValue; k>NewHomology[3][j]; k--) {
                NewHomology[0].insert(NewHomology[0].begin()+j+1, -1);
                NewHomology[1].insert(NewHomology[1].begin()+j+1, -1);
                NewHomology[2].insert(NewHomology[2].begin()+j+1, -1);
                NewHomology[3].insert(NewHomology[3].begin()+j+1,  k);
                NewHomology[4].insert(NewHomology[4].begin()+j+1, -1);
                NewHomology[5].insert(NewHomology[5].begin()+j+1, -1);
                CurValue--;
            }
            CurValue--;
        }

        for (int k=CurValue; k>0; k--) {
            NewHomology[0].insert(NewHomology[0].begin(), -1);
            NewHomology[1].insert(NewHomology[1].begin(), -1);
            NewHomology[2].insert(NewHomology[2].begin(), -1);
            NewHomology[3].insert(NewHomology[3].begin(),  k);
            NewHomology[4].insert(NewHomology[4].begin(), -1);
            NewHomology[5].insert(NewHomology[5].begin(), -1);
            CurValue--;
        }
        //}
    }
    break;
    case 0: { //Right Subtree
        //1. Convert the sampled Evolution History to the Corresponding Homology indices //{
        NewHomology.resize(6, IntVec(0,0));
        if (!HasExtraStartState) {
            int beforePos = TripleInit->GetLastIdxofSoanWhereRootIs(RegionStartIdxs[0]-1);

            for (int k=1; k<=beforePos; k++) {
                string tihl = TripleInit->GetTihlAtPosition(k, true);
                for (int j=0; j<tihl.size(); j++) {
                    int H_idx;
                    if      (j==0) {
                        H_idx = 5;
                    } else if (j==1) {
                        H_idx = 2;
                    } else if (j==2) {
                        H_idx = 3;
                    } else if (j==3) {
                        H_idx = 4;
                    }

                    if (tihl.at(j)=='E' || tihl.at(j)=='-') {
                        NewHomology[H_idx].push_back(-1);
                    } else {
                        NewHomology[H_idx].push_back(TripleInit->GetSoanIdxsAtPosition(k)[j]);
                    }
                }
            }
        }

        for (int i=0; i<NewEvolHistory.size(); i++) {
            for (int j=0; j<NewEvolHistory[i].size(); j++) {
                int H_idx;
                if      (i==0) {
                    H_idx = 5;
                } else if (i==1) {
                    H_idx = 2;
                } else if (i==2) {
                    H_idx = 3;
                } else if (i==3) {
                    H_idx = 4;
                }

                if (NewEvolHistory[i].at(j)=='E' || NewEvolHistory[i].at(j)=='-') {
                    NewHomology[H_idx].push_back(-1);
                } else {
                    NewHomology[H_idx].push_back(RegionStartIdxs[i]);
                    RegionStartIdxs[i]++;
                }
            }
        }

        if (!HasExtraEndState) {
            int afterPos = TripleInit->GetFirstIdxofSoanWhereRootIs(RegionStartIdxs[0]+1);

            for (int k=afterPos; k<TripleInit->Get_SoansIdxs().size()-1; k++) {
                string tihl = TripleInit->GetTihlAtPosition(k, true);
                for (int j=0; j<tihl.size(); j++) {
                    int H_idx;
                    if      (j==0) {
                        H_idx = 5;
                    } else if (j==1) {
                        H_idx = 2;
                    } else if (j==2) {
                        H_idx = 3;
                    } else if (j==3) {
                        H_idx = 4;
                    }

                    if (tihl.at(j)=='E' || tihl.at(j)=='-') {
                        NewHomology[H_idx].push_back(-1);
                    } else {
                        NewHomology[H_idx].push_back(TripleInit->GetSoanIdxsAtPosition(k)[j]);
                    }
                }
            }
        }

        //}

        //2. Form the "Skipped nodes" indices based on their parental node //{
        NewHomology[0].resize(NewHomology[4].size(), -1);
        NewHomology[1].resize(NewHomology[4].size(), -1);

        for (int j=0; j<NewHomology[4].size(); j++) {
            if (NewHomology[4][j]!=-1) {
                int pos = GetPosOfElement(TripleInit->Get_AlignmentHomologies()[4], NewHomology[4][j], 0);
                int val1 = TripleInit->Get_AlignmentHomologies()[0][pos];
                NewHomology[0][j] = val1;

                int val2 = TripleInit->Get_AlignmentHomologies()[1][pos];
                NewHomology[1][j] = val2;
            }
        }
        //}

        //3. Add the singleton values of the skipped nodes  //{
        //Node A Singletons
        //int CurValue = GetFirstNonNegativeValueReverse(TripleInit->Get_AlignmentHomologies()[0], TripleInit->Get_AlignmentHomologies()[0].size()-1);
        int CurValue = TripleInit->Get_Sequences()->GetBareSequenceAt(0).size();
        for (int j=NewHomology[0].size()-1; j>=0; j--) {
            if (NewHomology[0][j]==-1)       {
                continue;
            }
            if (NewHomology[0][j]==CurValue) {
                CurValue--;
                continue;
            }

            for (int k=CurValue; k>NewHomology[0][j]; k--) {
                NewHomology[0].insert(NewHomology[0].begin()+j+1,  k);
                NewHomology[1].insert(NewHomology[1].begin()+j+1, -1);
                NewHomology[2].insert(NewHomology[2].begin()+j+1, -1);
                NewHomology[3].insert(NewHomology[3].begin()+j+1, -1);
                NewHomology[4].insert(NewHomology[4].begin()+j+1, -1);
                NewHomology[5].insert(NewHomology[5].begin()+j+1, -1);
                CurValue--;
            }
            CurValue--;
        }

        for (int k=CurValue; k>0; k--) {
            NewHomology[0].insert(NewHomology[0].begin(),  k);
            NewHomology[1].insert(NewHomology[1].begin(), -1);
            NewHomology[2].insert(NewHomology[2].begin(), -1);
            NewHomology[3].insert(NewHomology[3].begin(), -1);
            NewHomology[4].insert(NewHomology[4].begin(), -1);
            NewHomology[5].insert(NewHomology[5].begin(), -1);
            CurValue--;
        }

        //Node B Singletons
        //CurValue = GetFirstNonNegativeValueReverse(TripleInit->Get_AlignmentHomologies()[1], TripleInit->Get_AlignmentHomologies()[1].size()-1);
        CurValue = TripleInit->Get_Sequences()->GetBareSequenceAt(1).size();
        for (int j=NewHomology[1].size()-1; j>=0; j--) {
            if (NewHomology[1][j]==-1)       {
                continue;
            }
            if (NewHomology[1][j]==CurValue) {
                CurValue--;
                continue;
            }

            for (int k=CurValue; k>NewHomology[1][j]; k--) {
                NewHomology[0].insert(NewHomology[0].begin()+j+1, -1);
                NewHomology[1].insert(NewHomology[1].begin()+j+1,  k);
                NewHomology[2].insert(NewHomology[2].begin()+j+1, -1);
                NewHomology[3].insert(NewHomology[3].begin()+j+1, -1);
                NewHomology[4].insert(NewHomology[4].begin()+j+1, -1);
                NewHomology[5].insert(NewHomology[5].begin()+j+1, -1);
                CurValue--;
            }
            CurValue--;
        }

        for (int k=CurValue; k>0; k--) {
            NewHomology[0].insert(NewHomology[0].begin(), -1);
            NewHomology[1].insert(NewHomology[1].begin(),  k);
            NewHomology[2].insert(NewHomology[2].begin(), -1);
            NewHomology[3].insert(NewHomology[3].begin(), -1);
            NewHomology[4].insert(NewHomology[4].begin(), -1);
            NewHomology[5].insert(NewHomology[5].begin(), -1);
            CurValue--;
        }
        //}
    }
    break;

    }
}

//Function to Estimate the new Homology Structure based on the re-sampled evolution history
void Menus::GetNewHomologyPair(int Branch, IntVecRef RegionStartIdxs, string RegionStartTihlLbl, bool HasExtraStartState, bool HasExtraEndState) {

    for (int i=0; i<RegionStartTihlLbl.size(); i++) {
        if (RegionStartTihlLbl.at(i)=='E' || RegionStartTihlLbl.at(i)=='-') {
            RegionStartIdxs[i]++;
        }
    }

    IntVec ResamplNodes; //The indices of the already Resampled Nodes

    switch(Branch) {
    case 1: {   //E-F (E is root)
        ResamplNodes.push_back(4);
        ResamplNodes.push_back(5);
        break;
    }
    case 2: {   //E-A (E is root)
        ResamplNodes.push_back(4);
        ResamplNodes.push_back(0);
        break;
    }
    case 3: {   //E-B (E is root)
        ResamplNodes.push_back(4);
        ResamplNodes.push_back(1);
        break;
    }
    case 4: {   //F-C (F is root)
        ResamplNodes.push_back(5);
        ResamplNodes.push_back(2);
        break;
    }
    case 5: {   //F-D (F is root)
        ResamplNodes.push_back(5);
        ResamplNodes.push_back(3);
        break;
    }
    case 6: {   //F-E (F is root)
        ResamplNodes.push_back(5);
        ResamplNodes.push_back(4);
        break;
    }
    }

    NewHomology.resize(6, IntVec(0,0));

    //1. Convert the sampled Evolution History to the Corresponding Homology indices //{
    if (!HasExtraStartState) {
        int beforePos = PairInit->GetLastIdxofSoanWhereRootIs(RegionStartIdxs[0]-1);

        for (int k=1; k<=beforePos; k++) {
            string tihl = PairInit->GetTihlAtPosition(k, true);
            for (int j=0; j<tihl.size(); j++) {
                int H_idx = ResamplNodes[j];
                if (tihl.at(j)=='E' || tihl.at(j)=='-') {
                    NewHomology[H_idx].push_back(-1);
                } else {
                    NewHomology[H_idx].push_back(PairInit->GetSoanIdxsAtPosition(k)[j]);
                }
            }
        }
    }

    for (int i=0; i<NewEvolHistory.size(); i++) {
        for (int j=0; j<NewEvolHistory[i].size(); j++) {
            int H_idx = ResamplNodes[i];
            if (NewEvolHistory[i].at(j)=='E' || NewEvolHistory[i].at(j)=='-') {
                NewHomology[H_idx].push_back(-1);
            } else {
                NewHomology[H_idx].push_back(RegionStartIdxs[i]);
                RegionStartIdxs[i]++;
            }
        }
    }

    if (!HasExtraEndState) {
        int afterPos = PairInit->GetFirstIdxofSoanWhereRootIs(RegionStartIdxs[0]+1);

        for (int k=afterPos; k<PairInit->Get_SoansIdxs().size()-1; k++) {
            string tihl = PairInit->GetTihlAtPosition(k, true);
            for (int j=0; j<tihl.size(); j++) {
                int H_idx = ResamplNodes[j];

                if (tihl.at(j)=='E' || tihl.at(j)=='-') {
                    NewHomology[H_idx].push_back(-1);
                } else {
                    NewHomology[H_idx].push_back(PairInit->GetSoanIdxsAtPosition(k)[j]);
                }
            }
        }
    }
    //}

    //2. Form the "Skipped nodes" indices based on their parental node //{
    //Case E-A or E-B (with E being the Root) -> Fix F Node //{
    if (ResamplNodes[0]==4 && ResamplNodes[1]!=5) {
        NewHomology[5].resize(NewHomology[4].size(), -1);
        for (int j=0; j<NewHomology[4].size(); j++) {
            if (NewHomology[4][j]!=-1) {
                int pos = GetPosOfElement(PairInit->Get_AlignmentHomologies()[4], NewHomology[4][j], 0);
                int val1 = PairInit->Get_AlignmentHomologies()[5][pos];
                NewHomology[5][j] = val1;
            }
        }

        //Now Add Node F Singletons (if any)
        int CurValue = GetFirstNonNegativeValueReverse(PairInit->Get_AlignmentHomologies()[5], PairInit->Get_AlignmentHomologies()[5].size()-1);
        for (int j=NewHomology[5].size()-1; j>=0; j--) {
            if (NewHomology[5][j]==-1)       {
                continue;
            }
            if (NewHomology[5][j]==CurValue) {
                CurValue--;
                continue;
            }

            for (int k=CurValue; k>NewHomology[5][j]; k--) {
                for (int l=0; l<NewHomology.size(); l++) {
                    if (NewHomology[l].size()>0 && l!=5) {
                        NewHomology[l].insert(NewHomology[l].begin()+j+1, -1);
                    } else if (NewHomology[l].size()>0 && l==5) {
                        NewHomology[l].insert(NewHomology[l].begin()+j+1,  k);
                    }
                }
                CurValue--;
            }
            CurValue--;
        }

        for (int k=CurValue; k>0; k--) {
            for (int l=0; l<NewHomology.size(); l++) {
                if (NewHomology[l].size()>0 && l!=5) {
                    NewHomology[l].insert(NewHomology[l].begin(), -1);
                } else if (NewHomology[l].size()>0 && l==5) {
                    NewHomology[l].insert(NewHomology[l].begin(),  k);
                }
            }
            CurValue--;
        }
        //}
    }

    //Case F-C or F-D (with F being the Root) -> Fix E Node //{
    if (ResamplNodes[0]==5 && ResamplNodes[1]!=4) {
        NewHomology[4].resize(NewHomology[5].size(), -1);
        for (int j=0; j<NewHomology[5].size(); j++) {
            if (NewHomology[5][j]!=-1) {
                int pos = GetPosOfElement(PairInit->Get_AlignmentHomologies()[5], NewHomology[5][j], 0);
                int val1 = PairInit->Get_AlignmentHomologies()[4][pos];
                NewHomology[4][j] = val1;
            }
        }

        //Now Add Node E Singletons (if any)
        int CurValue = GetFirstNonNegativeValueReverse(PairInit->Get_AlignmentHomologies()[4], PairInit->Get_AlignmentHomologies()[4].size()-1);
        for (int j=NewHomology[4].size()-1; j>=0; j--) {
            if (NewHomology[4][j]==-1)       {
                continue;
            }
            if (NewHomology[4][j]==CurValue) {
                CurValue--;
                continue;
            }

            for (int k=CurValue; k>NewHomology[4][j]; k--) {
                for (int l=0; l<NewHomology.size(); l++) {
                    if (NewHomology[l].size()>0 && l!=4) {
                        NewHomology[l].insert(NewHomology[l].begin()+j+1, -1);
                    } else if (NewHomology[l].size()>0 && l==4) {
                        NewHomology[l].insert(NewHomology[l].begin()+j+1,  k);
                    }
                }
                CurValue--;
            }
            CurValue--;
        }

        for (int k=CurValue; k>0; k--) {
            for (int l=0; l<NewHomology.size(); l++) {
                if (NewHomology[l].size()>0 && l!=4) {
                    NewHomology[l].insert(NewHomology[l].begin(), -1);
                } else if (NewHomology[l].size()>0 && l==4) {
                    NewHomology[l].insert(NewHomology[l].begin(),  k);
                }
            }
            CurValue--;
        }
    }
    //}

    //Now that Nodes E and F are formed, continue with nodes A,B, C and D //{
    for (int i=0; i<4; i++) { //Nodes A and B
        if (i==ResamplNodes[0] || i== ResamplNodes[1]) {
            continue;    //skip already sampled nodes
        }
        int BaseNode;
        if (i==0 || i==1) {
            BaseNode = 4;
        } else              {
            BaseNode = 5;
        }

        NewHomology[i].resize(NewHomology[BaseNode].size(), -1);
        for (int j=0; j<NewHomology[BaseNode].size(); j++) {
            if (NewHomology[BaseNode][j]!=-1) {
                int pos = GetPosOfElement(PairInit->Get_AlignmentHomologies()[BaseNode], NewHomology[BaseNode][j], 0);
                int val1 = PairInit->Get_AlignmentHomologies()[i][pos];
                NewHomology[i][j] = val1;
            }
        }
    }
    //}

    //Finally, Fix the Node A, B, C, D, Singletons
    for (int i=0; i<4; i++) { //Nodes A and B
        if (i==ResamplNodes[0] || i== ResamplNodes[1]) {
            continue;    //skip already sampled nodes
        }
        int CurValue = GetFirstNonNegativeValueReverse(PairInit->Get_AlignmentHomologies()[i], PairInit->Get_AlignmentHomologies()[i].size()-1);
        for (int j=NewHomology[i].size()-1; j>=0; j--) {
            if (NewHomology[i][j]==-1)       {
                continue;
            }
            if (NewHomology[i][j]==CurValue) {
                CurValue--;
                continue;
            }

            for (int k=CurValue; k>NewHomology[i][j]; k--) {
                for (int l=0; l<NewHomology.size(); l++) {
                    if (l!=i) {
                        NewHomology[l].insert(NewHomology[l].begin()+j+1, -1);
                    } else {
                        NewHomology[l].insert(NewHomology[l].begin()+j+1,  k);
                    }
                }
                CurValue--;
            }
            CurValue--;
        }

        for (int k=CurValue; k>0; k--) {
            for (int l=0; l<NewHomology.size(); l++) {
                if (l!=i) {
                    NewHomology[l].insert(NewHomology[l].begin(), -1);
                } else {
                    NewHomology[l].insert(NewHomology[l].begin(),  k);
                }
            }
            CurValue--;
        }

    }

    //}
}

//Function to Estimate the new Homology Structure based on the newly estimated Homology Structure
void Menus::GetNewAlignment() {
    for (int i=0; i<NewHomology.size()-2; i++) {
        stringstream ss;
        ss.str("");
        for (int j=0; j<NewHomology[i].size(); j++) {
            if (NewHomology[i][j]==-1) {
                ss<<"-";
            } else {
                ss<<TripleInit->Get_Sequences()->GetBareSequenceAt(i).at(NewHomology[i][j]-1);
            }
        }
        NewAlignment.push_back(ss.str());
    }
}

//Function to Print the Newly sampled Evolution History to the Standard Output Stream
void Menus::PrintNewEvolHistory() {
    cout<<"\nNew Evolution History:\n----------------------\n"<<NewEvolHistory<<endl;
}

//Function to Print the Newly estimated Homology Structure to the Standard Output Stream
void Menus::PrintNewHomology() {
    cout<<"\nNew Homology Structure:\n-----------------------\n";
    for (int i=0; i<NewHomology.size(); i++) {
        for (int j=0; j<NewHomology[i].size(); j++) {
            cout<<right<<setw(5)<<NewHomology[i][j];
        }
        cout<<endl;
    }
    cout<<endl;
}

//Function to Load a starting Homology from a file
void Menus::LoadNewHomology() {
    string filename = StartingHomologiesFile;
    string line;
    ifstream myfile (filename.c_str());
    if (myfile.is_open()) {
        while ( myfile.good() ) {
            getline (myfile,line);

            stringstream stream(line);
            string entry;
            IntVec Hom_line;
            while(getline(stream, entry, '\t') ){
                Hom_line.push_back( atoi(entry.c_str()) );
            }
            NewHomology.push_back(Hom_line);
        }
        myfile.close();
    } else cerr << "Unable to open file: "<<filename<<" for reading.";
}

//Function to Save the newly Estismated Homology to a file
void Menus::SaveNewHomology() {
    stringstream ss; ss.str("");

    for (int i=0; i<NewHomology.size(); i++) {
        for (int j=0; j<NewHomology[i].size(); j++) {
            ss<<NewHomology[i][j];
            if (j<NewHomology[i].size()-1){
                ss<<"\t";
            }
        }
        ss<<endl;
    }
    WriteToFile(StartingHomologiesFile, ss.str(), false, true);
}

//Function to Print the Newly estimated Alignment to the Standard Output Stream
void Menus::PrintNewAlignment() {
    cout<<"\nNew Alignment:\n--------------\n"<<NewAlignment<<endl;
}

//Function to Print the Created Intervals
void Menus::PrintIntervals(IntervalVecRef itvs) {
    for (int i=0; i<itvs.size(); i++) {
        cout<<"Interval "<<(i+1)<<": ["<<right<<setw(3)<<itvs[i].start<<" - "<<setw(3)<<itvs[i].end<<"]"<<endl;
    }
}

//Function to Print the specification for a given alignment region
void Menus::PrintAlignmentRegion(IntVecRef StartIdxs, IntVecRef EndIdxs, string StartTihlLbl, string StartSoanLbl, bool HasExtraStartState, bool HasExtraEndState, IntVecRef Dimensions, DoubleVecRef Lengths) {
    cout<<"\n StartIdxs: <";
    for (int i=0; i<StartIdxs.size(); i++) {
        cout<<setw(5)<<StartIdxs[i];
    }
    cout<<" >";
    cout<<"    StartTihlLbl: "<<StartTihlLbl<<"    StartSoanLbl: "<<StartSoanLbl<<endl;

    cout<<"   EndIdxs: <";
    for (int i=0; i<EndIdxs.size(); i++) {
        cout<<setw(5)<<EndIdxs[i];
    }
    cout<<" >"<<endl;

    cout<<"Dimensions: <";
    for (int i=0; i<Dimensions.size(); i++) {
        cout<<setw(5)<<Dimensions[i];
    }
    cout<<" >"<<endl;

    cout<<"   Lengths: <";
    for (int i=0; i<Lengths.size(); i++) {
        cout<<setw(5)<<Lengths[i];
    }
    cout<<" >"<<endl;

    if (HasExtraStartState) {
        cout<<"\t  HasExtraStartState: True \t";
    } else                    {
        cout<<"\t  HasExtraStartState: False\t";
    }

    if (HasExtraEndState)   {
        cout<<"  HasExtraEndState: True"<<endl;
    } else                    {
        cout<<"  HasExtraEndState: False"<<endl;
    }
}

//Function to save the Newly estimated Alignment back to the file
void Menus::UpdateNewAlignment() {
    stringstream ss;
    ss.str("");
    for (int i=0; i<NewAlignment.size(); i++) {
        ss<<TripleInit->Get_Sequences()->GetSequenceIdAt(i)<<endl;
        ss<<NewAlignment[i]<<endl<<endl;
    }

    WriteToFile(AlignedSequencesFile, ss.str(), false, true);
}

//Function to Find an "All Homologous" Column within the Overall Homologies Table starting from the given startPos and going left (if Derection==1) or right (if Direction==2).
int Menus::GetIdxOfAllHomologousColumn(IntVecRef SubTreeIdxs, IntMatRef OverallHomologies, int StartPos, int Direction) {
    int Hpos = -1;
    switch (Direction) {
    case 1: //Go to the left of StartPos until an Homology is found
        for (int j=StartPos; j>=0; j--) {
            bool HFound = true;
            for (int i=0; i<SubTreeIdxs.size(); i++) {
                if ( OverallHomologies[SubTreeIdxs[i]][j]==-1 ) {
                    HFound = false;
                    break;
                }
            }
            if (HFound) {
                Hpos = j;
                break;
            }
        }
        break;

    case 2: //Go to the right of StartPos until an Homology is found
        for (int j=StartPos; j<OverallHomologies[0].size(); j++) {
            bool HFound = true;
            for (int i=0; i<SubTreeIdxs.size(); i++) {
                if ( OverallHomologies[SubTreeIdxs[i]][j]==-1 ) {
                    HFound = false;
                    break;
                }
            }
            if (HFound) {
                Hpos = j;
                break;
            }
        }
        break;
    }

    //if      (Hpos==-1 && Direction==1) { Hpos = 0; }
    //else if (Hpos==-1 && Direction==2) { Hpos = Seqs->Get_IndicesFromAlignment()[0].size()-1; }

    return Hpos;
}

//Function to read the given file and store the alignment (ignoring the sequence headers) in the given vector of strings
void Menus::GetAlignmentFromFile(StringVecRef algn, string filename) {
    /** Read Alignment **/
    string line;
    ifstream myfile (filename.c_str());

    stringstream seqText;
    bool firstLine = true;

    if (myfile.is_open()) {
        while ( myfile.good() ) {
            getline (myfile,line);

            if (line.compare("")==0 || line.compare("\n")==0 || line.compare(" ")==0) continue; //skip empty lines

            if (StartsWith(line, ">") || StartsWith(line,";")) {
                if (!firstLine) {
                    algn.push_back(seqText.str());
                    seqText.str("");
                }
                firstLine = false;
            } else {
                seqText << line;
            }

        }
        algn.push_back(seqText.str());
    } else cerr << " Unable to open file: "<<filename<<" for reading.";
    myfile.close();

}
