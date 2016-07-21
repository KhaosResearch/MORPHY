#ifndef TYPEDEFINITIONS_H
#define TYPEDEFINITIONS_H

#include <iostream>
#include <vector>

using namespace std;

namespace TypeDefinitions {
    typedef vector<string>   StringVec;
    typedef vector<string>&  StringVecRef;
    typedef vector<string>*  StringVecPtr;

    typedef const vector<string>   StringCVec;
    typedef const vector<string>&  StringCVecRef;
    typedef const vector<string>*  StringCVecPtr;

    typedef vector<double>   DoubleVec;
    typedef vector<double>&  DoubleVecRef;
    typedef vector<double>*  DoubleVecPtr;

    typedef const vector<double>   DoubleCVec;
    typedef const vector<double>&  DoubleCVecRef;
    typedef const vector<double>*  DoubleCVecPtr;

    typedef vector<int>      IntVec;
    typedef vector<int>&     IntVecRef;
    typedef vector<int>*     IntVecPtr;

    typedef const vector<int>      IntCVec;
    typedef const vector<int>&     IntCVecRef;
    typedef const vector<int>*     IntCVecPtr;

    typedef vector <vector< vector<int > > >   IntHyp;
    typedef vector <vector< vector<int > > >&  IntHypRef;
    typedef vector <vector< vector<int > > >*  IntHypPtr;

    typedef const vector <vector< vector<int > > >   IntCHyp;
    typedef const vector <vector< vector<int > > >&  IntCHypRef;
    typedef const vector <vector< vector<int > > >*  IntCHypPtr;

    typedef vector< vector<string> >  StringMat;
    typedef vector< vector<string> >& StringMatRef;
    typedef vector< vector<string> >* StringMatPtr;

    typedef const vector< vector<string> >  StringCMat;
    typedef const vector< vector<string> >& StringCMatRef;
    typedef const vector< vector<string> >* StringCMatPtr;

    typedef vector< vector<double> >  DoubleMat;
    typedef vector< vector<double> >& DoubleMatRef;
    typedef vector< vector<double> >* DoubleMatPtr;

    typedef const vector< vector<double> >  DoubleCMat;
    typedef const vector< vector<double> >& DoubleCMatRef;
    typedef const vector< vector<double> >* DoubleCMatPtr;

    typedef vector <vector< vector<double > > >   DoubleHyp;
    typedef vector <vector< vector<double > > >&  DoubleHypRef;
    typedef vector <vector< vector<double > > >*  DoubleHypPtr;

    typedef const vector <vector< vector<double > > >   DoubleCHyp;
    typedef const vector <vector< vector<double > > >&  DoubleCHypRef;
    typedef const vector <vector< vector<double > > >*  DoubleCHypPtr;

    typedef vector< vector<int> >     IntMat;
    typedef vector< vector<int> >&    IntMatRef;
    typedef vector< vector<int> >*    IntMatPtr;

    typedef const vector< vector<int> >     IntCMat;
    typedef const vector< vector<int> >&    IntCMatRef;
    typedef const vector< vector<int> >*    IntCMatPtr;

    typedef long long int LongInt;

    typedef struct {
        string StartSoanLbl;
        IntVec StartSoanIdxs;
        string StartTihlLbl;
        string EndSoanLbl;
        IntVec EndSoanIdxs;
        bool   HasExtraStartState;
        bool   HasExtraEndState;
    } SoanBoundaries;

    typedef SoanBoundaries& SoanBoundariesRef;

    typedef struct {
        int start;
        int end;
    } Interval;

    typedef Interval& IntervalRef;
    typedef vector<Interval>  IntervalVec;
    typedef vector<Interval>& IntervalVecRef;

    const string ParametersFile            = "Parameter Files/Parameters.txt";
    const string BareSequencesFile         = "Sequences Files/Sequences.fasta";
    const string AlignedSequencesFile      = "Sequences Files/Alignment.fasta";
    const string ReferenceSequencesFile    = "Sequences Files/Reference.fasta";
    const string ExtraAlignSequencesFile   = "Sequences Files/ExtraAlign.fasta";
    const string TreeFile                  = "Sequences Files/Tree.nck";
    const string ExportTreeFile            = "Output Files/Tree.txt";
    const string LogHomologsFile           = "Output Files/LogH.txt";
    const string ProbabilitiesHomologsFile = "Output Files/ProbsH.txt";
    const string LogFile                   = "Output Files/LogN.txt";
    const string ProbabilitiesFile         = "Output Files/ProbsN.txt";
    const string AlignmentHomologiesFile   = "Output Files/AlignmentHomologies.txt";
    const string RegionHomologiesFile      = "Output Files/RegionHomologies.txt";
    const string FPTablesFile              = "Output Files/FPTables.txt";
    const string HomologiesFile            = "Output Files/Homologies.txt";
    const string StartingHomologiesFile    = "Parameter Files/StartingHomologies.txt";
    const string ReachableSoans2LeavedFile = "Parameter Files/ReachableSoans2Leaved.txt";
    const string ReachableSoans3LeavedFile = "Parameter Files/ReachableSoans3Leaved.txt";
    const string ReachableSoans4LeavedFile = "Parameter Files/ReachableSoans4Leaved.txt";
    const string AlignerStatusFile         = "Parameter Files/Status.txt";
};

#endif // TYPEDEFINITIONS_H
