#include "Parameters.h"

//Default (Empty) Constructor
Parameters::Parameters()  { }

//The Default Deconstructor
Parameters::~Parameters() { }

//Definition of The equality (==) operator
bool Parameters::operator==(const Parameters &params) {
    for (int i=0;i<4;i++){
        if (BaseFreq[i]!=params.Get_BaseFreq()[i]) return false;
    }
    if (Lambda!=params.Get_Lambda() || Gamma!=params.Get_Gamma() || SampleSize!= params.Get_SampleSize() || AlignFragSize!=params.Get_AlignFragSize() )
        return false;

    if (AlignmentOutput.compare(params.Get_AlignmentOutput())!=0) return false;

    if (SubModel.compare(params.Get_SubModel())!=0) return false;


    if (k!=params.Get_Kappa() || k1!=params.Get_AG_Rate() || k2!=params.Get_CT_Rate() || x1!=params.Get_X1_Param() ||
        x2!=params.Get_X2_Param() || x3!=params.Get_X3_Param() || x4!=params.Get_X4_Param() ||
        x5!=params.Get_X5_Param() || x6!=params.Get_X6_Param()
        ) return false;


    return true;

}

//Definition of The inequality (!=) operator
bool Parameters::operator!= (const Parameters& params) {
    return !(*(this)==params);
}

//Copy Constructor to Clone objects
Parameters::Parameters(const Parameters& other) {
    for (int i=0;i<4;i++){
        this->BaseFreq[i] = other.BaseFreq[i];
    }

    this->Lambda         = other.Get_Lambda();
    this->Gamma          = other.Get_Gamma();
    this->SampleSize     = other.Get_SampleSize();
    this->AlignFragSize  = other.Get_AlignFragSize();
    this->AlignmentOutput= other.Get_AlignmentOutput();
    this->SubModel       = other.Get_SubModel();
    this->k              = other.Get_Kappa();
    this->k1             = other.Get_AG_Rate();
    this->k1             = other.Get_AG_Rate();
    this->x1             = other.Get_X1_Param();
    this->x2             = other.Get_X6_Param();
    this->x3             = other.Get_X1_Param();
    this->x4             = other.Get_X1_Param();
    this->x5             = other.Get_X1_Param();
    this->x6             = other.Get_X1_Param();
}

//Function to Check whether the Parameters have been initiliazed
bool Parameters::HaveBeenInitialized() {
    double sum = 0;
    for (int i=0;i<4;i++){
        sum += BaseFreq[i];
    }
    if (sum!=1 || Lambda<0 || Gamma<1) return false;

    if (AlignmentOutput.compare("MP")!=0 && AlignmentOutput.compare("ML")!=0 && AlignmentOutput.compare("MS")!=0) return false;
    if (SubModel.compare("F84")!=0 && SubModel.compare("TN93")!=0 && SubModel.compare("GTR")!=0) return false;
    if (SampleSize<=0 || AlignFragSize<=0) return false;
    if (SubModel.compare("F84")==0) {
        if (k<0) return false;
    }else if (SubModel.compare("TN93")==0) {
        if (k1<0 || k2<0) return false;
    }else if (SubModel.compare("GTR")==0) {
        if (x1<0 || x2<0 || x3<0 || x4<0 || x5<0 || x6<0) return false;
    }

    return true;
}

//Function to Read Parameters from the CMD
void Parameters::ReadDataFromCmd() {
    cout << endl;
    cout << "Enter A frequency, p(A)  : ";
    cin >> BaseFreq[A];
    cout << "Enter C frequency, p(C)  : ";
    cin >> BaseFreq[C];
    cout << "Enter G frequency, p(G)  : ";
    cin >> BaseFreq[G];
    cout << "Enter T frequency, p(T)  : ";
    cin >> BaseFreq[T];
    cout << "Enter insertion rate (Lambda)                    : ";
    cin >> Lambda;
    cout << "Enter expected length of fragments (Gamma)       : ";
    cin >> Gamma;
    cout << "Enter desired length of simulated sequence       : ";
    cin >> SampleSize;
    cout << "Enter desired length of the alignment fragments  : ";
    cin >> AlignFragSize;
    cout << "Enter desired length of output type (ML, MS, MP) : ";
    cin >> AlignmentOutput;
    cout << "Enter desired Substitution Model(F84, TN93, GTR) : ";
    cin >> SubModel;

    if (SubModel.compare("F84")==0) {
        cout << "\tEnter (k) Transitions-Transversions Rate         : ";
        cin >> k;
        k1=k2=x1=x2=x3=x4=x5=x6=0;
    }else if (SubModel.compare("TN93")==0) {
        cout << "\tEnter the A<->G rate                             : ";
        cin >> k1;
        cout << "\tEnter the C<->T rate                             : ";
        cin >> k2;
        k=x1=x2=x3=x4=x5=x6=0;
    }else if (SubModel.compare("GTR")==0) {
        cout << "\tEnter x1 substitution rate parameter             : ";
        cin >> x1;
        cout << "\tEnter x2 substitution rate parameter             : ";
        cin >> x2;
        cout << "\tEnter x3 substitution rate parameter             : ";
        cin >> x3;
        cout << "\tEnter x4 substitution rate parameter             : ";
        cin >> x4;
        cout << "\tEnter x5 substitution rate parameter             : ";
        cin >> x5;
        cout << "\tEnter x6 substitution rate parameter             : ";
        cin >> x6;
        k=k1=k2=0;
    }

}

//Function to Read Parameters from an input file
void Parameters::ReadDataFromFile(string filename) {
    ifstream inputstr(filename.c_str());
    if (!inputstr) {
        cerr << "cannot open file " << filename.c_str() << endl;
        exit(1);
    }
    inputstr >> BaseFreq[A] >> BaseFreq[C] >> BaseFreq[G]>> BaseFreq[T] >> Lambda >> Gamma >> SampleSize >> AlignFragSize >> AlignmentOutput >> SubModel ;

    if (SubModel.compare("F84")==0) {
        inputstr >> k;
        k1=k2=x1=x2=x3=x4=x5=x6=0;
    }else if (SubModel.compare("TN93")==0) {
        inputstr >> k1 >> k2;
        k=x1=x2=x3=x4=x5=x6=0;
    }else if (SubModel.compare("GTR")==0) {
        inputstr >> x1 >> x2 >> x3 >> x4 >> x5 >> x6;
        k=k1=k2=0;
    }

    inputstr.close();
}

//Definition of the toString Operation
ostream& operator<<(ostream& ostr, const Parameters& p) {
    const double* freq = p.Get_BaseFreq();
    ostr << endl;
    ostr << "Initial Nucleotide Probabilities:\n";
    ostr << "\tp(A)  : ";
    ostr << freq[A] << endl;
    ostr << "\tp(C)  : ";
    ostr << freq[C] << endl;
    ostr << "\tp(G)  : ";
    ostr << freq[G] << endl;
    ostr << "\tp(T)  : ";
    ostr << freq[T] << endl;
    ostr << "Insertion rate (Lambda)              : ";
    ostr << p.Get_Lambda() << endl;
    ostr << "Expected length of fragments (Gamma) : ";
    ostr << p.Get_Gamma() << endl;
    ostr << "Desired size for sampling sequences  : ";
    ostr << p.Get_SampleSize() << endl;
    ostr << "Length of the alignment fragments    : ";
    ostr << p.Get_AlignFragSize() << endl;
    ostr << "Desired alignment Output Type        : ";
    if (p.Get_AlignmentOutput().compare("ML")==0){
        ostr << "Most Likely - Forward (ML)" << endl;
    }else if (p.Get_AlignmentOutput().compare("MP")==0){
        ostr << "Most Probable -Viterbi (MP)" << endl;
    }else if (p.Get_AlignmentOutput().compare("MS")==0){
        ostr << "Most Sampled - Forward (MS)" << endl;
    }

    ostr << "Substitution Model                   : ";
    ostr << p.Get_SubModel() << endl;

    if (p.Get_SubModel().compare("F84")==0) {
        ostr << "\t(k) Transitions-Transversions Rate: ";
        ostr << p.Get_Kappa() << endl;
    }else if (p.Get_SubModel().compare("TN93")==0) {
        ostr << "\t(k1) A<->G rate: ";
        ostr << p.Get_AG_Rate() << endl;
        ostr << "\t(k2) C<->T rate: ";
        ostr << p.Get_CT_Rate() << endl;
    }else if (p.Get_SubModel().compare("GTR")==0) {
        ostr << "\tx1 substitution rate parameter: ";
        ostr << p.Get_X1_Param() << endl;
        ostr << "\tx2 substitution rate parameter: ";
        ostr << p.Get_X2_Param() << endl;
        ostr << "\tx3 substitution rate parameter: ";
        ostr << p.Get_X3_Param() << endl;
        ostr << "\tx4 substitution rate parameter: ";
        ostr << p.Get_X4_Param() << endl;
        ostr << "\tx5 substitution rate parameter: ";
        ostr << p.Get_X5_Param() << endl;
        ostr << "\tx6 substitution rate parameter: ";
        ostr << p.Get_X6_Param() << endl;
    }


    return ostr;
}
