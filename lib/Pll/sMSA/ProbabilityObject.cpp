#include "ProbabilityObject.h"

//Default (empty Constructor)
ProbabilityObject::ProbabilityObject()  { OverallProbability = 0.0; }

//Virtual Decontructor
ProbabilityObject::~ProbabilityObject() {
    //for (int i=0;i<Contributors.size();i++){
    //    delete Contributors[i];
    // }
    vector<int>().swap(PosIn_ND);
    vector<string>().swap(ContributorTihls);
    vector<string>().swap(ContributorTihlsExt);
    vector<double>().swap(ContributorProbabilities);
    vector<double>().swap(ContibutorDirectTransition);
    vector<double>().swap(ContibutorSelfTransition);
    vector<bool>().swap(ContributorWingFolds);
}

//Definition of The equality (==) operator
bool ProbabilityObject::operator==(const ProbabilityObject &other){
    stringstream s1, s2; s1.str(""); s2.str(""); s1 << *(this); s2 << (other);
    if (s1.str().compare(s2.str())!= 0) return false;
    return true;
}

//Definition of The inequality (!=) operator
bool ProbabilityObject::operator!=(const ProbabilityObject &other){
    return !(*(this)==other);
}

//Overloading of the comparison operators
bool ProbabilityObject::operator<   (const ProbabilityObject& p) const {
    return PosIn_1D < p.Get_PosIn_1D();
}

//Overloading of the comparison operators
bool ProbabilityObject::operator<=  (const ProbabilityObject& p) const {
    return PosIn_1D <= p.Get_PosIn_1D();
}

//Overloading of the comparison operators
bool ProbabilityObject::operator>   (const ProbabilityObject& p) const {
    return PosIn_1D > p.Get_PosIn_1D();
}

//Overloading of the comparison operators
bool ProbabilityObject::operator>=  (const ProbabilityObject& p) const {
    return PosIn_1D >= p.Get_PosIn_1D();
}

//Function that adds a new contributor to the collection (vector) of contributos
void ProbabilityObject::AddContributor(ProbabilityObject* contr) {
    Contributors.push_back(contr);
}

//Function that returns the contributor whose position equals to the given -as argument- position
ProbabilityObject* ProbabilityObject::GetContributorAtPos(int pos) const{
    return Contributors[pos];
}

//Function that adds a new string corresponding to a contributor
void ProbabilityObject::AddContributorTihl(string tihl) {
    ContributorTihls.push_back(tihl);
}

//Function that returns the tihl of the contributor whose position equals to the given -as argument- position
string ProbabilityObject::GetContributorTihlAtPos(int pos) const {
    return ContributorTihls[pos];
}

//Function that adds a new string corresponding to a contributor
void ProbabilityObject::AddContributorTihlExt(string tihl) {
    ContributorTihlsExt.push_back(tihl);
}

//Function that returns the tihl of the contributor whose position equals to the given -as argument- position
string ProbabilityObject::GetContributorTihlExtAtPos(int pos) const {
    return ContributorTihlsExt[pos];
}

//Function that adds a new probability value corresponding to a contributor
void ProbabilityObject::AddContributorProbability(double val) {
    ContributorProbabilities.push_back(val);
}

//Function that returns the contributed probability of the contributor whose position equals to the given -as argument- position
double ProbabilityObject::GetContributorProbabilityAtPos(int pos) const {
    return ContributorProbabilities[pos];
}

//Function to increase the OverallProbability by the given value
void ProbabilityObject::IncreaseOverallProbabilityBy(double val) {
    OverallProbability = AddLog(OverallProbability, val);
}

//Function that adds true or false according to whether the given contributor is wingfolded or not
void ProbabilityObject::AddContributorWingFold(bool val){
    ContributorWingFolds.push_back(val);
}

//Function that adds the direct transition probability corresponding to a contributor
void ProbabilityObject::AddContributorDirectTransProb(double val) {
    ContibutorDirectTransition.push_back(val);
}

//Function that adds the indirect transition probability corresponding to a contributor
void ProbabilityObject::AddContributorInDirectTransProb(double val) {
    ContibutorInDirectTransition.push_back(val);
}

//Function that returns the direct transition probability corresponding to the given contributor
double ProbabilityObject::GetContributorDirectTransProbAtPos(int pos) const {
    return ContibutorDirectTransition[pos];
}

//Function that returns the indirect transition probability corresponding to the given contributor
double ProbabilityObject::GetContributorInDirectTransProbAtPos(int pos) const {
    return ContibutorInDirectTransition[pos];
}

//Function that adds the self transition probability (to the non emitting state) corresponding to a contributor
void ProbabilityObject::AddContributorSelfTransProb(double val) {
    ContibutorSelfTransition.push_back(val);
}

//Function that returns the self transition probability (to the non emitting state) corresponding to the given contributor
double ProbabilityObject::GetContributorSelfTransProbAtPos(int pos) const {
    return ContibutorSelfTransition[pos];
}

//Function that returns whether the contributor at the given position is wing folded or not
 bool ProbabilityObject::GetContributorWingFoldAtPos(int pos) const {
    return ContributorWingFolds[pos];
 }

//Function to return a brief description (i.e. without accoubting for its contributors) of current ProbabilityObject
string ProbabilityObject::GetBriefDescription() const {
    stringstream ss; ss.str("");
    ss << "State: " << this->Get_SoanLbl() << " ";
    ss << "[" << setw(3) << this->Get_PosIn_1D() << "] <=> [";

    for (unsigned int i=0;i<this->Get_PosIn_ND().size();i++){
        if (i<this->Get_PosIn_ND().size()-1){
            ss << setw(2) << this->Get_PosIn_ND()[i] << ", ";
        }else{
            ss << setw(2) << this->Get_PosIn_ND()[i];
        }
    }

    ss << "] ";

    if (this->Is_EndState()) { ss << " (End State) "; }
    ss << "  " << Contributors.size() << " Contributors. " << endl;

    return ss.str();
}

//Definition of the toString Operation
ostream& operator<<(ostream& ostr, const ProbabilityObject& prob) {
    ostr << "State: " << prob.Get_SoanLbl() << " ";
    ostr << "[" << setw(3) << prob.Get_PosIn_1D() << "] <=> [";

    for (unsigned int i=0;i<prob.Get_PosIn_ND().size();i++){
        if (i<prob.Get_PosIn_ND().size()-1){
            ostr << setw(2) << prob.Get_PosIn_ND()[i] << ", ";
        }else{
            ostr << setw(2) << prob.Get_PosIn_ND()[i];
        }
    }

    ostr << "] ";

    if (prob.Is_EndState()) { ostr << " (End State) "; }

    if (prob.Get_Contributors().size()>0) {
        ostr << "Prob: " << setw(9) << prob.Get_OverallProbability() << " Formed By: \n";

        for (unsigned int i=0;i<prob.Get_Contributors().size();i++){
            ostr << "\t" << "State: " << prob.GetContributorAtPos(i)->Get_SoanLbl() << " ";
            ostr << "[" << setw(3) << prob.GetContributorAtPos(i)->Get_PosIn_1D() << "] <=> [";

            for (unsigned int j=0;j<prob.GetContributorAtPos(i)->Get_PosIn_ND().size();j++){
                if (j<prob.GetContributorAtPos(i)->Get_PosIn_ND().size()-1){
                    ostr << setw(2) << prob.GetContributorAtPos(i)->Get_PosIn_ND()[j] << ", ";
                }else{
                    ostr << setw(2) << prob.GetContributorAtPos(i)->Get_PosIn_ND()[j];
                }
            }
            ostr << "] ";

            ostr << " " << " -> (" << prob.GetContributorTihlExtAtPos(i) <<")";
            ostr << " " << " Pr.: " << prob.GetContributorProbabilityAtPos(i) <<"";
            if (prob.GetContributorAtPos(i)->Is_EndState()) { ostr << " (End State) "; }
            if (prob.GetContributorWingFoldAtPos(i)) { ostr << " WF " << " DT:" << prob.GetContributorDirectTransProbAtPos(i) <<" ST:" <<prob.GetContributorSelfTransProbAtPos(i); }
            ostr << endl;
        }
    }else {
        ostr << "Prob: " << setw(9) << prob.Get_OverallProbability() <<"  (No Contributors)" << endl;
    }

    return ostr;
}
