#include "FIDModel.h"

//Default (Empty) Constructor
FIDModel::FIDModel()  { }

//Constructor that builds a FID object based on the given Lambda and Gamma parameters
FIDModel::FIDModel(uint GammaVal, double LambdaVal) {
    if (GammaVal < 1 || LambdaVal < 0) {
        cout<<"Gamma value should be >= 1 and Lambda value>=0!"<<endl;
        exit(1);
    }

    Gamma  = GammaVal;
    Lambda = LambdaVal;
}

//The Default Deconstructor
FIDModel::~FIDModel() {
}

//Definition of The equality (==) operator
bool FIDModel::operator==(const FIDModel &fid) {
    for (int i=0;i<3;i++){
        for (int j=0;j<3;j++){
            if (this->ProbabilitiesMatrix[i][j] != fid.ProbabilitiesMatrix[i][j]) return false;
        }
    }

    if (this->Lambda!=fid.Get_Lambda() || this->Gamma!=fid.Get_Gamma())
        return false;

    return true;

}

//Definition of The inequality (!=) operator
bool FIDModel::operator!= (const FIDModel& fid) {
    return !(*(this)==fid);
}

//Copy Constructor to Clone objects
FIDModel::FIDModel(const FIDModel& other) {
    for (int i=0;i<3;i++){
        for (int j=0;j<3;j++){
            this->ProbabilitiesMatrix[i][j] = other.ProbabilitiesMatrix[i][j];
        }
    }

    this->Lambda        = other.Get_Lambda();
    this->Gamma         = other.Get_Gamma();
}

//Function that evaluates all the transition probabilities for a given branch length
void FIDModel::EvaluateTransitionsMatrix(double length) {
    ProbabilitiesMatrix[BB][BB]       = (1 - ((1+(Lambda*length)-exp(-(Lambda*length)))/(Gamma*(1+(Lambda*length)))));
    ProbabilitiesMatrix[BB][B_]       = ((1-exp(-(Lambda*length)))/(Gamma*(1+(Lambda*length))));
    ProbabilitiesMatrix[BB][_B]       = ((Lambda*length)/(Gamma*(1+(Lambda*length))));

    ProbabilitiesMatrix[B_][BB]       = (((Lambda*length)*exp(-(Lambda*length)))/(Gamma*(1-exp(-(Lambda*length)))*(1+(Lambda*length))));
    ProbabilitiesMatrix[B_][B_]       = ((Gamma*(1+(Lambda*length))-1)/(Gamma*(1+(Lambda*length))));
    ProbabilitiesMatrix[B_][_B]       = ((1-exp(-(Lambda*length))*(1+(Lambda*length)))/(Gamma*(1-exp(-(Lambda*length)))*(1+(Lambda*length))));

    ProbabilitiesMatrix[_B][BB]       = ((exp(-(Lambda*length)))/(Gamma*(1+(Lambda*length))));
    ProbabilitiesMatrix[_B][B_]       = ((1-exp(-(Lambda*length)))/(Gamma*(1+(Lambda*length))));
    ProbabilitiesMatrix[_B][_B]       = ((Gamma*(1+(Lambda*length))-1)/(Gamma*(1+(Lambda*length))));
}

//Function that returns the transition probability of being at state prevEvent and transitioning to state nextEvent at time length
double FIDModel::GetProbability(EventType1 prevEvent, EventType1 nextEvent, double length) {
    EvaluateTransitionsMatrix(length);

    switch (prevEvent) {
        case BB:
            switch (nextEvent) {
                case BB:
                    return ProbabilitiesMatrix[BB][BB];
                    break;
                case _B:
                    return ProbabilitiesMatrix[BB][_B];
                    break;
                case B_:
                    return ProbabilitiesMatrix[BB][B_];
                    break;
            }
            break;

        case _B:
            switch (nextEvent) {
                case BB:
                    return ProbabilitiesMatrix[_B][BB];
                    break;
                case _B:
                    return ProbabilitiesMatrix[_B][_B];
                    break;
                case B_:
                    return ProbabilitiesMatrix[_B][B_];
                    break;
            }
            break;

        case B_:
            switch (nextEvent) {
                case BB:
                    return ProbabilitiesMatrix[B_][BB];
                    break;
                case _B:
                    return ProbabilitiesMatrix[B_][_B];
                    break;
                case B_:
                    return ProbabilitiesMatrix[B_][B_];
                    break;
            }
            break;
        }
}

//Function that returns the transition probability of being at state prevEvent and transitioning to state nextEvent at time length
double FIDModel::GetProbability(EventType2 prevEvent, EventType2 nextEvent, double length){
    EvaluateTransitionsMatrix(length);

    switch (prevEvent) {
        case H:
            switch (nextEvent) {
                case H:
                    return ProbabilitiesMatrix[BB][BB];
                    break;
                case B:
                    return ProbabilitiesMatrix[BB][_B];
                    break;
                case E:
                    return ProbabilitiesMatrix[BB][B_];
                    break;
            }
            break;

        case B:
            switch (nextEvent) {
                case H:
                    return ProbabilitiesMatrix[_B][BB];
                    break;
                case B:
                    return ProbabilitiesMatrix[_B][_B];
                    break;
                case E:
                    return ProbabilitiesMatrix[_B][B_];
                    break;
            }
            break;

        case E:
            switch (nextEvent) {
                case H:
                    return ProbabilitiesMatrix[B_][BB];
                    break;
                case B:
                    return ProbabilitiesMatrix[B_][_B];
                    break;
                case E:
                    return ProbabilitiesMatrix[B_][B_];
                    break;
            }
            break;
        }
}

//Function that returns the already estimated transition probabilities in a matrix form
const RealMatrix& FIDModel::Get_TransitionsMatrix() {
    RealMatrix probabilities;

    for (int i=0;i<3;i++){
        vector<double> values;
        for (int j=0;j<3;j++){
            values.push_back(ProbabilitiesMatrix[i][j]);
        }
        probabilities.push_back(values);
    }

    return probabilities;
}

//Function that estimates the transition probabilities for a given length and returns the already estimated transition probabilities in a matrix form
const RealMatrix& FIDModel::Get_TransitionsMatrix(double length) {
    EvaluateTransitionsMatrix(length);
    return Get_TransitionsMatrix();
}

//Definition of the toString Operation
ostream& operator<<(ostream& ostr, const FIDModel& model) {
    ostr << endl;
    ostr << "FID Model" << endl;
    ostr << "=========" << endl;
    ostr << "Expected Fragment Length (Gamma): " << model.Get_Gamma()  << endl;
    ostr << "Indel Rate Per Site (Lambda)    : " << model.Get_Lambda()  << endl;
    ostr << endl << endl;
}
