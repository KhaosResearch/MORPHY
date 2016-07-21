#include "TN93Model.h"

TN93Model::TN93Model(const double* freq, double kappa1, double kappa2) {
    k1 = kappa1;
    k2 = kappa2;

    //Populate Rate Matrix
    RateMatrix[AA] = -(freq[C] + k2*freq[G] + freq[T]);  //A->A
    RateMatrix[AC] = freq[C];                            //A->C
    RateMatrix[AG] = k2*freq[G];                         //A->G
    RateMatrix[AT] = freq[T];                            //A->T

    RateMatrix[CA] = freq[A];                            //C->A
    RateMatrix[CC] = -(freq[A] + k1*freq[T] + freq[G]);  //C->C
    RateMatrix[CG] = freq[G];                            //C->G
    RateMatrix[CT] = k1*freq[T];                         //C->T

    RateMatrix[GA] = k2*freq[A];                         //G->A
    RateMatrix[GC] = freq[C];                            //G->C
    RateMatrix[GG] = -(k2*freq[A] + freq[C] + freq[T]);  //G->G
    RateMatrix[GT] = freq[T];                            //G->T

    RateMatrix[TA] = freq[A];                            //T->A
    RateMatrix[TC] = k1*freq[C];                         //T->C
    RateMatrix[TG] = freq[G];                            //T->G
    RateMatrix[TT] = -(freq[A] + k1*freq[C] + freq[G]);  //T->T

}

TN93Model::~TN93Model()  {
}

//Function to set the length of the model and calculate substitution matrix
void TN93Model::SetLength(double length) {
    //Multiply Rate Matrix by time (length)
    double RateMatrixTemp[16];
    for (int i=0;i<16;i++){
        RateMatrixTemp[i] = RateMatrix[i]*length;
    }

    //Estimate Matrix Exponential of RateMatrixTemp
    double *MatrixExponential = expm2 ( 4, RateMatrixTemp );

    for (int i=0;i<16;i++){
        SubMatrix[i] = MatrixExponential[i];
    }
}

//Function to estimate the substitution probability between the given nucleotides at length time
double TN93Model::GetProbability(Alphabet a, Alphabet b, double length) const{
    //Multiply Rate Matrix by time (length)
    double RateMatrixTemp[16];
    for (int i=0;i<16;i++){
        RateMatrixTemp[i] = RateMatrix[i]*length;
    }

    //Estimate Matrix Exponential of RateMatrixTemp
    double *MatrixExponential = expm2 ( 4, RateMatrixTemp );

    return MatrixExponential[GetAlphabetPair(a, b)];
}

//Function to estimate the substitution probability between the given nucleotides
double TN93Model::GetProbability(Alphabet a, Alphabet b) const{
    return SubMatrix[GetAlphabetPair(a, b)];
}

//Function to return the calculate substitution matrix for the given length
double* TN93Model::GetSubstitutionMatrix(double length) {
    //Multiply Rate Matrix by time (length)
    double RateMatrixTemp[16];
    for (int i=0;i<16;i++){
        RateMatrixTemp[i] = RateMatrix[i]*length;
    }

    //Estimate Matrix Exponential of RateMatrixTemp
    double *MatrixExponential = expm2 ( 4, RateMatrixTemp );

    return MatrixExponential;
}
