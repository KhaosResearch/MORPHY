#include "GTRModel.h"

GTRModel::GTRModel(const double* freq, double x_1, double x_2, double x_3, double x_4, double x_5, double x_6) {
    x1 =  x_1;
    x2 =  x_2;
    x3 =  x_3;
    x4 =  x_4;
    x5 =  x_5;
    x6 =  x_6;


    //Populate Rate Matrix
    RateMatrix[AA] = -(x1+x2+x3);                                                       //A->A
    RateMatrix[AC] = (freq[A]*x1)/freq[C];                                              //A->C
    RateMatrix[AG] = (freq[A]*x2)/freq[G];                                              //A->G
    RateMatrix[AT] = (freq[A]*x3)/freq[T];                                              //A->T

    RateMatrix[CA] = x1;                                                                //C->A
    RateMatrix[CC] = -(((freq[A]*x1)/freq[C]) + x4 + x5);                               //C->C
    RateMatrix[CG] = (freq[C]*x4)/freq[G];                                              //C->G
    RateMatrix[CT] = (freq[C]*x5)/freq[T];                                              //C->T

    RateMatrix[GA] = x2;                                                                //G->A
    RateMatrix[GC] = x4;                                                                //G->C
    RateMatrix[GG] = -(((freq[A]*x2)/freq[G]) + ((freq[C]*x4)/freq[G]) + x6);           //G->G
    RateMatrix[GT] = (freq[G]*x6)/freq[T];                                              //G->T

    RateMatrix[TA] = x3;                                                                           //T->A
    RateMatrix[TC] = x5;                                                                           //T->C
    RateMatrix[TG] = x6;                                                                           //T->G
    RateMatrix[TT] = -(((freq[A]*x3)/freq[T]) + ((freq[C]*x5)/freq[T]) + ((freq[G]*x6)/freq[T]));  //T->T

}

GTRModel::~GTRModel()  {
}

//Function to set the length of the model and calculate substitution matrix
void GTRModel::SetLength(double length) {
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
double GTRModel::GetProbability(Alphabet a, Alphabet b, double length) const{
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
double GTRModel::GetProbability(Alphabet a, Alphabet b) const{
    return SubMatrix[GetAlphabetPair(a, b)];
}

//Function to return the calculate substitution matrix for the given length
double* GTRModel::GetSubstitutionMatrix(double length) {
    //Multiply Rate Matrix by time (length)
    double RateMatrixTemp[16];
    for (int i=0;i<16;i++){
        RateMatrixTemp[i] = RateMatrix[i]*length;
    }

    //Estimate Matrix Exponential of RateMatrixTemp
    double *MatrixExponential = expm2 ( 4, RateMatrixTemp );

    return MatrixExponential;
}
