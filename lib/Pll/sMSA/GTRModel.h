#ifndef GTRMODEL_H
#define GTRMODEL_H

#include "Mappings.h"
#include "matrix_exponential.hpp"

using namespace Mappings;

class GTRModel {
    public:
        GTRModel(const double* freq, double x_1, double x_2, double x_3, double x_4, double x_5, double x_6);  //Default Constructor
        virtual ~GTRModel();                                              //Virtual Deconstructor

        void   SetLength(double length);                                    //Function to set the length of the model and calculate substitution matrix
        double GetProbability(Alphabet a, Alphabet b, double length) const; //Function to estimate the substitution probability between the given nucleotides at length time
        double GetProbability(Alphabet a, Alphabet b)                const; //Function to estimate the substitution probability between the given nucleotides at length time
        double* GetSubstitutionMatrix(double length);                       //Function to return the calculate substitution matrix for the given length
    protected:
    private:
        double RateMatrix[4*4];             //the 4x4 Rate Matrix
        double SubMatrix[4*4];              //the 4x4 Substitution Matrix
        double x1;                          //The x1 substitution parameter
        double x2;                          //The x2 substitution parameter
        double x3;                          //The x3 substitution parameter
        double x4;                          //The x4 substitution parameter
        double x5;                          //The x5 substitution parameter
        double x6;                          //The x6 substitution parameter
};

#endif // GTRMODEL_H
