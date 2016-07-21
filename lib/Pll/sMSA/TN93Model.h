#ifndef TN93MODEL_H
#define TN93MODEL_H

#include "Mappings.h"
#include "matrix_exponential.hpp"

using namespace Mappings;

class TN93Model {
    public:
        TN93Model(const double* freq, double kappa1, double kappa2);             //Default Constructor
        virtual ~TN93Model();                                              //Virtual Deconstructor

        void   SetLength(double length);                                    //Function to set the length of the model and calculate substitution matrix
        double GetProbability(Alphabet a, Alphabet b, double length) const; //Function to estimate the substitution probability between the given nucleotides at length time
        double GetProbability(Alphabet a, Alphabet b)                const; //Function to estimate the substitution probability between the given nucleotides at length time
        double* GetSubstitutionMatrix(double length);                       //Function to return the calculate substitution matrix for the given length
    protected:
    private:
        double RateMatrix[4*4];             //the 4x4 Rate Matrix
        double SubMatrix[4*4];              //the 4x4 Substitution Matrix
        double k1;                          //The A<->G rate
        double k2;                          //The C<->T rate
};

#endif // TN93MODEL_H
