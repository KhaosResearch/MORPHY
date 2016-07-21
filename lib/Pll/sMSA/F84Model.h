#ifndef F84MODEL_H
#define F84MODEL_H

#include "Mappings.h"
#include "matrix_exponential.hpp"

using namespace Mappings;

class F84Model {
    public:
        F84Model(const double* freq, double kappa);             //Default Constructor
        virtual ~F84Model();                              //Virtual Deconstructor

        void   SetLength(double length);                                    //Function to set the length of the model and calculate substitution matrix
        double GetProbability(Alphabet a, Alphabet b, double length) const; //Function to estimate the substitution probability between the given nucleotides at length time
        double GetProbability(Alphabet a, Alphabet b) const;                //Function to estimate the substitution probability between the given nucleotides at length time
        double* GetSubstitutionMatrix(double length);                       //Function to return the calculate substitution matrix for the given length
    protected:
    private:
        double RateMatrix[4*4];             //the 4x4 Rate Matrix
        double SubMatrix[4*4];              //the 4x4 Substitution Matrix
        double k;                           //the kappa parameter that distinguishes between the rate of transitions and transversions
};

#endif // F84MODEL_H
