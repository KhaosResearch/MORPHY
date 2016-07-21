#ifndef BOOSTRANDOMGENERATOR_H
#define BOOSTRANDOMGENERATOR_H

/** Header file Inclusions **/
#include <iostream>
#include <ctime>
#include <random>

typedef std::mt19937 MyRNG;  // the Mersenne Twister with a popular choice of parameters

using namespace std;

class BoostRandomGenerator {
    public:
        BoostRandomGenerator() { };
        virtual ~BoostRandomGenerator() {};
        void SetSeed(int seed){
            if (seed==-1){
                rng.seed(std::time(0));
            }else{
                rng.seed(seed);
            }
        }

        int GetRandomInt(int min, int max){
            std::uniform_int_distribution<> dist(min, max);
            return dist(rng);
        }

        double GetRandomDouble(double min, double max){
            std::uniform_real_distribution<> dist(min, max);
            return dist(rng);
        }
        //Other possible distributions:
        //std::normal_distribution normal_dist(mean, stddeviation);  // N(mean, stddeviation)
    protected:
    private:
        MyRNG rng;                   // e.g. keep one global instance (per thread)
};

#endif // BOOSTRANDOMGENERATOR_H
