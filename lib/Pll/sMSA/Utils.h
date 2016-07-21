/****************************************************************************
 * This file is part of the program sMSA (statistical Multiple Sequence     *
 * Alignment) for joint sampling alignments and mutation parameter values   *
 * of n-given DNA-sequences.                                                *
 *                                                                          *
 *    Copyright (C) 2012  Dimitrios Lyras                                   *
 *    e-mail: dimlyras@gmail.com                                            *
 *                                                                          *
 *  This program is free software; you can redistribute it and/or modify    *
 *  it under the terms of the GNU General Public License as published by    *
 *  the Free Software Foundation; either version 2 of the License, or       *
 *  (at your option) any later version.                                     *
 *                                                                          *
 *  This program is distributed in the hope that it will be useful,         *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *  GNU General Public License for more details.                            *
 *                                                                          *
 *  You should have received a copy of the GNU General Public License       *
 *  along with this program; if not, write to the                           *
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330,             *
 *  Boston, MA  02111-1307  USA                                             *
 ****************************************************************************/

/***********************************************************************************/
/**                                 File Scopes:                                   */
/**  Aim: This file provides with set of useful functions in the form of a toolbox */
/**                              to be used at will.                               */
/***********************************************************************************/

#ifndef UTILS_H_INCLUDED
#define UTILS_H_INCLUDED

/** Header file Inclusions **/
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <map>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <unistd.h> //for the Jar Execution
#include <algorithm>
#include <functional>
#include <omp.h>    //For parallelizing the program


#include "BoostRandomGenerator.h"

using namespace std;


namespace Utils {

/***************************/
/** Templates for Vectors **/
/***************************/ //{
//To get the size (end) of a vector or array - call end(vectorName)
template<typename T, size_t N>
T * end(T (&ra)[N]) {
    return ra + N;
}

//Overload the toString Operator
template <typename T>
std::ostream& operator<<(std::ostream& ostr, const std::vector<T>& vect)     {
    for (int i=0; i<vect.size(); i++) {
        ostr << vect[i] << endl;
    }
    return ostr;
}

//Overloading the + operator
template <typename T>
std::vector<T> operator+(const std::vector<T>& a, const std::vector<T>& b) {
    assert(a.size() == b.size());

    std::vector<T> result;
    result.reserve(a.size());

    std::transform(a.begin(), a.end(), b.begin(),
                   std::back_inserter(result), std::plus<T>());
    return result;
}

//Overloading the - operator
template <typename T>
std::vector<T> operator-(const std::vector<T>& a, const std::vector<T>& b) {
    assert(a.size() == b.size());

    std::vector<T> result;
    result.reserve(a.size());

    std::transform(a.begin(), a.end(), b.begin(),
                   std::back_inserter(result), std::minus<T>());
    return result;
}

//Overloading the += operator
template <typename T>
std::vector<T>& operator+=(std::vector<T>& a, const std::vector<T>& b) {
    assert(a.size() == b.size());

    a = a + b;
    return a;
}

//Overloading the -= operator
template <typename T>
std::vector<T>& operator-=(std::vector<T>& a, const std::vector<T>& b) {
    assert(a.size() == b.size());

    a = a - b;
    return a;
}

//Overloading the += operator (scalar)
template <typename T>
std::vector<T>& operator+=(std::vector<T>& a, T b) {
    std::transform(a.begin(), a.end(), a.begin(), bind2nd(std::plus<T>(), b));
    //for multiplication change plus<T> to multiplies<T>
    return a;
}

//Overloading the -= operator (scalar)
template <typename T>
std::vector<T>& operator-=(std::vector<T>& a, T b) {
    std::transform(a.begin(), a.end(), a.begin(), bind2nd(std::minus<T>(), b));
    //for multiplication change plus<T> to multiplies<T>
    return a;
}

//Overloading the == operator
template <typename T>
bool operator==(const std::vector<T>& a, const std::vector<T>& b) {
    if (a.size() != b.size() ) return false;
    for (int i =0; i<a.size(); i++) {
        if (a[i] != b[i]) return false;
    }
    return true;
}

//Overloading the != operator
template <typename T>
bool operator!=(const std::vector<T>& a, const std::vector<T>& b) {
    if (a.size() == b.size() ) return false;
    for (int i =0; i<a.size(); i++) {
        if (a[i] == b[i]) return false;
    }
    return true;
}

//Function to Delete (Free memory) of a vector of pointers
template <class T>
void DeleteVectorOfPointers(std::vector<T>& vector) {
    typename T::iterator i;
    for (i = vector->begin(); i < vector->end(); ++i) {
        delete *i;
    }
    delete vector;
}


//Template to find a pair in a vector of pairs according to the .first() value
//Usage: it = std::find_if(vec.begin(), vec.begin(), match_first<int>(1));
//if (it!=vec.end()) -> means pair found
template <typename K>
struct match_first {
    const K _k;
    match_first(const K& k) : _k(k) {}
    template <typename V>
    bool operator()(const std::pair<K, V>& el) const {
        return _k == el.first;
    }
};

//Function to compare the values of 2 map pairs and return whether left map pair has a "second" value Larger to the right map pair
inline bool Larger(const std::pair<std::string, double>& lhs, const std::pair<std::string, double>& rhs) {
    return lhs.second < rhs.second;
}

//Function that returns the Largest "second" value of a map
inline double LargestValue(const std::map<std::string, double>& m) {
    return max_element(m.begin(), m.end(), Larger)->second;
}

//Function that returns the index of the Largest "second" value of a map
inline int IndexOfLargestValue(std::map<std::string, double>& m) {
    std::map<std::string, double>::iterator it = max_element(m.begin(), m.end(), Larger);
    return std::distance(m.begin(), it);
}

//Function to return the position of most probable evolutionary history based on the log of the overall probability of the evolutionary history
template <typename T>
int GetPosOfMaxElement(vector<T>& vec) {
    return max_element(vec.begin(), vec.end()) - vec.begin();
}

//Function to return -1 if the given vector does not contain the given value, or the pos (first appearance) of this value within the given vector, otherwise.
template <typename T>
int GetPosOfElement(vector<T>& vec, T value, int StartPos) {
    typename std::vector<T>::iterator it;
    it = std::find(vec.begin()+StartPos, vec.end(), value);

    if( it != vec.end()) {
        return it - vec.begin();
    } else {
        return -1;
    }
}

//Function to return the first non-negative value found within the given vector
template <typename T>
T GetFirstNonNegativeValue(vector<T>& vec, int StartPos) {
    for (int i=StartPos; i<vec.size(); i++) {
        if (vec[i]>=0) {
            return vec[i];
        }
    }
    return -1;
}

//Function to return the first non-negative value found within the given vector
template <typename T>
T GetFirstNonNegativeValue(const vector<T>& vec, int StartPos) {
    for (int i=StartPos; i<vec.size(); i++) {
        if (vec[i]>=0) {
            return vec[i];
        }
    }
    return -1;
}

//Function to return the position (index) of the first non-negative value found within the given vector (const version)
template <typename T>
T GetPosOfFirstNonNegativeValue(const vector<T>& vec, int StartPos) {
    for (int i=StartPos; i<vec.size(); i++) {
        if (vec[i]>=0) {
            return i;
        }
    }
    return -1;
}

//Function to return the first non-negative value found within the given vector but going backwards from the given start position
template <typename T>
T GetFirstNonNegativeValueReverse(vector<T>& vec, int StartPos) {
    for (int i=StartPos; i>=0; i--) {
        if (vec[i]>=0) {
            return vec[i];
        }
    }
    return -1;
}

//Function to return the first non-negative value found within the given vector but going backwards from the given start position
template <typename T>
T GetFirstNonNegativeValueReverse(const vector<T>& vec, int StartPos) {
    for (int i=StartPos; i>=0; i--) {
        if (vec[i]>=0) {
            return vec[i];
        }
    }
    return -1;
}

//Function to return the last non-negative value found within the given vector
template <typename T>
T GetLastNonNegativeValue(vector<T>& vec, int EndPos) {
    for (int i=vec.size()-1; i>=EndPos; i--) {
        if (vec[i]>=0) {
            return vec[i];
        }
    }
    return -1;
}

//Function to return the position (index) of the last non-negative value found within the given vector (const version)
template <typename T>
T GetPosOfLastNonNegativeValue(const vector<T>& vec, int StartPos) {
    for (int i=StartPos; i>=0; i--) {
        if (vec[i]>=0) {
            return i;
        }
    }
    return -1;
}

//Function to return the number of times the given value is found withing the the given vector - returns -1 if not found
template <typename T>
T GetValueCount(vector<T>& vec, T value) {
    int counter = 0;
    for (int i=0; i<vec.size(); i++) {
        if (vec[i]==value) {
            counter++;
        }
    }
    return counter;
}

//Function to Compare the Maps Elements based on their Keys
typedef std::map<int, std::string> M;
inline bool Map_Value_ComparerKey(M::value_type &i1, M::value_type &i2) {
    return i1.first<i2.first;
}

//Function to Compare the Maps Elements based on their Values
inline bool Map_Value_ComparerValue(M::value_type &i1, M::value_type &i2) {
    return i1.second<i2.second;
}

//Function to Get the Biggest Key Value of a Map
inline int GetBiggestKeyValue(M map){
    M::iterator itor = std::max_element(map.begin(), map.end(), Map_Value_ComparerKey);
    return (*itor).first;
}

//}

/*******************************/
/** File Management Functions **/
/*******************************/ //{
//Function to Save a string (with contents) to a File
inline bool WriteToFile(string filename, string contents, bool append, bool clearContents) {
    ofstream myfile;

    if (clearContents) {
        myfile.open(filename.c_str(),ios::out|ios::trunc);
        myfile.close();
    }

    if (append) {
        myfile.open(filename.c_str(), fstream::out | fstream::app);
    } else {
        myfile.open(filename.c_str(), fstream::out);
    }

    if (myfile.is_open()) {
        myfile << contents;
        myfile.close();

        return true;
    } else cerr << "Unable to open file: "<<filename<<" for writing.";

    return false;
}

//Function to Read Contents from File and get the returned result as a string
inline string ReadFromFile(string filename) {
    string contents = "";
    string line;
    ifstream myfile (filename.c_str());
    if (myfile.is_open()) {
        while ( myfile.good() ) {
            getline (myfile,line);
            contents+=line;
        }
        myfile.close();
    } else cerr << "Unable to open file: "<<filename<<" for reading.";

    return contents;
}

inline bool FileExists(string filename) {
    ifstream myfile (filename.c_str());
    if (myfile){
        return true;
    }else{
        return false;
    }
}

inline int DeleteFile(string filename) {
    if (remove(filename.c_str())==0){
        return 0;
    }else{
        return 1;
    }
}


//}

/************************************/
/** Runtime Intervention Functions **/
/************************************/ //{
//Function that Pauses the Command-line output Until a key is pressed by the user
inline void Pause() {
    fflush(stdin);
    cout<<endl<<"Press enter to continue..."<<endl;
    getchar();
}

//Function that Pauses (puts in sleep) the program for given milliseconds
inline void Sleep(int mseconds ) {
    clock_t goal = mseconds + clock();
    while (goal > clock());
}

//Function to Execute a CMD Application
inline void ExecuteProcess(string ProcessToExecute) {
    system(ProcessToExecute.c_str());
}

//Function to flush stdin
inline void ClearInputStream() {
    rewind(stdin);
}

//Function to clear CMD
inline void ClearScreen() {
    //cout << string( 100, '\n' ); //->safe code for all operating systems, though not very good solution
    system("cls");   //for windows
    //system("clear"); //for other systems
}

//Function to Run Specified Jar File
inline void ExecuteJar(string JarFileName) {
    execlp("java", "java", "-jar", JarFileName.c_str(), (char *)0);
}
//}

/****************************************/
/** Random Number Generation Functions **/
/****************************************/ //{
//Native c++: Function Function that Returns a double value withing the range [rangeMin, rangeMax]
inline double GetRandomRealNumber(double rangeMin, double rangeMax) {
    double f = (double)rand() / RAND_MAX;
    return rangeMin + f * (rangeMax - rangeMin);
}

//Function to Generate n-Random Numbers whose sum equals to the given argument value
inline void Get_N_RandomProbabilities(BoostRandomGenerator* rand, int n, vector<double>& numsVec) {
    double sum = 0;

    for (int i=0; i<n-1; i++) {
        int num = rand->GetRandomInt(1, 100 - (n-(i+1)) - sum -1);
        numsVec.push_back(num);
        sum += num;
    }
    numsVec.push_back(100-sum);

    for (int i=0; i<numsVec.size(); i++) {
        numsVec[i] /= 100.0;
        //cout<<numsVec[i]<<endl;
    }
}

//Function to Generate n-Random Sequences in the form of Probabilities Vectors and save results to File
inline void Get_N_RandomSequences(int n, vector<int>& seqLens) {
    //Init the random Generator
    BoostRandomGenerator rand;
    rand.SetSeed(-1); //to seed with random values

    WriteToFile("SeqsVecs.txt", "", false, true);
    for (int l=0; l<n; l++) {
        vector<vector<double> > randVec;
        for (int i=0; i<seqLens[l]; i++) {
            vector<double> tempRand;
            Get_N_RandomProbabilities(&rand, 4, tempRand);
            randVec.push_back(tempRand);
        }

        for (int i=0; i<randVec[0].size(); i++) {
            for (int j=0; j<randVec.size(); j++) {
                stringstream ss;
                ss.str("");
                ss << randVec[j][i]<<"\t";
                WriteToFile("SeqsVecs.txt", ss.str(), true, false);
            }
            WriteToFile("SeqsVecs.txt", "\n", true, false);
        }

        WriteToFile("SeqsVecs.txt", "\n", true, false);
    }
}
//}

/*********************************/
/** Sequence Creation Functions **/
/*********************************/ //{
//Function that Returns a Random Nucleotide based on the given array of Base Frequencies and a given Generator
inline string GetRandomNucleotide(const double* freq, BoostRandomGenerator& numberGenerator) {
    double rand_num = numberGenerator.GetRandomDouble(0, 1);
    int j;
    for (j=0; rand_num>0; j++) {
        rand_num-= *(freq+j);
    }

    switch(j-1) {
    case 0:
        return "A";
    case 1:
        return "C";
    case 2:
        return "G";
    case 3:
        return "T";
    }
}

//Function that Creates and Returns a Random Sequence of seqLen Nucleotides
inline string GetRandomNucleotideSequence(int seqLen, const double* freq, bool UseRandomSeed) {
    const int rangeMin = 0;
    const int rangeMax = 1;

    BoostRandomGenerator rand;
    if (UseRandomSeed) {
        rand.SetSeed(std::time(0)); // seed with the current time
    } else {
        rand.SetSeed(59u); // seed with the specific value to repeat experiment at will
    }

    stringstream rand_seq;
    rand_seq.str("");
    for (int i=0; i<seqLen; i++) {
        rand_seq<<GetRandomNucleotide(freq, rand);
    }

    return rand_seq.str();
}
//}

/****************************/
/** Mathematical Functions **/
/****************************/ //{
//Function that Calculates the differnece between 2 double values and returns whether they differ less than a very small value or not
inline bool PracticallyEquals(const double & v1, const double & v2, const double & epsilon =1e-9) {
    const double diff = v1 - v2;
    return diff < epsilon && diff > -epsilon;
}

//Function that returns a random exponential value of the given double value
inline double RandExp(double x) {
    double R=rand()/(RAND_MAX+1.0);
    if ( R>0.99 ) return -log(R)/x;
    else return -log(0.99)/x+RandExp(x);
}

//Function that returns a random Quotient value of the given double value
inline double RandQuot(double q) {
    double x=q/(1+q), r=x+rand()/(RAND_MAX*10.0)-1.0/20.0;
    x = (r<0 ? 1+r : (r>=1 ? r-1 : r));
    return  x/(1-x);
}

//Function that returns the sum of 2 logarithmetic input values
inline double AddLog(double a, double b) {
    if (a>b) return a+log(1+exp(b-a));
    else return b+log(1+exp(a-b));
}

//Function that returns the difference of 2 logarithmetic input values
inline double RemoveLog(double a, double b) {
    if (a>b) return a-log(1+exp(b-a));
    else return b-log(1+exp(a-b));
}

//Function that returns the sum of 3 logarithmetic input values
inline double AddLog(double a, double b, double c) {
    return AddLog(AddLog(a,b),c);
}

//Function that returns the sum of 3 logarithmetic input values
inline double RemoveLog(double a, double b, double c) {
    return RemoveLog(RemoveLog(a,b),c);
}

//Function that returns the maximumn between 2 input values
inline double maximum(double a, double b) {
    return ( a>b ? a : b );
}

//Function that returns the maximumn between 3 input values
inline double maximum(double a, double b, double c) {
    return maximum(maximum(a,b),c);
}

//Function that returns the maximumn between 2 input values
inline int maximum(int a, int b) {
    return ( a>b ? a : b );
}

//Function that returns the maximumn between 3 input values
inline int maximum(int a, int b, int c) {
    return maximum(maximum(a,b),c);
}

//Function that returns the minimum between 2 input values
inline double minimum(double a, double b) {
    return ( a<b ? a : b );
}

//Function that returns the minimum between 3 input values
inline double minimum(double a, double b, double c) {
    return minimum(minimum(a,b),c);
}

//Function that returns the minimum between 2 input values
inline int minimum(int a, int b) {
    return ( a<b ? a : b );
}

//Function that returns the minimum between 3 input values
inline int minimum(int a, int b, int c) {
    return minimum(minimum(a,b),c);
}
//}

/***********************************/
/** String Manipulation Functions **/
/***********************************/ //{
//Function that counts how many times a string is contained in another string
inline size_t StringCount(const std::string& referenceString, const std::string& subString) {
    const size_t step = subString.size();

    size_t count(0);
    size_t pos(0) ;

    while( (pos=referenceString.find(subString, pos)) !=std::string::npos) {
        pos +=step;
        ++count ;
    }

    return count;

    /*  EDIT: This function does not allow for overlaps, i.e. searching for sub-string "AA" in string "AAAAAAAA" results in a count of 4. To allow for overlap, this line

        pos += step
        should be replaced by

        ++pos
        This will result in a count of 7.
    */

}

//Function that replaces every appearance of part of a string with another string
inline void ReplaceAll(std::string& str, const std::string& from, const std::string& to) {
    size_t start_pos = 0;
    while((start_pos = str.find(from, start_pos)) != std::string::npos) {
        str.replace(start_pos, from.length(), to);
        start_pos += to.length(); // In case 'to' contains 'from', like replacing 'x' with 'yx'
    }
}

//Function that finds whether a string ends with another string or not
inline bool EndsWith (std::string const &fullString, std::string const &ending) {
    if (fullString.length() >= ending.length()) {
        return (0 == fullString.compare (fullString.length() - ending.length(), ending.length(), ending));
    } else {
        return false;
    }
}

//Function that returns the reversed version of a given string
inline string GetReversedString(string source) {
    return string ( source.rbegin(), source.rend() );
}

//Function that Returns the first position (index) where a string (target) is found in another string (source) or returns -1 if not found
inline int FirstIndexOf(string source, string target) {
    string::size_type loc = source.find( target, 0 );
    if( loc != string::npos ) {
        return loc;
    } else {
        return -1;
    }
}

//Overloaded version of the previous Function. Acts as previously but starting after the given starting position (stratPos)
inline int FirstIndexOf(string source, string target, int startPos) {
    string::size_type loc = source.find( target, startPos );
    if( loc != string::npos ) {
        return loc;
    } else {
        return -1;
    }
}

//Function that Returns the last position (index) where a string (target) is found in another string (source) or returns -1 if not found
inline int LastIndexOf(string source, string target) {
    int pos = FirstIndexOf(GetReversedString(source), target);
    if (pos==-1) return -1;
    else return source.size() -1 - pos;
}

//Function that Returns whether a given string contains another string or not
inline bool Contains(string source, string target) {
    size_t found = source.find(target);
    if (found!=string::npos) return true;
    else                     return false;
}

//Function that Returns whether a string starts with the given preffix of not
inline bool StartsWith (std::string const &fullString, std::string const &prefix) {
    return prefix.length() <= fullString.length() && equal(prefix.begin(), prefix.end(), fullString.begin());
}

//Function that returns a vector with the indices where a given subString is found in a given referenceString
inline void GetIndicesOf(vector<int>& indices, const std::string& referenceString, const std::string& subString) {
    const size_t step = subString.size();
    size_t pos(0) ;

    while( (pos=referenceString.find(subString, pos)) !=std::string::npos) {
        indices.push_back(pos);
        pos +=step;
    }
}

//Function that Splits a Given String According to a delimiter and returns an already constructed container also given as argument
inline std::vector<std::string> &Split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while(std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

//Function that Splits a Given String According to a delimiter and returns an already constructed container also given as argument
inline void SplitValues(const std::string &s, char delim, std::vector<double> &elems) {
    std::stringstream ss(s);
    std::string item;
    int i=0;
    while(std::getline(ss, item, delim)) {
        elems[i++] = (atof(item.c_str()));
    }
}

//Function that Splits a Given String According to a delimiter and returns a new container with the results
inline std::vector<std::string> Split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    return Split(s, delim, elems);
}
//}

/*****************************/
/** Combinatorics Functions **/
/*****************************/ //{
inline void GetAllCombinations(std::string s, size_t comboLength, size_t numOfValuesToUse, const vector<string>& labels, string TreeType, bool clearVector, vector<string>& combinations ) {
    static int seqLength = comboLength;

    if (clearVector) {
        seqLength = comboLength;
    }

    if( comboLength )
        for( size_t digit = 0; digit < numOfValuesToUse; ++digit )
            GetAllCombinations( s + static_cast<char>( digit + '0' ), comboLength - 1, numOfValuesToUse, labels, TreeType, false, combinations);
    else {
        string seq(s);

        for (unsigned int i=0; i<labels.size(); i++) {
            std::stringstream ss;
            ss << i;
            ReplaceAll(seq, ss.str(), labels[i]);
        }

        if (TreeType.compare("TihlTree")==0) {
            if (StringCount(seq, "B")<=1) { //No more than one B(irth) is allowed per tihl!
                if (StringCount(seq, "-")==0 && StringCount(seq, "B")==0 ) {   //To account for the B--, -B-, and --B cases
                    if      (seqLength==2 && StartsWith(seq, "H") ) { combinations.push_back(seq); }
                    else if (seqLength>2)                           { combinations.push_back(seq); }
                } else if (StringCount(seq, "B")==1 && StringCount(seq, "-")==seqLength-1 ) {
                    if      (seqLength==2 && !StartsWith(seq, "B")) { combinations.push_back(seq); }
                    else if (seqLength>2)                           { combinations.push_back(seq); }
                }
            }
        } else if (TreeType.compare("SoanTree")==0) {
            combinations.push_back(seq);
        }
    }
}
//}

/*****************************************/
/** For Dynamic Multidimensional Arrays **/
/*****************************************/ //{
//Function that returns the product of the Dimensions of a vector starting from the Given starting position
inline long long int GetDimensionsProduct(int StartPos, vector<int> Dimensions) {
    if (StartPos ==  Dimensions.size()) return 1;
    long long int product = 1;

    for (unsigned int i=StartPos; i<Dimensions.size(); i++) {
        product *= Dimensions[i];
    }
    return product;
}

//Function that returns the equivalent 1d position in a vector for the given Indexes of an N-d (N=Dimesions.size()) array or vector
inline long long int GetPosForIndexes(vector<int> Indexes, vector<int> Dimensions) {
    long long int pos = 0;
    for (int i=0; i<Indexes.size(); i++) {
        pos += Indexes[i] * GetDimensionsProduct(i+1, Dimensions);
    }

    return pos;
}

//Function that returns a vector of indices for an N-d (N=Dimesions.size()) array or vector based on the given index of the equivalent 1-d array or vector
inline vector<int> GetIndexesForPos(int Pos, const vector<int>& Dimensions) {
    vector<int> Indexes (Dimensions.size(), 0 );
    int i_val = 0, r_val = 0, idx = 1;
    int DimensionsProduct = (int)GetDimensionsProduct(idx, Dimensions);
    i_val = (int)Pos / DimensionsProduct;
    r_val = (int)Pos % DimensionsProduct;

    if(i_val == r_val) {
        for (int i=0; i<Dimensions.size(); i++) {
            Indexes[i] = i_val;
        }
    }

    //while(i_val != r_val){
    for (int i=0; i<Dimensions.size(); i++) {
        Indexes[idx-1] = i_val;
        idx++;
        DimensionsProduct = (int)GetDimensionsProduct(idx, Dimensions);
        i_val = (int)r_val / DimensionsProduct;
        r_val = (int)r_val % DimensionsProduct;
    }

    return Indexes;
}

//Function that returns a vector of indices for an N-d (N=Dimesions.size()) array or vector based on the given index of the equivalent 1-d array or vector
inline const vector<int>&GetIndexesForPos(int Pos, const vector<int>& Dimensions, vector<int>& Indexes) {
    Indexes.resize(Dimensions.size(), 0 );
    int i_val = 0, r_val = 0, idx = 1;
    int DimensionsProduct = (int)GetDimensionsProduct(idx, Dimensions);
    i_val = (int)Pos / DimensionsProduct;
    r_val = (int)Pos % DimensionsProduct;

    if(i_val == r_val) {
        for (int i=0; i<Dimensions.size(); i++) {
            Indexes[i] = i_val;
        }
    }

    //while(i_val != r_val){
    for (int i=0; i<Dimensions.size(); i++) {
        Indexes[idx-1] = i_val;
        idx++;
        DimensionsProduct = (int)GetDimensionsProduct(idx, Dimensions);
        i_val = (int)r_val / DimensionsProduct;
        r_val = (int)r_val % DimensionsProduct;
    }

    return Indexes;
}
//}

}



#endif // UTILS_H_INCLUDED
