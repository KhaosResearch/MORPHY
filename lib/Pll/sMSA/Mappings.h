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
/**  Aim: This files defines with the proper structires to efficiently represent   */
/**   a nucleotides->integer mapping as well as the mapping between tihl events    */
/***********************************************************************************/

#ifndef MAPPINGS_H_INCLUDED
#define MAPPINGS_H_INCLUDED

/** Header file Inclusions **/
#include <string>
#include <iostream>

using namespace std;

/** Global Definition of some constants needed for the Nucleotides Mapping **/
const unsigned int ALPHABET_SIZE = 4;
const unsigned int PURINE = 1;

namespace Mappings{
    //Nucleotides Mapping
    enum Alphabet   { A=0, C=1, G=2, T=3 };

    //AminoAcids Mapping
    enum AlphabetAA   { AmA=0, AmR=1, AmN=2, AmD=3, AmC=4, AmQ=5, AmE=6, AmG=7, AmH=8, AmI=9, AmL=10, AmK=11, AmM=12, AmF=13, AmP=14, AmS=15, AmT=16, AmW=17, AmY=18, AmV=19, AmB=20, AmJ=21, AmZ=22, AmX=23, AmGAP=24 };

    //Nucleotides Substitutions Mapping
    enum AlphabetPair   { AA=0, AC=1, AG=2, AT=3, CA=4, CC=5, CG=6, CT=7, GA=8, GC=9, GG=10, GT=11, TA=12, TC=13, TG=14, TT=15  };

    //Events Mapping in BB, B_, _B form
    enum EventType1 { BB = 0, B_ = 1, _B = 2 };

    //Events Mapping in H, B, E, N form
    enum EventType2 { H = 0, B = 1, N = 2, E = 3} ;

    /** For Alphabet **/
    //Function that returns the proper Alphabet according to the given character
    inline Alphabet GetAlphabet(const char& c) {
        switch (c) {
        case 'A' :
        case 'a' :
            return A;
        case 'C' :
        case 'c' :
            return C;
        case 'G' :
        case 'g' :
            return G;
        case 'T' :
        case 't' :
            return T;
        default :
            std::cout << "Error! - Invalid nucleotide " << c << "\n";
        }
    }

    //Function that returns the proper Alphabet according to the given string (overloaded Version)
    inline Alphabet GetAlphabet(const string& s) {
        if      (s.compare("A")==0 || s.compare("a")==0 ) return A;
        else if (s.compare("C")==0 || s.compare("c")==0 ) return C;
        else if (s.compare("G")==0 || s.compare("g")==0 ) return G;
        else if (s.compare("T")==0 || s.compare("t")==0 ) return T;
        else { std::cerr << "Error! - Invalid nucleotide " << s << "\n"; }
    }

    //Function that returns the proper Alphabet according to the given integer (overloaded Version)
    inline Alphabet GetAlphabet(const int& s) {
        if      ( s==0 ) return A;
        else if ( s==1 ) return C;
        else if ( s==2 ) return G;
        else if ( s==3 ) return T;
        else { std::cerr << "Error! - Invalid nucleotide " << s << "\n"; }
    }

    //Function that returns the proper Alphabet in string format according to the given Alphabet (nucleotide)
    inline string GetAlphabet_String(const Alphabet& s) {
        if      ( s==A ) return "A";
        else if ( s==C ) return "C";
        else if ( s==G ) return "G";
        else if ( s==T ) return "T";
        else { std::cerr << "Error! - Invalid nucleotide " << s << "\n"; }
    }

    //Function that returns the proper AlphabetAA according to the given integer (overloaded Version)
    inline AlphabetAA GetAlphabetAA(const int& s) {
		if      (s==0 ) return AmA;
        else if (s==1 ) return AmR;
        else if (s==2 ) return AmN;
        else if (s==3 ) return AmD;
        else if (s==4 ) return AmC;
        else if (s==5 ) return AmQ;
        else if (s==6 ) return AmE;
        else if (s==7 ) return AmG;
        else if (s==8 ) return AmH;
        else if (s==9 ) return AmI;
        else if (s==10) return AmL;
        else if (s==11) return AmK;
        else if (s==12) return AmM;
        else if (s==13) return AmF;
        else if (s==14) return AmP;
        else if (s==15) return AmS;
        else if (s==16) return AmT;
        else if (s==17) return AmW;
        else if (s==18) return AmY;
        else if (s==19) return AmV;
        else if (s==20) return AmB;
        else if (s==21) return AmJ;
        else if (s==22) return AmZ;
        else if (s==23) return AmX;
        else if (s==24) return AmGAP;
        else { std::cerr << "Error! - Invalid nucleotide " << s << "\n"; }
    }

    //Function that returns the proper AlphabetAA according to the given string (overloaded Version)
    inline int GetAlphabetAA_Int(const string& s) {
        if      (s.compare("A")==0 || s.compare("a")==0 ) return 0;
        else if (s.compare("R")==0 || s.compare("r")==0 ) return 1;
        else if (s.compare("N")==0 || s.compare("n")==0 ) return 2;
        else if (s.compare("D")==0 || s.compare("d")==0 ) return 3;
        else if (s.compare("C")==0 || s.compare("c")==0 ) return 4;
        else if (s.compare("Q")==0 || s.compare("q")==0 ) return 5;
        else if (s.compare("E")==0 || s.compare("e")==0 ) return 6;
        else if (s.compare("G")==0 || s.compare("g")==0 ) return 7;
        else if (s.compare("H")==0 || s.compare("h")==0 ) return 8;
        else if (s.compare("I")==0 || s.compare("i")==0 ) return 9;
        else if (s.compare("L")==0 || s.compare("l")==0 ) return 10;
        else if (s.compare("K")==0 || s.compare("k")==0 ) return 11;
        else if (s.compare("M")==0 || s.compare("m")==0 ) return 12;
        else if (s.compare("F")==0 || s.compare("f")==0 ) return 13;
        else if (s.compare("P")==0 || s.compare("p")==0 ) return 14;
        else if (s.compare("S")==0 || s.compare("s")==0 ) return 15;
        else if (s.compare("T")==0 || s.compare("t")==0 ) return 16;
        else if (s.compare("W")==0 || s.compare("w")==0 ) return 17;
        else if (s.compare("Y")==0 || s.compare("y")==0 ) return 18;
        else if (s.compare("V")==0 || s.compare("v")==0 ) return 19;
        else if (s.compare("B")==0 || s.compare("b")==0 ) return 20;
        else if (s.compare("J")==0 || s.compare("j")==0 ) return 21;
        else if (s.compare("Z")==0 || s.compare("z")==0 ) return 22;
        else if (s.compare("X")==0 || s.compare("x")==0 ) return 23;
        else if (s.compare("-")==0                      ) return 24;
        else { std::cerr << "Error! - Invalid nucleotide " << s << "\n"; }
    }

    //Function that returns the proper AlphabetAA in string format according to the given Alphabet (nucleotide)
    inline string GetAlphabetAA_String(const Alphabet& s) {
        if      ( s==AmA   ) return "A";
        else if ( s==AmR   ) return "R";
        else if ( s==AmN   ) return "N";
        else if ( s==AmD   ) return "D";
        else if ( s==AmC   ) return "C";
        else if ( s==AmQ   ) return "Q";
        else if ( s==AmE   ) return "E";
        else if ( s==AmG   ) return "G";
        else if ( s==AmH   ) return "H";
        else if ( s==AmI   ) return "I";
        else if ( s==AmL   ) return "L";
        else if ( s==AmK   ) return "K";
        else if ( s==AmM   ) return "M";
        else if ( s==AmF   ) return "F";
        else if ( s==AmP   ) return "P";
        else if ( s==AmS   ) return "S";
        else if ( s==AmT   ) return "T";
        else if ( s==AmW   ) return "W";
        else if ( s==AmY   ) return "Y";
        else if ( s==AmV   ) return "V";
        else if ( s==AmB   ) return "B";
        else if ( s==AmJ   ) return "J";
        else if ( s==AmZ   ) return "Z";
        else if ( s==AmX   ) return "X";
        else if ( s==AmGAP ) return "-";
        else { std::cerr << "Error! - Invalid aminoacid " << s << "\n"; }
    }

     //Function that returns the proper AlphabetAA in string format according to the given Alphabet (nucleotide)
    inline string GetAlphabetAA_String(const int& s) {
        if      ( s==0 ) return "A";
        else if ( s==1 ) return "R";
        else if ( s==2 ) return "N";
        else if ( s==3 ) return "D";
        else if ( s==4 ) return "C";
        else if ( s==5 ) return "Q";
        else if ( s==6 ) return "E";
        else if ( s==7 ) return "G";
        else if ( s==8 ) return "H";
        else if ( s==9 ) return "I";
        else if ( s==10) return "L";
        else if ( s==11) return "K";
        else if ( s==12) return "M";
        else if ( s==13) return "F";
        else if ( s==14) return "P";
        else if ( s==15) return "S";
        else if ( s==16) return "T";
        else if ( s==17) return "W";
        else if ( s==18) return "Y";
        else if ( s==19) return "V";
        else if ( s==20) return "B";
        else if ( s==21) return "J";
        else if ( s==22) return "Z";
        else if ( s==23) return "X";
        else if ( s==24) return "-";
        else { std::cerr << "Error! - Invalid aminoacid " << s << "\n"; }
    }

    //Function that returns the proper Alphabet in string format according to the given integer (overloaded Version)
    inline string GetAlphabet_String(const int& s) {
        if      ( s==0 ) return "A";
        else if ( s==1 ) return "C";
        else if ( s==2 ) return "G";
        else if ( s==3 ) return "T";
        else { std::cerr << "Error! - Invalid nucleotide " << s << "\n"; }
    }

    //Function that returns the proper AlphabetAA according to the given character
    inline AlphabetAA GetAlphabetAA(const char& c) {
        switch (c) {
        case 'A' :
        case 'a' :
            return AmA;
        case 'R' :
        case 'r' :
            return AmR;
        case 'N' :
        case 'n' :
            return AmN;
        case 'D' :
        case 'd' :
            return AmD;
        case 'C' :
        case 'c' :
            return AmC;
        case 'Q' :
        case 'q' :
            return AmQ;
        case 'E' :
        case 'e' :
            return AmE;
        case 'G' :
        case 'g' :
            return AmG;
        case 'H' :
        case 'h' :
            return AmH;
        case 'I' :
        case 'i' :
            return AmI;
        case 'L' :
        case 'l' :
            return AmL;
        case 'K' :
        case 'k' :
            return AmK;
        case 'M' :
        case 'm' :
            return AmM;
        case 'F' :
        case 'f' :
            return AmF;
        case 'P' :
        case 'p' :
            return AmP;
        case 'S' :
        case 's' :
            return AmS;
        case 'T' :
        case 't' :
            return AmT;
        case 'W' :
        case 'w' :
            return AmW;
        case 'Y' :
        case 'y' :
            return AmY;
        case 'V' :
        case 'v' :
            return AmV;
        case 'B' :
        case 'b' :
            return AmB;
        case 'J' :
        case 'j' :
            return AmJ;
        case 'Z' :
        case 'z' :
            return AmZ;
        case 'X' :
        case 'x' :
            return AmX;
        case '-' :
            return AmGAP;
        default :
            std::cout << "Error! - Invalid nucleotide " << c << "\n";
        }
    }

    //Function that returns the proper AlphabetAA according to the given string (overloaded Version)
    inline AlphabetAA GetAlphabetAA(const string& s) {
        if      (s.compare("A")==0 || s.compare("a")==0 ) return AmA;
        else if (s.compare("R")==0 || s.compare("r")==0 ) return AmR;
        else if (s.compare("N")==0 || s.compare("n")==0 ) return AmN;
        else if (s.compare("D")==0 || s.compare("d")==0 ) return AmD;
        else if (s.compare("C")==0 || s.compare("c")==0 ) return AmC;
        else if (s.compare("Q")==0 || s.compare("q")==0 ) return AmQ;
        else if (s.compare("E")==0 || s.compare("e")==0 ) return AmE;
        else if (s.compare("G")==0 || s.compare("g")==0 ) return AmG;
        else if (s.compare("H")==0 || s.compare("h")==0 ) return AmH;
        else if (s.compare("I")==0 || s.compare("i")==0 ) return AmI;
        else if (s.compare("L")==0 || s.compare("l")==0 ) return AmL;
        else if (s.compare("K")==0 || s.compare("k")==0 ) return AmK;
        else if (s.compare("M")==0 || s.compare("m")==0 ) return AmM;
        else if (s.compare("F")==0 || s.compare("f")==0 ) return AmF;
        else if (s.compare("P")==0 || s.compare("p")==0 ) return AmP;
        else if (s.compare("S")==0 || s.compare("s")==0 ) return AmS;
        else if (s.compare("T")==0 || s.compare("t")==0 ) return AmT;
        else if (s.compare("W")==0 || s.compare("w")==0 ) return AmW;
        else if (s.compare("Y")==0 || s.compare("y")==0 ) return AmY;
        else if (s.compare("V")==0 || s.compare("v")==0 ) return AmV;
        else if (s.compare("B")==0 || s.compare("b")==0 ) return AmB;
        else if (s.compare("J")==0 || s.compare("j")==0 ) return AmJ;
        else if (s.compare("Z")==0 || s.compare("z")==0 ) return AmZ;
        else if (s.compare("X")==0 || s.compare("x")==0 ) return AmX;
        else if (s.compare("-")==0                      ) return AmGAP;
        else { std::cerr << "Error! - Invalid nucleotide " << s << "\n"; }
    }

    /** For Alphabet **/
    //Function that returns the proper AlphabetPair according to the given string
    inline AlphabetPair GetAlphabetPair(const string& s) {
        if      (s.compare("AA")==0 || s.compare("aa")==0 ) return AA;
        else if (s.compare("AC")==0 || s.compare("ac")==0 ) return AC;
        else if (s.compare("AG")==0 || s.compare("ag")==0 ) return AG;
        else if (s.compare("AT")==0 || s.compare("at")==0 ) return AT;
        else if (s.compare("CA")==0 || s.compare("ca")==0 ) return CA;
        else if (s.compare("CC")==0 || s.compare("cc")==0 ) return CC;
        else if (s.compare("CG")==0 || s.compare("cg")==0 ) return CG;
        else if (s.compare("CT")==0 || s.compare("ct")==0 ) return CT;
        else if (s.compare("GA")==0 || s.compare("ga")==0 ) return GA;
        else if (s.compare("GC")==0 || s.compare("gc")==0 ) return GC;
        else if (s.compare("GG")==0 || s.compare("gg")==0 ) return GG;
        else if (s.compare("GT")==0 || s.compare("gt")==0 ) return GT;
        else if (s.compare("TA")==0 || s.compare("ta")==0 ) return TA;
        else if (s.compare("TC")==0 || s.compare("tc")==0 ) return TC;
        else if (s.compare("TG")==0 || s.compare("tg")==0 ) return TG;
        else if (s.compare("TT")==0 || s.compare("tt")==0 ) return TT;
        else { std::cerr << "Error! - Invalid nucleotide pair" << s << "\n"; }
    }

    //Function that returns the proper AlphabetPair according to the given strings
    inline AlphabetPair GetAlphabetPair(const string& s1, const string& s2) {
        if      (s1.compare("A")==0 || s1.compare("a")==0 ) {
            if      (s2.compare("A")==0 || s2.compare("a")==0 ) return AA;
            else if (s2.compare("C")==0 || s2.compare("c")==0 ) return AC;
            else if (s2.compare("G")==0 || s2.compare("g")==0 ) return AG;
            else if (s2.compare("T")==0 || s2.compare("t")==0 ) return AT;
            else { std::cerr << "Error! - Invalid nucleotide " << s2 << "\n"; }
        }
        else if (s1.compare("C")==0 || s1.compare("c")==0 ) {
            if      (s2.compare("A")==0 || s2.compare("a")==0 ) return CA;
            else if (s2.compare("C")==0 || s2.compare("c")==0 ) return CC;
            else if (s2.compare("G")==0 || s2.compare("g")==0 ) return CG;
            else if (s2.compare("T")==0 || s2.compare("t")==0 ) return CT;
            else { std::cerr << "Error! - Invalid nucleotide " << s2 << "\n"; }
        }
        else if (s1.compare("G")==0 || s1.compare("g")==0 ) {
            if      (s2.compare("A")==0 || s2.compare("a")==0 ) return GA;
            else if (s2.compare("C")==0 || s2.compare("c")==0 ) return GC;
            else if (s2.compare("G")==0 || s2.compare("g")==0 ) return GG;
            else if (s2.compare("T")==0 || s2.compare("t")==0 ) return GT;
            else { std::cerr << "Error! - Invalid nucleotide " << s2 << "\n"; }
        }
        else if (s1.compare("T")==0 || s1.compare("t")==0 ) {
            if      (s2.compare("A")==0 || s2.compare("a")==0 ) return TA;
            else if (s2.compare("C")==0 || s2.compare("c")==0 ) return TC;
            else if (s2.compare("G")==0 || s2.compare("g")==0 ) return TG;
            else if (s2.compare("T")==0 || s2.compare("t")==0 ) return TT;
            else { std::cerr << "Error! - Invalid nucleotide " << s2 << "\n"; }
        }
        else { std::cerr << "Error! - Invalid nucleotide " << s1 << "\n"; }
    }

    //Function that returns the proper AlphabetPair according to the given characters
	inline AlphabetPair GetAlphabetPair(const char& s1, const char& s2) {
        if      (s1 == 'A' || s1 == 'a' ) {
            if      (s2 == 'A' || s2 == 'a' ) return AA;
            else if (s2 == 'C' || s2 == 'c' ) return AC;
            else if (s2 == 'G' || s2 == 'g' ) return AG;
            else if (s2 == 'T' || s2 == 't' ) return AT;
            else { std::cerr << "Error! - Invalid nucleotide " << s2 << "\n"; }
        }
        else if (s1 == 'C' || s1 == 'c' ) {
            if      (s2 == 'A' || s2 == 'a' ) return CA;
            else if (s2 == 'C' || s2 == 'c' ) return CC;
            else if (s2 == 'G' || s2 == 'g' ) return CG;
            else if (s2 == 'T' || s2 == 't' ) return CT;
            else { std::cerr << "Error! - Invalid nucleotide " << s2 << "\n"; }
        }
        else if (s1 == 'G' || s1 == 'g' ) {
            if      (s2 == 'A' || s2 == 'a' ) return GA;
            else if (s2 == 'C' || s2 == 'c' ) return GC;
            else if (s2 == 'G' || s2 == 'g' ) return GG;
            else if (s2 == 'T' || s2 == 't' ) return GT;
            else { std::cerr << "Error! - Invalid nucleotide " << s2 << "\n"; }
        }
        else if (s1 == 'T' || s1 == 't' ) {
            if      (s2 == 'A' || s2 == 'a' ) return TA;
            else if (s2 == 'C' || s2 == 'c' ) return TC;
            else if (s2 == 'G' || s2 == 'g' ) return TG;
            else if (s2 == 'T' || s2 == 't' ) return TT;
            else { std::cerr << "Error! - Invalid nucleotide " << s2 << "\n"; }
        }
        else { std::cerr << "Error! - Invalid nucleotide " << s1 << "\n"; }
    }

    //Function that returns the proper AlphabetPair according to the given integer and character
	inline AlphabetPair GetAlphabetPair(const int& s1, const char& s2) {
        if      (s1 == 0 ) {
            if      (s2 == 'A' || s2 == 'a' ) return AA;
            else if (s2 == 'C' || s2 == 'c' ) return AC;
            else if (s2 == 'G' || s2 == 'g' ) return AG;
            else if (s2 == 'T' || s2 == 't' ) return AT;
            else { std::cerr << "Error! - Invalid nucleotide " << s2 << "\n"; }
        }
        else if (s1 == 1 ) {
            if      (s2 == 'A' || s2 == 'a' ) return CA;
            else if (s2 == 'C' || s2 == 'c' ) return CC;
            else if (s2 == 'G' || s2 == 'g' ) return CG;
            else if (s2 == 'T' || s2 == 't' ) return CT;
            else { std::cerr << "Error! - Invalid nucleotide " << s2 << "\n"; }
        }
        else if (s1 == 2 ) {
            if      (s2 == 'A' || s2 == 'a' ) return GA;
            else if (s2 == 'C' || s2 == 'c' ) return GC;
            else if (s2 == 'G' || s2 == 'g' ) return GG;
            else if (s2 == 'T' || s2 == 't' ) return GT;
            else { std::cerr << "Error! - Invalid nucleotide " << s2 << "\n"; }
        }
        else if (s1 == 3 ) {
            if      (s2 == 'A' || s2 == 'a' ) return TA;
            else if (s2 == 'C' || s2 == 'c' ) return TC;
            else if (s2 == 'G' || s2 == 'g' ) return TG;
            else if (s2 == 'T' || s2 == 't' ) return TT;
            else { std::cerr << "Error! - Invalid nucleotide " << s2 << "\n"; }
        }
        else { std::cerr << "Error! - Invalid nucleotide " << s1 << "\n"; }
    }

    //Function that returns the proper AlphabetPair according to the given integers
	inline AlphabetPair GetAlphabetPair(const int& s1, const int& s2) {
        if      (s1 == 0) {
            if      (s2 == 0) return AA;
            else if (s2 == 1) return AC;
            else if (s2 == 2) return AG;
            else if (s2 == 3) return AT;
            else { std::cerr << "Error! - Invalid nucleotide " << s2 << "\n"; }
        }
        else if (s1 == 1) {
            if      (s2 == 0) return CA;
            else if (s2 == 1) return CC;
            else if (s2 == 2) return CG;
            else if (s2 == 3) return CT;
            else { std::cerr << "Error! - Invalid nucleotide " << s2 << "\n"; }
        }
        else if (s1 == 2) {
            if      (s2 == 0) return GA;
            else if (s2 == 1) return GC;
            else if (s2 == 2) return GG;
            else if (s2 == 3) return GT;
            else { std::cerr << "Error! - Invalid nucleotide " << s2 << "\n"; }
        }
        else if (s1 == 3) {
            if      (s2 == 0) return TA;
            else if (s2 == 1) return TC;
            else if (s2 == 2) return TG;
            else if (s2 == 3) return TT;
            else { std::cerr << "Error! - Invalid nucleotide " << s2 << "\n"; }
        }
        else { std::cerr << "Error! - Invalid nucleotide " << s1 << "\n"; }
    }

    //Function that returns the proper AlphabetPair according to the given integer (overloaded Version)
    inline AlphabetPair GetAlphabetPair(const int& s) {
        if      ( s==0 ) return AA;
        else if ( s==1 ) return AC;
        else if ( s==2 ) return AG;
        else if ( s==3 ) return AT;
        else if ( s==4 ) return CA;
        else if ( s==5 ) return CC;
        else if ( s==6 ) return CG;
        else if ( s==7 ) return CT;
        else if ( s==8 ) return GA;
        else if ( s==9 ) return GC;
        else if ( s==10) return GG;
        else if ( s==11) return GT;
        else if ( s==12) return TA;
        else if ( s==13) return TC;
        else if ( s==14) return TG;
        else if ( s==15) return TT;
        else { std::cerr << "Error! - Invalid nucleotide pair " << s << "\n"; }
    }

    //Function that returns the proper AlphabetPair in string format according to the given Alphabet (nucleotide)
    inline string GetAlphabetPair_String(const AlphabetPair& s) {
        if      ( s==AA ) return "AA";
        else if ( s==AC ) return "AC";
        else if ( s==AG ) return "AG";
        else if ( s==AT ) return "AT";
        else if ( s==CA ) return "CA";
        else if ( s==CC ) return "CC";
        else if ( s==CG ) return "CG";
        else if ( s==CT ) return "CT";
        else if ( s==GA ) return "GA";
        else if ( s==GC ) return "GC";
        else if ( s==GG ) return "GG";
        else if ( s==GT ) return "GT";
        else if ( s==TA ) return "TA";
        else if ( s==TC ) return "TC";
        else if ( s==TG ) return "TG";
        else if ( s==TT ) return "TT";
        else { std::cerr << "Error! - Invalid nucleotide pair" << s << "\n"; }
    }

    //Function that returns the proper AlphabetPair in string format according to the given integer (overloaded Version)
    inline string GetAlphabetPair_String(const int& s) {
        if      ( s==0 ) return "AA";
        else if ( s==1 ) return "AC";
        else if ( s==2 ) return "AG";
        else if ( s==3 ) return "AT";
        else if ( s==4 ) return "CA";
        else if ( s==5 ) return "CC";
        else if ( s==6 ) return "CG";
        else if ( s==7 ) return "CT";
        else if ( s==8 ) return "GA";
        else if ( s==9 ) return "GC";
        else if ( s==10) return "GG";
        else if ( s==11) return "GT";
        else if ( s==12) return "TA";
        else if ( s==13) return "TC";
        else if ( s==14) return "TG";
        else if ( s==15) return "TT";
        else { std::cerr << "Error! - Invalid nucleotide pair" << s << "\n"; }
    }

    //Function to return the Corresponding AlphabetPair for the 2 given Alphabet
    inline AlphabetPair GetAlphabetPair(const Alphabet& a, const Alphabet& b) {
        if (a==A){
            if      (b==A) return AA;
            else if (b==C) return AC;
            else if (b==G) return AG;
            else if (b==T) return AT;
        }else if (a==C){
            if      (b==A) return CA;
            else if (b==C) return CC;
            else if (b==G) return CG;
            else if (b==T) return CT;
        }else if (a==G){
            if      (b==A) return GA;
            else if (b==C) return GC;
            else if (b==G) return GG;
            else if (b==T) return GT;
        }else if (a==T){
            if      (b==A) return TA;
            else if (b==C) return TC;
            else if (b==G) return TG;
            else if (b==T) return TT;
        }
    }

    /** For EventType 1 **/
    //Function that returns the proper EvenType1 according to the given character
    inline EventType1 GetEventType1(const string& s) {
        if      (s.compare("BB")==0 || s.compare("bb")==0 || s.compare("bB")==0 || s.compare("Bb")==0 ) return BB;
        else if (s.compare("B_")==0 || s.compare("b_")==0 ) return B_;
        else if (s.compare("_B")==0 || s.compare("_b")==0 ) return _B;
        else { std::cerr << "Error! - Invalid event! " << s << "\n"; }
    }

    //Function that returns the proper EventType1 according to the given integer (overloaded Version)
    inline EventType1 GetEventType1(const int& s) {
        if      ( s==0 ) return BB;
        else if ( s==1 ) return B_;
        else if ( s==2 ) return _B;
        else { std::cerr << "Error! - Invalid event! " << s << "\n"; }
    }

    //Function that returns the proper EventType1 in string format according to the given EventType1 (event)
    inline string GetEventType1_String(const EventType1& s) {
        if      ( s==BB ) return "BB";
        else if ( s==B_ ) return "B_";
        else if ( s==_B ) return "_B";
        else { std::cerr << "Error! - Invalid event! " << s << "\n"; }
    }

    //Function that returns the proper EventType1 in string format according to the given integer (overloaded Version)
    inline string GetEventType1_String(const int& s) {
        if      ( s==0 ) return "BB";
        else if ( s==1 ) return "B_";
        else if ( s==2 ) return "_B";
        else { std::cerr << "Error! - Invalid event! " << s << "\n"; }
    }

    /** For EventType2 **/
    //Function that returns the proper EvenType2 according to the given character
    inline EventType2 GetEventType2(const char& c) {
        switch (c) {
            case 'H' :
            case 'h' :
                return H;
            case 'B' :
            case 'b' :
                return B;
            case 'N' :
            case 'n' :
                return N;
            case 'E' :
            case 'e' :
                return E;
            default :
                std::cout << "Error! - Invalid event! " << c << "\n";
        }
    }

    //Function that returns the proper EvenType2 according to the given string (overloaded Version)
    inline EventType2 GetEventType2(const string& s) {
        if      (s.compare("H")==0 || s.compare("h")==0 ) return H;
        else if (s.compare("B")==0 || s.compare("b")==0 ) return B;
        else if (s.compare("N")==0 || s.compare("n")==0 ) return N;
        else if (s.compare("E")==0 || s.compare("e")==0 ) return E;
        else { std::cerr << "Error! - Invalid event! " << s << "\n"; }
    }

    //Function that returns the proper EventType2 according to the given integer (overloaded Version)
    inline EventType2 GetEventType2(const int& s) {
        if      ( s==0 ) return H;
        else if ( s==1 ) return B;
        else if ( s==2 ) return N;
        else if ( s==3 ) return E;
        else { std::cerr << "Error! - Invalid event! " << s << "\n"; }
    }

    //Function that returns the proper EventType2 in string format according to the given EventType2 (event)
    inline string GetEventType2_String(const EventType2& s) {
        if      ( s==H ) return "H";
        else if ( s==B ) return "B";
        else if ( s==N ) return "N";
        else if ( s==E ) return "E";
        else { std::cerr << "Error! - Invalid event! " << s << "\n"; }
    }

    //Function that returns the proper EventType2 in string format according to the given integer (overloaded Version)
    inline string GetEventType2_String(const int& s) {
        if      ( s==0 ) return "H";
        else if ( s==1 ) return "H";
        else if ( s==2 ) return "N";
        else if ( s==3 ) return "E";
        else { std::cerr << "Error! - Invalid event! " << s << "\n"; }
    }
}

#endif // MAPPINGS_H_INCLUDED
