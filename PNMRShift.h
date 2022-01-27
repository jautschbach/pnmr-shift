#ifndef _PNMRSHIFT2_H
#define _PNMRSHIFT2_H

#include <string>
#include <cctype>
#include "XYZFile.h"
#include "SimpleMatrix.h"
#include "ThrowStream.h"

using namespace std;

//Structures and classes

struct gTensor
{
    SimpleMatrixD dia,para,total;
    SimpleMatrixD gemat, deltagisomat, deltagtilde;

    gTensor()
    {
        dia = para = total = gemat = deltagisomat = deltagtilde = SimpleMatrixD(3,3);
    }
};

struct ATensor
{
    string el;
    unsigned long id;

    SimpleMatrixD fcsd,psoso,total;

    //nr=nonrelativistic, so=spinorbit, fc=fermi contact, sd=spin dipole, as=asymmetric
    SimpleMatrixD nrfc, nrsd, sofc, sosd, as;

    ATensor()
    {
        fcsd = psoso = total = nrfc = nrsd = sofc = sosd = as = SimpleMatrixD(3,3);
        el = "xx";
        id = 0;
    }
};


struct OrbShieldTensor
{
    string el;
    unsigned long id;

    SimpleMatrixD dia,para,total;

    OrbShieldTensor()
    {
        dia = para = total = SimpleMatrixD(3,3);
        el = "xx";
        id = 0;
    }
};

typedef SimpleMatrixD DTensor;


//Structure to hold different contributions
struct Contributions
{
    SimpleMatrixD fermiContact, pseudoContact; //Moon
    SimpleMatrixD detailed[15]; //hrobarik, etc;

    Contributions(const Contributions & c)
    {
        fermiContact = c.fermiContact;
        pseudoContact = c.pseudoContact;

        for(int i = 0; i < 15; i++)
            detailed[i] = c.detailed[i];
    }
    Contributions()
    {
        fermiContact = pseudoContact = SimpleMatrixD(3,3);
        for(int i = 0; i < 15; i++)
            detailed[i] = SimpleMatrixD(3,3);

    }
};

//Isotropic value of a tensor
template<typename T>
T IsotropicValue(const SimpleMatrix<T> & sm)
{
    T avg = 0;
    if(sm.nrows() != sm.ncols())
        THROWSTREAM << "Can't get an isotropic value of a " << sm.nrows() << "x" << sm.ncols() << " matrix!";
    for(unsigned long i = 0; i < sm.nrows(); i++)
        avg += sm(i,i);
    return avg/sm.nrows();
}


//Convert a string to lower case
inline string StringToLowerCopy(const string & str)
{
    string tmp;
    for (size_t i = 0; i < str.length(); i++)
        tmp.append(1, tolower(str[i]));
    return tmp;
}

inline void StringToLower(string & str)
{
    str = StringToLowerCopy(str);
}

inline string ElementCaseCopy(const string & el)
{
    string s;

    if (el.length() == 0)
        THROWSTREAM << "Trying to convert an empty element string!";

    s.append(1, toupper(el[0]));
    for (unsigned int i = 1; i < el.length(); i++)
        s.append(1, tolower(el[i]));

    return s;
}

inline void ElementCase(string & el)
{
    el = ElementCaseCopy(el);
}

inline void Trim(string & str)
{
    if(str.length() == 0)
        return;
    size_t begin = str.find_first_not_of(" \t");
    size_t end = str.find_last_not_of(" \t");
    if(begin == str.npos || end == str.npos)
        str = "";
    else if(end < begin)
        THROWSTREAM << "Trimming error - begin is " << begin << " and end is " << end;
    else
        str = str.substr(begin,end-begin+1);
}

//Function prototypes
//Output.cpp
void PrintTensor(const string & title, const SimpleMatrixD & mat);
void PrintHyperfine(const string & title, vector<ATensor> atensors);
void PrintOrbitalShield(const string & title, vector<OrbShieldTensor> orbtensors);
void PrintCoords(const XYZFile & coords);
void PrintResults_PNMRShift(const vector<Contributions> & cont,
                  const vector<ATensor> atensors,
                  const vector<OrbShieldTensor> & orbshield, const XYZFile & xyz,
                  const vector<string> & avgs, bool detail);
//Help.cpp
void PrintHelp_PNMRShift(void);
void PrintHelp_gFind(void);

//ConvertSelectionString.cpp
vector<long> ConvertSelectionString(const string & str, bool allowdupe);

//ReadFile.cpp
void ReadFile_PNMRShift(string filename, gTensor & gtens, vector<ATensor> & atensvec,
              vector<OrbShieldTensor> & orbtensvec, DTensor & dtens, vector<string> & averaging);

// simpler input version
void ReadFile_PNMRShift_simple(string filename, gTensor & gtens, vector<ATensor> & atensvec,
                               vector<OrbShieldTensor> & orbtensvec, DTensor & dtens, vector<string> & averaging,
                               bool fullg);
void ReadFile_gFind(string filename, SimpleMatrixD  & gtens, vector<double> & shifts, XYZFile & xyz);


#endif
// Copyright 2012-2015 Ben Pritchard, Bob Martin, and Jochen Autschbach
// This file is part of PNMRShift.
//
// PNMRShift is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
//  PNMRShift is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with PNMRShift.  If not, see <http://www.gnu.org/licenses/>.
