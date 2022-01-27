#ifndef CONSTANTS_H_
#define CONSTANTS_H_

#define CONSTANT_GE 2.00231930436163l
#define CONSTANT_BE 927.400968e-26
#define CONSTANT_BN 5.05078353e-27
#define CONSTANT_H 6.62606957e-34
#define CONSTANT_K 1.3806488e-23
#define CONSTANT_C 299792458
#define CONSTANT_PI 3.14159265358979323846l
#define CONSTANT_MU0 4*CONSTANT_PI*1E-7

#include <string>
#include <climits>

using namespace std;

inline double NuclearG(const string & str)
{
    string s = StringToLowerCopy(str);
    if (s == "h")
        return 5.5857;
    else if (s == "he")
        return -4.25524;
    else if (s == "li")
        return 2.17096;
    else if (s == "be")
        return -0.785066667;
    else if (s == "b")
        return 1.7924;
    else if (s == "c")
        return 1.40482;
    else if (s == "n")
        return 0.40376;
    else if (s == "o")
        return -0.75752;
    else if (s == "f")
        return 5.25774;
    else if (s == "ne")
        return -0.4412;
    else if (s == "na")
        return 1.4783466667;
    else if (s == "mg")
        return -0.34218;
    else if (s == "al")
        return 1.456604;
    else if (s == "si")
        return -1.1106;
    else if (s == "p")
        return 2.2632;
    else if (s == "s")
        return 0.429213333;
    else if (s == "cl")
        return 0.547913333;
    else if (s == "ar")
        return -0.4542857143;
    else if (s == "k")
        return 0.260973333;
    else if (s == "ca")
        return -0.3763714286;
    else if (s == "sc")
        return 1.35899714286;
    else if (s == "ti")
        return -0.315392;
    else if (s == "v")
        return 1.47105885714;
    else if (s == "cr")
        return -0.31636;
    else if (s == "mn")
        return 1.38748;
    else if (s == "fe")
        return 0.1812;
    else if (s == "co")
        return 1.32285714286;
    else if (s == "ni")
        return -0.500013333;
    else if (s == "cu")
        return 1.4822;
    else if (s == "zn")
        return 0.3502;
    else if (s == "ga")
        return 1.344393333;
    else if (s == "ge")
        return -0.1954371111;
    else if (s == "as")
        return 0.959646666;
    else if (s == "se")
        return 1.07012;
    else if (s == "br")
        return 1.404266667;
    else if (s == "kr")
        return -0.2157108889;
    else if (s == "rb")
        return 0.5412;
    else if (s == "sr")
        return -0.24288889;
    else if (s == "y")
        return -0.27484;
    else if (s == "zr")
        return -0.521448;
    else if (s == "nb")
        return 1.371222222;
    else if (s == "mo")
        return -0.36568;
    else if (s == "tc")
        return 1.263266667;
    else if (s == "ru")
        return -0.28752;
    else if (s == "rh")
        return -0.1768;
    else if (s == "pd")
        return -0.2568;
    else if (s == "ag")
        return -0.22714;
    else if (s == "cd")
        return -1.1897722;
    else if (s == "in")
        return 1.23133333;
    else if (s == "sn")
        return -2.0946;
    else if (s == "sb")
        return 1.3452;
    else if (s == "te")
        return -1.777;
    else if (s == "i")
        return 1.12532;
    else if (s == "xe")
        return -1.556;
    else if (s == "cs")
        return 0.73771428571;
    else if (s == "ba")
        return 0.6249333333;
    else if (s == "la")
        return 0.79514285714;
    else if (s == "ce")
        return 0.31428571429;
    else if (s == "pr")
        return 1.71;
    else if (s == "nd")
        return -0.3057142857;
    else if (s == "pm")
        return 0.74285714286;
    else if (s == "sm")
        return -0.2328571429;
    else if (s == "eu")
        return 0.6132;
    else if (s == "gd")
        return -0.226666667;
    else if (s == "tb")
        return 1.342666667;
    else if (s == "dy")
        return 0.2692;
    else if (s == "ho")
        return 1.19142857143;
    else if (s == "er")
        return -0.1611142857;
    else if (s == "tm")
        return -0.464;
    else if (s == "yb")
        return -0.271956;
    else if (s == "lu")
        return 0.63791428571;
    else if (s == "hf")
        return 0.22671428571;
    else if (s == "ta")
        return 0.67714285714;
    else if (s == "w")
        return 0.2355696;
    else if (s == "re")
        return 1.28788;
    else if (s == "os")
        return 0.439953333;
    else if (s == "ir")
        return 0.10933333;
    else if (s == "pt")
        return 1.219;
    else if (s == "au")
        return 0.097166667;
    else if (s == "hg")
        return 1.01177;
    else if (s == "tl")
        return 3.27643;
    else if (s == "pb")
        return 1.18516;
    else if (s == "bi")
        return 0.913555556;
    else if (s == "po")
        return 1.54;
    else if (s == "at")
        return 0;
    else if (s == "rn")
        return 1.2;
    else if (s == "fr")
        return 0.78;
    else if (s == "ra")
        return 0.180666667;
    else if (s == "ac")
        return 0.733333333;
    else if (s == "th")
        return 0.184;
    else if (s == "pa")
        return 1.34;
    else if (s == "u")
        return -0.1085714286;
    else if (s == "np")
        return 1.256;
    else if (s == "pu")
        return 0.406;
    else if (s == "am")
        return 0.6;
    else if (s == "cm")
        return 0.082222222;
    else if (s == "bk")
        return 0.5714285714;
    else if (s == "es")
        return 1.1714285714;
    else if(s == "xx") //dummy atom
        return 0; //silently
    else
    {
        std::cout << "=== Warning - I don't have a nuclear g factor for element \"" << s << "\"\n";
        return 0;
    }
}

#endif /* CONSTANTS_H_ */
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
