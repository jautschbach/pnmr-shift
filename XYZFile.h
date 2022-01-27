#ifndef __IG_BPLIB_XYZFILE_H__
#define __IG_BPLIB_XYZFILE_H__

#include <string>
#include <vector>
#include <fstream>

#include "ThrowStream.h"

using namespace std;

struct AtomCoord
{
    string el;
    double x,y,z;
    unsigned long id;

    AtomCoord(const string & line);
    AtomCoord()
    {
        x = y = z = id = 0;
    }

    AtomCoord(const string & El, unsigned long Id, double X, double Y, double Z)
    {
        x = X;
        y = Y;
        z = Z;
        id = Id;
        el = El;
    }
};

struct XYZFile
{
    unsigned long atomcount;
    string title;
    vector<AtomCoord> atoms;

    void ReadXYZFile(const string & filename);

    XYZFile(const string & filename)
    {
        ReadXYZFile(filename);
    }

    XYZFile() { }
};

#endif /* XYZFILE_H_ */
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
