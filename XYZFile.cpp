#include "XYZFile.h"
#include "PNMRShift.h"

AtomCoord::AtomCoord(const string & line)
{
    id = 0;

    stringstream ss(line);
    ss.exceptions(stringstream::failbit | stringstream::badbit);

    try
    {
        ss >> el >> x >> y >> z;
        if(!ss.eof())
            THROWSTREAM << "Error - extra stuff on line.";
    }
    catch(std::exception & ex)
    {
        THROWSTREAM << "Error parsing XYZ file line.";
    }

    ElementCase(el);

}

void XYZFile::ReadXYZFile(const string & filename)
{
    ifstream file(filename.c_str());

    if(!file.is_open())
        THROWSTREAM << "Unable to open file " << filename;

    string line;
    getline(file, line);

    stringstream ss(line);
    ss >> atomcount;
    if(ss.fail() | ss.bad())
        THROWSTREAM << "Error reading number of atoms";

    getline(file, title);

    unsigned long linenumber = 3;

    try
    {
        while(getline(file, line).good())
        {
            Trim(line);
            if(line.length() > 0)
            {
                atoms.push_back(line);
                atoms[atoms.size()-1].id = linenumber-2;
            }
            linenumber++;
        }
    }
    catch(std::exception & ex)
    {
        file.close(); // needed?
        THROWSTREAMAPPEND(ex) << "Error on line " << linenumber << " of xyz file " << filename;
    }

    if(!file.eof())
        THROWSTREAM << "Error reading file " << filename << " - expecting EOF";

    if(atoms.size() != atomcount)
        THROWSTREAM << "File " << filename << " is supposed to have " << atomcount
                    << " atoms, but I only count " << atoms.size();
}

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
