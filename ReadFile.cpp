
#include <fstream>

#include "PNMRShift.h"
#include "Constants.h"


// Note: Will not clear or zero std::vectors or matrices!
void ReadFile_PNMRShift(string filename, gTensor & gtens, vector<ATensor> & atensvec,
                        vector<OrbShieldTensor> & orbtensvec, DTensor & dtens, vector<string> & averaging)
{
    try
    {
        ifstream file(filename.c_str());

        if(!file.is_open())
            THROWSTREAM << "Unable to open file " << filename;

        string wholeline,line;
        unsigned long lineno = 0;

        while(getline(file, line).good())
        {
            lineno++;

            Trim(line);
            if(line.length() == 0)
                continue;

            StringToLower(line);

            wholeline = line;
            stringstream ss(line);
            line.clear();
            ss >> line;

            //cout << file.good() << file.bad() << file.eof() << file.fail() << " "
            //     << ss.good() << ss.bad() << ss.eof() << ss.fail()
            //     << " READ FILE: Line " << lineno << " line is \"" << line << "\"\n";


            if(line == "gtensor")
            {
                if(!ss.eof())
                    THROWSTREAM << "Extra characters on line " << lineno << " after 'gtensor'";

                for(unsigned int i = 0; i < 3; i++)
                {
                    getline(file, line);
                    lineno++;
                    Trim(line);

                    stringstream ssgtens(line);
                    ssgtens >> gtens.dia(i,0) >> gtens.dia(i,1) >> gtens.dia(i,2);

                    if(ssgtens.bad() || ssgtens.fail() || !ssgtens.eof())
                        THROWSTREAM << "Extra characters or parsing error on line " << lineno
                                    << " - good,bad,eof,fail = "
                                    << ssgtens.good() << ssgtens.bad() << ssgtens.eof() << ssgtens.fail();

                }
                for(unsigned int i = 0; i < 3; i++)
                {
                    getline(file, line);
                    lineno++;
                    Trim(line);

                    stringstream ssgtens(line);
                    ssgtens >> gtens.para(i,0) >> gtens.para(i,1) >> gtens.para(i,2);

                    if(ssgtens.bad() || ssgtens.fail() || !ssgtens.eof())
                        THROWSTREAM << "Extra characters or parsing error on line " << lineno
                                    << " - good,bad,eof,fail = "
                                    << ssgtens.good() << ssgtens.bad() << ssgtens.eof() << ssgtens.fail();
                }

                gtens.dia *= 0.001;
                gtens.para *= 0.001;
                gtens.total = gtens.dia + gtens.para;

                gtens.gemat = CONSTANT_GE * IdentityMatrix<double>(3);

                gtens.deltagisomat = IdentityMatrix<double>(3) * IsotropicValue(gtens.total);

                gtens.deltagtilde = gtens.total - gtens.deltagisomat;
                gtens.total += gtens.gemat;
            }

            else if(line == "zfstensor")
            {
                if(!ss.eof())
                    THROWSTREAM << "Extra characters on line " << lineno << " after 'zfstensor'";

                for(unsigned int i = 0; i < 3; i++)
                {
                    getline(file, line);
                    lineno++;
                    Trim(line);

                    stringstream ssdtens(line);
                    ssdtens >> dtens(i,0) >> dtens(i,1) >> dtens(i,2);

                    if(ssdtens.bad() || ssdtens.fail() || !ssdtens.eof())
                        THROWSTREAM << "Extra characters or parsing error on line " << lineno
                                    << " - good,bad,eof,fail = "
                                    << ssdtens.good() << ssdtens.bad() << ssdtens.eof() << ssdtens.fail();
                }
            }

            else if(line == "atensor")
            {
                ATensor atens;
                ss >> atens.id >> atens.el;

                if(ss.fail() || ss.bad())
                    THROWSTREAM << "Error reading atensor atom number and element line " << lineno;

                if(!ss.eof())
                    THROWSTREAM << "Extra characters on line " << lineno << " after '" << atens.el << "'";


                ElementCase(atens.el);

                for(unsigned int i = 0; i < 3; i++)
                {
                    getline(file, line);
                    lineno++;
                    Trim(line);

                    stringstream ssatens(line);
                    ssatens >> atens.fcsd(i,0) >> atens.fcsd(i,1) >> atens.fcsd(i,2);

                    if(ssatens.bad() || ssatens.fail() || !ssatens.eof())
                        THROWSTREAM << "Extra characters or parsing error on line " << lineno
                                    << " - good,bad,eof,fail = "
                                    << ssatens.good() << ssatens.bad() << ssatens.eof() << ssatens.fail();
                }
                for(unsigned int i = 0; i < 3; i++)
                {
                    getline(file, line);
                    lineno++;
                    Trim(line);

                    stringstream ssatens(line);
                    ssatens >> atens.psoso(i,0) >> atens.psoso(i,1) >> atens.psoso(i,2);

                    if(ssatens.bad() || ssatens.fail() || !ssatens.eof())
                        THROWSTREAM << "Extra characters or parsing error on line " << lineno
                                    << " - good,bad,eof,fail = "
                                    << ssatens.good() << ssatens.bad() << ssatens.eof() << ssatens.fail();
                }
                atens.total = atens.fcsd + atens.psoso;

                //decompose further
                atens.nrfc = IdentityMatrix<double>(3) * IsotropicValue(atens.fcsd);
                atens.nrsd = atens.fcsd - atens.nrfc;

                // Symmetric part
                SimpleMatrixD sym = 0.5 * (atens.psoso + Transpose(atens.psoso));
                atens.sofc = IdentityMatrix<double>(3) * IsotropicValue(sym);
                atens.sosd = sym - atens.sofc;

                // Asymmetric part
                atens.as = 0.5 * (atens.psoso - Transpose(atens.psoso));


                // Replace the vector element
                if(atensvec.size() < atens.id)
                    THROWSTREAM << "ATensor ID " << atens.id << " is out of range of the coordinates (" << atensvec.size() << " coords)";
                if(atens.id == 0)
                    THROWSTREAM << "How do I have an A-Tensor for atom zero? They start at 1 you know...";

                atensvec[atens.id - 1] = atens;
            }
            else if(line == "orbtensor")
            {
                OrbShieldTensor orb;
                ss >> orb.id >> orb.el;

                if(ss.fail() || ss.bad())
                    THROWSTREAM << "Error reading orbtensor atom number and element line " << lineno;

                if(!ss.eof())
                    THROWSTREAM << "Extra characters on line " << lineno << " after '" << orb.el << "'";

                ElementCase(orb.el);

                for(unsigned int i = 0; i < 3; i++)
                {
                    getline(file, line);
                    lineno++;
                    Trim(line);

                    stringstream ssotens(line);
                    ssotens >> orb.dia(i,0) >> orb.dia(i,1) >> orb.dia(i,2);

                    if(ssotens.bad() || ssotens.fail() || !ssotens.eof())
                        THROWSTREAM << "Extra characters or parsing error on line " << lineno
                                    << " - good,bad,eof,fail = "
                                    << ssotens.good() << ssotens.bad() << ssotens.eof() << ssotens.fail();
                }
                for(unsigned int i = 0; i < 3; i++)
                {
                    getline(file, line);
                    lineno++;
                    Trim(line);

                    stringstream ssotens(line);
                    ssotens >> orb.para(i,0) >> orb.para(i,1) >> orb.para(i,2);

                    if(ssotens.bad() || ssotens.fail() || !ssotens.eof())
                        THROWSTREAM << "Extra characters or parsing error on line " << lineno
                                    << " - good,bad,eof,fail = "
                                    << ssotens.good() << ssotens.bad() << ssotens.eof() << ssotens.fail();
                }

                orb.total = orb.dia + orb.para;


                if(orbtensvec.size() < orb.id)
                    THROWSTREAM << "Orb Tensor ID " << orb.id << " is out of range of the coordinates (" << orbtensvec.size() << " coords)";
                if(orb.id == 0)
                    THROWSTREAM << "How do I have an Orbital shielding tensor for atom zero? They start at 1 you know...";

                orbtensvec[orb.id-1] = orb;
            }
            else if(line == "averaging")
            {
                if(!ss.eof())
                    THROWSTREAM << "Extra characters on line " << lineno << " after 'averaging'";


                while(getline(file, line).good())
                {
                    lineno++;
                    Trim(line);
                    if(line.length() == 0)
                        break;

                    averaging.push_back(line);
                }

            }
            else if(!(line.size() > 0 && line[0] == '#')) //line is not a comment
                THROWSTREAM << "Unknown command. What does \"" << wholeline << "\" mean?";

            //averaging will read past the last line given for that section
            // so this will prevent crashing and burning as getline tries to read
            // past the end of the file
            if(!file.good() && file.eof())
                break;

        }

        if(!file.eof())
            THROWSTREAM << "Error reading file " << filename << " - state (good,bad,eof,fail) = ("
                        << file.good() << file.bad() << file.eof() << file.fail() << ")";

        file.close();
    }
    catch(std::exception & ex)
    {
        THROWSTREAMAPPEND(ex) << "Error parsing file \"" << filename << "\"";
    }

}

// Note: Will not clear or zero std::vectors or matrices!
void ReadFile_PNMRShift_simple(string filename, gTensor & gtens, vector<ATensor> & atensvec,
                               vector<OrbShieldTensor> & orbtensvec, DTensor & dtens, vector<string> & averaging,
                               bool fullg)
{
    try
    {
        ifstream file(filename.c_str());

        if(!file.is_open())
            THROWSTREAM << "Unable to open file " << filename;

        string wholeline,line;
        unsigned long lineno = 0;

        while(getline(file, line).good())
        {
            lineno++;

            Trim(line);
            if(line.length() == 0)
                continue;

            StringToLower(line);

            wholeline = line;
            stringstream ss(line);
            line.clear();
            ss >> line;

            //cout << file.good() << file.bad() << file.eof() << file.fail() << " "
            //     << ss.good() << ss.bad() << ss.eof() << ss.fail()
            //     << " READ FILE: Line " << lineno << " line is \"" << line << "\"\n";


            if(line == "gtensor")
            {
                if(!ss.eof())
                    THROWSTREAM << "Extra characters on line " << lineno << " after 'gtensor'";

                for(unsigned int i = 0; i < 3; i++)
                {
                    getline(file, line);
                    lineno++;
                    Trim(line);

                    stringstream ssgtens(line);
                    ssgtens >> gtens.total(i,0) >> gtens.total(i,1) >> gtens.total(i,2);

                    if(ssgtens.bad() || ssgtens.fail() || !ssgtens.eof())
                        THROWSTREAM << "Extra characters or parsing error on line " << lineno
                                    << " - good,bad,eof,fail = "
                                    << ssgtens.good() << ssgtens.bad() << ssgtens.eof() << ssgtens.fail();

                }

                //clear these
                gtens.dia *= 0.0;
                gtens.para *= 0.0;

                gtens.gemat = CONSTANT_GE * IdentityMatrix<double>(3);

                if(!fullg)
                  gtens.total = 0.001*gtens.total + gtens.gemat;
    
                gtens.deltagisomat = (IdentityMatrix<double>(3) * IsotropicValue(gtens.total - gtens.gemat));
                gtens.deltagtilde = gtens.total - gtens.deltagisomat - gtens.gemat;
            }

            else if(line == "zfstensor")
            {
                if(!ss.eof())
                    THROWSTREAM << "Extra characters on line " << lineno << " after 'zfstensor'";

                for(unsigned int i = 0; i < 3; i++)
                {
                    getline(file, line);
                    lineno++;
                    Trim(line);

                    stringstream ssdtens(line);
                    ssdtens >> dtens(i,0) >> dtens(i,1) >> dtens(i,2);

                    if(ssdtens.bad() || ssdtens.fail() || !ssdtens.eof())
                        THROWSTREAM << "Extra characters or parsing error on line " << lineno
                                    << " - good,bad,eof,fail = "
                                    << ssdtens.good() << ssdtens.bad() << ssdtens.eof() << ssdtens.fail();
                }
            }

            else if(line == "atensor")
            {
                ATensor atens;
                ss >> atens.id >> atens.el;

                if(ss.fail() || ss.bad())
                    THROWSTREAM << "Error reading atensor atom number and element line " << lineno;

                if(!ss.eof())
                    THROWSTREAM << "Extra characters on line " << lineno << " after '" << atens.el << "'";


                ElementCase(atens.el);

                // store the tensor in fscd, ignore psoso
                for(unsigned int i = 0; i < 3; i++)
                {
                    getline(file, line);
                    lineno++;
                    Trim(line);

                    stringstream ssatens(line);
                    ssatens >> atens.fcsd(i,0) >> atens.fcsd(i,1) >> atens.fcsd(i,2);

                    if(ssatens.bad() || ssatens.fail() || !ssatens.eof())
                        THROWSTREAM << "Extra characters or parsing error on line " << lineno
                                    << " - good,bad,eof,fail = "
                                    << ssatens.good() << ssatens.bad() << ssatens.eof() << ssatens.fail();
                }

                // clear these matrices
                atens.psoso *= 0.0;
                atens.sofc *= 0.0;
                atens.sosd *= 0.0;
                atens.as *= 0.0;

                atens.total = atens.fcsd;

                //decompose further
                atens.nrfc = IdentityMatrix<double>(3) * IsotropicValue(atens.fcsd);
                atens.nrsd = atens.fcsd - atens.nrfc;

                // Replace the vector element
                if(atensvec.size() < atens.id)
                    THROWSTREAM << "ATensor ID " << atens.id << " is out of range of the coordinates (" << atensvec.size() << " coords)";
                if(atens.id == 0)
                    THROWSTREAM << "How do I have an A-Tensor for atom zero? They start at 1 you know...";

                atensvec[atens.id - 1] = atens;
            }
            else if(line == "orbtensor")
            {
                OrbShieldTensor orb;
                ss >> orb.id >> orb.el;

                if(ss.fail() || ss.bad())
                    THROWSTREAM << "Error reading orbtensor atom number and element line " << lineno;

                if(!ss.eof())
                    THROWSTREAM << "Extra characters on line " << lineno << " after '" << orb.el << "'";

                ElementCase(orb.el);

                for(unsigned int i = 0; i < 3; i++)
                {
                    getline(file, line);
                    lineno++;
                    Trim(line);

                    stringstream ssotens(line);
                    ssotens >> orb.total(i,0) >> orb.total(i,1) >> orb.total(i,2);

                    if(ssotens.bad() || ssotens.fail() || !ssotens.eof())
                        THROWSTREAM << "Extra characters or parsing error on line " << lineno
                                    << " - good,bad,eof,fail = "
                                    << ssotens.good() << ssotens.bad() << ssotens.eof() << ssotens.fail();
                }

                orb.dia *= 0.0;
                orb.para *= 0.0;

                if(orbtensvec.size() < orb.id)
                    THROWSTREAM << "Orb Tensor ID " << orb.id << " is out of range of the coordinates (" << orbtensvec.size() << " coords)";
                if(orb.id == 0)
                    THROWSTREAM << "How do I have an Orbital shielding tensor for atom zero? They start at 1 you know...";

                orbtensvec[orb.id-1] = orb;
            }
            else if(line == "averaging")
            {
                if(!ss.eof())
                    THROWSTREAM << "Extra characters on line " << lineno << " after 'averaging'";


                while(getline(file, line).good())
                {
                    lineno++;
                    Trim(line);
                    if(line.length() == 0)
                        break;

                    averaging.push_back(line);
                }

            }
            else if(!(line.size() > 0 && line[0] == '#')) //line is not a comment
                THROWSTREAM << "Unknown command. What does \"" << wholeline << "\" mean?";

            //averaging will read past the last line given for that section
            // so this will prevent crashing and burning as getline tries to read
            // past the end of the file
            if(!file.good() && file.eof())
                break;

        }

        if(!file.eof())
            THROWSTREAM << "Error reading file " << filename << " - state (good,bad,eof,fail) = ("
                        << file.good() << file.bad() << file.eof() << file.fail() << ")";

        file.close();
    }
    catch(std::exception & ex)
    {
        THROWSTREAMAPPEND(ex) << "Error parsing file \"" << filename << "\"";
    }

}

void ReadFile_gFind(string filename, SimpleMatrixD  & gtens, vector<double> & shifts, XYZFile & xyz)
{
    try
    {
        gtens = SimpleMatrixD(3,3);
        XYZFile xyzshifts;

        ifstream file(filename.c_str());

        if(!file.is_open())
            THROWSTREAM << "Unable to open file " << filename;

        string wholeline,line;
        unsigned long lineno = 0;

        while(getline(file, line).good())
        {
            lineno++;

            Trim(line);
            if(line.length() == 0)
                continue;

            StringToLower(line);

            wholeline = line;
            stringstream ss(line);
            line.clear();
            ss >> line;

            //cout << file.good() << file.bad() << file.eof() << file.fail() << " "
            //     << ss.good() << ss.bad() << ss.eof() << ss.fail()
            //     << " READ FILE: Line " << lineno << " line is \"" << line << "\"\n";


            if(line == "gguess")
            {
                if(!ss.eof())
                    THROWSTREAM << "Extra characters on line " << lineno << " after 'gtensor'";

                getline(file, line);
                lineno++;
                Trim(line);

                stringstream ssgtens(line);
                ssgtens >> gtens(0,0) >> gtens(1,1) >> gtens(2,2);

                if(ssgtens.bad() || ssgtens.fail() || !ssgtens.eof())
                    THROWSTREAM << "Extra characters or parsing error on line " << lineno
                                << " - good,bad,eof,fail = "
                                << ssgtens.good() << ssgtens.bad() << ssgtens.eof() << ssgtens.fail();

            }
            else if(line == "shifts")
            {
                if(!ss.eof())
                    THROWSTREAM << "Extra characters on line " << lineno << " after 'shifts'";


                while(getline(file, line).good())
                {
                    lineno++;
                    Trim(line);
                    if(line.length() == 0)
                        break;

                    stringstream sssh(line);
                    double shift = 0;
                    unsigned long id = 0;
                    string el;
                    sssh >> id >> el >> shift;

                    if(sssh.bad() || sssh.fail() || !sssh.eof())
                        THROWSTREAM << "Extra characters or parsing error on line " << lineno
                                    << " - good,bad,eof,fail = "
                                    << sssh.good() << sssh.bad() << sssh.eof() << sssh.fail();

                    ElementCase(el);

                    //Get the info from the XYZ file

                    if(xyz.atoms.size() < id || id == 0)
                        THROWSTREAM << "Atom id " << id << " for FC guess is out of range of "
                                    << xyz.atoms.size() << " atoms in the XYZ file.";

                    if(xyz.atoms[id-1].el != el)
                        THROWSTREAM << "Shift for atom id " << id << " says it should be element '" << el
                                    << "' but the XYZ file says it should be '" << xyz.atoms[id-1].el << "'";

                    xyzshifts.atoms.push_back(xyz.atoms[id-1]);
                    shifts.push_back(shift);

                }
            }
            else if(!(line.size() > 0 && line[0] == '#')) //line is not a comment
                THROWSTREAM << "Unknown command. What does \"" << wholeline << "\" mean?";

            //some sections will read past the last line given for that section
            // so this will prevent crashing and burning as getline tries to read
            // past the end of the file
            if(!file.good() && file.eof())
                break;

        }

        if(!file.eof())
            THROWSTREAM << "Error reading file " << filename << " - state (good,bad,eof,fail) = ("
                        << file.good() << file.bad() << file.eof() << file.fail() << ")";

        file.close();

        xyzshifts.atomcount = xyzshifts.atoms.size();
        xyz = xyzshifts;

    }
    catch(std::exception & ex)
    {
        THROWSTREAMAPPEND(ex) << "Error parsing file \"" << filename << "\"";
    }
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
