#include <iostream>
#include <cmath>
#include <complex>

#include "PNMRShift.h"
#include "Constants.h"

using namespace std;



string GetNextArg(int & i, int argc, char ** argv)
{
    i++;
    if(i >= argc)
        THROWSTREAM << "This is the last argument; I was expecting more...";

    return argv[i];
}

vector<SimpleMatrixD> CalculateATensor(const XYZFile & oldcoords, const XYZFile & coords, SimpleMatrixD gtens, unsigned long centerAtom)
{
    vector<SimpleMatrixD> atensvec;
    atensvec.resize(coords.atomcount);

    //check for the center atom
    if(centerAtom <= 0 || centerAtom > oldcoords.atomcount)
        THROWSTREAM << "Error - center atom out of range. It must be in the "
                    << "range [1," << oldcoords.atomcount << "]. I got " << centerAtom;

    //Get the coordinates of the center atom, converting to meters from angstroms
    double centerx = oldcoords.atoms[centerAtom-1].x * 1E-10;
    double centery = oldcoords.atoms[centerAtom-1].y * 1E-10;
    double centerz = oldcoords.atoms[centerAtom-1].z * 1E-10;

    // this is missing GN and the powers of r
    // this also include a conversion from J to MHz
    double prefac = (1E-6 * CONSTANT_MU0 * CONSTANT_BE * CONSTANT_BN) / (4*CONSTANT_PI * CONSTANT_H);

    for(unsigned long i = 0; i < coords.atomcount; i++)
    {
        //get an atom coordinates as a vector, converting to meters
        SimpleMatrixD rvec(3,1);
        rvec(0,0) = 1E-10*coords.atoms[i].x - centerx;
        rvec(1,0) = 1E-10*coords.atoms[i].y - centery;
        rvec(2,0) = 1E-10*coords.atoms[i].z - centerz;

        double r = sqrt(pow(rvec(0,0),2) + pow(rvec(1,0),2) + pow(rvec(2,0),2));
        SimpleMatrixD r2I = pow(r,2)*IdentityMatrix<double>(3);
        double r5 = pow(r,5);

        double prefac2 = prefac * NuclearG(coords.atoms[i].el) / r5;
        //rvec is a column vector, stored as a matrix
        // so we can do an outer product this way
        // (3x1) * (1x3) = (3x3)
        atensvec[i] = prefac2 * ((3*rvec*Transpose(rvec)) - r2I) * gtens;
    }

    return atensvec;


}


vector<double> CalculateShifts(const SimpleMatrixD & gtens, const XYZFile & oldcoords, const XYZFile & coords, unsigned long centerAtom, double temperature, double spin)
{
    vector<double> results;
    results.resize(coords.atomcount);

    vector<SimpleMatrixD> atensvec = CalculateATensor(oldcoords, coords, gtens, centerAtom);
    //PrintHyperfine("Hyperfine tensor guess", atensvec);

    //We put a positive sign here to represent these in terms of shifts
    double prefac = (+1.0e12 * CONSTANT_H * CONSTANT_BE * spin * (spin+1.0) )/
                    (3.0*CONSTANT_BN * CONSTANT_K * temperature);

    for(size_t i = 0; i < atensvec.size(); i++)
    {
        double prefac2 = prefac / NuclearG(coords.atoms[i].el);

        SimpleMatrixD afc = IsotropicValue(atensvec[i])*IdentityMatrix<double>(3);

        results[i] = IsotropicValue(prefac2 * gtens * Transpose(atensvec[i])-afc);
    }

    return results;
}

int main(int argc, char ** argv)
{

    try
    {
        //Parse the command line
        if(argc == 1)
        {
            PrintHelp_gFind();
            return 0;
        }

        SimpleMatrixD gtens;

        vector<double> shifts;

        double temperature = 0;
        double spin = 0;
        unsigned long centerAtom = 0;

        string inpfile;
        string coordfile;

        string argvstr;
        for(int i = 1; i < argc; i++)
        {
            argvstr = argv[i];
            if(argvstr == "-h" || argvstr == "--help")
            {
                PrintHelp_gFind();
                return 0;
            }
            else if (argvstr == "-t")
            {
                stringstream ss(GetNextArg(i,argc,argv));
                ss >> temperature;
                if(ss.fail() || ss.bad() || !ss.eof())
                    THROWSTREAM << "Error reading temperature from the command line";
            }
            else if (argvstr == "-s")
            {
                stringstream ss(GetNextArg(i,argc,argv));
                ss >> spin;
                if(ss.fail() || ss.bad() || !ss.eof())
                    THROWSTREAM << "Error reading spin from the command line";
            }

            else if (argvstr == "-f")
                inpfile = GetNextArg(i,argc,argv);
            else if (argvstr == "-c")
                coordfile = GetNextArg(i,argc,argv);
            else if (argvstr == "--center")
            {
                stringstream ss(GetNextArg(i,argc,argv));
                ss >> centerAtom;
                if(ss.fail() || ss.bad() || !ss.eof())
                    THROWSTREAM << "Error reading the center atom from the command line";
            }
            else
                THROWSTREAM << "Error reading command line: I don't know what \"" << argvstr << "\" means!";
        }

        if(inpfile == "")
            THROWSTREAM << "Error: I don't have an input file (-f option)";

        if(coordfile == "")
            THROWSTREAM << "Error: I don't have an xyz coordinate file (-c option)";

        if(temperature == 0)
            THROWSTREAM << "I need a temperature (-t option)";

        if(spin <= 0)
            THROWSTREAM << "I need a positive spin (-s option)";

        if(centerAtom <= 0)
            THROWSTREAM << "I need the paramagnetic center atom (--center option)";

        cout << "Temperature: " << temperature << "\n";
        cout << "Spin: " << spin << "\n";

        XYZFile coords(coordfile);


        XYZFile newcoords = coords;
        ReadFile_gFind(inpfile, gtens, shifts, newcoords);


        PrintTensor("Initial guess at g-tensor", gtens);
        cout << "\n\n";

        cout << "****************************************************\n"
             << "* Calculating the hyperfine tensor SD term using   *\n"
             << "*  purely geometric/dipolar formula                *\n"
             << "****************************************************\n";

        PrintCoords(newcoords);
        cout << "\n\n";



        for(int i = 0; i < 200; i++)
        {
            cout << "=====================================\n"
                 << "Step " << i << "\n"
                 << "=====================================\n";

            PrintTensor("Guess at g-tensor", gtens);

            vector<double> calcshifts = CalculateShifts(gtens, coords, newcoords, centerAtom, temperature, spin);

            double sumsq = 0;
            for(size_t i = 0; i < calcshifts.size(); i++)
                sumsq += pow(calcshifts[i] - shifts[i], 2);
            cout << "Sum of squares: " << sumsq << "\n";


            SimpleMatrixD gtensmod(3,3);

            // derivative w.r.t. gxx
            gtensmod = gtens;
            gtensmod(0,0) += 0.0001;
            calcshifts = CalculateShifts(gtensmod, coords, newcoords, centerAtom, temperature, spin);
            double sumsqpxx = 0;
            for(size_t i = 0; i < calcshifts.size(); i++)
                sumsqpxx += pow(calcshifts[i] - shifts[i], 2);
            double dgpxx = (sumsqpxx - sumsq)/0.0001;
            cout << "dgpxx: Sum of squares: " << sumsqpxx << " deriv = " << dgpxx << "\n";
            
            gtensmod = gtens;
            gtensmod(0,0) -= 0.0001;
            calcshifts = CalculateShifts(gtensmod, coords, newcoords, centerAtom, temperature, spin);
            double sumsqmxx = 0;
            for(size_t i = 0; i < calcshifts.size(); i++)
                sumsqmxx += pow(calcshifts[i] - shifts[i], 2);
            double dgmxx = (sumsq - sumsqmxx)/0.0001;
            cout << "dgmxx: Sum of squares: " << sumsqmxx << " deriv = " << dgmxx << "\n";

            // derivative w.r.t. gyy
            gtensmod = gtens;
            gtensmod(1,1) += 0.0001;
            calcshifts = CalculateShifts(gtensmod, coords, newcoords, centerAtom, temperature, spin);
            double sumsqpyy = 0;
            for(size_t i = 0; i < calcshifts.size(); i++)
                sumsqpyy += pow(calcshifts[i] - shifts[i], 2);
            double dgpyy = (sumsqpyy - sumsq)/0.0001;
            cout << "dgpyy: Sum of squares: " << sumsqpyy << " deriv = " << dgpyy << "\n";

            gtensmod = gtens;
            gtensmod(1,1) -= 0.0001;
            calcshifts = CalculateShifts(gtensmod, coords, newcoords, centerAtom, temperature, spin);
            double sumsqmyy = 0;
            for(size_t i = 0; i < calcshifts.size(); i++)
                sumsqmyy += pow(calcshifts[i] - shifts[i], 2);
            double dgmyy = (sumsq - sumsqmyy)/0.0001;
            cout << "dgmyy: Sum of squares: " << sumsqmyy << " deriv = " << dgmyy << "\n";

            // derivative w.r.t. gzz
            gtensmod = gtens;
            gtensmod(2,2) += 0.0001;
            calcshifts = CalculateShifts(gtensmod, coords, newcoords, centerAtom, temperature, spin);
            double sumsqpzz = 0;
            for(size_t i = 0; i < calcshifts.size(); i++)
                sumsqpzz += pow(calcshifts[i] - shifts[i], 2);
            double dgpzz = (sumsqpzz - sumsq)/0.0001;
            cout << "dgpzz: Sum of squares: " << sumsqpzz << " deriv = " << dgpzz << "\n";

            gtensmod = gtens;
            gtensmod(2,2) -= 0.0001;
            calcshifts = CalculateShifts(gtensmod, coords, newcoords, centerAtom, temperature, spin);
            double sumsqmzz = 0;
            for(size_t i = 0; i < calcshifts.size(); i++)
                sumsqmzz += pow(calcshifts[i] - shifts[i], 2);
            double dgmzz = (sumsq - sumsqmzz)/0.0001;
            cout << "dgmzz: Sum of squares: " << sumsqmzz << " deriv = " << dgmzz << "\n";


            //first derivatives
            double dgxx = (dgpxx+dgmxx)/2.0;
            double dgyy = (dgpyy+dgmyy)/2.0;
            double dgzz = (dgpzz+dgmzz)/2.0;

            //second derivatives
            double dg2xx = (dgpxx-dgmxx)/(2.0*0.0001);
            double dg2yy = (dgpyy-dgmyy)/(2.0*0.0001);
            double dg2zz = (dgpzz-dgmzz)/(2.0*0.0001);


            cout << "First Derivatives:\n" << dgxx << "\n" << dgyy << "\n" << dgzz << "\n\n";
            cout << "Second Derivatives:\n" << dg2xx << "\n" << dg2yy << "\n" << dg2zz << "\n\n";


            //Newton's method
            gtens(0,0) -= 0.4*(dgxx/dg2xx);
            gtens(1,1) -= 0.4*(dgyy/dg2yy);
            gtens(2,2) -= 0.4*(dgzz/dg2zz);
            
        }


    }
    catch(std::exception & ex)
    {
        cout << "\nEXCEPTION!\n";
        cout << ex.what();
        cout << "\n\n";
        return 1;
    }

    return 0;

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
