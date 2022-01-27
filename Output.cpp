#include <iostream>
#include <boost/format.hpp>

#include "PNMRShift.h"

using namespace std;

const char * contribNames[15] = {"GE.ANRFC", "GE.ANRSD", "GE.ASOFC", "GE.ASOSD", "GE.AS",
                                 "dGiso.ANRFC", "dGiso.ANRSD", "dGiso.ASOFC", "dGiso.ASOSD", "dGiso.AS",
                                 "dGtilde.ANRFC", "dGtilde.ANRSD", "dGtilde.ASOFC", "dGtilde.ASOSD", "dGtilde.AS"
                                };

void PrintTensor(const string & title, const SimpleMatrixD & mat)
{
    if(mat.nrows() != 3 || mat.ncols() != 3)
        THROWSTREAM << "Can't print a not 3x3 tensor - this one is " << mat.nrows() << "x" << mat.ncols();

    boost::format matrixrowfmt("%1$12.5f %2$12.5f %3$12.5f");

    cout << title << "\n"
         << "=============================================\n";
    for(unsigned short i = 0; i < 3; i++)
        cout << matrixrowfmt % mat(i,0) % mat(i,1) % mat(i,2) << "\n";

    cout << "---------------------------------------------\n"
         << "  Isotropic value = " << IsotropicValue(mat) << "\n";
}


void PrintHyperfine(const string & title, vector<ATensor> atensors)
{
    boost::format atensfmt(" %1$5s %2$5s   %3$9.5f   %4$9.5f %5$9.5f %6$9.5f   %7$9.5f %8$9.5f %9$9.5f   %10$9.5f %11$9.5f %12$9.5f");
    boost::format atensheadfmt(" %1$=5s %2$=5s   %3$=9.5f   %4$=9.5f %5$=9.5f %6$=9.5f   %7$=9.5f %8$=9.5f %9$=9.5f   %10$=9.5f %11$=9.5f %12$=9.5f");

    cout << "\n"
         << title << "\n"
         << atensheadfmt % "#" % "El" % "Iso. Value" % "FCSD x" % "FCSD y" % "FCSD z" % "PSOSO x" % "PSOSO y" % "PSOSO z" % "Total x" % "Total y" % "Total z" << "\n"
         << "========================================================================================================================\n";

    for(size_t i = 0; i < atensors.size(); i++)
    {
        cout << atensfmt % (i+1) % atensors[i].el % IsotropicValue(atensors[i].total) % atensors[i].fcsd(0,0) % atensors[i].fcsd(0,1) % atensors[i].fcsd(0,2)
             % atensors[i].psoso(0,0) % atensors[i].psoso(0,1) % atensors[i].psoso(0,2)
             % atensors[i].total(0,0) % atensors[i].total(0,1) % atensors[i].total(0,2) << "\n";
        cout << atensfmt % ""    % ""             % ""                                % atensors[i].fcsd(1,0) % atensors[i].fcsd(1,1) % atensors[i].fcsd(1,2)
             % atensors[i].psoso(1,0) % atensors[i].psoso(1,1) % atensors[i].psoso(1,2)
             % atensors[i].total(1,0) % atensors[i].total(1,1) % atensors[i].total(1,2) << "\n";
        cout << atensfmt % ""    % ""             % ""                                % atensors[i].fcsd(2,0) % atensors[i].fcsd(2,1) % atensors[i].fcsd(2,2)
             % atensors[i].psoso(2,0) % atensors[i].psoso(2,1) % atensors[i].psoso(2,2)
             % atensors[i].total(2,0) % atensors[i].total(2,1) % atensors[i].total(2,2) << "\n";
        cout << "\n";
    }

    cout << "\n\n";
}


void PrintOrbitalShield(const string & title, vector<OrbShieldTensor> orbtensors)
{
    boost::format atensfmt(" %1$5s %2$5s   %3$9.5f   %4$9.5f %5$9.5f %6$9.5f   %7$9.5f %8$9.5f %9$9.5f   %10$9.5f %11$9.5f %12$9.5f");
    boost::format atensheadfmt(" %1$=5s %2$=5s   %3$=9.5f   %4$=9.5f %5$=9.5f %6$=9.5f   %7$=9.5f %8$=9.5f %9$=9.5f   %10$=9.5f %11$=9.5f %12$=9.5f");

    cout << "\n"
         << title << "\n"
         << atensheadfmt % "#" % "El" % "Iso. Value" % "DIA x" % "DIA y" % "DIA z" % "PARA x" % "PARA y" % "PARA z" % "Total x" % "Total y" % "Total z" << "\n"
         << "========================================================================================================================\n";

    for(size_t i = 0; i < orbtensors.size(); i++)
    {
        cout << atensfmt % (i+1) % orbtensors[i].el % IsotropicValue(orbtensors[i].total) % orbtensors[i].dia(0,0) % orbtensors[i].dia(0,1) % orbtensors[i].dia(0,2)
             % orbtensors[i].para(0,0) % orbtensors[i].para(0,1) % orbtensors[i].para(0,2)
             % orbtensors[i].total(0,0) % orbtensors[i].total(0,1) % orbtensors[i].total(0,2) << "\n";
        cout << atensfmt % ""    % ""             % ""                                % orbtensors[i].dia(1,0) % orbtensors[i].dia(1,1) % orbtensors[i].dia(1,2)
             % orbtensors[i].para(1,0) % orbtensors[i].para(1,1) % orbtensors[i].para(1,2)
             % orbtensors[i].total(1,0) % orbtensors[i].total(1,1) % orbtensors[i].total(1,2) << "\n";
        cout << atensfmt % ""    % ""             % ""                                % orbtensors[i].dia(2,0) % orbtensors[i].dia(2,1) % orbtensors[i].dia(2,2)
             % orbtensors[i].para(2,0) % orbtensors[i].para(2,1) % orbtensors[i].para(2,2)
             % orbtensors[i].total(2,0) % orbtensors[i].total(2,1) % orbtensors[i].total(2,2) << "\n";
        cout << "\n";
    }

    cout << "\n\n";
}

void PrintCoords(const XYZFile & coords)
{
    boost::format coordfmt(" %1$5s %2$5s   %3$12.5f %4$12.5f %5$12.5f");
    boost::format coordheadfmt(" %1$5s %2$5s   %3$=12s %4$=12s %5$=12s");

    std::cout << "\nAtomic Coordinates (Angstroms)\n"
              << coordheadfmt % "#" % "El" % "x" % "y" % "z" << "\n"
              << "=============================================================\n";

    for(size_t i = 0; i < coords.atoms.size(); i++)
        cout << coordfmt % (i+1) % coords.atoms[i].el % coords.atoms[i].x % coords.atoms[i].y % coords.atoms[i].z << "\n";
}






void PrintResults_PNMRShift(const vector<Contributions> & cont,
                  const vector<ATensor> atensors,
                  const vector<OrbShieldTensor> & orbshield, const XYZFile & xyz,
                  const vector<string> & avgs, bool detail)
{
    boost::format     fmt(" %1$7s %2$7s   %3$10.4f %4$14.5f %5$14.5f %6$14.5f %7$14.5f %8$14.5f");
    boost::format headfmt(" %1$7s %2$7s   %3$10s %4$=14s %5$=14s %6$=14s %7$=14s %8$=14s");

    if(cont.size() != xyz.atoms.size() || orbshield.size() != xyz.atoms.size() ||  orbshield.size() != atensors.size())
        THROWSTREAM << "Final check error! cont,orb,xyz,atens = "
                    << cont.size() << ","
                    << orbshield.size() << ","
                    << xyz.atoms.size() << ","
                    << atensors.size();

    cout << "\nFinal results.  Shifts and Shieldings in ppm.\n"
         << "Elements for which I don't have the nuclear g factor may show \"inf\" or \"nan\"\n\n"
         << headfmt % "#" % "El" % "Aiso (MHz)" % "Orb Shielding" % "FC Shield" % "PC Shield" % "FC+PC Shield" % "Total**" << "\n"
         << "================================================================================================================\n";

    for(size_t i = 0; i < xyz.atoms.size(); i++)
    {
        if(orbshield[i].el != xyz.atoms[i].el && orbshield[i].el != "xx")
            THROWSTREAM << "BAD ATOM ORDERING IN ORBITAL SHIELDING. CHECK ABOVE";
        if(atensors[i].el != xyz.atoms[i].el && atensors[i].el != "xx")
            THROWSTREAM << "BAD ATOM ORDERING IN ORBITAL SHIELDING. CHECK ABOVE";

        double orbiso = IsotropicValue(orbshield[i].total);
        double fciso = IsotropicValue(cont[i].fermiContact);
        double pciso = IsotropicValue(cont[i].pseudoContact);

        cout << fmt % (i+1) % xyz.atoms[i].el % IsotropicValue(atensors[i].total)  % orbiso % fciso % pciso
             % (fciso + pciso) % (orbiso + fciso + pciso) << "\n";
    }

    cout << "\n";

    // Now do averages
    boost::format avgfmt(" %1$15s   %2$10.4f %3$14.5f %4$14.5f %5$14.5f %6$14.5f %7$14.5f");

    for(size_t i = 0; i < avgs.size(); i++)
    {
        string label, selstr;
        size_t colon = avgs[i].find_last_of(':');
        if(colon == avgs[i].npos)
            THROWSTREAM << "Error - bad formatting of averaging parameter: " << avgs[i];

        label = avgs[i].substr(0,colon);
        selstr = avgs[i].substr(colon+1); //gets to the end

        vector<long> sel = ConvertSelectionString(selstr, true);

        double avgorb = 0, avgfc = 0, avgpc = 0, avgsum = 0, avgtot = 0, avgatens = 0;

        for(size_t j = 0; j < sel.size(); j++)
        {
            if(sel[j] > (long)xyz.atoms.size() || sel[j] <= 0)
                THROWSTREAM << "Atom selection " << sel[j] << " for average '" << label << "' is not valid.";

            avgorb += IsotropicValue(orbshield[sel[j]-1].total);
            avgfc += IsotropicValue(cont[sel[j]-1].fermiContact);
            avgpc += IsotropicValue(cont[sel[j]-1].pseudoContact);
            avgsum += (IsotropicValue(cont[sel[j]-1].fermiContact) + IsotropicValue(cont[sel[j]-1].pseudoContact));
            avgtot += (IsotropicValue(orbshield[sel[j]-1].total) + IsotropicValue(cont[sel[j]-1].fermiContact) + IsotropicValue(cont[sel[j]-1].pseudoContact));
            avgatens += IsotropicValue(atensors[sel[j]-1].total);
        }

        avgorb /= sel.size();
        avgfc /= sel.size();
        avgpc /= sel.size();
        avgsum /= sel.size();
        avgtot /= sel.size();
        avgatens /= sel.size();

        cout << avgfmt % label % avgatens % avgorb % avgfc % avgpc % avgsum % avgtot << "\n";
    }

    cout << "\n";

    cout << "================================================================================================================\n"
         << "** Subtract this from the reference shielding to get total chemical shift.\n\n";


    //detailed results
    boost::format detailedatomhead(" Atom %1$4f %2$4s ( %3$9.5f %4$9.5f %5$9.5f )");
    boost::format detailedheadfmt(" %1$16s %2$15s     %3$=15s %4$=15s %5$=15s");
    boost::format     detailedfmt(" %1$16s %2$15.5f   %3$15.5f %4$15.5f %5$15.5f");


    if(detail)
    {
        cout << "\n\n\nDETAILED BREAKDOWN OF CONTRIBUTIONS\n";

        for(size_t i = 0; i < xyz.atoms.size(); i++)
        {

            SimpleMatrixD total(3,3);

            cout << "================================================================================================================\n";
            cout << detailedatomhead % (i+1) % xyz.atoms[i].el % xyz.atoms[i].x % xyz.atoms[i].y % xyz.atoms[i].z << "\n\n";
            cout << detailedheadfmt % "Contrib" % "Iso (ppm)" % "X" % "Y" % "Z" << "\n";
            cout << "================================================================================================================\n";

            for(int j = 0; j < 15; j++)
            {
                cout << detailedfmt % contribNames[j] % IsotropicValue(cont[i].detailed[j]) % cont[i].detailed[j](0,0) %  cont[i].detailed[j](0,1) %  cont[i].detailed[j](0,2) << "\n";
                cout <<                                               detailedfmt % "" % "" % cont[i].detailed[j](1,0) %  cont[i].detailed[j](1,1) %  cont[i].detailed[j](1,2) << "\n";
                cout <<                                               detailedfmt % "" % "" % cont[i].detailed[j](2,0) %  cont[i].detailed[j](2,1) %  cont[i].detailed[j](2,2) << "\n\n";

                total = total + cont[i].detailed[j];

            }

            cout << detailedfmt % "SUM"    % IsotropicValue(total) % total(0,0) %  total(0,1) %  total(0,2) << "\n";
            cout <<                          detailedfmt % "" % "" % total(1,0) %  total(1,1) %  total(1,2) << "\n";
            cout <<                          detailedfmt % "" % "" % total(2,0) %  total(2,1) %  total(2,2) << "\n\n";

            cout << detailedfmt % "ORB"    % IsotropicValue(orbshield[i].total) % orbshield[i].total(0,0) %  orbshield[i].total(0,1) %  orbshield[i].total(0,2) << "\n";
            cout <<                          detailedfmt % "" % "" % orbshield[i].total(1,0) %  orbshield[i].total(1,1) %  orbshield[i].total(1,2) << "\n";
            cout <<                          detailedfmt % "" % "" % orbshield[i].total(2,0) %  orbshield[i].total(2,1) %  orbshield[i].total(2,2) << "\n\n";

            total = total + orbshield[i].total;

            cout << detailedfmt % "TOTAL SHIELD"  % IsotropicValue(total) % total(0,0) %  total(0,1) %  total(0,2) << "\n";
            cout <<                          detailedfmt % "" % "" % total(1,0) %  total(1,1) %  total(1,2) << "\n";
            cout <<                          detailedfmt % "" % "" % total(2,0) %  total(2,1) %  total(2,2) << "\n\n";
            cout << "================================================================================================================\n";
            cout << "\n\n";

        }

        //Do averages
        for(size_t i = 0; i < avgs.size(); i++)
        {
            SimpleMatrixD total(3,3);

            string label, selstr;
            size_t colon = avgs[i].find_last_of(':');
            if(colon == avgs[i].npos)
                THROWSTREAM << "Error - bad formatting of averaging parameter: " << avgs[i];

            label = avgs[i].substr(0,colon);
            selstr = avgs[i].substr(colon+1); //gets to the end

            vector<long> sel = ConvertSelectionString(selstr, true);

            SimpleMatrixD avgcont[15], avgshield(3,3);
            for(size_t j = 0; j < 15; j++)
                avgcont[j] = SimpleMatrixD(3,3);

            for(size_t j = 0; j < sel.size(); j++)
            {
                if(sel[j] > (long)xyz.atoms.size() || sel[j] <= 0)
                    THROWSTREAM << "Atom selection " << sel[j] << " for average '" << label << "' is not valid.";

                for(unsigned short k = 0; k < 15; k++)
                    avgcont[k] += cont[sel[j]-1].detailed[k];

                avgshield += orbshield[sel[j]-1].total;

            }

            for(size_t j = 0; j < 15; j++)
                avgcont[j]/= sel.size();

            avgshield /= sel.size();


            cout << "================================================================================================================\n";
            cout << "Average: " << label << "\n\n";
            cout << detailedheadfmt % "Contrib" % "Iso (ppm)" % "X" % "Y" % "Z" << "\n";
            cout << "================================================================================================================\n";

            for(int j = 0; j < 15; j++)
            {
                cout << detailedfmt % contribNames[j] % IsotropicValue(avgcont[j]) % avgcont[j](0,0) %  avgcont[j](0,1) %  avgcont[j](0,2) << "\n";
                cout <<                                               detailedfmt % "" % "" % avgcont[j](1,0) %  avgcont[j](1,1) %  avgcont[j](1,2) << "\n";
                cout <<                                               detailedfmt % "" % "" % avgcont[j](2,0) %  avgcont[j](2,1) %  avgcont[j](2,2) << "\n\n";

                total = total + avgcont[j];

            }

            cout << detailedfmt % "SUM"    % IsotropicValue(total) % total(0,0) %  total(0,1) %  total(0,2) << "\n";
            cout <<                          detailedfmt % "" % "" % total(1,0) %  total(1,1) %  total(1,2) << "\n";
            cout <<                          detailedfmt % "" % "" % total(2,0) %  total(2,1) %  total(2,2) << "\n\n";

            cout << detailedfmt % "ORB"    % IsotropicValue(avgshield) % avgshield(0,0) %  avgshield(0,1) %  avgshield(0,2) << "\n";
            cout <<                          detailedfmt % "" % "" % avgshield(1,0) %  avgshield(1,1) %  avgshield(1,2) << "\n";
            cout <<                          detailedfmt % "" % "" % avgshield(2,0) %  avgshield(2,1) %  avgshield(2,2) << "\n\n";

            total = total + avgshield;

            cout << detailedfmt % "TOTAL SHIELD"  % IsotropicValue(total) % total(0,0) %  total(0,1) %  total(0,2) << "\n";
            cout <<                          detailedfmt % "" % "" % total(1,0) %  total(1,1) %  total(1,2) << "\n";
            cout <<                          detailedfmt % "" % "" % total(2,0) %  total(2,1) %  total(2,2) << "\n\n";
            cout << "================================================================================================================\n";
            cout << "\n\n";

        }
    }
}




void PrintHelp_PNMRShift(void)
{
    boost::format helpfmt("   %1$-15s   %2$-50s");

    cout << "\nUsage: PNMRShift [options]\n\n\n";

    cout << "General Parameters\n";
    cout << "------------------\n";
    cout << helpfmt % "-t" % "Temperature (in K)" << "\n";
    cout << helpfmt % "-s" % "Overall Spin (doublet = 0.5, triplet = 1.0, etc)" << "\n";

    cout << "\n\n";
    cout << helpfmt % "-f" % "Tensor input file (see manual for format)" << "\n";
    cout << helpfmt % "-c" % "Coordinates (standard XYZ file, including" << "\n";
    cout << helpfmt % "" %  "  lines with number of atoms and title (title" << "\n";
    cout << helpfmt % "" %  "  may be empty)" << "\n";


    cout << "\n\n";
    cout << "Optional Arguments\n";
    cout << "------------------\n";
    cout << helpfmt % "--geosd"      % "Use a purely-dipolar hyperfine tensor based solely on geometry" << "\n";
    cout << helpfmt % "--keepfc"     % "When used with --geosd, also keep the FC portions of the" << "\n";
    cout << helpfmt % ""             % "hyperfine tensor. Default is to discard them" << "\n";
    cout << helpfmt % "--correctA"   % "Correct the hyperfine tensor with the g-tensor. Default" << "\n";
    cout << helpfmt % ""             % "  is not to correct." << "\n";
    cout << helpfmt % "--detail"     % "Show detailed contributions (15 terms; see Hrobarik, Liimatainen," << "\n";
    cout << helpfmt % ""             % "  Pennanen, etc)." << "\n";
    cout << helpfmt % "-a"           % "Also show averaging of these nuclei (see manual for more info)\n";
    cout << helpfmt % "--deltag"     % "Input contains delta-g (in ppt) rather than the full g-tensor\n";
    cout << helpfmt % "--splitinp"   % "Input contains tensors split into components (g-tensor into\n";
    cout << helpfmt % ""             % "  dia- and paramagnetic, hyperfine into FCSD and PSOSO). Implies --deltag\n";
    cout << helpfmt % "--pvzfs"      % "Calculates the ZFS correction using the method of Pennanen and Vaara.\n";
    cout << helpfmt % "" 	     % "  Default is to use method of Soncini and Van den Heuvel.\n";
    cout << "\n\n";
}

void PrintHelp_gFind(void)
{
    cout << "\n\nHelp! I'm trapped in a program factory!\n\n";
}
// Copyright 2012-2015 Ben Pritchard, Bob Martin, and Jochen Autschbach
// This file is part of PNMRShift.
//
// PNMRShift is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// PNMRShift is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with PNMRShift.  If not, see <http://www.gnu.org/licenses/>.
