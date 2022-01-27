#include <iostream>
#include <cmath>
#include <complex>
#include <vector>
#include <algorithm>

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




int main(int argc, char ** argv)
{

    try
    {
        //Parse the command line
        if(argc == 1)
        {
            PrintHelp_PNMRShift();
            return 0;
        }

        gTensor gtens;
        DTensor dtens = IdentityMatrix<double>(3);

        vector<ATensor> atensvec;
        vector<OrbShieldTensor> orbshieldvec;

        double temperature = 0;
        double spin = 0;
        bool correctA = false;
        bool detail = false;
        bool geometricSD = false;
        bool keepfc = false;
        unsigned long centerAtom = 0;
        bool simpleinput = true;
        bool fullg = true;
	bool pvsso = false;

        vector<string> averaging;

        string inpfile;
        string coordfile;

        string argvstr;
        for(int i = 1; i < argc; i++)
        {
            argvstr = argv[i];
            if(argvstr == "-h" || argvstr == "--help")
            {
                PrintHelp_PNMRShift();
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
            else if (argvstr == "-a")
            {
                //averaging
                averaging.push_back(GetNextArg(i,argc, argv));
            }

            else if (argvstr == "-f")
                inpfile = GetNextArg(i,argc,argv);
            else if (argvstr == "-c")
                coordfile = GetNextArg(i,argc,argv);
            else if (argvstr == "--correctA")
                correctA = true;
            else if (argvstr == "--detail")
                detail = true;
            else if (argvstr == "--keepfc")
                keepfc = true;
            else if (argvstr == "--geosd")
            {
                geometricSD = true;
                stringstream ss(GetNextArg(i,argc,argv));
                ss >> centerAtom;
                if(ss.fail() || ss.bad() || !ss.eof())
                    THROWSTREAM << "Error reading the center atom from the command line";
            }
            else if (argvstr == "--splitinp")
            {
                simpleinput = false;
                fullg = false;
            }
            else if (argvstr == "--deltag")
                fullg = false;
	    else if (argvstr == "--pvzfs")
		pvsso = true;
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

        if(keepfc && !geometricSD)
            THROWSTREAM << "--keepfc can only be used with --geosd";

        if(detail && simpleinput)
            THROWSTREAM << "Detail shouldn't be used without --splitinp";

	cout << "    ###################" << "\n";
	cout << "    pNMRShift program" << "\n";
	cout << "    ###################" << "\n" << "\n";
	cout << "developed by B. Pritchard and B. Martin" << "\n" << "\n";
	cout << "research group of Jochen Autschbach, SUNY Buffalo" << "\n" << "\n";
	cout << "recommended citations:" << "\n" << "\n";
	cout << "[1] Autschbach, J.; Patchkovskii, S.; Pritchard, B.," << "\n";
	cout << "   J. Chem. Theory Comput. 2011, 7, 2175 - 2188." << "\n";
	cout << "   DOI:10.1021/ct2000143w" << "\n" << "\n";
	cout << "[2] Martin, B.; Autschbach, J.," << "\n";
	cout << "   J. Chem. Phys. 2015, 142, 054108." << "\n";
	cout << "   DOI: 10.1063/1.4906318" << "\n" << "\n";
	cout << "=============================================================" << "\n" << "\n";


	//Print Input Parameters
        cout << "Temperature: " << temperature << "\n";
        cout << "Spin: " << spin << "\n";

        XYZFile coords(coordfile);

        //fill these vectors with empty tensors
        atensvec.resize(coords.atomcount);
        orbshieldvec.resize(coords.atomcount);

        if(simpleinput)
          ReadFile_PNMRShift_simple(inpfile, gtens, atensvec, orbshieldvec, dtens, averaging, fullg);
        else
          ReadFile_PNMRShift(inpfile, gtens, atensvec, orbshieldvec, dtens, averaging);

        PrintCoords(coords);
        cout << "\n\n";

        if(!simpleinput)
        {
          PrintTensor("Diamagnetic Delta-g tensor (ppt)", 1000.0*gtens.dia);
          cout << "\n\n";

          PrintTensor("Paramagnetic Delta-g tensor (ppt)", 1000.0*gtens.para);
          cout << "\n\n";
        }

        PrintTensor("Total Delta-g tensor (ppt)", 1000.0*(gtens.deltagisomat + gtens.deltagtilde));
        cout << "\n\n";

        PrintTensor("Total g-tensor", gtens.total);
        cout << "\n\n";

        PrintHyperfine("Hyperfine tensors from input file (MHz)", atensvec);

        if(geometricSD)
        {
            cout << "****************************************************\n"
                 << "* Calculating the hyperfine tensor SD term using   *\n"
                 << "*  purely geometric/dipolar formula                *\n"
                 << "****************************************************\n";

            //check for the center atom
            if(centerAtom <= 0 || centerAtom > coords.atomcount)
                THROWSTREAM << "Error - center atom out of range. It must be in the "
                            << "range [1," << coords.atomcount << "]. I got " << centerAtom;

            //Get the coordinates of the center atom, converting to meters from angstroms
            double centerx = coords.atoms[centerAtom-1].x * 1E-10;
            double centery = coords.atoms[centerAtom-1].y * 1E-10;
            double centerz = coords.atoms[centerAtom-1].z * 1E-10;

            // this is missing GN and the powers of r
            // this also include a conversion from J to MHz
            double prefac = (1E-6 * CONSTANT_MU0 * CONSTANT_GE * CONSTANT_BE * CONSTANT_BN) / (4*CONSTANT_PI * CONSTANT_H);

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
                atensvec[i].nrsd = prefac2 * ((3*rvec*Transpose(rvec)) - r2I);

                //Leave the contact terms alone, if available. Then recalculate the new total
                atensvec[i].sosd = atensvec[i].as = SimpleMatrixD(3,3); //zeroes the matrices

                atensvec[i].el = coords.atoms[i].el;

                //for debugging, we can disble this
                if(!keepfc)
                {
                    atensvec[i].nrfc.Zero();
                    atensvec[i].sofc.Zero();
                }

                atensvec[i].fcsd = atensvec[i].nrfc + atensvec[i].nrsd;
                atensvec[i].psoso = atensvec[i].sofc;
                atensvec[i].total = atensvec[i].fcsd + atensvec[i].psoso;
            }

            PrintHyperfine("Hyperfine tensors with purely dipolar SD terms (MHz)", atensvec);
        }


        if(correctA)
        {
            cout << "****************************************************\n"
                 << "* Correcting the hyperfine tensor                  *\n"
                 << "****************************************************\n";

            SimpleMatrixD correction = gtens.total / CONSTANT_GE;
            for(vector<ATensor>::iterator it = atensvec.begin(); it != atensvec.end(); ++it)
            {
                it->fcsd *= correction;
                it->psoso *= correction;
                it->total *= correction;
                it->nrfc *= correction;
                it->nrsd *= correction;
                it->sofc *= correction;
                it->sosd *= correction;
                it->as *= correction;
            }
            PrintHyperfine("Corrected Hyperfine tensors (MHz)", atensvec);
        }

        PrintOrbitalShield("Orbital shielding tensor (ppm)", orbshieldvec);
        cout << "\n\n";

        PrintTensor("ZFS Tensor (cm-1)", dtens);
        cout << "\n\n";


        // Do some sanity checks
        // note that the 'id' members contain the atom number (1,2,3,4...) so it is not zero-based
        for(size_t i = 0; i < coords.atoms.size(); i++)
        {
            if(i != (coords.atoms[i].id - 1))
                THROWSTREAM << "Error - bad indexing in coords. Index " << i << " has atom " << coords.atoms[i].id
                            << " but should have " << i+1;

            if(orbshieldvec[i].id != 0 && i != (orbshieldvec[i].id - 1))
                THROWSTREAM << "Error - bad indexing in orbshieldvec. Index " << i << " has atom " << orbshieldvec[i].id
                            << " but should have " << i+1;

            if(atensvec[i].id != 0 && i != (atensvec[i].id - 1))
                THROWSTREAM << "Error - bad indexing in atensvec. Index " << i << " has atom " << atensvec[i].id
                            << " but should have " << i+1;



            if(orbshieldvec[i].el != coords.atoms[i].el && orbshieldvec[i].el != "xx")
                THROWSTREAM << "Bad atom in orbital shielding. The shielding for atom " << orbshieldvec[i].id
                            << " states that it is '" << orbshieldvec[i].el << "' but it should be '" << coords.atoms[i].el
                            << "' according to the XYZ file";


            if(atensvec[i].el != coords.atoms[i].el && atensvec[i].el != "xx")
                THROWSTREAM << "Bad atom in atensvec. The atensor for atom " << atensvec[i].id
                            << " states that it is '" << atensvec[i].el << "' but it should be '" << coords.atoms[i].el
                            << "' according to the XYZ file";


        }





        //calculate the contributions
        SimpleMatrixD ss0(3,3);

        /********************************************************
          Calculate Z Matrix
         ********************************************************/
        //check to make sure that the spin is a multiple of 1/2
        double fractpart,intpart;
        fractpart = modf(spin * 2.0, &intpart);

        if(fractpart != 0)
            THROWSTREAM << "Error - spin must be a multiple of 1/2";

        unsigned int nstates = 2.0*spin+1;

        //Smat[0] -> Sx
        //Smat[1] -> Sy
        //Smat[2] -> Sz
        // Sp, Sm = S+, S-
        SimpleMatrix<complex<double> > Smat[3], Sp(nstates,nstates), Sm(nstates,nstates);

        //S+
        for(unsigned int r = 0; r < (nstates-1); r++)
        {
            double m = (spin-1) - r;
            Sp(r,r+1) = sqrt((spin*(spin+1))-(m*(m+1)));
        }

        //S-
        for(unsigned int r = 1; r < nstates; r++)
        {
            double m = spin - (r-1);
            Sm(r,r-1) = sqrt((spin*(spin+1))-(m*(m-1)));
        }
        //Sz
        Smat[2] = SimpleMatrix<complex<double> >(nstates,nstates); //needs to be allocated first
        for(unsigned int i = 0; i <  nstates; i++)
            Smat[2](i,i) = spin - i;

        //Sx
        Smat[0] = 0.5*(Sp + Sm);

        //Sy
        Smat[1] = -0.5*complex<double>(0,1)*(Sp - Sm);

        //cout << "nstates: " << nstates << "\n";
        //cout << "S+:\n" << Sp << "\n\n";
        //cout << "S-:\n" << Sm << "\n\n";
        //cout << "Sx:\n" << Smat[0] << "\n\n";
        //cout << "Sy:\n" << Smat[1] << "\n\n";
        //cout << "Sz:\n" << Smat[2] << "\n\n";
        //cout << "S2:\n" << Smat[0]*Smat[0] + Smat[1]*Smat[1] + Smat[2]*Smat[2] << "\n\n";

        //first, convert zfs to SI units (joules)
        SimpleMatrixD zfstensorjr = dtens * (100.0 * CONSTANT_C * CONSTANT_H);
        //convert to a complex matrix so we can multiply it with other complex matrices
        SimpleMatrix<complex<double> > zfstensorj = zfstensorjr.Convert<complex<double> >();


        //cout << "zfstensor: " << dtens << "\n\n";
        //cout << "zfstensorj: " << zfstensorj << "\n\n";

        //Calculate SDS
        SimpleMatrix<complex<double> > sds(nstates,nstates);

        for(unsigned short i = 0; i <= 2; i++)
            for(unsigned short j = 0; j <= 2; j++)
                sds += zfstensorj(i,j) * (Smat[i]*Smat[j]);

        //cout << "sds: " << sds << "\n\n";

        SimpleMatrix<complex<double> > eigval = sds.ZHEEV();
        //cout << "eigval: " << eigval << "\n\n";
        //cout << "eigvec: " << sds << "\n\n";

        //check orthonormality of eigvectors
        //cout << "Checking orthonormality\n";
        //cout << ConjugateTranspose(sds)*sds << "\n\n";

        //actually calculate Z Matrix
	SimpleMatrix<complex<double> > ss0complex(3,3);
	double denom = 0;
	double beta = CONSTANT_K * temperature;
	if (pvsso){

	for(unsigned int i = 0; i < nstates; i++)
	    denom += exp(-1.0*real(eigval(i,0))/beta);

	for(int a = 0; a < 3; a++)
	    for(int b = 0; b < 3; b++)
	    {

	        for(unsigned int i = 0; i < nstates; i++)
	        {
	            //grab a vector
	            SimpleMatrix<complex<double> > vec(nstates,1);
	            for(unsigned int j = 0; j < nstates; j++)
			 vec(j,0) = sds(j,i);

                    //cout << "vec: " << vec << "\n";
                    //cout << "cvec: " << ConjugateTranspose(vec) << "\n";
                    SimpleMatrix<complex<double> > scalar = (ConjugateTranspose(vec) * (Smat[a] * (Smat[b] * vec)));
                    if(scalar.nrows() != 1 || scalar.ncols() != 1)
                        THROWSTREAM << "Expecting a single value. Instead, got a " << scalar.nrows() << "x" << scalar.ncols() << " matrix.";
                        ss0complex(a,b) += exp(-1.0 * eigval(i,0)/beta) * scalar(0,0);
                    }
                 
                 ss0complex(a,b) ;
	   }
	}
	else {

	std::vector<std::vector<int> > tmpvec;
	tmpvec.clear();
        std::vector<int> myvector;
	for (unsigned int i = 0; i < nstates; i++){
		myvector.clear();
		for (unsigned int j = 0; j < nstates; j++){
			if ( i == j ) {
			myvector.push_back(i);
			}
			else {
				if ( abs(real(eigval(i,0))-real(eigval(j,0))) < (CONSTANT_K/70000.00) ) {
				myvector.push_back(j);
				}
			}
			}
	
		
	tmpvec.push_back(myvector);
	}	

	tmpvec.erase(std::unique(tmpvec.begin(), tmpvec.end()), tmpvec.end());


	SimpleMatrix<complex <double> > tmpss0(1, 1);
        for(unsigned int d = 0; d < nstates; d++)
            denom += exp(-1.0*real(eigval(d,0))/beta);

            for(unsigned int i = 0; i < tmpvec.size(); i++)
                for(unsigned int j = 0; j < tmpvec[i].size(); j++)
          		for(unsigned int k = 0; k < tmpvec.size(); k++)
                    		for(unsigned int l = 0; l < tmpvec[k].size(); l++)
				        for(int a = 0; a < 3; a++)
				            for(int b = 0; b < 3; b++)
            {
                //for(unsigned int ii = 0; ii < nstates; ii++)
                {
                    //grab a vector
                    SimpleMatrix<complex<double> > vecj(nstates,1);
                    SimpleMatrix<complex<double> > vecl(nstates,1);
                    //unsigned int j;
                    for(unsigned int m = 0; m < nstates; m++){
                        vecj(m,0) = sds(m,tmpvec[i][j]);
			vecl(m,0) = sds(m,tmpvec[k][l]);
		    }

                    //cout << "vecj: " << vecj << << "\n";
                    //cout << "cvec: " << ConjugateTranspose(vecj) << "\n";
                    SimpleMatrix<complex<double> > scalar = ((ConjugateTranspose(vecj) * Smat[a] * vecl) * (ConjugateTranspose(vecl) * Smat[b] * vecj));
                    if(scalar.nrows() != 1 || scalar.ncols() != 1)
                        THROWSTREAM << "Expecting a single value. Instead, got a " << scalar.nrows() << "x" << scalar.ncols() << " matrix.";
                    if ( i==k ) {
			tmpss0(0, 0) = exp(-1.0 * eigval(tmpvec[i][j],0)/beta) * scalar(0,0);
                        ss0complex(a,b) += tmpss0(0, 0);
                    }
                    else {
		    	tmpss0(0, 0) = ((2*beta/(real(eigval(tmpvec[k][l],0))-real(eigval(tmpvec[i][j],0))))*(exp(-1.0 * eigval(tmpvec[i][j],0)/beta) ) ) * scalar(0,0);
                    	ss0complex(a,b) += tmpss0(0, 0);
                    }
                }
                ss0complex(a,b) ;
                //ss0complex(a,b) /= denom;

            }
	}
        //cout << "\nss0complex:\n" << ss0complex;
        for(unsigned int i = 0; i < nstates; i++)
        {
            if( abs(imag(eigval(i,0))/real(eigval(i,0))) > 1e-10 )
            {
                cout << "===============================================\n";
                cout << "WARNING WARNING WARNING WARNING WARNING WARNING\n";
                cout << "  Found element of eigval that has a significant\n";
                cout << "  imaginary part. Element [" << i << "]\n";
                cout << "       Real part: " << real(eigval(i,0)) << "\n";
                cout << "  Imaginary part: " << imag(eigval(i,0)) << "\n";
                cout << "WARNING WARNING WARNING WARNING WARNING WARNING\n";
                cout << "===============================================\n";
            }
        }
        for(unsigned int i = 0; i < 3; i++)
        {
            for(int j = 0; j < 3; j++)
            {
                //check to see if the imaginary parts are negligible
                if( abs(imag(ss0complex(i,j))/real(ss0complex(i,j))) > 1e-6 )
                {
                    cout << "===============================================\n";
                    cout << "WARNING WARNING WARNING WARNING WARNING WARNING\n";
                    cout << "  Found element of Z Matrix that has a significant\n";
                    cout << "  imaginary part. Element [" << i << "," << j << "]\n";
                    cout << "       Real part: " << real(ss0complex(i,j)) << "\n";
                    cout << "  Imaginary part: " << imag(ss0complex(i,j)) << "\n";
                    cout << "WARNING WARNING WARNING WARNING WARNING WARNING\n";
                    cout << "===============================================\n";
                }

                ss0(i,j) = real(ss0complex(i,j))/denom;
            }
        }
        /********************************************************************************
         Done with Z Matrix
         ********************************************************************************/

        PrintTensor("Z Matrix", ss0);
        cout << "\n\n";

        //We put a negative sign here to represent these in terms of shielding
        double prefac = (-1.0e12 * CONSTANT_H * CONSTANT_BE)/
                        (CONSTANT_BN * CONSTANT_K * temperature);


        //g . ss0 . A
        // apply ss0 to g. Product is associative
        gtens.dia *= ss0;
        gtens.para *= ss0;
        gtens.total *= ss0;

        gtens.gemat *= ss0;
        gtens.deltagisomat *= ss0;
        gtens.deltagtilde *= ss0;

        vector<Contributions> cont;
        cont.resize(atensvec.size());

        for(size_t i = 0; i < atensvec.size(); i++)
        {
            SimpleMatrixD acon(atensvec[i].nrfc + atensvec[i].sofc);
            SimpleMatrixD adip(atensvec[i].nrsd + atensvec[i].sosd + atensvec[i].as);

            double prefac2 = prefac / NuclearG(coords.atoms[i].el);

            cont[i].fermiContact = prefac2 * gtens.total * Transpose(acon);
            cont[i].pseudoContact = prefac2 * gtens.total * Transpose(adip);

            //calculate the detailed contributions
            cont[i].detailed[0] = prefac2 * gtens.gemat * Transpose(atensvec[i].nrfc);
            cont[i].detailed[1] = prefac2 * gtens.gemat * Transpose(atensvec[i].nrsd);
            cont[i].detailed[2] = prefac2 * gtens.gemat * Transpose(atensvec[i].sofc);
            cont[i].detailed[3] = prefac2 * gtens.gemat * Transpose(atensvec[i].sosd);
            cont[i].detailed[4] = prefac2 * gtens.gemat * Transpose(atensvec[i].as);

            cont[i].detailed[5] = prefac2 * gtens.deltagisomat * Transpose(atensvec[i].nrfc);
            cont[i].detailed[6] = prefac2 * gtens.deltagisomat * Transpose(atensvec[i].nrsd);
            cont[i].detailed[7] = prefac2 * gtens.deltagisomat * Transpose(atensvec[i].sofc);
            cont[i].detailed[8] = prefac2 * gtens.deltagisomat * Transpose(atensvec[i].sosd);
            cont[i].detailed[9] = prefac2 * gtens.deltagisomat * Transpose(atensvec[i].as);

            cont[i].detailed[10] = prefac2 * gtens.deltagtilde * Transpose(atensvec[i].nrfc);
            cont[i].detailed[11] = prefac2 * gtens.deltagtilde * Transpose(atensvec[i].nrsd);
            cont[i].detailed[12] = prefac2 * gtens.deltagtilde * Transpose(atensvec[i].sofc);
            cont[i].detailed[13] = prefac2 * gtens.deltagtilde * Transpose(atensvec[i].sosd);
            cont[i].detailed[14] = prefac2 * gtens.deltagtilde * Transpose(atensvec[i].as);
        }


        PrintResults_PNMRShift(cont, atensvec, orbshieldvec, coords, averaging, detail);


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
