#include <string>
#include <vector>
#include <algorithm>
#include "ThrowStream.h"

using namespace std;

vector<long> ConvertSelectionString(const string & str, bool allowdupe)
{
    vector<long> ret;
    stringstream temp;
    unsigned int num;

    for(unsigned int i = 0; i < str.length(); i++)
    {
        if( isdigit(str[i]) )
            temp << str[i];

        else if ( str[i] == ',' )
        {
            if(temp.str().length() != 0)
            {
                temp >> num;
                ret.push_back(num);
                temp.clear();
            }
            else
                THROWSTREAM << "Bad formatting of selection string!";
        }

        else if( str[i] == '-' )
        {
            unsigned int start, end;

            temp >> start;
            temp.clear();
            end = 0;
            i++;

            while( isdigit(str[i]) && i < str.length() )
            {
                temp << str[i];
                i++;
            }

            temp >> end;
            temp.clear();

            if(start > end)
                THROWSTREAM << "Bad formatting of selection string (start of a range is larger than the end)";

            if(end == 0)
                THROWSTREAM << "Bad formatting of selection string!";


            for(unsigned int p = start; p <= end; p++)
                ret.push_back(p);

            if ( i != str.length() && str[i] != ',' )
                THROWSTREAM << "Bad formatting of selection string!";
        }

        else
            THROWSTREAM << "Bad formatting of selection string!";

        if ( i == (str.length() - 1) )
        {
            if(isdigit(str[i]))
            {
                // was taken care above
                //temp << str[i];

                temp >> num;
                ret.push_back(num);
                temp.clear();
            }
            else
                THROWSTREAM << "Bad formatting of selection string!";
        }
    }

    //remove duplicates
    if(allowdupe)
    {
        sort(ret.begin(), ret.end());
        return ret;
    }

    vector<long> actualret;

    for(unsigned int m = 0; m < ret.size(); m++)
    {
        bool dupe = false;

        for(unsigned int n = 0; n < actualret.size(); n++)
        {
            if(ret[m] == actualret[n])
                dupe = true;
        }

        if(!dupe)
            actualret.push_back(ret[m]);
    }

    //sort

    sort(actualret.begin(), actualret.end());

    return actualret;
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
