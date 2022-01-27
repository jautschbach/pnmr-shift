#ifndef __IG_BPLIB_THROWSTREAM_H__
#define __IG_BPLIB_THROWSTREAM_H__

#include <exception>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>

class ThrowStream;
std::ostream & operator<<(std::ostream & os, const ThrowStream & ts);

struct ThrowStreamLocation
{
    unsigned long line;
    std::string file;
    std::string function;

    ThrowStreamLocation(unsigned long iline, const std::string & ifile, const std::string & ifunction)
    {
        line = iline;
        file = ifile;
        function = ifunction;
    }
};

#define TSHERE ThrowStreamLocation(__LINE__,__FILE__,__FUNCTION__)

class ThrowStream : public std::exception
{
private:
    std::string _desc;
public:
    ThrowStream(const unsigned long line, const std::string & file, const std::string & function);
    ThrowStream(const std::exception & ex, unsigned long line, const std::string & file, const std::string & function);
    ~ThrowStream() throw() { };

    ThrowStream & Append(const unsigned long line, const std::string & file, const std::string & function);
    ThrowStream & Append(const std::exception & ex, const unsigned long line,
                         const std::string & file, const std::string & function);

    char const* what() const throw();

    template<typename T>
    ThrowStream & operator<<(const T & rhs);

    friend std::ostream & operator<<(std::ostream & os, const ThrowStream & ts);
};


#define THROWSTREAM throw ThrowStream(__LINE__, __FILE__, __FUNCTION__)
#define THROWSTREAMAPPEND(ex) throw ThrowStream( (ex), __LINE__, __FILE__, __FUNCTION__)
#define THROWSTREAMOBJ(ex) ThrowStream (ex)(__LINE__, __FILE__, __FUNCTION__); (ex)
#define THROWSTREAMOBJAPPEND(ex) (ex).Append(__LINE__, __FILE__, __FUNCTION__)
#define THROWSTREAMOBJAPPENDCOPY(ex,ey) ex.Append((ey), __LINE__, __FILE__, __FUNCTION__)
#define THROWSTREAMOBJCOPY(ex,ey) ThrowStream (ex)((ey),__LINE__, __FILE__, __FUNCTION__ ); (ex)

template<typename T>
ThrowStream &
ThrowStream::operator<<(const T & rhs)
{
    std::stringstream ss;
    ss << rhs;
    _desc.append(ss.str());
    return *this;
}

inline ThrowStream::ThrowStream(const unsigned long line, const std::string & file, const std::string & function)
{
    Append(line, file, function);
}

inline ThrowStream::ThrowStream(const std::exception & ex, unsigned long line, const std::string & file, const std::string & function)
{
    Append(ex, line, file, function);
}

inline ThrowStream &
ThrowStream::Append(const unsigned long line, const std::string & file, const std::string & function)
{
    std::stringstream ss;

    ss << '\n';

#ifdef THROWSTREAM_EXCEPTIONSOURCE
    ss << "( " << file << ":" << line << " , in " << function << "() )    ->  ";
#endif

    _desc.append(ss.str());

    return *this;
}

inline ThrowStream &
ThrowStream::Append(const std::exception & ex, const unsigned long line,
                    const std::string & file, const std::string & function)
{
    //depends on if this is actually a throwstream
    const ThrowStream * pts;
    if((pts = dynamic_cast<const ThrowStream *>(&ex)))
    {
        _desc.append(pts->_desc);
    }
    else
    {
        //will end up a double append, but oh well, it's the
        //best I can do with only an std::exception
        Append(line,file,function) << ex.what();
    }

    Append(line,file,function);

    return *this;
}

inline char const* ThrowStream::what() const throw()
{
    return _desc.c_str();
}

inline std::ostream & operator<<(std::ostream & os, const ThrowStream & ts)
{
    os << ts._desc;
    return os;
}


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
