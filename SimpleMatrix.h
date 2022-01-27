//Simple matrix class
#include <typeinfo>
#include <complex>
#include <cstdlib>
#include <ostream>
#include "ThrowStream.h"

extern "C" {
    /* ZHEEV prototype */
    void zheev_( char* jobz, char* uplo, int* n, complex<double>* a, int* lda,
                 double* w, complex<double>* work, int* lwork, double* rwork, int* info );
}

template<typename T>
class SimpleMatrix
{
private:
    T * _data;
    unsigned long _nrows, _ncols;

    void _Allocate(unsigned long nrows, unsigned long ncols)
    {
        if(_data != NULL)
            delete [] _data;

        _data = new T[nrows*ncols];
        _nrows = nrows;
        _ncols = ncols;
        Zero();
    }

    void _Deallocate(void)
    {
        _nrows = _ncols = 0;
        if(_data != NULL)
            delete [] _data;
        _data = NULL;
    }

public:
    SimpleMatrix()
    {
        _data = NULL;
        _nrows = _ncols = 0;
    }

    SimpleMatrix(unsigned long nrows, unsigned long ncols)
    {
        _data = NULL;
        _nrows = _ncols = 0;
        _Allocate(nrows, ncols);
    }

    SimpleMatrix(const SimpleMatrix & s)
    {
        _data = NULL;
        _nrows = _ncols = 0;
        if(s._data != NULL)
        {
            _Allocate(s._nrows, s._ncols);
            for(unsigned long i = 0; i < _nrows*_ncols; i++)
                _data[i] = s._data[i];
        }
        else
        {
            _data = NULL;
            _nrows = _ncols = 0;
        }
    }

    SimpleMatrix(T * const array, unsigned long nrows, unsigned long ncols)
    {
        _data = NULL;
        _nrows = _ncols = 0;
        _Allocate(nrows, ncols);
        for(unsigned long i = 0; i < nrows*ncols; i++)
            _data[i] = array[i];
    }

    SimpleMatrix(T ** const array, unsigned long nrows, unsigned long ncols)
    {
        _data = NULL;
        _nrows = _ncols = 0;
        _Allocate(nrows, ncols);
        for(unsigned long i = 0; i < nrows; i++)
            for(unsigned long j = 0; j < ncols; i++)
                (*this)(i,j) = array[i][j];
    }


    ~SimpleMatrix()
    {
        delete [] _data;
    }

    unsigned long nrows(void) const
    {
        return _nrows;
    }
    unsigned long ncols(void) const
    {
        return _ncols;
    }

    void Zero(void)
    {
        for(unsigned long i = 0; i < _nrows*_ncols; i++)
            _data[i] = 0;
    }

    T & operator() (unsigned long row, unsigned long col)
    {
        if(row >= _nrows)
            THROWSTREAM << "Error: attempted to get row " << row << " from matrix with " << _nrows << " rows.";
        if(col >= _ncols)
            THROWSTREAM << "Error: attempted to get row " << col << " from matrix with " << _ncols << " rows.";
        return _data[col+row*_ncols];
    }

    const T & operator() (unsigned long row, unsigned long col) const
    {
        if(row >= _nrows)
            THROWSTREAM << "Error: attempted to get row " << row << " from matrix with " << _nrows << " rows.";
        if(col >= _ncols)
            THROWSTREAM << "Error: attempted to get row " << col << " from matrix with " << _ncols << " rows.";
        return _data[col+row*_ncols];
    }

    void IdentityMatrix(void)
    {
        if(_nrows != _ncols)
            THROWSTREAM << "Can't make a non-square matrix an identity matrix";

        for(unsigned long row = 0; row < _nrows; row++)
            for(unsigned long col = 0; col < _ncols; col++)
            {
                if(row == col)
                    (*this)(row,col) = 1;
                else
                    (*this)(row,col) = 0;
            }
    }

    template<typename V>
    SimpleMatrix<V> Convert(void) const
    {
        //attempt to convert one data type to another
        SimpleMatrix<V> result(_nrows,_ncols);
        for(unsigned long i = 0; i < _nrows; i++)
            for(unsigned long j = 0; j < _ncols; j++)
                result(i,j) = (*this)(i,j);
        return result;
    }


    SimpleMatrix & operator=(const SimpleMatrix & rhs)
    {
        if(this == &rhs) //self assignment
            return *this;
        if(_nrows != rhs._nrows || _ncols != rhs._ncols)
            _Allocate(rhs._nrows, rhs._ncols);
//           THROWSTREAM << "Cannot assign matrices of two different sizes! Trying to assign a " << rhs._nrows << "x" << rhs._ncols << " to a " << _nrows << "x" << _ncols;

        for(unsigned long i = 0; i < _nrows*_ncols; i++)
            _data[i] = rhs._data[i];
        return (*this);
    }

    SimpleMatrix & operator+=(const SimpleMatrix & rhs)
    {
        if(_nrows != rhs._nrows || _ncols != rhs._ncols)
            THROWSTREAM << "Cannot add matrices of two different sizes! Trying to add a " << rhs._nrows << "x" << rhs._ncols << " to a " << _nrows << "x" << _ncols;
        for(unsigned long i = 0; i < _nrows*_ncols; i++)
            _data[i] += rhs._data[i];
        return (*this);
    }

    SimpleMatrix & operator-=(const SimpleMatrix & rhs)
    {
        if(_nrows != rhs._nrows || _ncols != rhs._ncols)
            THROWSTREAM << "Cannot subtract two matrices of two different sizes! Trying to add a " << rhs._nrows << "x" << rhs._ncols << " to a " << _nrows << "x" << _ncols;
        for(unsigned long i = 0; i < _nrows*_ncols; i++)
            _data[i] -= rhs._data[i];
        return (*this);
    }


    const SimpleMatrix operator+(const SimpleMatrix & rhs) const
    {
        if(_nrows != rhs._nrows || _ncols != rhs._ncols)
            THROWSTREAM << "Cannot add matrices of two different sizes! Trying to add a " << rhs._nrows << "x" << rhs._ncols << " to a " << _nrows << "x" << _ncols;
        SimpleMatrix<T> result(*this);
        result += rhs;
        return result;
    }

    const SimpleMatrix operator-(const SimpleMatrix & rhs) const
    {
        if(_nrows != rhs._nrows || _ncols != rhs._ncols)
            THROWSTREAM << "Cannot subtract matrices of two different sizes! Trying to add a " << rhs._nrows << "x" << rhs._ncols << " to a " << _nrows << "x" << _ncols;
        SimpleMatrix result(*this);
        result -= rhs;
        return result;
    }

    const SimpleMatrix operator*(const SimpleMatrix & rhs) const
    {
        if(_ncols != rhs._nrows)
            THROWSTREAM << "Cannot multiply matrices! Trying to multiply a " << rhs._nrows << "x" << rhs._ncols << " to a " << _nrows << "x" << _ncols;
        SimpleMatrix result(_nrows,rhs._ncols);

        //could be sped up using raw data
        for(unsigned long i = 0; i < _nrows; i++)
            for(unsigned long j = 0; j < rhs._ncols; j++)
            {
                result(i,j) = 0; // in case I change the class to NOT zero on initialization
                for(unsigned long k = 0; k < _ncols; k++)
                    result(i,j) += (*this)(i,k)*rhs(k,j);
            }

        return result;
    }

    SimpleMatrix & operator*=(const SimpleMatrix & rhs)
    {
        if(_ncols != rhs._nrows)
            THROWSTREAM << "Cannot multiply matrices! Trying to multiply a " << rhs._nrows << "x" << rhs._ncols << " to a " << _nrows << "x" << _ncols;

        SimpleMatrix lhs(*this);
        if(this == &rhs)
        {
            //multiplying a matrix my itself
            SimpleMatrix rhs2(*this);
            *this = lhs*rhs2;
        }
        else
        {
            _Allocate(lhs._nrows,rhs._ncols);
            *this = lhs*rhs;
        }

        return *this;
    }


    template<typename S>
    SimpleMatrix & operator*=(const S & s)
    {
        for(unsigned long i = 0; i < (_nrows*_ncols); i++)
            _data[i] *= s;
        return *this;
    }

    template<typename S>
    SimpleMatrix & operator/=(const S & s)
    {
        for(unsigned long i = 0; i < (_nrows*_ncols); i++)
            _data[i] /= s;
        return *this;
    }

    template<typename S>
    const SimpleMatrix operator*(const S & s) const
    {
        SimpleMatrix result(*this);
        result *= s;
        return result;
    }

    template<typename S>
    const SimpleMatrix operator/(const S & s) const
    {
        SimpleMatrix result(*this);
        result /= s;
        return result;
    }

    SimpleMatrix ZHEEV(void)
    {
        //Check if this is a complex double matrix and error if it isn't
        if(typeid(T) != typeid(complex<double>))
            THROWSTREAM << "Error - SVD only works with matrices of type complex<double> for now";

        //The fortran COMPLEX(16) type is equivalent to std::complex<double>
        // (http://software.intel.com/sites/products/documentation/hpc/mkl/mkl_userguide_lnx/GUID-9A806754-1453-4433-9254-6FEA2C66B21B.htm)
        if(_nrows != _ncols)
            THROWSTREAM << "Error - matrix for ZHEEV must be a square matrix";
        if(_nrows == 0) //therefore _ncols = 0
            THROWSTREAM << "Error - empty matrix";


        //We need a triangular part. We will use the upper triangular part
        // However, note that fortran uses different array layout so
        // we will have to store the lower part in the upper part
        for(unsigned long i = 0; i < _ncols; i++)
            for(unsigned long j = 0; j <= i; j++)
                (*this)(j,i) = (*this)(i,j);
        for(unsigned long i = 0; i < _nrows; i++)
            for(unsigned long j = 0; j < i; j++)
                (*this)(i,j) = 0;

        complex<double> wkopt;
        int info, lwork;
        int n = (int)_nrows; //could be dangerous for large matrices! Should maybe check for value in _nrows (if greater than max int)
        int lda = n;
        double * w = new double[n];
        double * rwork = new double[3*n-2];

        //query for workspace
        lwork = -1;
        char v = 'V';
        char l = 'L';
        zheev_(&v, &l, &n, _data, &lda, w, &wkopt, &lwork, rwork, &info);
        lwork = (int)real(wkopt);
        complex<double> * work = new complex<double>[lwork];

        //actually do it
        info = 1;
        zheev_( &v, &l, &n, _data, &lda, w, work, &lwork, rwork, &info );

        // this matrix object will now contain the eigenvectors (stored in rows)
        // convert to columns
        T val;
        for(unsigned long i = 0; i < _nrows; i++)
            for(unsigned long j = 0; j < i; j++)
            {
                val = (*this)(i,j);
                (*this)(i,j) = (*this)(j,i);
                (*this)(j,i) = val;
            }


        // the w element contains the eigenvalues
        SimpleMatrix<T> eigval(n,1);
        for(int i = 0; i < n; i++)
            eigval(i,0) = w[i];

        //free up space
        delete [] w;
        delete [] rwork;
        delete [] work;

        if(info != 0)
            THROWSTREAM << "Error - zheev_ failed to converge. Info is " << info;

        return eigval;
    }


};


template<typename T, typename S>
const SimpleMatrix<T> operator*(const S & s, const SimpleMatrix<T> & mat)
{
    return mat * s;
}

template<typename T, typename S>
const SimpleMatrix<T> operator/(const S & s, const SimpleMatrix<T> & mat)
{
    return mat / s;
}

template<typename T>
const SimpleMatrix<T> IdentityMatrix(unsigned long size)
{
    SimpleMatrix<T> result(size,size);
    for(unsigned long i = 0; i < size; i++)
        result(i,i) = 1;
    return result;
}

template<typename T>
ostream& operator<<(ostream & os, const SimpleMatrix<T> & m)
{
    unsigned long nrows = m.nrows();
    unsigned long ncols = m.ncols();

    for(unsigned long i = 0; i < nrows; i++)
    {
        for(unsigned long j = 0; j < ncols; j++)
            os << m(i,j) << "  ";
        os << "\n";
    }

    return os;
}

template<typename T>
const SimpleMatrix<T> Transpose(const SimpleMatrix<T> & sm)
{
    SimpleMatrix<T> result(sm.ncols(), sm.nrows());

    for(unsigned long i = 0; i < sm.ncols(); i++)
        for(unsigned long j = 0; j < sm.nrows(); j++)
            result(i,j) = sm(j,i);

    return result;
}

template<typename T>
const SimpleMatrix<complex<T> > Conjugate(const SimpleMatrix<complex<T> > & sm) 
{
    SimpleMatrix<complex<T> > result(sm.nrows(),sm.ncols());
    for(unsigned long i = 0; i < sm.nrows(); i++)
    for(unsigned long j = 0; j < sm.ncols(); j++)
    {
        complex<T> val = sm(i,j);
        result(i,j) = complex<T>(val.real(),-1*val.imag());
    }
    return result;
}

//for non-complex matrices
template<typename T>
const SimpleMatrix<T> Conjugate(const SimpleMatrix<T> & sm) 
{
    return sm;
}

template<typename T>
const SimpleMatrix<T> ConjugateTranspose(const SimpleMatrix<T> & sm)
{
    //could probably be improved, speed wide
    return Conjugate(Transpose(sm));
}

typedef SimpleMatrix<float> SimpleMatrixF;
typedef SimpleMatrix<double> SimpleMatrixD;

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
