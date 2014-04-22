//**************************************************************************************//
//     AUTHOR: Malik Kirchner "malik.kirchner@gmx.net"                                  //
//                                                                                      //
//     This program is free software: you can redistribute it and/or modify             //
//     it under the terms of the GNU General Public License as published by             //
//     the Free Software Foundation, either version 3 of the License, or                //
//     (at your option) any later version.                                              //
//                                                                                      //
//     This program is distributed in the hope that it will be useful,                  //
//     but WITHOUT ANY WARRANTY; without even the implied warranty of                   //
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                    //
//     GNU General Public License for more details.                                     //
//                                                                                      //
//     You should have received a copy of the GNU General Public License                //
//     along with this program.  If not, see <http://www.gnu.org/licenses/>.            //
//                                                                                      //
//     Dieses Programm ist Freie Software: Sie können es unter den Bedingungen          //
//     der GNU General Public License, wie von der Free Software Foundation,            //
//     Version 3 der Lizenz oder (nach Ihrer Option) jeder späteren                     //
//     veröffentlichten Version, weiterverbreiten und/oder modifizieren.                //
//                                                                                      //
//     Dieses Programm wird in der Hoffnung, dass es nützlich sein wird, aber           //
//     OHNE JEDE GEWÄHRLEISTUNG, bereitgestellt; sogar ohne die implizite               //
//     Gewährleistung der MARKTFÄHIGKEIT oder EIGNUNG FÜR EINEN BESTIMMTEN ZWECK.       //
//     Siehe die GNU General Public License für weitere Details.                        //
//                                                                                      //
//     Sie sollten eine Kopie der GNU General Public License zusammen mit diesem        //
//     Programm erhalten haben. Wenn nicht, siehe <http://www.gnu.org/licenses/>.       //
//                                                                                      //
//**************************************************************************************//
/*! \file */ 
#pragma once

#include <vector>
#include <error/baseerror.hpp>

template< typename T, unsigned N >
struct TupleA {
    T data[N];

    TupleA() { for (unsigned k = 0; k < N; k++ ) data[k] = 0; }

    template< typename ... _T >
    TupleA( const _T ... in ) { 
        static_assert( (N == sizeof...(in)) || (1 == sizeof...(in)), "Invalid number of arguments in constructor");

        if ( N == sizeof...(in) ) {
            T val [] = { static_cast<T>( in ) ... };
            for ( unsigned k = 0; k < N; k++ )
                data[k] = val[k];
        } else {
            T val [] = { static_cast<T>( in ) ... };
            for ( unsigned k = 0; k < N; k++ )
                data[k] = val[0];
        }

    }

    inline T& operator()( const unsigned k ) { return data[k]; }
    inline const T& operator()( const unsigned k ) const { return data[k]; }
};

template< typename T, unsigned N >
inline std::ostream& operator<< ( std::ostream& out, const TupleA<T, N>& v ) {
    out << "( ";
    for ( unsigned k = 0; k < N; k++ )     
        out << v.data[k] << ((k < N - 1 ) ? ", " : "");
    out << " )";
    return out;
}




template< typename T >
struct TupleB {
    typedef T* _Tp;
    
    const size_t size;
    const _Tp    data;

    TupleB( const size_t s ) : size(s), data( new T [s] ) { for (unsigned k = 0; k < s; k++ ) data[k] = 0; }

    template< typename ... _T >
    TupleB( const _T ... in ) : size(sizeof...(in)), data( new T [size] )  { 
            T val [] = { static_cast<T>( in ) ... };
            for ( unsigned k = 0; k < size; k++ )
                data[k] = val[k];
    }
    
    template< typename ... _T >
    TupleB( const unsigned N, const _T ... in ) : size(sizeof...(in)), data( new T [size] ) { 
        if ( (N != sizeof...(in)) && (1 != sizeof...(in)) ) throw BaseError(/*"Invalid number of arguments in constructor"*/__ERROR_INFO__);

        if ( N == sizeof...(in) ) {
            T val [] = { static_cast<T>( in ) ... };
            for ( unsigned k = 0; k < N; k++ )
                data[k] = val[k];
        } else {
            T val [] = { static_cast<T>( in ) ... };
            for ( unsigned k = 0; k < N; k++ )
                data[k] = val[0];
        }

    }
    
    inline T& operator()( const unsigned k ) { return data[k]; }
    inline const T& operator()( const unsigned k ) const { return data[k]; }
};

template< typename T >
inline std::ostream& operator<< ( std::ostream& out, const TupleB<T>& v ) {
    out << "( ";
    for ( unsigned k = 0; k < v.size; k++ )     
        out << v.data[k] << ((k < v.size - 1 ) ? ", " : "");
    out << " )";
    return out;
}
