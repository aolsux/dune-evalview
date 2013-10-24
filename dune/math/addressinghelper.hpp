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

#include <utils/tuple.hpp>

namespace math {

/** @addtogroup AddressingHelper
 *  
 *  @{
 */
 
//! row-major matrix index (C/C++)
template< unsigned M, unsigned N >
inline const unsigned rmat_idx( const unsigned m, const unsigned n ) { return N*m + n; }

//! column-major matrix index (FOLTRAN)
template< unsigned M, unsigned N >
inline const unsigned cmat_idx( const unsigned m, const unsigned n ) { return m + M*n; }

//! row-major matrix index (C/C++)
inline const unsigned rmat_idx( const unsigned m, const unsigned n, const unsigned ColN ) { return ColN*m + n; }

//! column-major matrix index (FOLTRAN)
inline const unsigned cmat_idx( const unsigned m, const unsigned n, const unsigned RowM ) { return m + RowM*n; }



template< int ... N, typename T >
inline const int rten_idx( const TupleA< T, sizeof...(N) >& idx ) {
    const int len[] = {N...};

    int l = idx(0);
    for ( int k = 1; k < (int)sizeof...(N); k++ )
        l = len[k]*l + idx(k);

    return l;
}

template< int N0, typename T > inline const int rten_idx( const TupleA< T, 1 >& idx ) 
{ return idx(0); }

template< int N0, int N1, typename T > inline const int rten_idx( const TupleA< T, 2 >& idx ) 
{ return idx(1) + N1*idx(0); }

template< int N0, int N1, int N2, typename T > inline const int rten_idx( const TupleA< T, 3 >& idx ) 
{ return idx(2) + N2*(idx(1) + N1*idx(0)); }

template< int N0, int N1, int N2, int N3, typename T > inline const int rten_idx( const TupleA< T, 4 >& idx ) 
{ return idx(3) + N3*(idx(2) + N2*(idx(1) + N1*idx(0))); }

template< int N0, int N1, int N2, int N3, int N4, typename T > inline const int rten_idx( const TupleA< T, 5 >& idx ) 
{ return idx(4) + N4*(idx(3) + N3*(idx(2) + N2*(idx(1) + N1*idx(0)))); }


template< typename T >
inline const int rten_idx( const TupleB< T >& idx, const int N ) {
    int l = idx(0);
    for ( int k = 1; k < (int)idx.size; k++ )
        l = N*l + idx(k);

    return l;
}


template< int ... N, typename T >
inline const int cten_idx( const TupleA< T, sizeof...(N) >& idx ) {
    const int len[] = {N...};

    int l = idx(sizeof...(N)-1);
    for ( int k = (int)sizeof...(N)-2; k >= 0; k--)
        l = len[k]*l + idx(k);

    return l;
}

template< int N0, typename T > inline const int cten_idx( const TupleA< T, 1 >& idx ) 
{ return idx(0); }

template< int N0, int N1, typename T > inline const int cten_idx( const TupleA< T, 2 >& idx ) 
{ return idx(0) + N0*idx(1); }

template< int N0, int N1, int N2, typename T > inline const int cten_idx( const TupleA< T, 3 >& idx ) 
{ return idx(0) + N0*(idx(1) + N1*idx(2)); }

template< int N0, int N1, int N2, int N3, typename T > inline const int cten_idx( const TupleA< T, 4 >& idx ) 
{ return idx(0) + N0*(idx(1) + N1*(idx(2) + N2*idx(3))); }

template< int N0, int N1, int N2, int N3, int N4, typename T > inline const int cten_idx( const TupleA< T, 5 >& idx ) 
{ return idx(0) + N0*(idx(1) + N1*(idx(2) + N2*(idx(3) + N3*idx(4)))); }


template< int ... N, typename ... types >
inline const int rten_idx( const types ... i ) {
    static_assert( sizeof...(N) == sizeof...(types), "Invalid number of indices.");

    const int idx[] = {static_cast<int>(i)...};
    const int len[] = {static_cast<int>(N)...};

    int l = idx[0];
    for ( int k = 1; k < (int)sizeof...(i); k++ )
        l = len[k]*l + idx[k];

    return l;
}

template< int N0 > inline const int rten_idx( const int i0 ) 
{ return i0; }

template< int N0, int N1 > inline const int rten_idx( const int i0, const int i1 ) 
{ return i1 + N1*i0; }

template< int N0, int N1, int N2 > inline const int rten_idx( const int i0, const int i1, const int i2 ) 
{ return i2 + N2*(i1 + N1*i0); }

template< int N0, int N1, int N2, int N3 > inline const int rten_idx( const int i0, const int i1, const int i2, const int i3 ) 
{ return i3 + N3*(i2 + N2*(i1 + N1*i0)); }

template< int N0, int N1, int N2, int N3, int N4 > inline const int rten_idx( const int i0, const int i1, const int i2, const int i3, const int i4 ) 
{ return i4 + N4*(i3 + N3*(i2 + N2*(i1 + N1*i0))); }



template< int ... N, typename ... types >
inline const int cten_idx( const types ... i ) {
    static_assert( sizeof...(N) == sizeof...(types), "Invalid number of indices.");
    
    const int idx[] = {static_cast<int>(i)...};
    const int len[] = {static_cast<int>(N)...};
    
    int l = idx[sizeof...(i)-1];
    for ( int k = (int)sizeof...(i)-2; k >= 0; k--)
        l = len[k]*l + idx[k];
    
    return l;
}

template< int N0 > inline const int cten_idx( const int i0 ) 
{ return i0; }

template< int N0, int N1 > inline const int cten_idx( const int i0, const int i1 ) 
{ return i0 + N0*i1; }

template< int N0, int N1, int N2 > inline const int cten_idx( const int i0, const int i1, const int i2 ) 
{ return i0 + N0*(i1 + N1*i0); }

template< int N0, int N1, int N2, int N3 > inline const int cten_idx( const int i0, const int i1, const int i2, const int i3 ) 
{ return i0 + N0*(i1 + N1*(i2 + N2*i3)); }

template< int N0, int N1, int N2, int N3, int N4 > inline const int cten_idx( const int i0, const int i1, const int i2, const int i3, const int i4 ) 
{ return i0 + N0*(i1 + N1*(i2 + N2*(i3 + N3*i4))); }

/** @} */
}
