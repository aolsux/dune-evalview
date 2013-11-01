//**************************************************************************************//
//     AUTHOR: Malik Kirchner "malik.kirchner@gmx.net"                                  //
//             Martin Rückl "martin.rueckl@physik.hu-berlin.de"                         //
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

#include <cmath>
#include <cstring>
#include <vector>

#include <error/matherror.hpp>
#include <math/smallmatrix.hpp>
#include <utils/fastmemcpy.h>

extern "C" {
#include <emmintrin.h>
#include <pmmintrin.h>
#include <immintrin.h>   // (Meta-header, for GCC only)
}
#include <smmintrin.h>

#define SHORT_VECTOR_INIT_ZERO

namespace math {

/** @addtogroup ShortVector
 *
 *  @{
 */

/*!***************************************************************************************
 * @class ShortVector
 * @brief Plain representation of a math vector
 *
 * @url http://en.wikipedia.org/wiki/Euclidean_vector
 *****************************************************************************************/
template< typename T, unsigned N >
struct ShortVector {
    T data[N];

    /*!***************************************************************************************
    * Contruct and initialize with zero iff SHORT_VECTOR_INIT_ZERO is defined.
    *****************************************************************************************/
    ShortVector() {
        #ifdef SHORT_VECTOR_INIT_ZERO
        memset( data, 0, sizeof(T)*N );
        #endif
    }

    /*!***************************************************************************************
    * Copy/Move-Contruct.
    *****************************************************************************************/
    ShortVector( const ShortVector< T, N >&  rhs ) = default;
    ShortVector( ShortVector< T, N >&& rhs ) = default;

    /*!***************************************************************************************
    * Contruct from one or N scalars.
    *****************************************************************************************/
    template< typename ... _T >
    ShortVector( const _T ... in ) {
        static_assert( (N == sizeof...(in)) || (1 == sizeof...(in)), "Invalid number of arguments in constructor");

        T val [] = { static_cast<T>( in ) ... };

        if ( N == sizeof...(in) )
            for ( unsigned k = 0; k < N; k++ )
                data[k] = val[k];
         else
            for ( unsigned k = 0; k < N; k++ )
                data[k] = val[0];
    }

    /*!***************************************************************************************
    * Access k-th component by reference.
    *****************************************************************************************/
    inline T& operator () ( const unsigned k ) { return data[k]; }
    /*!***************************************************************************************
    * Access k-th component by constant reference.
    *****************************************************************************************/
    inline const T& operator()( const unsigned k ) const { return data[k]; }

    /*!***************************************************************************************
    * Assign the same scalar right hand side to each component.
    *****************************************************************************************/
     ShortVector< T, N >& operator = ( const T rhs ) {
        for ( unsigned k = 0; k < N; k++ )
            data[k] = rhs;
        return *this;
    }

    /*!***************************************************************************************
    * Assign an other ShortVector.
    *****************************************************************************************/
    ShortVector< T, N >& operator = ( const ShortVector< T, N >& rhs ) = default;

//     /*!***************************************************************************************
//     * Assign a \f$N\times 1\f$ Matrix.
//     *****************************************************************************************/
//     inline const ShortVector< T, N >& operator = ( const SmallMatrix< T, N, 1 >& rhs ) {
//         fast_memcpy( data, rhs.data, sizeof(T)*N );
//         return *this;
//     }

    /*!***************************************************************************************
    * Vector sum in place.
    *****************************************************************************************/
    inline const ShortVector< T, N >& operator += ( const ShortVector< T, N >& rhs ) {
        for ( unsigned k = 0; k < N; k++ )
            data[k] += rhs.data[k];
        return *this;
    }

    /*!***************************************************************************************
    * Vector difference in place.
    *****************************************************************************************/
    inline const ShortVector< T, N >& operator -= ( const ShortVector< T, N >& rhs ) {
        for ( unsigned k = 0; k < N; k++ )
            data[k] += rhs.data[k];
        return *this;
    }

    /*!***************************************************************************************
    * Multiply every component in place by a scalar right hand side.
    *****************************************************************************************/
    inline const ShortVector< T, N >& operator *= ( const T rhs ) {
        for ( unsigned k = 0; k < N; k++ )
            data[k] *= rhs;
        return *this;
    }

    /*!***************************************************************************************
    * Devide every component in place by a scalar right hand side.
    *****************************************************************************************/
    inline const ShortVector< T, N >& operator /= ( const T rhs ) {
        const T aux = static_cast<T>(1.)/rhs;
        for ( unsigned k = 0; k < N; k++ )
            data[k] *= aux;
        return *this;
    }
};

// template< unsigned N >
// struct ShortVector<float, N> {
//     const unsigned size;
//     __m128 data[((N%4) ? (N + 4 - N%4) : (N))/4];
//
//     /*!***************************************************************************************
//     * Contruct and initialize with zero iff SHORT_VECTOR_INIT_ZERO is defined.
//     *****************************************************************************************/
//     ShortVector() : size(((N%4) ? (N + 4 - N%4) : (N))/4) {
//         const __m128 aux = _mm_set_ps(0.f, 0.f, 0.f, 0.f);
//         for ( unsigned k = 0; k < size; k++ )
//             data[k] = aux;
//     }
//
//     /*!***************************************************************************************
//     * Copy-Contruct.
//     *****************************************************************************************/
//     ShortVector( const ShortVector< float, N >& rhs ) : size(((N%4) ? (N + 4 - N%4) : (N))/4) {
//         for ( unsigned k = 0; k < size; k++ )
//             data[k] = rhs.data[k];
//     }
//
//     /*!***************************************************************************************
//     * Contruct from one or N scalars.
//     *****************************************************************************************/
//     template< typename ... _T >
//     ShortVector( const _T ... in ) : size(((N%4) ? (N + 4 - N%4) : (N))/4) {
//         static_assert( (N == sizeof...(in)) || (1 == sizeof...(in)), "Invalid number of arguments in constructor");
//
//         const __m128 aux = _mm_set_ps(0.f, 0.f, 0.f, 0.f);
//         for ( unsigned k = 0; k < size; k++ )
//             data[k] = aux;
//
//         float val [] = { static_cast<float>( in ) ... };
//
//         if ( N == sizeof...(in) )
//             for ( unsigned k = 0; k < N; k++ )
//                 reinterpret_cast<float*>(data)[k] = val[k];
//         else
//             for ( unsigned k = 0; k < N; k++ )
//                 reinterpret_cast<float*>(data)[k] = val[0];
//
//     }
//
//     /*!***************************************************************************************
//     * Access k-th component by reference.
//     *****************************************************************************************/
//     inline float& operator () ( const unsigned k ) { return reinterpret_cast<float*>(data)[k]; }
//     /*!***************************************************************************************
//     * Access k-th component by constant reference.
//     *****************************************************************************************/
//     inline const float& operator()( const unsigned k ) const { return reinterpret_cast<float*>(data)[k]; }
//
//     /*!***************************************************************************************
//     * Assign the same scalar right hand side to each component.
//     *****************************************************************************************/
//     inline const ShortVector< float, N >& operator = ( const float rhs ) {
//         const __m128 aux = _mm_set_ps(rhs, rhs, rhs, rhs);
//         for ( unsigned k = 0; k < size; k++ )
//             data[k] = aux;
//         return *this;
//     }
//
//     /*!***************************************************************************************
//     * Assign an other ShortVector.
//     *****************************************************************************************/
//     inline const ShortVector< float, N >& operator = ( const ShortVector< float, N >& rhs ) {
//         fast_memcpy( data, rhs.data, sizeof(__m128)*size );
//         return *this;
//     }
//
//     /*!***************************************************************************************
//     * Devide every component in place by a scalar right hand side.
//     *****************************************************************************************/
//     inline const ShortVector< float, N >& operator /= ( const float rhs ) {
//         const float  aux0 = 1.f/rhs;
//         const __m128 aux1 = _mm_set_ps(aux0, aux0, aux0, aux0);
//         for ( unsigned k = 0; k < size; k++ )
//             data[k] = _mm_mul_ps(data[k], aux1);
//         return *this;
//     }
// };

/*!***************************************************************************************
 * Vector sum \f$A+B\f$.
 *****************************************************************************************/
template< typename T, unsigned N >
inline const ShortVector< T, N > operator + ( const ShortVector< T, N >& A, const ShortVector< T, N >& B ) {
    ShortVector< T, N > C;
    for ( unsigned k = 0; k < N; k++ )
        C.data[k] = A.data[k] + B.data[k];
    return C;
}

/*!***************************************************************************************
 * Vector difference \f$A-B\f$.
 *****************************************************************************************/
template< typename T, unsigned N >
inline const ShortVector< T, N > operator - ( const ShortVector< T, N >& A, const ShortVector< T, N >& B ) {
    ShortVector< T, N > C;
    for ( unsigned k = 0; k < N; k++ )
        C.data[k] = A.data[k] - B.data[k];
    return C;
}

/*!***************************************************************************************
 * Multiply A by \f$-1\f$.
 *****************************************************************************************/
template< typename T, unsigned N >
inline const ShortVector< T, N > operator - ( const ShortVector< T, N >& A ) {
    ShortVector< T, N > C;
    for ( unsigned k = 0; k < N; k++ )
        C.data[k] = -A.data[k];
    return C;
}

/*!***************************************************************************************
 * Multiply every component of B by a scalar A.
 *****************************************************************************************/
template< typename T, unsigned N >
inline const ShortVector< T, N > operator * ( const T A, const ShortVector< T, N >& B ) {
    ShortVector< T, N > C;
    for ( unsigned k = 0; k < N; k++ )
        C.data[k] = A*B.data[k];
    return C;
}

/*!***************************************************************************************
 * Multiply every component of A by a scalar B.
 *****************************************************************************************/
template< typename T, unsigned N >
inline const ShortVector< T, N > operator * ( const ShortVector< T, N >& A, const T B ) {
    ShortVector< T, N > C;
    for ( unsigned k = 0; k < N; k++ )
        C.data[k] = A.data[k]*B;
    return C;
}

/*!***************************************************************************************
 * Devide every component of A by a scalar B.
 *****************************************************************************************/
template< typename T, unsigned N >
inline const ShortVector< T, N > operator / ( const ShortVector< T, N >& A, const T B ) {
    ShortVector< T, N > C;
    const T aux = static_cast<T>(1.)/B;
    for ( unsigned k = 0; k < N; k++ )
        C.data[k] = A.data[k]*aux;
    return C;
}

/*!***************************************************************************************
 * Scalar product \f$ \left<\cdot,\cdot\right> : \{a, b\}\in R^n \times R^n \rightarrow c\in R, \quad c = \sum_k\, a_k b_k \f$
 *****************************************************************************************/
template< typename T, unsigned N >
inline const T dot( const ShortVector< T, N >& A, const ShortVector< T, N >& B ) {
    T C = A.data[0]*B.data[0];
    for ( unsigned k = 1; k < N; k++ )
        C += A.data[k]*B.data[k];
    return C;
}

///*!***************************************************************************************
// * Scalar product (float) \f$ \left<\cdot,\cdot\right> : \{a, b\}\in R^n \times R^n \rightarrow c\in R, \quad c = \sum_k\, a_k b_k \f$
// *****************************************************************************************/
// template< unsigned N >
// inline const float dot( const ShortVector< float, N >& A, const ShortVector< float, N >& B ) {
//     const __m128 aux = _mm_set_ps(0.f, 0.f, 0.f, 0.f);
//
//
//     for ( unsigned k = 0; k < A.size; k++ )
//         aux = _mm_add_ps( aux,  _mm_add_ps( A.data[k], B.data[k] ));
//
//     const float* C = reinterpret_cast<const float*>(aux);
//     return C[0]+C[0]+C[0]+C[0];
// }

/*!***************************************************************************************
 * Cross product \f$ \times : \{a, b\}\in R^3 \times R^3 \rightarrow c\in R^3, \quad c_i = \sum_{j,k}\,  a_j b_k \, \varepsilon^{ijk} \f$
 *****************************************************************************************/
template< typename T >
inline const ShortVector< T, 3 > cross( const ShortVector< T, 3 >& A, const ShortVector< T, 3 >& B ) {
    ShortVector< T, 3 > C;

    C.data[0] = A.data[1]*B.data[2];
    C.data[0]-= A.data[2]*B.data[1];
    C.data[1] = A.data[2]*B.data[0];
    C.data[1]-= A.data[0]*B.data[2];
    C.data[2] = A.data[0]*B.data[1];
    C.data[2]-= A.data[1]*B.data[0];

    return C;
}

/*!***************************************************************************************
 * Triple product \f$ \left[\cdot,\cdot,\cdot\right] : \{a, b, c\}\in R^3 \times R^3 \times R^3 \rightarrow d\in R, \quad d = a \cdot (b\times c)= (a \times b)\cdot c \f$
 *****************************************************************************************/
template< typename T >
inline const T triple( const ShortVector< T, 3 >& A, const ShortVector< T, 3 >& B, const ShortVector< T, 3 >& C ) {
    return dot( A, cross(B,C) );
}

/*!***************************************************************************************
 * Norm squared \f$ \left< A,A \right> \f$
 *****************************************************************************************/
template< typename T, unsigned N >
inline const T norm2( const ShortVector< T, N >& A ) {
    return dot(A,A);
}

/*!***************************************************************************************
 * Norm \f$ \sqrt{\left< A,A \right>} \f$
 *****************************************************************************************/
template< typename T, unsigned N >
inline const T norm( const ShortVector< T, N >& A ) {
    return sqrt(dot(A,A));
}

/*!***************************************************************************************
 * Normalized vector \f$ A/\sqrt{\left< A,A \right>} \f$
 *****************************************************************************************/
template< typename T, unsigned N >
inline const ShortVector< T, N > normalized( const ShortVector< T, N >& A ) {
    ShortVector< T, N > C = A;
    const T aux           = static_cast<T>(1.)/norm(C);
    C                    *= aux;
    return C;
}

/*!***************************************************************************************
 * Normalized vector \f$ A/\sqrt{\left< A,A \right>} \f$
 *****************************************************************************************/
template< typename T, unsigned N >
inline void normalize( ShortVector< T, N >& A ) {
    const T aux  = static_cast<T>(1.)/norm(A);
    A           *= aux;
}

/*!***************************************************************************************
 * Angle between vectors \f$ \text{acos}\left(\frac{\left< A,B \right>}{\|A\|\cdot\|B\|}\right) \f$
 *****************************************************************************************/
template< typename T, unsigned N >
inline const T angle( const ShortVector< T, N >& A, const ShortVector< T, N >& B ) {
    return acos( dot(A,B)/sqrt( norm2(A) * norm2(B) ) );
}

/*!***************************************************************************************
 * Set all components to zero.
 *****************************************************************************************/
template< typename T, unsigned N >
inline void zero( ShortVector< T, N >& C ) {
    memset( C.data, 0, sizeof(T)*N );
}

/*!***************************************************************************************
 * Stream formated.
 *****************************************************************************************/
template< typename T, unsigned N >
inline std::ostream& operator<< ( std::ostream& out, const ShortVector<T, N>& v ) {
    out << "( ";
    for ( unsigned k = 0; k < N; k++ )
        out << v.data[k] << ((k < N - 1 ) ? ", " : "");
    out << " )";
    return out;
}


/** @} */
}
