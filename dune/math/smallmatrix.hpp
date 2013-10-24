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
#include <math/addressinghelper.hpp>
#include <utils/fastmemcpy.h>

#define SMALL_MATRIX_INIT_ZERO

namespace math {

/** @addtogroup SmallMatrix
 *
 *  @{
 */

template< typename T, unsigned M, unsigned N >
struct SmallMatrix {
    T               data [M*N];
    const unsigned  MxN;

    SmallMatrix() : MxN(M*N) {
        #ifdef SMALL_MATRIX_INIT_ZERO
        memset( data, 0, sizeof(T)*MxN );
        #endif
    }

    SmallMatrix( const SmallMatrix< T, M, N >& rhs ) : MxN(M*N) {
        fast_memcpy( data, rhs.data, sizeof(T)*MxN );
    }

    template< typename ... _T >
    SmallMatrix( const _T ... in ) : MxN(M*N) {
        static_assert( (M*N == sizeof...(in)) || ((1 == sizeof...(in)) && (M==N)), "Invalid number of arguments in constructor" );

        T val [] = { static_cast<T>( in ) ... };

        if ( N == sizeof...(in) )
            for ( unsigned k = 0; k < MxN; k++ )
                data[k] = val[k];
        else
            for ( unsigned m = 0; m < M; m++ ) {
                data[rmat_idx<M,N>(m,m)] = val[0];
                for ( unsigned n = 0; n < m; n++ ) {
                    data[rmat_idx<M,N>(m,n)] = static_cast<T>(0.);
                    data[rmat_idx<M,N>(n,m)] = static_cast<T>(0.);
                }
            }
    }

    inline T& operator () ( const unsigned m, const unsigned n ) { return data[rmat_idx<M,N>(m,n)]; }
    inline const T& operator () ( const unsigned m, const unsigned n ) const { return data[rmat_idx<M,N>(m,n)]; }

    inline const SmallMatrix< T, M, N >& operator = ( const SmallMatrix< T, M, N >& rhs ) {
        fast_memcpy( data, rhs.data, sizeof(T)*MxN );
        return *this;
    }

    inline const SmallMatrix< T, M, N >& operator = ( const T rhs ) {
        for ( unsigned m = 0; m < M; m++ ) {
            data[rmat_idx<M,N>(m,m)] = rhs;
            for ( unsigned n = 0; n < m; n++ ) {
                data[rmat_idx<M,N>(m,n)] = static_cast<T>(0.);
                data[rmat_idx<M,N>(n,m)] = static_cast<T>(0.);
            }
        }
        return *this;
    }

    inline const SmallMatrix< T, M, N >& operator += ( const SmallMatrix< T, M, N >& rhs ) {
        for ( unsigned k = 0; k < MxN; k++ )
            data[k] += rhs.data[k];
        return *this;
    }

    inline const SmallMatrix< T, M, N >& operator -= ( const SmallMatrix< T, M, N >& rhs ) {
        for ( unsigned k = 0; k < MxN; k++ )
            data[k] += rhs.data[k];
        return *this;
    }

    inline const SmallMatrix< T, M, N >& operator *= ( const T rhs ) {
        for ( unsigned k = 0; k < MxN; k++ )
            data[k] *= rhs;
        return *this;
    }

    inline const SmallMatrix< T, M, N >& operator /= ( const T rhs ) {
        const T aux = static_cast<T>(1.)/rhs;
        for ( unsigned k = 0; k < MxN; k++ )
            data[k] *= aux;
        return *this;
    }
};

template< typename T, unsigned M, unsigned N >
inline const SmallMatrix< T, M, N > operator + ( const SmallMatrix< T, M, N >& A, const SmallMatrix< T, M, N >& B ) {
    SmallMatrix< T, M, N > C;
    for ( unsigned k = 0; k < C.MxN; k++ )
        C.data[k] = A.data[k] + B.data[k];
    return C;
}

template< typename T, unsigned M, unsigned N >
inline const SmallMatrix< T, M, N > operator - ( const SmallMatrix< T, M, N >& A, const SmallMatrix< T, M, N >& B ) {
    SmallMatrix< T, M, N > C;
    for ( unsigned k = 0; k < C.MxN; k++ )
        C.data[k] = A.data[k] - B.data[k];
    return C;
}

template< typename T, unsigned M, unsigned N >
inline const SmallMatrix< T, M, N > operator - ( const SmallMatrix< T, M, N >& A ) {
    SmallMatrix< T, M, N > C;
    for ( unsigned k = 0; k < C.MxN; k++ )
        C.data[k] = -A.data[k];
    return C;
}

template< typename T, unsigned L, unsigned M, unsigned N >
inline const SmallMatrix< T, L, N > operator * ( const SmallMatrix< T, L, M >& A, const SmallMatrix< T, M, N >& B ) {
    SmallMatrix< T, L, N > C;
    for ( unsigned l = 0; l < L; l++ )
        for ( unsigned n = 0; n < N; n++ ) {
            C.data[rmat_idx<L,N>(l,n)] = A.data[rmat_idx<L,M>(l,0)]*B.data[rmat_idx<M,N>(0,n)];
            for ( unsigned m = 1; m < M; m++ )
                C.data[rmat_idx<L,N>(l,n)] += A.data[rmat_idx<L,M>(l,m)]*B.data[rmat_idx<M,N>(m,n)];
        }
    return C;
}

template< typename T, unsigned M, unsigned N >
inline const SmallMatrix< T, M, N > operator * ( const T A, const SmallMatrix< T, M, N >& B ) {
    SmallMatrix< T, M, N > C;
    for ( unsigned k = 0; k < C.MxN; k++ )
        C.data[k] = A*B.data[k];
    return C;
}

template< typename T, unsigned M, unsigned N >
inline const SmallMatrix< T, M, N > operator * ( const SmallMatrix< T, M, N >& A, const T B ) {
    SmallMatrix< T, M, N > C;
    for ( unsigned k = 0; k < C.MxN; k++ )
        C.data[k] = A.data[k]*B;
    return C;
}

template< typename T, unsigned M, unsigned N >
inline const SmallMatrix< T, M, N > operator / ( const SmallMatrix< T, M, N >& A, const T B ) {
    SmallMatrix< T, M, N > C;
    const T aux = 1./B;
    for ( unsigned k = 0; k < C.MxN; k++ )
        C.data[k] = A.data[k]*aux;
    return C;
}

//======= determinant ===================================================================
// --> Laplacescher Entwicklungssatz:
//     http://de.wikipedia.org/wiki/Determinante#Laplacescher_Entwicklungssatz

template< typename T, unsigned N >
struct __functor_SmallMatrix_det {
    inline const T operator()( const SmallMatrix< T, N, N > &rhs ) const {
        __functor_SmallMatrix_det< T, N - 1 > det;
        SmallMatrix  < T, N - 1, N - 1 > minorant;
        T D = 0.;

        for ( unsigned j1 = 0; j1 < N; j1++ ) {
            for ( unsigned i = 1; i < N; i++ ) {
                unsigned j2 = 0;
                for ( unsigned j = 0; j < N; j++ ) {
                    if ( j == j1 ) continue;
                    minorant(i-1,j2) = rhs(i,j);
                    j2++;
                }
            }
            if (j1 & 1)
                D -= rhs(0,j1) * det( minorant );
            else
                D += rhs(0,j1) * det( minorant );
        }

        return D;
    }
};

template< typename T > struct __functor_SmallMatrix_det<T,1> {
    inline const T operator()( const SmallMatrix< T, 1, 1 > &rhs ) const {
        return   rhs.data[0];
    }
};

template< typename T > struct __functor_SmallMatrix_det<T,2> {
    inline const T operator()( const SmallMatrix< T, 2, 2 > &rhs ) const {
        return   rhs(0,0)*rhs(1,1) - rhs(0,1)*rhs(1,0);
    }
};

template< typename T > struct __functor_SmallMatrix_det<T,3> {
    inline const T operator()( const SmallMatrix< T, 3, 3 > &rhs ) const {
        return   rhs(0,0)*(rhs(1,1)*rhs(2,2) - rhs(2,1)*rhs(1,2))
               - rhs(1,0)*(rhs(0,1)*rhs(2,2) - rhs(2,1)*rhs(0,2))
               + rhs(2,0)*(rhs(0,1)*rhs(1,2) - rhs(1,1)*rhs(0,2));
    }
};

template< typename T, unsigned N >
inline const T det( const SmallMatrix< T, N, N > &rhs ) {
    __functor_SmallMatrix_det< T, N > func_det;
    return func_det( rhs );
}

//======= permanent =====================================================================
// --> Laplacescher Entwicklungssatz:
//     http://de.wikipedia.org/wiki/Determinante#Laplacescher_Entwicklungssatz

template< typename T, unsigned N >
struct __functor_SmallMatrix_perm {
    inline const T operator()( const SmallMatrix< T, N, N > &rhs ) const {
        __functor_SmallMatrix_det< T, N - 1 > det;
        SmallMatrix  < T, N - 1, N - 1 > minorant;
        T P = 0.;

        for ( unsigned j1 = 0; j1 < N; j1++ ) {
            for ( unsigned i = 1; i < N; i++ ) {
                unsigned j2 = 0;
                for ( unsigned j = 0; j < N; j++ ) {
                    if ( j == j1 ) continue;
                    minorant(i-1,j2) = rhs(i,j);
                    j2++;
                }
            }
            P += rhs(0,j1) * perm( minorant );
        }

        return P;
    }
};

template< typename T > struct __functor_SmallMatrix_perm<T,1> {
    inline const T operator()( const SmallMatrix< T, 1, 1 > &rhs ) const {
        return   rhs.data[0];
    }
};

template< typename T > struct __functor_SmallMatrix_perm<T,2> {
    inline const T operator()( const SmallMatrix< T, 2, 2 > &rhs ) const {
        return   rhs(0,0)*rhs(1,1) + rhs(0,1)*rhs(1,0);
    }
};

template< typename T > struct __functor_SmallMatrix_perm<T,3> {
    inline const T operator()( const SmallMatrix< T, 3, 3 > &rhs ) const {
        return   rhs(0,0)*(rhs(1,1)*rhs(2,2) + rhs(2,1)*rhs(1,2))
               + rhs(1,0)*(rhs(0,1)*rhs(2,2) + rhs(2,1)*rhs(0,2))
               + rhs(2,0)*(rhs(0,1)*rhs(1,2) + rhs(1,1)*rhs(0,2));
    }
};

template< typename T, unsigned N >
inline const T perm( const SmallMatrix< T, N, N > &rhs ) {
    __functor_SmallMatrix_perm< T, N > func_perm;
    return func_perm( rhs );
}

//======= functions on matrices =========================================================

template< typename T, unsigned M, unsigned N >
inline const SmallMatrix< T, N, M > transpose ( const SmallMatrix< T, M, N >& A ) {
    SmallMatrix< T, N, M > C;
    for ( unsigned m = 0; m < M; m++ )
        for ( unsigned n = 0; n < N; n++ )
            C.data[rmat_idx<N,M>(n,m)] = A.data[rmat_idx<M,N>(m,n)];
    return C;
}

template< typename T, unsigned N >
inline const T trace( const SmallMatrix< T, N, N >& A ) {
    T C = 0.;
    for ( unsigned n = 0; n < N; n++ )
        C += A.data[rmat_idx<N,N>(n,n)];
    return C;
}

template< typename T, unsigned N >
inline const SmallMatrix< T, N, 1 > diag( const SmallMatrix< T, N, N >& A ) {
    SmallMatrix< T, N, 1 > C;
    for ( unsigned n = 0; n < N; n++ )
        C.data[n] = A.data[rmat_idx<N,N>(n,n)];
    return C;
}

template< typename T, unsigned N >
inline const SmallMatrix< T, N, N > diag( const SmallMatrix< T, N, 1 >& A ) {
    SmallMatrix< T, N, N > C;
    for ( unsigned m = 0; m < N; m++ ) {
        C.data[rmat_idx<N,N>(m,m)] = A.data[m];
        #ifndef SMALL_MATRIX_INIT_ZERO
        for ( unsigned n = 0; n < m; n++ ) {
            C.data[rmat_idx<N,N>(m,n)] = static_cast<T>(0.);
            C.data[rmat_idx<N,N>(n,m)] = static_cast<T>(0.);
        }
        #endif
    }
    return C;
}

template< typename T, unsigned N >
inline const SmallMatrix< T, N, N > diag( const SmallMatrix< T, 1, N >& A ) {
    SmallMatrix< T, N, N > C;
    for ( unsigned m = 0; m < N; m++ ) {
        C.data[rmat_idx<N,N>(m,m)] = A.data[m];
        #ifndef SMALL_MATRIX_INIT_ZERO
        for ( unsigned n = 0; n < m; n++ ) {
            C.data[rmat_idx<N,N>(m,n)] = static_cast<T>(0.);
            C.data[rmat_idx<N,N>(n,m)] = static_cast<T>(0.);
        }
        #endif
    }
    return C;
}

template< typename T, unsigned N >
inline void identity( SmallMatrix< T, N, N >& C ) {
    for ( unsigned m = 0; m < N; m++ ) {
        C.data[rmat_idx<N,N>(m,m)] = static_cast<T>(1.);
        for ( unsigned n = 0; n < m; n++ ) {
            C.data[rmat_idx<N,N>(m,n)] = static_cast<T>(0.);
            C.data[rmat_idx<N,N>(n,m)] = static_cast<T>(0.);
        }
    }
}

template< typename T, unsigned M, unsigned N >
inline void zero( SmallMatrix< T, M, N >& C ) {
    memset( C.data, 0, sizeof(T)*C.MxN );
}


/** @} */
}
