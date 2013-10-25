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

#include <cmath>
#include <error/matherror.hpp>

namespace math {

/** @addtogroup MathHelper
 *  
 *  @{
 */

// inline const unsigned biggest2( const unsigned a ) {
//     if ( (a == 0) || (a == 1) ) return 0;
//     unsigned k = 0;
//     unsigned b = a;
//     while ( b > 0 ) {
//         b >> 1;
//         k++;
//     }
//     return k;
// }
 
template < typename T > 
inline const T pow( const T b, const unsigned e ) {
    T aux = 1.;
    for ( unsigned k = 0; k < e; k++ )
        aux *= b;
    return aux;
}

//Compute product of all parameters
template < int ...P, typename T = int >
inline const T prod() {
    T res         = static_cast<T>(1);
    T elemArray[] = {P...};
    for ( unsigned k = 0; k < sizeof...(P); k++ )
        res*= static_cast<T>(elemArray[k]);
    return res;
}

//All parameters equal?
template < int ...P >
inline const bool equal() {
    bool res = true;
    int elemArray[] = {P...};
    for ( unsigned k = 1; k < sizeof...(P); k++ )
        res &= elemArray[k] == elemArray[0];
    return res;
}


/** @} */
}
