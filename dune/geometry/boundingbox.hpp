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
#include <math/helper.hpp>
#include <math/shortvector.hpp>

namespace geometry {
 
template< typename T, unsigned dim >
class BoundingBox {
public:
    math::ShortVector< T, dim > dimension;
    math::ShortVector< T, dim > corner;
    
    BoundingBox() {            
        dimension = 1.;
        corner    = 0.;
    }
    
    BoundingBox( const math::ShortVector< T, dim >& c0, const math::ShortVector< T, dim >& d ) :
        dimension(d), corner(c0)
    {
    }
    
    BoundingBox( const BoundingBox< T, dim >& bb ) : dimension(bb.dimension), corner(bb.corner) {}
    
    const bool isInside( const math::ShortVector< T, dim >& p ) const {
        const math::ShortVector< T, dim > aux = p - corner;
        for ( unsigned k = 0; k < dim; k++ ) 
            if ( aux[k] > dimension[k] ) return false;
        return true;
    }
    
    const BoundingBox<T, dim> split( const unsigned orientation, const bool left ) const {
        BoundingBox<T, dim> bb( *this );
        
        bb.dimension(orientation) *= .5;
        if ( !left )
            bb.corner(orientation) += bb.dimension(orientation);
        
        return bb;
    }
};
    
}
