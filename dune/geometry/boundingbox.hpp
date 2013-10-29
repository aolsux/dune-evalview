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

#include <utils/utils.hpp>
#include <math/shortvector.hpp>
#include <limits>
#include <string.h>

#include <boost/serialization/serialization.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>

namespace geometry {

template< typename BT, unsigned dim >
class BoundingBox {
private:
    friend class boost::serialization::access;
    
    template< class Archive >
    void serialize(Archive & ar, const unsigned int version )    {
        using namespace boost::serialization;
        for ( unsigned k = 0; k < dim; k++) {
            ar & make_nvp( ("corner_"   + asString(k)).c_str() , _corner(k));
        }
        for ( unsigned k = 0; k < dim; k++) {
            ar & make_nvp( ("dimension_" + asString(k)).c_str() , _dimension(k));
        }
        
        update();
    }
    
    bool _origin_unknown;
    
protected:
    math::ShortVector< BT, dim > _corner;
    math::ShortVector< BT, dim > _dimension;
    math::ShortVector< BT, dim > _center;
    math::ShortVector< BT, dim > _upper;
    
    
    inline void update() {
        _center = _corner + (BT).5*_dimension;
        _upper  = _corner + _dimension;
    }
    
public:
    
    BoundingBox()  : 
        _origin_unknown(true), 
        _corner(.0), 
        _dimension(.0), 
        _center(.0), 
        _upper(.0) {}

    BoundingBox( const BoundingBox<BT, dim>& bb )  :
        _origin_unknown(bb._origin_unknown), 
        _corner(bb._corner), 
        _dimension(bb._dimension), 
        _center(bb._center), 
        _upper(bb._upper) {}
    
    BoundingBox( const std::vector< BT >&    bbc, const std::vector< BT >&    bbd ) :
        _origin_unknown(false), 
        _corner( bbc ), 
        _dimension( bbd )
    {
        update();
    }

    inline const void                include(const math::ShortVector< BT, dim >& p);
    inline const void                include(const BoundingBox< BT, dim>& bb);
    inline const void                enlarge( const BT alpha );
    inline const BT                  maxDimension() const;
    inline const bool                checkInside( const math::ShortVector< BT, dim >& p ) const;
    
    
    inline math::ShortVector< BT, dim > dimension() { return _dimension; }
    inline math::ShortVector< BT, dim > center   () { return _center;    }
    inline math::ShortVector< BT, dim > corner   () { return _corner;    }
    inline math::ShortVector< BT, dim > lower    () { return _corner;    }
    inline math::ShortVector< BT, dim > upper    () { return _upper;     }
    inline BT                           dimension( const unsigned d ) { return _dimension(d); }
    inline BT                           center   ( const unsigned d ) { return _center(d);    }    
    inline BT                           corner   ( const unsigned d ) { return _corner(d);    }
    inline BT                           lower    ( const unsigned d ) { return _corner(d);    }
    inline BT                           upper    ( const unsigned d ) { return _upper(d);     }
    
    inline const math::ShortVector< BT, dim > dimension() const { return _dimension; }
    inline const math::ShortVector< BT, dim > center   () const { return _center;    }
    inline const math::ShortVector< BT, dim > corner   () const { return _corner;    }
    inline const math::ShortVector< BT, dim > lower    () const { return _corner;    }
    inline const math::ShortVector< BT, dim > upper    () const { return _upper;     }
    inline const BT                           dimension( const unsigned d ) const { return _dimension(d); }
    inline const BT                           center   ( const unsigned d ) const { return _center(d);    }    
    inline const BT                           corner   ( const unsigned d ) const { return _corner(d);    }
    inline const BT                           lower    ( const unsigned d ) const { return _corner(d);    }
    inline const BT                           upper    ( const unsigned d ) const { return _upper(d);     }
    
    inline void set( math::ShortVector< BT, dim > c, math::ShortVector< BT, dim > d ) {
        _origin_unknown = false;
        _corner     = c;
        _dimension  = d;
        update();
    }
    
    inline void set( const std::vector< BT >&    bbc, const std::vector< BT >&    bbd ) {
        _origin_unknown = false;
        _corner     = math::ShortVector<BT, dim>( bbc );
        _dimension  = math::ShortVector<BT, dim>( bbd );
        update();
    }
    
    inline const BoundingBox<BT, dim> split( const unsigned orientatio, const BT ratio, const bool left );
    
    // compound assignment operators ===========================================
    inline BoundingBox< BT, dim >& operator+=( const math::ShortVector< BT, dim > &rhs ) {
        _corner += rhs;
        update();
        return *this;
    }

    inline BoundingBox< BT, dim >& operator-=( const math::ShortVector< BT, dim > &rhs ) {
        _corner -= rhs;
        update();
        return *this;
    }

    inline BoundingBox< BT, dim >& operator*=( const BT rhs ) {
        _center     = _corner + (BT)(.5)*_dimension;
        _dimension *= rhs;
        _corner     = _center - (BT)(.5)*_dimension;
        update();
        return *this;
    }
    
    inline BoundingBox< BT, dim >& operator=( const math::ShortVector< BT, dim > &rhs) {
        _corner      = rhs;
        _dimension   = 0.;
        update();
        return *this;
    }
    
    const void print() const;
};

template< typename BT, unsigned dim >
const void BoundingBox< BT, dim>::include(const math::ShortVector< BT, dim>& p) {
    if ( !_origin_unknown ) {
        const math::ShortVector< BT, dim > aux = p - _corner;
        for ( unsigned k = 0; k < dim; k++ ) {
            if ( aux(k) < 0. ) {
                _dimension(k) -= aux(k);
                _corner(k)    += aux(k);
            } else {
                _dimension(k) = std::max( aux(k), _dimension(k));
            }
        }
    } else {
        _corner      = p;
        _dimension   = 0.;
    }
    
    update();
}


template< typename BT, unsigned dim >
const void BoundingBox< BT, dim>::include(const BoundingBox< BT, dim>& bb) {
    if ( !_origin_unknown )
        for (unsigned i=0; i<dim; i++) {
            _corner(i)   = std::min( _corner(i)    , bb._corner(i)    );
            _dimension(i)= std::max( _dimension(i) , bb._dimension(i) );
        }
    else {
        _origin_unknown  = bb._origin_unknown;        
        _corner          = bb._corner;
        _dimension       = bb._dimension;
    }
    
    update();
}

template< typename BT, unsigned dim >
const bool BoundingBox< BT, dim>::checkInside( const math::ShortVector< BT, dim >& p ) const {
    const math::ShortVector< BT, dim > aux = p - _corner;
    for ( unsigned k = 0; k < dim; k++ )
        if ( aux(k) > _dimension(k) ) return false;
    return true;
}

template< typename BT, unsigned dim >
const void BoundingBox< BT, dim>::enlarge( const BT alpha ) {
    (*this) *= 1.+alpha;
    
    update();
}

template< typename BT, unsigned dim >
inline const BoundingBox<BT, dim> BoundingBox< BT, dim>::split( const unsigned orientatio, const BT ratio, const bool left ) {
    BoundingBox< BT, dim> bb( *this );
    
    if ( left ) {
        bb._dimension(orientatio) *=    ratio;
    } else {
        bb._corner(orientatio)    += ratio*bb._dimension(orientatio);
        bb._dimension(orientatio) *= 1.-ratio;
    }
    
    bb.update();    
    return bb;
}

template< typename BT, unsigned dim >
const BT BoundingBox< BT, dim>::maxDimension() const {
    BT res  = _dimension(0);
    for (unsigned i=1; i<dim; i++)
        res = std::max(res, _dimension(i));
    return res;
}

template< typename BT, unsigned dim >
const void BoundingBox< BT, dim>::print() const {
    std::cout << "Bounding Box: c" << _corner << "; d" << _dimension;
}

template< class Archive, typename BT, unsigned dm, unsigned dim >
void serialize(Archive & ar, BoundingBox< BT, dim >& g, const unsigned int version )
{
    using namespace boost::serialization;
    for ( unsigned k = 0; k < dim; k++) {
        ar & make_nvp( ("corner_"   + asString(k)).c_str() , g.corner(k));
    }
    for ( unsigned k = 0; k < dim; k++) {
        ar & make_nvp( ("dimension_" + asString(k)).c_str() , g.dimension(k));
    }
    
    g.update();
}


} //namespace geometry

template< typename T, unsigned N >
inline std::ostream& operator<< ( std::ostream& out, const geometry::BoundingBox<T, N>& bb ) {
    out << "lower corner    " << bb.corner()                  << std::endl;
    out << "higher corner   " << bb.corner() + bb.dimension() << std::endl;
    out << "dimension       " << bb.dimension()               << std::endl;
    out << "center          " << bb.center()                  << std::endl;

    return out;
}
