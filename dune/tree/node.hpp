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

#pragma once

#include <limits>
#include <iostream>
#include <utils/utils.hpp>
#include <geometry/boundingbox.hpp>
#include <assert.h>
#include <boost/iterator/iterator_concepts.hpp>

namespace tree {



template< class GV >
class Node {
public:
    struct Traits {
        typedef GV                                  GridView;
        typedef typename GridView::Grid             GridType;
        
        static constexpr unsigned                   dim         = GridView::dimension;
        static constexpr unsigned                   dimw        = GridView::dimensionworld;               
                
        typedef typename GridView::ctype            Real;
        typedef math::ShortVector< Real, dim >      LinaVector;
        typedef geometry::BoundingBox< Real, dim >  BoundingBox;        
        typedef typename GridType::template Codim<0>::EntitySeed    EntitySeed;
        typedef typename GridType::template Codim<0>::Entity        Entity;
    };

protected:
    typedef typename Traits::Real           Real;
    typedef typename Traits::LinaVector     LinaVector;
    typedef typename Traits::BoundingBox    BoundingBox;
    typedef typename Traits::EntitySeed     EntitySeed;
    typedef typename Traits::Entity         Entity;
    typedef typename Traits::GridType       GridType;
    typedef typename Traits::GridView       GridView;
    
    
    static constexpr unsigned dim     = Traits::dim;
    static constexpr unsigned dimw    = Traits::dimw;
    

    struct Vertex
    {
        std::vector<const EntitySeed*>   _entity_seed;
        LinaVector                       _global;
        unsigned                         _id;
        unsigned                         _idx;
        
        Vertex() :
            _global( 0. ), 
            _id    ( 0  ), 
            _idx   ( 0  ) {}
            
        
        Vertex( const Vertex& v ) :
            _entity_seed (v._entity_seed ), 
            _global      ( v._global ), 
            _id          ( v._id  ), 
            _idx         ( v._idx  ) {}

    };

protected:
    const Node<GridView>*           _parent;
    Node<GridView>*                 _child[2];
    std::vector< Vertex* >          _vertex;
    const GridView&                 _gridView;
    BoundingBox                     _bounding_box;
    LinaVector                      _normal;            //!> the normal of the plane that splits this node    
    unsigned                        _orientation;       //!> the dimension that is split by this node
    unsigned                        _level;             //!> the depth of the node in the tree
    bool                            _isLeaf;
    
protected:
    Node() = delete;

    //Only needed for Root!
    Node( const Node<GridView>* parent, const GridView& gv) :
        _parent(parent), 
        _child({NULL, NULL}), 
        _gridView(gv), 
        _orientation(0), 
        _normal(0.), 
        _level(0), 
        _isLeaf(false) 
    {
        _normal(_orientation) = 1.;
    }

    bool left(const LinaVector& p)  const { return p(_orientation) < _bounding_box.center(_orientation); }
    bool right(const LinaVector& p) const { return !left( p ); }

    // vid contain indices of all vertices that reside in the space defined by boundingbox
    // p array of global space coordinates corresponding to the indices stored int vid
    // s size of p and vid TODO: remove s we dont need it
    template< class Iterator >
    void put( Iterator it_begin, Iterator it_end ) {
        _vertex.clear();

        for ( auto p = it_begin; p!= it_end; ++p) 
            _vertex.push_back(*p);

        // abort the recursion if there is only one vertex left within this node
        if ( _vertex.size() <= 1 ) {
            _isLeaf = true;
            return;
        }

        split();

        std::vector< Vertex* > l,r;
        for ( auto vec : _vertex )
        {
            if( left(vec->_global) )
                l.push_back( vec );
            else
                r.push_back( vec );
        }

        if ( l.size() > 0 ) 
            _child[0]->put( l.begin(), l.end() );
        else 
            _child[0]->_isLeaf = true;
            
        if ( r.size() > 0 ) 
            _child[1]->put( r.begin(), r.end() );
        else 
            _child[1]->_isLeaf = true;

    }

    void split() {
        assert( _child[0] == NULL );
        assert( _child[1] == NULL );
        // construct the two childs
        _child[0] = new Node( this, _bounding_box.split(_orientation,true) , _level+1 );
        _child[1] = new Node( this, _bounding_box.split(_orientation,false), _level+1 );
        _isLeaf   = false;
    }

public:

    Node( const Node<GridView>& node ) = delete;
    Node& operator = ( const Node<GridView>& node ) = delete;

    Node( const Node<GridView>* parent, const BoundingBox& box, const unsigned level) :
        _parent(parent),
        _gridView(parent->_gridView), 
        _level(level),
        _bounding_box(box),
        _normal(0.),
        _orientation(level%dim),
        _child( {NULL, NULL} )
    {
        _normal( _orientation ) = 1.;
//         std::cout << "level         "   << _level << std::endl;
//         std::cout << "orientation   "   << _orientation << std::endl;
//         std::cout << "boundingbox\n"; _bounding_box.operator<<(std::cout);
//         std::cout << "normal        "   << _normal << std::endl;
        
        if ( level > 1000 ) throw;
    }

    virtual ~Node() {
        safe_delete( _child[0] );
        safe_delete( _child[1] );
    }

    const Node*         child(const unsigned i) const { assert( i < 2 ); return _child[i];   }
    const bool          isLeaf()                const { return _isLeaf;     }
    const unsigned      level()                 const { return _level;      }
    const unsigned      orientation()           const { return _orientation;}
    const LinaVector    normal()                const { return _normal;     }


    // iterate over all entities of the node
    std::vector<const Entity&> entities() const {}


public:
    struct TreeStats {
        unsigned numNodes;
        unsigned numLeafs;
        unsigned numVertices;
                
        unsigned minLevel;
        Real     aveLevel;
        unsigned maxLevel;
        
        unsigned minVertices;
        Real     aveVertices;
        unsigned maxVertices;
        
        unsigned minEntitiesPerLeaf;
        Real     aveEntitiesPerLeaf;
        unsigned maxEntitiesPerLeaf;
        
        TreeStats() : 
            numNodes( 0 ), 
            numLeafs( 0 ),         
            numVertices( 0 ), 
            minLevel( std::numeric_limits<unsigned>::max() ), 
            maxLevel( std::numeric_limits<unsigned>::min() ), 
            aveLevel( 0. ), 
            minVertices( std::numeric_limits<unsigned>::max() ), 
            maxVertices( std::numeric_limits<unsigned>::min() ), 
            aveVertices( 0. ), 
            minEntitiesPerLeaf( std::numeric_limits<unsigned>::max() ), 
            maxEntitiesPerLeaf( std::numeric_limits<unsigned>::min() ), 
            aveEntitiesPerLeaf( 0. ) {}
            
        std::ostream& operator<< ( std::ostream& out ) const {
            
            out << "Number of Nodes                     " << numNodes           << std::endl;
            out << "Number of Leafs                     " << numLeafs           << std::endl;
            out << "Number of Vertices                  " << numVertices        << std::endl << std::endl;
                    
            out << "Minimum Level                       " << minLevel           << std::endl;
            out << "Average Level                       " << aveLevel           << std::endl;
            out << "Maximum Level                       " << maxLevel           << std::endl << std::endl;
            
            out << "Minimum number of Vertices per Node " << minVertices        << std::endl;
            out << "Average number of Vertices per Node " << aveVertices        << std::endl;
            out << "Maximum number of Vertices per Node " << maxVertices        << std::endl << std::endl;
            
            out << "Minimum number of Entities per Leaf " << minEntitiesPerLeaf << std::endl;
            out << "Average number of Entities per Leaf " << aveEntitiesPerLeaf << std::endl;
            out << "Maximum number of Entities per Leaf " << maxEntitiesPerLeaf << std::endl;
            
            return out;
        }
    };
    
protected:
    
    virtual void fillTreeStats( TreeStats& ts ) const {
        ts.minLevel = std::min( ts.minLevel , _level );
        ts.maxLevel = std::max( ts.maxLevel , _level );
        ts.aveLevel += static_cast<Real>(_level);
        
        const unsigned vs = _vertex.size();
        ts.minVertices  = std::min( ts.minVertices , vs );
        ts.maxVertices  = std::max( ts.maxVertices , vs );
        ts.aveVertices += static_cast<Real>( vs );
        
        ts.numNodes++;
        if ( _isLeaf ) {
            ts.numLeafs++;
            
            if ( vs > 0 ) {
                assert( vs == 1 );
                const unsigned    vss  = _vertex[0]->_entity_seed.size();
                ts.minEntitiesPerLeaf  = std::min( ts.minEntitiesPerLeaf , vss );
                ts.maxEntitiesPerLeaf  = std::max( ts.maxEntitiesPerLeaf , vss );
                ts.aveEntitiesPerLeaf += static_cast<Real>( vss );
            }
        } else {
            assert( _child[0] != NULL );
            assert( _child[1] != NULL );
            _child[0]->fillTreeStats( ts );
            _child[1]->fillTreeStats( ts );
        }
    }
    

};


}
