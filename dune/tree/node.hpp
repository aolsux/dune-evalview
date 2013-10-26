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
#include <unordered_map>
#include <error/duneerror.hpp>
#include <utils/utils.hpp>
#include <geometry/boundingbox.hpp>
#include <assert.h>
#include <boost/iterator/iterator_concepts.hpp>
#include <fem/dune.h>

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
        typedef Dune::FieldVector< Real, dim >      FieldVector;
        typedef geometry::BoundingBox< Real, dim >  BoundingBox;        
        typedef typename GridType::template Codim<dim>::EntitySeed          VertexSeed;
        typedef typename GridType::template Codim<dim>::EntityPointer       VertexPointer;
        typedef typename GridType::template Codim<0>::EntitySeed            EntitySeed;
        typedef typename GridType::template Codim<0>::EntityPointer         EntityPointer;
        typedef typename GridType::template Codim<0>::Entity                Entity;
        
    };

protected:
    typedef typename Traits::Real           Real;
    typedef typename Traits::LinaVector     LinaVector;
    typedef typename Traits::FieldVector    FieldVector;
    typedef typename Traits::BoundingBox    BoundingBox;
    typedef typename Traits::VertexSeed     VertexSeed;
    typedef typename Traits::VertexPointer  VertexPointer;
    typedef typename Traits::EntitySeed     EntitySeed;
    typedef typename Traits::EntityPointer  EntityPointer;
    typedef typename Traits::Entity         Entity;
    typedef typename Traits::GridType       GridType;
    typedef typename Traits::GridView       GridView;
    
    
    static constexpr unsigned dim     = Traits::dim;
    static constexpr unsigned dimw    = Traits::dimw;
    

    struct VertexContainer
    {
        std::vector<unsigned>   _entity_seeds;
        LinaVector              _global;
        VertexSeed              _seed;
        unsigned                _id;
        
        VertexContainer() :
            _entity_seeds    (  ), 
            _global          (0.), 
            _id              (0 )
        {}
          
        VertexContainer( const VertexSeed& seed ) :
            _entity_seeds    (    ), 
            _global          (0.  ), 
            _seed            (seed), 
            _id              (0   )
        {}
        
        VertexContainer( const VertexContainer& v ) :
            _entity_seeds    (v._entity_seeds     ), 
            _global          (v._global           ), 
            _seed            (v._seed             ), 
            _id              (v._id               )
        {}

        void remove_duplicates() {
            // remove duplicate entities
            std::sort( _entity_seeds.begin(), _entity_seeds.end()); 
            auto lastE = std::unique(_entity_seeds.begin(), _entity_seeds.end());
            _entity_seeds.erase(lastE, _entity_seeds.end());
        }        
    };

public: //protected:
    const Node<GridView>*           _parent;
    Node<GridView>*                 _child[2];
    std::vector< VertexContainer* > _vertices;
    const GridView&                 _gridView;
    const GridType&                 _grid;
    BoundingBox                     _bounding_box;
    LinaVector                      _normal;            //!> the normal of the plane that splits this node    
    unsigned                        _orientation;       //!> the dimension that is split by this node
    unsigned                        _level;             //!> the depth of the node in the tree
    bool                            _isLeaf;
    bool                            _isEmpty;
        
protected:
    Node() = delete;

    //Only needed for Root!
    Node( const Node<GridView>* parent, const GridView& gv) :
        _parent(parent), 
        _child({NULL, NULL}), 
        _gridView(gv), 
        _grid(_gridView.grid()), 
        _orientation(0), 
        _normal(0.), 
        _level(0), 
        _isLeaf(false), 
        _isEmpty(true)
    {
        _normal(_orientation) = 1.;
    }

    bool left (const LinaVector& p) const { return p(_orientation) < _bounding_box.center(_orientation); }
    bool right(const LinaVector& p) const { return !left( p ); }

    // vid contain indices of all vertices that reside in the space defined by boundingbox
    // p array of global space coordinates corresponding to the indices stored int vid
    // s size of p and vid TODO: remove s we dont need it
    template< class Iterator >
    void put( Iterator it_begin, Iterator it_end ) {
        _vertices.clear();
        _vertices.reserve( it_end - it_begin );
        for ( auto p = it_begin; p!= it_end; ++p) 
            _vertices.push_back(*p);
        _vertices.shrink_to_fit();
        _isEmpty = _vertices.size() < 1;
        
        // abort the recursion if there is only one vertex left within this node
        if ( _vertices.size() <= 1 ) {
            _isLeaf     = true;
            return;
        }

        split();

        std::vector< VertexContainer* > l,r;
        for ( auto vec : _vertices ) {
            if( left(vec->_global) )
                l.push_back( vec );
            else
                r.push_back( vec );
        }
        
        _child[0]->put( l.begin(), l.end() );
        _child[1]->put( r.begin(), r.end() );
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
        _grid(_gridView.grid()), 
        _level(level),
        _bounding_box(box),
        _normal(0.),
        _orientation(level%dim),
        _child( {NULL, NULL} ), 
        _isLeaf(false), 
        _isEmpty(true)
    {
        _normal( _orientation ) = 1.;        
        if ( level > 1000 ) throw GridError( "Tree depth > 1000!", __ERROR_INFO__ );
    }

    virtual ~Node() {
        release();
    }
    
    virtual void release() {
        safe_delete( _child[0] );
        safe_delete( _child[1] );
    }

    const Node*             child (const unsigned i)    const { return _child[i]; }
    const VertexContainer*  vertex(const unsigned i)    const { return _vertices[i]; }
    const unsigned          vertex_size()               const { return _vertices.size(); }
    const bool              isLeaf()                    const { return _isLeaf;     }
    const bool              isEmpty()                   const { return _isEmpty;    }
    const unsigned          level()                     const { return _level;      }
    const unsigned          orientation()               const { return _orientation;}
    const LinaVector        normal()                    const { return _normal;     }

public:
    struct TreeStats {
        unsigned numNodes;
        unsigned numLeafs;
        unsigned numVertices;
                
        unsigned minLevel;
        Real     aveLevel;
        unsigned maxLevel;
        
        unsigned minLeafLevel;
        Real     aveLeafLevel;
        unsigned maxLeafLevel;
        
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
            minLeafLevel( std::numeric_limits<unsigned>::max() ), 
            maxLeafLevel( std::numeric_limits<unsigned>::min() ), 
            aveLeafLevel( 0. ), 
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
            
            out << "Minimum Leaf Level                  " << minLeafLevel       << std::endl;
            out << "Average Leaf Level                  " << aveLeafLevel       << std::endl;
            out << "Maximum Leaf Level                  " << maxLeafLevel       << std::endl << std::endl;
            
            out << "Minimum number of Vertices per Node " << minVertices        << std::endl;
            out << "Average number of Vertices per Node " << aveVertices        << std::endl;
            out << "Maximum number of Vertices per Node " << maxVertices        << std::endl << std::endl;
            
            out << "Minimum number of Entities per Leaf " << minEntitiesPerLeaf << std::endl;
            out << "Average number of Entities per Leaf " << aveEntitiesPerLeaf << std::endl;
            out << "Maximum number of Entities per Leaf " << maxEntitiesPerLeaf << std::endl;
            
            return out;
        }
    };
    
public:
    
    virtual void fillTreeStats( TreeStats& ts ) const {
        ts.minLevel = std::min( ts.minLevel , _level );
        ts.maxLevel = std::max( ts.maxLevel , _level );
        ts.aveLevel += static_cast<Real>(_level);
        
        const unsigned vs = _vertices.size();
        ts.minVertices  = std::min( ts.minVertices , vs );
        ts.maxVertices  = std::max( ts.maxVertices , vs );
        ts.aveVertices += static_cast<Real>( vs );
        
        ts.numNodes++;
        if ( _isLeaf ) {
            ts.numLeafs++;
            
            ts.minLeafLevel = std::min( ts.minLeafLevel , _level );
            ts.maxLeafLevel = std::max( ts.maxLeafLevel , _level );
            ts.aveLeafLevel += static_cast<Real>(_level);
            
            if ( vs > 0 ) {
                assert( vs == 1 );
                const unsigned    vss  = _vertices[0]->_entity_seeds.size();
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
    
    const Node* searchDown( const LinaVector& x ) const {
        if ( _isLeaf ) return this;
            
        if ( left(x) ) return _child[0]->searchDown(x);    
        return _child[1]->searchDown(x);
    }
    
    
    struct DepthFirstResult {
        const EntitySeed es;
        const bool       found;
        
        DepthFirstResult() : es(), found(false) {}
        DepthFirstResult( const EntitySeed& es_ ) : es(es_), found(true) {}
        DepthFirstResult( const DepthFirstResult& r ) : es(r.es), found(r.found) {}
    };
    
    struct EntityContainer {
        EntitySeed  _seed;
        LinaVector  _global;
        unsigned    _id;
        
        EntityContainer( const EntitySeed& seed ) : _seed(seed), _global(0.), _id(0) {}
    };
    
    const DepthFirstResult searchUp( const FieldVector& xg, const std::vector<EntityContainer*>& _entities, const Node* caller = NULL ) const {
        const auto res = searchDown( xg, _entities, caller );
        if ( res.found ) return res;
        
        if ( _parent != NULL )
            return _parent->searchUp( xg, _entities, this );        
        
        return DepthFirstResult( );
    } 
    
        
    const DepthFirstResult searchDown( const FieldVector& xg, const std::vector<EntityContainer*>& _entities, const Node* caller = NULL ) const {
        if ( _isEmpty ) return DepthFirstResult( );
           
        if ( _isLeaf  ) {
           
            for ( auto es = vertex(0)->_entity_seeds.begin(); es != vertex(0)->_entity_seeds.end(); ++es ) {
                const EntityPointer ep( _grid.entityPointer( _entities[*es]->_seed ) );
                const Entity&   e   = *ep; 
                const auto&     geo = e.geometry();   
                const auto&     gre = Dune::GenericReferenceElements< Real, dim >::general(geo.type());
                const auto      xl  = geo.local( xg );
                if ( gre.checkInside( xl ) ) {
                    return DepthFirstResult( e.seed() );
                }
            }
           
        } else {
            if ( caller != _child[0] ) {
                const auto res0 = _child[0]->searchDown( xg, _entities, this );
                if ( res0.found ) return res0;
            }
            if ( caller != _child[1] ) {
                const auto res1 = _child[1]->searchDown( xg, _entities, this );
                if ( res1.found ) return res1;
            }
        }        
        
        return DepthFirstResult( );
    }
};


}
