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
#include <algorithm>
#include <error/duneerror.hpp>
#include <utils/utils.hpp>
#include <geometry/boundingbox.hpp>
#include <fem/helper.hpp>
#include <assert.h>

namespace tree {


template< class GV >
class Node {
//=======================================================================================================
// public traits
//=======================================================================================================
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

    
//=======================================================================================================
// protected data
//=======================================================================================================
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

    struct EntityContainer {
        EntitySeed                      _seed;
        std::vector<unsigned>           _neighbour_seeds;
        geometry::BoundingBox<Real,dim> _bb;
        LinaVector                      _global;
        unsigned                        _id;
        unsigned                        _idx;

        EntityContainer() :
            _seed           (),
            _neighbour_seeds(),
            _bb             (), 
            _global         (0.),
            _id             (0), 
            _id             (0)
        {}
        
        EntityContainer( const EntitySeed& seed ) : _seed(seed), _neighbour_seeds(), _bb(), _global(0.), _id(0), _idx(0) {}
        
        EntityContainer( const EntityContainer& ec ) :
            _seed           (ec._seed             ),
            _neighbour_seeds(ec._neighbour_seeds  ),
            _bb             (ec._bb               ), 
            _global         (ec._global           ),
            _id             (ec._id               ), 
            _idx            (ec._idx               )
        {}
        
        void remove_duplicates( const std::vector< EntityContainer* >& entities ) {
            // remove duplicate entities
            struct Comp {
                const std::vector< EntityContainer* >& entities;
                const LinaVector&                      g;
                const bool operator ()(const unsigned a, const unsigned b) const {
                    return math::norm2( entities[a]->_global - g ) < math::norm2( entities[b]->_global - g ); 
                };
                Comp( const std::vector< EntityContainer* >& entities_, const LinaVector& g_ ) : entities(entities_), g(g_) {}
            };            
            Comp comp( entities, _global );
            
            std::sort( _neighbour_seeds.begin(), _neighbour_seeds.end(), comp );
            auto lastE = std::unique(_neighbour_seeds.begin(), _neighbour_seeds.end());
            _neighbour_seeds.erase(lastE, _neighbour_seeds.end());
            for ( auto it = _neighbour_seeds.begin(); it != _neighbour_seeds.end(); ++it )
                if ( *it == _idx ) { _neighbour_seeds.erase(it); break; }
        }
    };
    
    Node<GridView>*                 _parent;
    Node<GridView>*                 _child[2];
    Real                            _median;
    std::vector< EntityContainer* > _entities;
    const GridView&                 _gridView;
    const GridType&                 _grid;
    BoundingBox                     _bounding_box;
    unsigned                        _orientation;       //!< the dimension that is split by this node
    unsigned                        _level;             //!< the depth of the node in the tree
    bool                            _isLeaf;
    bool                            _isEmpty;
    bool                            _balanced;
    int                             _balance_factor;

    
//=======================================================================================================
// protected methods
//=======================================================================================================
protected:    
    //== constructot / destructor =======================================================================
    //Only needed for Root!
    Node( Node<GridView>* parent, const GridView& gv, const bool bal = false ) :
        _parent(parent),
        _child{NULL, NULL},
        _median(0.),
        _entities(), 
        _gridView(gv),
        _grid(_gridView.grid()),
        _bounding_box(),
        _orientation(0),
        _level(0),        
        _isLeaf(false),
        _isEmpty(true),
        _balanced( bal ), 
        _balance_factor(0)
    { }
    
    //Only needed for Nodes themselfs!
    Node( Node<GridView>* parent, const BoundingBox& box, const unsigned level, const unsigned ori, const bool bal ) :
        _parent(parent),
        _child{NULL, NULL},
        _median(0.),
        _entities(), 
        _gridView(parent->_gridView),
        _grid(_gridView.grid()),
        _bounding_box(box),
        _orientation(ori%dim),
        _level(level),
        _isLeaf(false),
        _isEmpty(true),
        _balanced(bal), 
        _balance_factor(0)
    {
        if ( level > 1000 ) throw GridError( "Tree depth > 1000!", __ERROR_INFO__ );
    }
    
    //== build tree =====================================================================================
    inline const bool left (const LinaVector& p) const { return p(_orientation) < _median; }
    inline const bool right(const LinaVector& p) const { return !left( p ); }
    
    inline const bool left (const LinaVector& p, const unsigned ori ) const { return p(ori) < _median; }
    inline const bool right(const LinaVector& p, const unsigned ori ) const { return !left( p, ori ); }

    // vid contain indices of all vertices that reside in the space defined by boundingbox
    // p array of global space coordinates corresponding to the indices stored int vid
    // s size of p and vid TODO: remove s we dont need it
    template< class Iterator >
    void put( Iterator it_begin, Iterator it_end ) {
        _entities.clear();
        _entities.reserve( it_end - it_begin );
        for ( auto p = it_begin; p!= it_end; ++p)
            _entities.push_back(*p);
        _entities.shrink_to_fit();

        _isEmpty = _entities.size() <  1;
        _isLeaf  = _entities.size() == 1;
        // abort the recursion if there is only one vertex left within this node
        if ( _isLeaf || _isEmpty ) return;

        _median = _bounding_box.corner(_orientation) + .5*_bounding_box.dimension(_orientation);
        std::vector< EntityContainer* > l,r;
        for ( auto vec : _entities ) {
            if( left(vec->_global) )
                l.push_back( vec );
            else
                r.push_back( vec );
        }

        split( .5 );
        _child[0]->put( l.begin(), l.end() );
        _child[1]->put( r.begin(), r.end() );
    }
    
    template< class Iterator >
    void reput( Iterator it_begin, Iterator it_end ) {
        _entities.clear();
        _entities.reserve( it_end - it_begin );
        for ( auto p = it_begin; p!= it_end; ++p)
            _entities.push_back(*p);
        _entities.shrink_to_fit();
        
        _isEmpty = _entities.size() <  1;
        _isLeaf  = _entities.size() == 1;
        // abort the recursion if there is only one vertex left within this node
        if ( _isLeaf || _isEmpty ) return;
            
        _median = _bounding_box.corner(_orientation) + .5*_bounding_box.dimension(_orientation);
        std::vector< VertexContainer* > l,r;
        for ( auto vec : _entities ) {
            if( left(vec->_global) )
                l.push_back( vec );
            else
                r.push_back( vec );
        }

        if (_child[0]) _child[0]->reput( l.begin(), l.end() );
        if (_child[1]) _child[1]->reput( r.begin(), r.end() );
    }
    
    void reput( ) {
        _entities.shrink_to_fit();
        
        _isEmpty = _entities.size() <  1;
        _isLeaf  = _entities.size() == 1;
        // abort the recursion if there is only one vertex left within this node
        if ( _isLeaf || _isEmpty ) return;
            
        _median = _bounding_box.corner(_orientation) + .5*_bounding_box.dimension(_orientation);
        std::vector< VertexContainer* > l,r;
        for ( auto vec : _entities ) {
            if( left(vec->_global) )
                l.push_back( vec );
            else
                r.push_back( vec );
        }

        if (_child[0]) _child[0]->reput( l.begin(), l.end() );
        if (_child[1]) _child[1]->reput( r.begin(), r.end() );
    }

    void split( const Real ratio ) {
        assert( _child[0] == NULL );
        assert( _child[1] == NULL );
        // construct the two childs
        _child[0] = new Node( this, _bounding_box.split(_orientation, ratio, true) , _level+1, _orientation+1, _balanced );
        _child[1] = new Node( this, _bounding_box.split(_orientation, ratio, false), _level+1, _orientation+1, _balanced );
        _isLeaf   = false;
    }
    
    //== balance tree ===================================================================================
    void leftRotate() {
        throw NotImplemented(__ERROR_INFO__);
        if ( _parent   == NULL ) return;                                    // can't rotate root
        if ( _child[0] == NULL ) return;                                    // can't rotate leafs
        if ( _child[1] == NULL ) return;                                    // can't rotate leafs
        
        const unsigned c    = ( _parent->_child[0] == this ) ? 0 : 1;       // which child am I? 
        Node* anchor        = _parent;                                      // backup parent
        _parent             = _child[1];                                    // right child to parent
        anchor->_child[c]   = _child[1];                                    // 
        _child[1]           = _parent->_child[0];                           // take over new parents left child
        _child[1]->_parent  = this;
        _parent->_child[0]  = this;                                         // make me new parents left child
        
        anchor->update();                                                   // update sub-tree
    }
    
    void rightRotate() {
        throw NotImplemented(__ERROR_INFO__);
        if ( _parent   == NULL ) return;                                    // can't rotate root
        if ( _child[0] == NULL ) return;                                    // can't rotate leafs
        if ( _child[1] == NULL ) return;                                    // can't rotate leafs
        
        const unsigned c    = ( _parent->_child[0] == this ) ? 0 : 1;       // which child am I? 
        Node* anchor        = _parent;                                      // backup parent
        _parent             = _child[0];                                    // left child to parent
        anchor->_child[c]   = _child[0];                                    // 
        _child[0]           = _parent->_child[1];                           // take over new parents left child
        _child[0]->_parent  = this;
        _parent->_child[1]  = this;                                         // make me new parents left child
        
        anchor->update();                                                   // update sub-tree
    }
    
    void balance() {
        throw NotImplemented(__ERROR_INFO__);
        if (_child[0]) _child[0]->balance();
        if (_child[1]) _child[1]->balance();
        
        if (_isLeaf) return;
        
        for ( unsigned c = 0; c < 2; c++ )
        if (_child[c]) {
            while ( (_child[c]->_balance_factor*_child[c]->_balance_factor > 1) ) {
                if ( _child[c]->_balance_factor < -1 ) {
                    _child[c]->leftRotate(); 
                    std::cout << "level a" << _level << ",   c " << c << ",   bf " <<  _child[c]->_balance_factor << std::endl;
                    return;
                }
                
                if ( _child[c]->_balance_factor > +1 ) {
                    _child[c]->rightRotate(); 
                    std::cout << "level b" << _level << ",   c " << c << ",   bf " <<  _child[c]->_balance_factor << std::endl;
                    return;
                }
            }
        }
    }
    
    //== update state machine ===========================================================================
    const int updateBalanceFactor() {
        const int l = _child[0] ? _child[0]->updateBalanceFactor() : 0;
        const int r = _child[1] ? _child[1]->updateBalanceFactor() : 0;
        
        _balance_factor = l-r;
        
//         std::cout << "level " << _level << ",   l " <<  l << ",   r " <<  r << ",   bf " <<  _balance_factor << std::endl;
        
        return std::max( l, r ) + 1;
    }
    
    void updateState( const unsigned lv = 0 ) {
        _isEmpty = _entities.size() < 1;
        _isLeaf  = ((_child[0] == NULL) && (_child[1] == NULL));
        _level   = lv;
        
        if (_child[0]) _child[0]->updateState( lv + 1 );
        if (_child[1]) _child[1]->updateState( lv + 1 );
    }
    
    void updateBoundingBox() {
        if ( _isLeaf ) return;
            
        const Real ratio = (_median - _bounding_box.corner(_orientation))/_bounding_box.dimension(_orientation);
        _child[0]->_bounding_box = _bounding_box.split(_orientation, ratio, true);
        _child[1]->_bounding_box = _bounding_box.split(_orientation, ratio, false);
        
        _child[0]->updateBoundingBox();
        _child[1]->updateBoundingBox();
    }
    
    //== optimize for size ==============================================================================
    void deleteEmpty() {
        if ( _child[0] ) {
            _child[0]->deleteEmpty();
            if ( _child[0]->isEmpty() ) 
                safe_delete(_child[0]);
        }
        
        if ( _child[1] ) {
            _child[1]->deleteEmpty();
            if ( _child[1]->isEmpty() ) 
                safe_delete(_child[1]);
        }
    }
    
    void removeSingles() {
        if ( _child[0] ) _child[0]->removeSingles();
        if ( _child[1] ) _child[1]->removeSingles();
        
        if ( _child[0] ) {
            if ( (!(_child[0]->_child[0])) && (_child[0]->_child[1]) ) {
                auto aux    = _child[0];
                _child[0]   = aux->_child[1];
                aux->_child[1] = NULL;
                safe_delete( aux );
            } 
            if ( (!(_child[0]->_child[1])) && (_child[0]->_child[0]) ) {
                auto aux    = _child[0];
                _child[0]   = aux->_child[0];
                aux->_child[0] = NULL;
                safe_delete( aux );
            }
        }
        
        if ( _child[1] ) {
            if ( (!(_child[1]->_child[0])) && (_child[1]->_child[1]) ) {
                auto aux    = _child[1];
                _child[1]   = aux->_child[1];
                aux->_child[1] = NULL;
                safe_delete( aux );
            } 
            if ( (!(_child[1]->_child[1])) && (_child[1]->_child[0]) ) {
                auto aux    = _child[1];
                _child[1]   = aux->_child[0];
                aux->_child[0] = NULL;
                safe_delete( aux );
            }
        }
    }
    
//=======================================================================================================
// public data
//=======================================================================================================
public:
    struct TreeStats {
        unsigned depth;
        
        unsigned numNodes;
        unsigned numLeafs;
        unsigned numVertices;
        unsigned numEmpty;
        unsigned numBadChildren;

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
            depth( 0 ),
            numNodes( 0 ),
            numLeafs( 0 ),
            numVertices( 0 ),
            numEmpty( 0 ),
            numBadChildren( 0 ),
            minLevel( std::numeric_limits<unsigned>::max() ),
            aveLevel( 0. ),
            maxLevel( std::numeric_limits<unsigned>::min() ),
            minLeafLevel( std::numeric_limits<unsigned>::max() ),
            aveLeafLevel( 0. ),
            maxLeafLevel( std::numeric_limits<unsigned>::min() ),
            minVertices( std::numeric_limits<unsigned>::max() ),
            aveVertices( 0. ),
            maxVertices( std::numeric_limits<unsigned>::min() ),
            minEntitiesPerLeaf( std::numeric_limits<unsigned>::max() ),
            aveEntitiesPerLeaf( 0. ), 
            maxEntitiesPerLeaf( std::numeric_limits<unsigned>::min() )
	{}

        std::ostream& operator<< ( std::ostream& out ) const {
            out << "Depth                               " << depth              << std::endl << std::endl;
        
            out << "Number of Nodes                     " << numNodes           << std::endl;
            out << "Number of Leafs                     " << numLeafs           << std::endl;
            out << "Number of Vertices                  " << numVertices        << std::endl ;
            out << "Number of empty Nodes               " << numEmpty           << std::endl ;
            out << "Number of bad Children              " << numBadChildren     << std::endl << std::endl;

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
    
    struct DepthFirstResult {
        const EntitySeed    es;
        const FieldVector   xl;
        const bool          found;

        DepthFirstResult( const DepthFirstResult&  ) = default;
        DepthFirstResult(       DepthFirstResult&& ) = default;        
        DepthFirstResult() : es(), xl(0.), found(false) {}
        DepthFirstResult( const EntitySeed& es_, const FieldVector& xl_ ) : es(es_), xl(xl_), found(true) {}
    };

    inline const Node*             child (const unsigned i)    const { return _child[i];    }
    inline const EntityContainer*  entity(const unsigned i)    const { return _entities[i]; }
    inline const unsigned          entity_size()               const { return _entities.size(); }
    inline const bool              isLeaf()                    const { return _isLeaf;      }
    inline const bool              isEmpty()                   const { return _isEmpty;     }
    inline const bool              balanced()                  const { return _balanced;    }
    inline const unsigned          level()                     const { return _level;       }
    inline const unsigned          orientation()               const { return _orientation; }
    
//=======================================================================================================
// public methods
//=======================================================================================================
public:
    //== constructot / destructor =======================================================================
    Node() = delete;
    Node( const Node<GridView>& node ) = delete;
    Node& operator = ( const Node<GridView>& node ) = delete;
    
    virtual ~Node() {
        release();
    }

    virtual void release() {
        safe_delete( _child[0] );
        safe_delete( _child[1] );
    }
    
    //== update state machine  ==========================================================================
    void update() {
        updateState();
        updateBoundingBox();
        updateBalanceFactor();
    }
    
    //== information on tree ============================================================================
    virtual void fillTreeStats( TreeStats& ts ) const {
        ts.minLevel = std::min( ts.minLevel , _level );
        ts.maxLevel = std::max( ts.maxLevel , _level );
        ts.aveLevel += static_cast<Real>(_level);

        const unsigned vs = _entities.size();
        ts.minVertices  = std::min( ts.minVertices , vs );
        ts.maxVertices  = std::max( ts.maxVertices , vs );
        ts.aveVertices += static_cast<Real>( vs );

        ts.numNodes++;
        
        if (_isEmpty ) ts.numEmpty++;
            
        if ( (_child[0]==NULL) && (_child[1]==NULL) && !_isLeaf ) ts.numBadChildren++;
        if ( (_child[0]==NULL) ^  (_child[1]==NULL)             ) ts.numBadChildren++;
        
        if ( _isLeaf ) {
            ts.numLeafs++;

            ts.minLeafLevel = std::min( ts.minLeafLevel , _level );
            ts.maxLeafLevel = std::max( ts.maxLeafLevel , _level );
            ts.aveLeafLevel += static_cast<Real>(_level);

            if ( vs > 0 ) {
//                 assert( vs == 1 );
                const unsigned    vss  = _entities[0]->_neighbour_seeds.size();
                ts.minEntitiesPerLeaf  = std::min( ts.minEntitiesPerLeaf , vss );
                ts.maxEntitiesPerLeaf  = std::max( ts.maxEntitiesPerLeaf , vss );
                ts.aveEntitiesPerLeaf += static_cast<Real>( vss );
            }
        } else {
            if (_child[0]) _child[0]->fillTreeStats( ts );
            if (_child[1]) _child[1]->fillTreeStats( ts );
        }
    }

    
    //== search / iterate tree  =========================================================================
    const Node* searchDown( const LinaVector& x ) const {
        if ( _isLeaf ) return this;

        if ( left(x) ) return _child[0]->searchDown(x);
        return _child[1]->searchDown(x);
    }

    const DepthFirstResult  searchUp( const FieldVector& xg, const std::vector<EntityContainer*>& entities, const Node* caller = NULL ) const {
        const auto res = searchDown( xg, entities, caller );
        if ( res.found ) return res;

        if ( _parent != NULL )
            return _parent->searchUp( xg, entities, this );

        return DepthFirstResult( );
    }

    const DepthFirstResult  searchDown( const FieldVector& xg, const std::vector<EntityContainer*>& entities, const Node* caller = NULL ) const {
        if ( _isEmpty ) return DepthFirstResult( );

        if ( _isLeaf  ) {
            const LinaVector x = fem::asShortVector<Real, dim>(xg);

            {
                const EntityPointer ep( _grid.entityPointer( entity(0)->_seed ) );
                const Entity&   e   = *ep;
                const auto&     geo = e.geometry();
                const auto&     gre = Dune::GenericReferenceElements< Real, dim >::general(geo.type());
                const auto      xl  = geo.local( xg );
                if ( gre.checkInside( xl ) ) {
                    return DepthFirstResult( e.seed(), xl );
                }
            }                
                
                
            {
                for ( auto es = entity(0)->_neighbour_seeds.begin(); es != entity(0)->_neighbour_seeds.end(); ++es ) {
                    if ( !entities[*es]->_bb.inside(x) ) continue;
                    const EntityPointer ep( _grid.entityPointer( entities[*es]->_seed ) );
                    const Entity&   e   = *ep;
                    const auto&     geo = e.geometry();
                    const auto&     gre = Dune::GenericReferenceElements< Real, dim >::general(geo.type());
                    const auto      xl  = geo.local( xg );
                    if ( gre.checkInside( xl ) ) {
                        return DepthFirstResult( e.seed(), xl );
                    }
                }
            }

        } else {
            if ( (caller != _child[0]) && _child[0] ) {
                const auto res0 = _child[0]->searchDown( xg, entities, this );
                if ( res0.found ) return res0;
            }
            if ( (caller != _child[1]) && _child[1] ) {
                const auto res1 = _child[1]->searchDown( xg, entities, this );
                if ( res1.found ) return res1;
            }
        }

        return DepthFirstResult( );
    }
    
    
};


}
