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
#include <vector>
#include <unordered_map>

#include <fem/helper.hpp>
#include <tree/node.hpp>
#include <tree/leafview.hpp>
#include <tree/levelview.hpp>
#include <error/duneerror.hpp>



namespace tree {


template< class GV >
class Root : public Node<GV> {
public:
    typedef typename Node<GV>::Traits Traits;

protected:
    using Node<GV>::_parent;
    using Node<GV>::_gridView;
    using Node<GV>::_grid;
    using Node<GV>::_vertices;
    using Node<GV>::_bounding_box;
    using Node<GV>::split;
    using Node<GV>::put;
    using Node<GV>::searchDown;
    using Node<GV>::searchUp;


    typedef typename Node<GV>::EntityContainer  EntityContainer;
    typedef typename Node<GV>::VertexContainer  VertexContainer;
    typedef typename Node<GV>::DepthFirstResult DepthFirstResult;
    typedef typename Traits::Real               Real;
    typedef typename Traits::Entity             Entity;
    typedef typename Traits::EntitySeed         EntitySeed;
    typedef typename Traits::EntityPointer      EntityPointer;
    typedef typename Traits::VertexSeed         VertexSeed;
    typedef typename Traits::VertexPointer      VertexPointer;
    typedef typename Traits::GridView           GridView;
    typedef typename Traits::GridType           GridType;
    typedef typename Traits::LinaVector         LinaVector;
    typedef typename Traits::FieldVector        FieldVector;

    static constexpr unsigned dim     = Traits::dim;
    static constexpr unsigned dimw    = Traits::dimw;

public:
    std::map< unsigned, unsigned > id2idxEntity;
    std::map< unsigned, unsigned > id2idxVertex;

    std::vector<EntityContainer*>                                       _entities;
//     Dune::HierarchicSearch< GridType, typename GridType::LeafIndexSet > _hr_locator;

public:
    Root( const Root<GridView>& root ) {};

    virtual ~Root( ) {
        release();
    };

    Root( const GridView& gridview, const bool bal = false ) :
        Node<GV>(NULL,gridview, bal)
//         _hr_locator( gridview.grid(), gridview.grid().leafIndexSet() )
    {
        build();
    }


    virtual void release() {
        Node<GV>::release();
        for ( auto e : _entities )
            safe_delete( e );
        for ( auto v : _vertices )
            safe_delete( v );
        _entities.clear();
        _vertices.clear();
    }

    void build() {
        std::vector< VertexContainer* > _l_vertices;

        const auto& idSet = _grid.globalIdSet();

        // collect cells on leaf view
        for( auto e = _gridView.template begin<0>(); e != _gridView.template end<0>(); ++e ) {
            _entities.push_back( new EntityContainer(e->seed()) );
            _entities.back()->_id = idSet.id(*e);
            id2idxEntity[idSet.id(*e)] = _entities.size()-1;
        }

        // collect vertices on leaf view
        for( auto e = _gridView.template begin<dim>(); e != _gridView.template end<dim>(); ++e ) {
            _l_vertices.push_back( new VertexContainer(e->seed()) );
            _l_vertices.back()->_id = idSet.id(*e);
            id2idxVertex[idSet.id(*e)] = _l_vertices.size()-1;
        }

        // fill container of all entity seeds
        for( auto e = _gridView.template begin<0>(); e != _gridView.template end<0>(); ++e ) {
            const unsigned idx = id2idxEntity[ idSet.id(*e) ];
            const auto&    geo = e->geometry();
            const auto&    gre = Dune::GenericReferenceElements< Real, dim >::general(geo.type());

            const unsigned v_size = (unsigned)gre.size(dim);

            for ( unsigned k = 0; k < v_size; k++ ) {
                const auto& pc = e->template subEntity<dim>(k);
                const auto& c  = *pc;
                typename Traits::LinaVector gl = fem::asShortVector<Real, dim>( geo.global( gre.position(k,dim) ) ) ;

                VertexContainer* _v = _l_vertices[ id2idxVertex[ idSet.id(c) ] ];

                // store global coordinates of all vertices
                _bounding_box.append(gl);
                _entities[idx]->_bb.append(gl);
                _v->_global = gl;
                _v->_entity_seeds.push_back( idx );
            }
        }

        // generate list of vertices
        this->put( _l_vertices.begin(), _l_vertices.end() );
    }

    void rebuild() {
        release();
        build();
    }

     // iterate over all leafs of the node
    LeafView<GridView> leafView() const {
        return LeafView<GridView>( *this );
    }

     // iterate over all leafs of the node
    LevelView<GridView> levelView(unsigned level) const {
        return LevelView<GridView>( *this, level );
    }


    virtual void fillTreeStats( typename Node<GridView>::TreeStats& ts ) const {
        Node<GV>::fillTreeStats(ts);

        ts.numVertices         = _vertices.size();
        ts.aveLevel           /= static_cast<Real>( ts.numNodes );
        ts.aveLeafLevel       /= static_cast<Real>( ts.numLeafs );
        ts.aveVertices        /= static_cast<Real>( ts.numNodes );
        ts.aveEntitiesPerLeaf /= static_cast<Real>( ts.numLeafs );
    }

    void printTreeStats( std::ostream& out ) const {
        typename Node<GridView>::TreeStats ts;
        fillTreeStats(ts);
        ts.operator<<(out) << std::endl;
    }


    struct EntityData {
        const EntityPointer                 pointer;
        const Entity&                       entity;
        const FieldVector                   xl;

        EntityData( const EntityPointer pointer_,
                    const Entity&       entity_,
                    const FieldVector   xl_  ) : pointer(pointer_),  entity(entity_), xl(xl_) {}
    };

    const EntityData findEntity( const LinaVector& x )  {
        // find node containing all possible cells
        const Node<GridView>* node = searchDown( x );
        const auto fx  = fem::asFieldVector(x);
        const auto res = node->searchUp( fx, _entities, node );

        if ( res.found ) {
            const auto      ep  = _grid.entityPointer( res.es );
            const Entity&   e   = *ep;
            const auto&     geo = e.geometry();
            const auto&     gre = Dune::GenericReferenceElements< Real, dim >::general(geo.type());
            const auto      xl  = geo.local( fem::asFieldVector(x));
            return EntityData( ep, e, xl );
        }

//         const EntityPointer ep( _hr_locator.findEntity( fem::asFieldVector(x) ) );
//         const Entity&   e   = *ep;
//         const auto&     geo = e.geometry();
//         const auto&     gre = Dune::GenericReferenceElements< Real, dim >::general(geo.type());
//         const auto      xl  = geo.local( fem::asFieldVector(x));
//         return EntityData( ep, e, xl );
        // eff              100%
        // kd-tree          9.8000e-01 111.8x speed-up
        // hr-loc           1.0960e+02
        // pnt-loc          1.0906e+02

        throw GridError( "Global coordinates are outside the grid!", __ERROR_INFO__ );
    }
};


}
