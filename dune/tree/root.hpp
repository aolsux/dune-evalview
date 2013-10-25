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

#include <vector>
#include <map>

#include <tree/node.hpp>
#include <tree/leafview.hpp>
#include <tree/levelview.hpp>



namespace tree {


template< class GV >
class Root : public Node<GV> {
public:
    typedef typename Node<GV>::Traits Traits;
    
protected:
    using Node<GV>::_parent;
    using Node<GV>::_gridView;
    using Node<GV>::_vertex;
    using Node<GV>::_bounding_box;
    using Node<GV>::split;
    using Node<GV>::put;

    
    typedef typename Node<GV>::Vertex           Vertex;
    typedef typename Traits::Real               Real;
    typedef typename Traits::Entity             Entity;
    typedef typename Traits::EntitySeed         EntitySeed;
    typedef typename Traits::GridView           GridView;
    
    static constexpr unsigned dim     = Traits::dim;
    static constexpr unsigned dimw    = Traits::dimw;
    
protected:
    std::vector<EntitySeed>     _entities;

    // map each vertex id to its corresponding entity index. where entity index is the
    // index of the entity seed for the container above.
    std::map< unsigned, std::vector<unsigned> > _mapping;

public:
    Root( const Root<GridView>& root ) {};
    
    virtual ~Root( ) {
        for ( auto v : _vertex )
            safe_delete( v );
        _vertex.clear();
    };

    Root( const GridView& gridview ) :
        Node<GV>(NULL,gridview)
    {
        std::vector< Vertex* > _l_vertex;
        // create container of all entity seeds
        for( auto e = gridview.template begin<0>(); e != gridview.template end<0>(); ++e ) {
            _entities.push_back( e->seed() );
            auto geo = e->geometry();
            
            auto& gidSet    = _gridView.grid().globalIdSet();
            auto& gidxSet   = _gridView.grid().leafIndexSet();
            
            for ( unsigned k = 0; k < geo.corners(); k++ ) {
                Vertex* _v = new Vertex;

                // store global coordinates of all vertices
                auto g = geo.corner(k);
                for(unsigned u = 0; u < dim; u++)
                    _v->_global(u) = g[u];

//                 _v._id  = v.id();
//                 _v._idx = v.idx();
                _v->_entity_seed.push_back( &(_entities.back()) );

                _bounding_box.append(_v->_global);                
                _l_vertex.push_back( _v );                     // TODO: use mapping to avoid duplication!!!!
            }
        }

        std::cout << "Bounding box\n" ; _bounding_box.operator<<(std::cout) << std::endl;
        std::cout << "Number of vertices " << _l_vertex.size() << std::endl;
        
        // generate list of vertices
//         this->put( _l_vertex.begin(), _l_vertex.end() );
    }

     // iterate over all leafs of the node
    LeafView<GridView> leafView() const
    {
        return LeafView<GridView>( *this );
    }

     // iterate over all leafs of the node
    LevelView<GridView> levelView(unsigned level) const
    {
        return LevelView<GridView>( *this, level );
    }

    
    virtual void fillTreeStats( typename Node<GridView>::TreeStats& ts ) const {
        Node<GV>::fillTreeStats(ts);
        
        ts.numVertices         = _vertex.size();
        ts.aveLevel           /= static_cast<Real>( ts.numNodes );
        ts.aveVertices        /= static_cast<Real>( ts.numNodes );
        ts.aveEntitiesPerLeaf /= static_cast<Real>( ts.numLeafs );    
    }
    
    void printTreeStats( std::ostream& out ) const {
        typename Node<GridView>::TreeStats ts;
        fillTreeStats(ts);
        ts.operator<<(out) << std::endl;
    }
};


}
