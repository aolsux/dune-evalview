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

#include <math/boundingbox.hpp>
#include <assert.h>

namespace evalview {



template< class GV >
class Node {
protected:
    struct Traits {
        typedef GV                               GridView;
        typedef typename GridView::GridType      GridType;
        typedef typename GridView::ctype         Real;
        enum { dim = GridType::dimension, dimw = GridType::dimensionworld};
    };

    typedef math::BoundingBox<Traits::Real,Traits::dim> BoundingBox;

    enum {dim = Traits::dim};

    struct Vertex
    {
        std::vector<EntitySeed&> _element_ids;
        unsigned _id;
        unsigned _idx;
        shortvector<real,dim> _global;

        template<class V>
        Vertex(const V& v)
        {
            auto g = v.position(0);
            for(unsigned u = 0; u < dim; u++)
                _global(u) = g[u];
            _id = v.id();
            _idx = v.idx();
        }

    };

protected:
    Node<GV>*                    _parent;
    Node* _child[2];

    std::vector<Vertex& >        _vertex;

    typename Traits::GridView&   _gridView;

    BoundingBox                  _bounding_box;

    // the normal of the plane that splits this node
    shortvector<real,dim>        _normal;

    // the dimension that is split by this node
    unsigned                     _orientation;

    // the depth of the node in the tree
    unsigned                     _level;

    bool                         _isLeaf;


protected:
    Node() = delete;

    //TODO: we probably dont need it.
    Node(Node* parent, const typename Traits::GridView& gv) :
        _parent(parent), _gridView(gv) {}

    bool left(const shortvector<real,dim>& p)  const { return dot((p-_bounding_box.center),_normal)<0;     }
    bool right(const shortvector<real,dim>& p) const { return !left(p);}

    // vid contain indices of all vertices that reside in the space defined by boundingbox
    // p array of global space coordinates corresponding to the indices stored int vid
    // s size of p and vid TODO: remove s we dont need it
    void put(beginIt, endIt)
    {
        vertex.clear();
        std::copy( beginIt, endIt, _vertex.begin());

        // abort the recursion if there is only one vertex left within this node
        if(_vertex_ids.size()<=1)
            return;

        split();

        std::vector<unsigned> l,r;
        for(unsigned k : _vertex_ids)
        {
            assert(k < s,"index out of bounds");
            if(left(p[k]))
                l.push_back(k);
            else
                r.push_back(k);
        }

        _child[0]->put(l,p,s);
        _child[1]->put(r,p,s);

    }

    void split()
    {
        assert(_child[0] == NULL, "Child 0 exists!" );
        assert(_child[1] == NULL, "Child 1 exists!" );
        // construct the two childs
        _child[0] = new Node(this, _bounding_box.split(_orientation,true), _level+1);
        _child[1] = new Node(this, _bounding_box.split(_orientation,false), _level+1);
    }

public:

    Node( const Node& node ) = delete;

    Node( const Node* parent, const BoundingBox& box, unsigned level) :
        _parent(parent),
        _level(level),
        _bounding_box(box),
        _normal(0.),
        _orientation(level%dim),
        _child({NULL,NULL})
    {
        _normal(_orientation) = 1.;
    }

    virtual ~Node()
    {
        safe_delete(_child[0]);
        safe_delete(_child[1]);
    }

    const Node* child(const unsigned i) const    {        return _child[i];    }



    // iterate over all entities of the node
    std::vector<const Entity&> entities() const {}



};


}
