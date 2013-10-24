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

public:
    enum Orientation {
        co_left     = -1,
        co_right    = +1
    };



protected:
    Node<GV>*                    _father;
    typename Traits::GridView&   _gridView;
    Orientation                  _orientation;
    BoundingBox                  _bounding_box;

protected:
    Node() = delete;

    Node(Node* father, const typename Traits::GridView& gv) :
        _father(father), _gridView(gv) {}

public:

    Node( const Node& node ) = delete;

    Node( const Node* father, const Orientation orientation ) :
        _father(father), _orientation(orientation)
    {

    }

   

    // iterate over all entities of the node
    std::vector<const Entity&> entities() const {}



};


}
