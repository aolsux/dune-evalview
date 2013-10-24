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

#include <fem/dune.h>

template< class GridType_, class LOP_ >
struct PointLocatorTraits {
    enum {dim   = GridType_::dimension,
          dimw  = GridType_::dimensionworld};

    typedef GridType_                               GridType;
    typedef LOP_                                    LocalOperator;
    typedef typename GridType::LeafGridView::ctype  Coord;
    typedef typename LOP_::Result                   Result;

    struct ReturnType {
        Result      res;
        bool        found;

        ReturnType() : res(),  found(false) {}
        ReturnType( const ReturnType& rt ) : res(rt.res),  found(rt.found) {}
    };
};

template< class GridType_, class LOP_, class Coord_ >
class PointLocator {
public:
    typedef PointLocatorTraits< GridType_, LOP_ > Traits;

private:

protected:
    typename Traits::GridType&        grid;
    typename Traits::LocalOperator&   lop;

public:
    PointLocator( typename Traits::GridType& grid_, typename Traits::LocalOperator& lop_ ) : grid(grid_), lop (lop_) {}

    template<class Iterator, class FieldU >
    const typename Traits::ReturnType eval( Iterator& it, const typename Dune::FieldVector<typename Traits::Coord, Traits::dimw>& x, const FieldU& field  ) {
        const auto         xl  = it->geometry().local(x);
        Dune::GeometryType gt  = it->geometry().type();
        const auto&        gre = Dune::GenericReferenceElements<typename Traits::Coord, Traits::dim>::general(gt);
        typename Traits::ReturnType dx;

        if ( gre.checkInside(xl) ) {
            if ( it->isLeaf() ) {
                dx.res      = lop.eval( it, xl, field );
                dx.found    = true;
            } else {
                for ( auto son = it->hbegin(it->level()+1); son != it->hend(it->level()+1); ++son ) {
                    dx = eval( son, x, field );
                    if ( dx.found ) break;
                }
            }
        }

        return dx;
    }

    template< class FieldU >
    const typename Traits::ReturnType eval( const typename Dune::FieldVector<typename Traits::Coord, Traits::dimw>& x, const FieldU& field  ) {
        // start tree search at grid level 0
        typename Traits::GridType::LevelGridView lv = grid.levelView(0);

        // iterate grid level 0, depth-first search, return if found
        typename Traits::ReturnType dx;
        for ( auto it = lv.template begin<0>(); it != lv.template end<0>();  ++it ) {
            dx = eval( it, x, field );
            if ( dx.found ) break;
        }

        return dx;
    }
};

//
// template< class GridType_, class LOP_, class Coord_ >
// class ODEWalker {
// public:
//     typedef PointLocatorTraits< GridType_, LOP_ > Traits;
//
// private:
//     Traits::Real dt;                                          // time step
//     Traits::Real fr;                                           // friction
//     ShortVector<Real, Traits::dimw> xo;
//     ShortVector<Real, Traits::dimw> xn;
//     ShortVector<Real, Traits::dimw> vo;
//     ShortVector<Real, Traits::dimw> vn;
//
// protected:
//     typename Traits::GridType&        grid;
//     typename Traits::LocalOperator&   lop;
//
// public:
//     ODEWalker( typename Traits::GridType& grid_, typename Traits::LocalOperator& lop_ ) : grid(grid_), lop (lop_) {
//         dt = .004;                                          // time step
//         fr = .02;                                           // friction
//         xo = .80;
//         xn = .80;
//         vo = .04;
//         vn = .04;
//         vo(0) = .0;
//         vn(0) = .0;
//     }
//
//     /*
//     *  (i)     ODE lösen,  solange in gleicher Zelle -> (i)/(ii)
//     *  (ii)    wird Zelle verlassen, Iterator ändern, zurück zu PointLocator -> (i)/(iii)/(iv)
//     *  (iii)   verlässt Trajektorie das Gebiet? -> throw
//     *  (iv)    ist nächster Aufpunkt nicht in nächsten Nachbarn, Schrittweite verringern -> (i)
//     *
//     */
//
//     template<class Iterator, class FieldU >
//     const typename Traits::ReturnType eval( Iterator& it, const typename Dune::FieldVector<typename Traits::Coord, Traits::dimw>& x, const FieldU& field  ) {
//         bool next_cell = true;
//
//         for ( Real t = 0.; t < 200. + .1*dt; t+=dt ) {
//             const auto         xl  = it->geometry().local(x);
//             Dune::GeometryType gt  = it->geometry().type();
//             const auto&        gre = Dune::GenericReferenceElements<typename Traits::Coord, Traits::dim>::general(gt);
//             typename Traits::ReturnType dx;
//
//             if ( gre.checkInside(xl) ) {
//                 dx.res      = lop.eval( it, xl, field );
//             } else {
//
//             }
//
//                 du0 = rhs( xo, fieldH );
//
//                 vo = (1.-fr*dt)*vn - .1*dt*du0.du;
//                 xo = xn +    dt*vn;
//
//                 du1 = rhs( xo, fieldH );
//
//                 vn = (1.-.5*fr*dt)*vn - .5*.1*dt*(du0.du+du1.du);
//                 xn = xn + .5*   dt*(vo+vn);
//
// //                 traj.push_back( XT<Real, Traits::dimw>(xn, t) );
//         }
//
//
//
//         return dx;
//     }
// };
