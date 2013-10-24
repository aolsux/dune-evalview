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

// constraints parameter class for selecting boundary condition type
class BCTypeParam
: public Dune::PDELab::DirichletConstraintsParameters /*@\label{bcp:base}@*/
{
public:

    template<typename I>
    bool isDirichlet(
        const I & intersection,   /*@\label{bcp:name}@*/
        const Dune::FieldVector<typename I::ctype, I::dimension-1> & coord
      ) const
    {
//         Dune::FieldVector<typename I::ctype, I::dimension>
//         x = intersection.geometry().global( coord );
/*
        if ( (x[0]<-1+1E-6) || (x[0]>1-1E-6) || (x[1]<-1+1E-6) || (x[1]>1-1E-6) )
            return true;
        */
        return true;  // Dirichlet b.c. on all other boundaries
                                            // Dirichlet b.c. on all other boundaries
    }
};

template<typename GV, typename RF>
class BCExtension : public Dune::PDELab::GridFunctionBase<Dune::PDELab::
                                 GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> >, BCExtension<GV,RF> >
{
    const GV& gv;
public:
    typedef Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> > Traits;

    //! construct from grid view
    BCExtension (const GV& gv_) : gv(gv_) {}

    //! evaluate extended function on element
    inline void evaluate (const typename Traits::ElementType& e,
                          const typename Traits::DomainType&  xlocal,
                                typename Traits::RangeType&   y) const
    {
//         const int dim = Traits::GridViewType::Grid::dimension;
//         typedef typename Traits::GridViewType::Grid::ctype ctype;
//         Dune::FieldVector<ctype,dim> x = e.geometry().global(xlocal);

//         if ( (x[0]<-1+1E-6) || (x[0]>1-1E-6) || (x[1]<-1+1E-6) || (x[1]>1-1E-6) )
//             y = 0.0;
//         else
            y = 0.0;

        return;
    }

    //! get a reference to the grid view
    inline const GV& getGridView () {return gv;}
};
