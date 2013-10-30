//**************************************************************************************//
//     AUTHOR: Malik Kirchner "malik.kirchner@gmx.net"                                  //
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

#include <dune/pdelab/common/function.hh>


namespace Dune {
namespace PDELab {


template<typename GF, class Locator>
class LocatingGridFunctionToFunctionAdapter : public FunctionInterface<FunctionTraits<  typename GF::Traits::GridViewType::ctype,
                                                                                        GF::Traits::GridViewType::dimensionworld,
                                                                                        Dune::FieldVector<typename GF::Traits::GridViewType::ctype,
                                                                                        GF::Traits::GridViewType::dimensionworld>,                                                                                    
                                                                                        typename GF::Traits::RangeFieldType,
                                                                                        GF::Traits::dimRange,
                                                                                        Dune::FieldVector<typename GF::Traits::RangeFieldType,
                                                                                        GF::Traits::dimRange> >,
                                                                                        GridFunctionToFunctionAdapter<GF> >
{
public:
    //! \brief Export type traits
    typedef FunctionTraits< typename GF::Traits::GridViewType::ctype,
                            GF::Traits::GridViewType::dimensionworld,
                            Dune::FieldVector<typename GF::Traits::GridViewType::ctype,
                            GF::Traits::GridViewType::dimensionworld >,
                            typename GF::Traits::RangeFieldType,
                            GF::Traits::dimRange,
                            Dune::FieldVector<typename GF::Traits::RangeFieldType,
                            GF::Traits::dimRange>
    > Traits;

    //! make a GridFunctionToFunctionAdapter
    LocatingGridFunctionToFunctionAdapter( const GF &gf_, const Locator& locator_) : gf(gf_) , locator(locator_) { }

    /** \brief Evaluate all basis function at given position

    Evaluates all shape functions at the given position and returns
    these values in a vector.
     */
    static_assert(std::is_same<typename Traits::DomainType,  Dune::FieldVector<double, 3> >::value, "domain type mismatch");
    static_assert(std::is_same<typename Traits::RangeType,  Dune::FieldVector<double, 1> >::value, "range type mismatch");
    inline void evaluate (const typename Traits::DomainType& x,    typename Traits::RangeType& y) const
    {
        typename GF::Traits::GridViewType::Grid::Traits::template Codim<0>::EntityPointer
        ep = locator.findEntity(x);
        gf.evaluate(*ep, ep->geometry().local(x), y);
    }

private:
    const GF&         gf;
    const Locator&    locator;
};


}
}
