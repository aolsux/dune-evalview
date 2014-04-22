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

#include <fstream>

#define ENABLE_ALBERTA 0
#define ALBERTA_DIM    2
#define ENABLE_ALUGRID 1

#include <dune/config.h>

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/alugrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/alugrid/2d/alugrid.hh>
#include <dune/grid/alugrid/2d/alu2dgridfactory.hh>
#include <dune/grid/alugrid/3d/alugrid.hh>
#include <dune/grid/alugrid/3d/alu3dgridfactory.hh>
// #include <dune/grid/albertagrid.hh>
// #include <dune/grid/albertagrid/gridfactory.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/common/entity.hh>
#include <dune/grid/common/entitypointer.hh>
#include <dune/grid/io/file/dgfparser/dgfwriter.hh>

#include <dune/geometry/generalvertexorder.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dune/istl/bvector.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/io.hh>
#include <dune/istl/superlu.hh>

#include <dune/pdelab/finiteelementmap/p0fem.hh>
#include <dune/pdelab/finiteelementmap/p1fem.hh>
#include <dune/pdelab/finiteelementmap/q1fem.hh>
#include <dune/pdelab/finiteelementmap/rannacher_turek2dfem.hh>
#include <dune/pdelab/finiteelementmap/hangingnodeconstraints.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include <dune/pdelab/gridfunctionspace/genericdatahandle.hh>
#include <dune/pdelab/gridfunctionspace/interpolate.hh>
#include <dune/pdelab/constraints/constraints.hh>
#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/common/vtkexport.hh>
#include <dune/pdelab/common/functionutilities.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>
#include <dune/pdelab/backend/backendselector.hh>
#include <dune/pdelab/backend/eigenvectorbackend.hh>
#include <dune/pdelab/backend/eigenmatrixbackend.hh>
#include <dune/pdelab/backend/eigensolverbackend.hh>
#include <dune/pdelab/backend/istlvectorbackend.hh>
#include <dune/pdelab/backend/istlmatrixbackend.hh>
#include <dune/pdelab/backend/istlsolverbackend.hh>
#include <dune/pdelab/stationary/linearproblem.hh>
#include <dune/pdelab/adaptivity/adapt.hh>
#include <dune/pdelab/common/geometrywrapper.hh>
#include <dune/pdelab/localoperator/defaultimp.hh>
#include <dune/pdelab/localoperator/pattern.hh>
#include <dune/pdelab/localoperator/flags.hh>
#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/constraints/constraintsparameters.hh>

