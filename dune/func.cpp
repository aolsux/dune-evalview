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

#define ENABLE_ALUGRID 1

#include <dune/grid-config.h>
#include <dune/common/fvector.hh>

#include <dune/grid/alugrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/utility/hierarchicsearch.hh>
#include <tree/pointlocator.hpp>

#include <random>

#include <utils/utils.hpp>
#include <fem/helper.hpp>
#include <fem/kdtreegridfunction.hpp>
#include <dune/grid/geometrygrid/entity.hh>
#include <dune/grid/io/file/dgfparser/dgfgridfactory.hh>
#include <dune/pdelab/gridfunctionspace/interpolate.hh>

#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/common/function.hh>

#include <test.h>

template< typename GV, typename F>
struct FuncTraits
{
    enum {dim   = GV::dimension,
          dimw  = GV::dimensionworld};
    typedef GV             GridView;
    typedef typename GridView::ctype                    Real;

    typedef F                                                                              FEM;

    typedef Dune::PDELab::NoConstraints                                                    Constraints;
    typedef Dune::PDELab::ISTLVectorBackend<1>                                             VectorBackend;
    typedef Dune::PDELab::GridFunctionSpace<GridView,FEM,Constraints,VectorBackend >       GridFunctionSpace;
    typedef Dune::PDELab::ISTLBlockVectorContainer<GridFunctionSpace, Real, 1>             VectorContainer;
    typedef Dune::PDELab::DiscreteGridFunction<GridFunctionSpace,VectorContainer>          GridFunction;
    
    class AnalyticFunction  
      : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,double,1>,AnalyticFunction>
    {
        public:
        typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,double,1> Traits;
        typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,AnalyticFunction> BaseType;

        AnalyticFunction (const GV& gv): BaseType(gv)
        {}

        inline void evaluateGlobal (const typename Traits::DomainType& x, typename Traits::RangeType& u) const
        {
            double abs = x.two_norm();
            u = std::pow(abs, 3)-std::pow(abs, 2);
        }

    }; // AnalyticFunction

};

template <typename Traits, typename Grid,  typename Locator, typename FieldVector>
double interpolate(Traits traits, const Locator& locator, const Grid& grid, const std::vector<FieldVector>& coordinates,  unsigned loops)
{
    
//     typedef typename Traits::GridView          GridView;
    typedef typename Traits::GridFunction      GridFunction;
    typedef typename Traits::AnalyticFunction  AnalyticFunction;
//     typedef typename Traits::VectorBackend     VectorBackend;
    typedef typename Traits::VectorContainer   VectorContainer;
    typedef typename Traits::GridFunctionSpace GridFunctionSpace;
    typedef typename Traits::FEM               FEM;
    typedef typename Traits::Constraints       Constraints;
    
    typedef Dune::PDELab::LocatingGridFunctionToFunctionAdapter<GridFunction, Locator> InterpolatedFunction;
    
    FEM               fem;
    Constraints       ce;
    GridFunctionSpace gfs(grid.leafView(), fem, ce);
    VectorContainer   vb(gfs);    
    GridFunction gridfunction(gfs, vb);
    
    AnalyticFunction     analytic(grid.leafView());
    InterpolatedFunction interpolation(gridfunction, locator);
    
    Dune::PDELab::interpolate(analytic,gfs,vb);
    
    double error = 0;
    
    for (unsigned u = 0; u < loops; u++)
        for (const auto& v :  coordinates)
        {
            Dune::FieldVector<double, 1> iresult,  aresult;
            interpolation.evaluate(v, iresult);
            analytic.evaluateGlobal(v,  aresult);
            error += std::abs((iresult-aresult)/aresult);
        }    
    
    return error;
}

// compare dune hierarchic search with new kd-tree search
template < class Grid >
void benchmark(const Grid& grid) {  
    static constexpr  unsigned              dim  = Grid::dimension;
    static constexpr  unsigned              dimw = Grid::dimensionworld;
    typedef typename  Grid::LeafGridView    GridView;
    typedef typename  GridView::ctype       Real;
    typedef Dune::FieldVector<Real, dimw>   FieldVector;
    
    const unsigned nV = 1000;
    const unsigned nL = 2000;
    const unsigned sp = 200;

    // create a list of random coordinates in the unitcube
    std::vector<FieldVector> fv(nV);
    for ( auto& v : fv) randomize<Real, dimw>(v);
    
    // search for the entities containing the coordinates
    Timer t;
    
    typedef typename Dune::PDELab::P1LocalFiniteElementMap< Real, Real, dim >                                      FEM;
    typedef FuncTraits<GridView,  FEM> Traits;
    
   std::cout << CE_STATUS <<  "building k-d-Tree ..."<< CE_RESET <<  std::endl;
   tree::PointLocator< GridView > kd_locator(grid.leafView(),false);
   std::cout << CE_STATUS <<  "k-d-Tree statistics"<< CE_RESET <<  std::endl;
   kd_locator.printTreeStats( std::cout );
    
    t.tic();
    std::cout << CE_STATUS << "kd-tree " << CE_RESET;
    const double resKD = interpolate( Traits(),kd_locator,  grid, fv, nL)/nL/nV;
    const Real ta = t.toc();
    std::cout << CE_STATUS << "time: " << ta << ",  rate: " << (double)(nV*nL)/(ta) << ",  error: " <<  resKD << CE_RESET <<  std::endl;
    
    // search for the entities containing the coordinates
    Dune::HierarchicSearch< typename GridView::Grid, typename GridView::IndexSet > hr_locator( grid,grid.leafIndexSet() );
    std::cout << CE_STATUS << "hr-tree " << CE_RESET;
    t.tic();
    const double resHR = interpolate(Traits(), hr_locator, grid, fv, nL/sp)*sp/nL/nV;
    const Real tb = sp*t.toc();
    std::cout << CE_STATUS << "time: " << tb << ",  rate: " << (double)(nV*nL)/(tb) << ",  error: " <<  resHR << CE_RESET <<  std::endl;
    std::cout << CE_STATUS << "SPEED-UP  " << CE_RESET << tb/ta << "x" << CE_RESET  << std::endl;
}



int main ( int argc, char **argv ) {
    Dune::MPIHelper::instance( argc, argv );

    srand((unsigned)std::time(NULL));

    std::cout.setf( std::ios::scientific );
    std::cout.precision( 4 );

    try {
       
        Dune::FieldVector<double, 3> ll,ur;
        ll[0] = ll[1] = ll[2] = 0.;
        ur[0] = ur[1] = ur[2] = 1.;
         
        Dune::array<unsigned,3> elements;
        elements[0] = elements[1] = elements[2] = 3;
       
//        {
//          // use structured cube grid 
//          typedef Dune::ALUCubeGrid< 3, 3 > Grid;
//          typedef Dune::StructuredGridFactory<Grid> GridFactory;
//          
//          Dune::shared_ptr<Grid> pgrid = GridFactory::createCubeGrid(ll,ur,elements);
//          
//          benchmark(*pgrid);
//        }
       
       {
         // use structured simplex grid 
         typedef Dune::ALUSimplexGrid< 3, 3 > Grid;
         typedef Dune::StructuredGridFactory<Grid> GridFactory;
         
         Dune::shared_ptr<Grid> pgrid = GridFactory::createSimplexGrid(ll,ur,elements);
         
         benchmark(*pgrid);
       }
      
      
           

    } catch ( std::exception & e) {
        std::cout << " STL ERROR : " << e.what () << std::endl;
        return 1;
    } catch ( Dune::Exception & e ) {
        std::cout << " DUNE ERROR : " << e.what () << std::endl;
        return 1;
    } catch (...) {
        std::cout << " Unknown ERROR " << std::endl;
        return 1;
    }

    return 0;
}
