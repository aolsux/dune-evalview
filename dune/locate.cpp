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

#include <test.h>

template < class Grid >
void random_refine(Grid& grid, double fraction = 0.25)
{
   std::mt19937 rng;
   // produce random numbers 0 or 1 with prbabilities (1-fraction) and fraction
   std::discrete_distribution<int> distribution({1-fraction,fraction});
   auto view = grid.leafView();
   
   // prepare the grid for refinement
   grid.preAdapt();
   
   // adapt the grid
   grid.adapt();
   for ( auto cell = view.template begin<0>(); cell != view.template end<0>(); ++cell )
      grid.mark( distribution(rng), *cell );
   
   // clean up
   grid.postAdapt();
}



// the new search algorithm
template< class GridView , class Locator>
unsigned locate( const Locator& locator, const GridView& gridview, const std::vector<Dune::FieldVector<typename GridView::ctype, GridView::dimension> >&  coordinates, unsigned loops )
{
    typedef typename GridView::Grid::template Codim<0>::EntityPointer EntityPointer;
    typedef typename GridView::Grid::template Codim<0>::Entity        Entity;
    typedef typename  GridView::Grid        GridType;
    static constexpr  unsigned              dim  = GridType::dimension;
    static constexpr  unsigned              dimw = GridType::dimensionworld;
    typedef typename  GridView::ctype       Real;

    auto& gids = gridview.grid().globalIdSet();
    unsigned res = 0;    
    for (unsigned u = 0; u < loops; u++)
        for(const auto& x : coordinates) {
            const EntityPointer ep = locator.findEntity( x );
            res += gids.id(*ep);
    }
   
   return res;
}

// compare dune hierarchic search with new kd-tree search
template < class Grid >
bool benchmark(const Grid& grid) {  
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
    
   std::cout << CE_STATUS <<  "building k-d-Tree ..."<< CE_RESET <<  std::endl;
//    ProfilerStart("treebuild.prof");
   tree::PointLocator< GridView > kd_locator(grid.leafView(),false);
//    ProfilerStop();
   std::cout << CE_STATUS <<  "k-d-Tree statistics"<< CE_RESET <<  std::endl;
   kd_locator.printTreeStats( std::cout );
    
    t.tic();
    std::cout << CE_STATUS << "kd-tree " << CE_RESET;
    const unsigned resKD = locate( kd_locator, grid.leafView(), fv, nL);
    const Real ta = t.toc();
    std::cout << "seconds " << ta << ",  rate " << (double)(nV*nL)/(ta) << std::endl;
    
    // search for the entities containing the coordinates
    Dune::HierarchicSearch< typename GridView::Grid, typename GridView::IndexSet > hr_locator( grid,grid.leafIndexSet() );
    std::cout << CE_STATUS << "hr-tree " << CE_RESET;
    t.tic();
    const unsigned resHR = locate(hr_locator, grid.leafView(), fv, nL/sp)*sp;
    const Real tb = sp*t.toc();
    std::cout << "seconds " << tb << ",  rate " << (double)(nV*nL)/(tb) << std::endl;
    std::cout << CE_STATUS << "SPEED-UP  " << CE_RESET << tb/ta << "x" << std::endl;
    
    return resKD == resHR;
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
       
       {
         // use structured cube grid 
         typedef Dune::ALUCubeGrid< 3, 3 > Grid;
         typedef Dune::StructuredGridFactory<Grid> GridFactory;
         
         Dune::shared_ptr<Grid> pgrid = GridFactory::createCubeGrid(ll,ur,elements);
         
         bool pass = benchmark(*pgrid);
         
         if ( !pass ) std::cout <<  CE_WARNING << "kd-locator and hr-locator have different result!" << CE_RESET << std::endl;
       }
       
       {
         // use structured simplex grid 
         typedef Dune::ALUSimplexGrid< 3, 3 > Grid;
         typedef Dune::StructuredGridFactory<Grid> GridFactory;
         
         Dune::shared_ptr<Grid> pgrid = GridFactory::createSimplexGrid(ll,ur,elements);
         
         bool pass = benchmark(*pgrid);
         
         if ( !pass ) std::cout <<  CE_WARNING << "kd-locator and hr-locator have different result!" << CE_RESET << std::endl;
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
