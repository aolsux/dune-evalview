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

#include <dune/config.h>
#include <dune/common/fvector.hh>

#include <dune/grid/alugrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <tree/pointlocator.hpp>

#include <random>

#include <utils/utils.hpp>


// #include <gperftools/profiler.h>

template<unsigned D, class T>
void randomize(Dune::FieldVector<T,D>& v)
{
   static std::mt19937 rng;
   static std::uniform_real_distribution<double> distribution(0.0,1.0);
    for (unsigned d = 0 ; d < D ; d++ )
         v[d] = distribution(rng);
}


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
template< class GridView >
void locate_kdtree( const GridView& gridview, const std::vector<Dune::FieldVector<typename GridView::ctype, GridView::dimension> >&  coordinates, unsigned loops )
{
   std::cout << CE_STATUS <<  "building k-d-Tree ..."<< CE_RESET <<  std::endl;
//    ProfilerStart("treebuild.prof");
   tree::PointLocator< GridView > locator(gridview,false);
//    ProfilerStop();
   std::cout << CE_STATUS <<  "k-d-Tree statistics"<< CE_RESET <<  std::endl;
   locator.printTreeStats( std::cout );
  
   for (unsigned u = 0; u < loops; u++)
      for(const auto& x : coordinates)
//          try
//          {
            auto info = locator.findEntity( x );
//          }
//          catch (GridError& e)
//          {
            // the coordinate is not within the grid
//          }
   
}


// the old search algorithm doing lots of coordinate transformations
template < class GridView >
void locate_hierarchic(const GridView& gridview, const std::vector<Dune::FieldVector<typename GridView::ctype, GridView::dimension> >&  coordinates, unsigned loops )
{
   Dune::HierarchicSearch< typename GridView::GridType, typename GridView::LeafIndexSet > locator( gridview.grid(), gridview.grid().leafIndexSet() );
   for (unsigned u = 0; u < loops; u++)
      for(const auto& x : coordinates)
//          try
//          {
              auto info = locator.findEntity( x );
//          }
//          catch (GridError& e)
//          {
            // the coordinate is not within the grid
//          }
}


// compare dune hierarchic search with new kd-tree search
template < class Grid >
void benchmark(const Grid& grid) 
{
   typedef Dune::FieldVector<Grid::ctype, Grid::dimw> FieldVector;
   
   const unsigned nV = 100;
   const unsigned nL = 2000;

   // create a list of random coordinates in the unitcube
   std::vector<FieldVector> fv; fv.reserve(nV);
   for ( auto v : fv)
      randomize(v);
  
   // search for the entities containing the coordinates
   Timer t;
   t.tic();
   std::cout << CE_STATUS << "kd-tree " << CE_RESET;
   locate_kdtree(grid.leafView(), fv, nL);
   const Real ta = t.toc();
   std::cout << ta << std::endl;
   
   // search for the entities containing the coordinates
   std::cout << CE_STATUS << "hr-tree " << CE_RESET;
   t.tic();
   locate_hierarchic(grid.leafView(), fv, nL);
   const Real tb = t.toc();
   std::cout << ta << std::endl;
   std::cout << CE_STATUS << "SPEED-UP  " << CE_RESET << tb/ta << "x" << std::endl;
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
        elements[0] = elements[1] = elements[2] = 10;
       
       {
         // use structured cube grid 
         typedef Dune::ALUCubeGrid< 3, 3 > Grid;
         typedef Dune::StructuredGridFactory<Grid> GridFactory;
         
         Grid* grid = GridFactory::createCubeGrid(ll,ur,elements);
         
         benchmark(*grid);
       }
       
       {
         // use structured simplex grid 
         typedef Dune::ALUSimplexGrid< 3, 3 > Grid;
         typedef Dune::StructuredGridFactory<Grid> GridFactory;
         
         Grid* grid = GridFactory::createSimplexGrid(ll,ur,elements);
         
         benchmark(*grid);
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
