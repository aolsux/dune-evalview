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

#include "test.h"
#include <omp.h>
#include <gperftools/profiler.h>



class FemLocalOperator :
    public Dune::PDELab::NumericalJacobianApplyVolume< FemLocalOperator >,
    public Dune::PDELab::NumericalJacobianVolume< FemLocalOperator >,
    public Dune::PDELab::FullVolumePattern,
    public Dune::PDELab::LocalOperatorDefaultFlags
{
public:
    // pattern assembly flags
    enum { doPatternVolume = true };

    // residual assembly flags
    enum { doAlphaVolume = true };

    FemLocalOperator (unsigned int intorder_=2) : intorder(intorder_) {}

    // volume integral depending on test and ansatz functions
    template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
    void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
    {
        // extract some types
        typedef typename LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::DomainFieldType   DF;
        typedef typename LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType    RF;
        typedef typename LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::JacobianType      JacobianType;
        typedef typename LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType         RangeType;
        typedef typename LFSU::Traits::SizeType                                                             size_type;


        // dimensions
        const int dim   = EG::Geometry::dimension;
        const int dimw  = EG::Geometry::dimensionworld;

        // select quadrature rule
        Dune::GeometryType gt = eg.geometry().type();
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);

        // loop over quadrature points
        for (typename Dune::QuadratureRule<DF,dim>::const_iterator
            it=rule.begin(); it!=rule.end(); ++it)
        {
            // evaluate basis functions on reference element
            std::vector<RangeType> phi(lfsu.size());
            lfsu.finiteElement().localBasis().evaluateFunction(it->position(),phi);

            // compute u at integration point
            RF u = 0.0;
            for ( size_type i = 0; i < lfsu.size(); i++ )
                u += x[i]*phi[i];

            // evaluate gradient of basis functions on reference element
            std::vector<JacobianType> js(lfsu.size());
            lfsu.finiteElement().localBasis().evaluateJacobian(it->position(),js);

            // transform gradients from reference element to real element
            const Dune::FieldMatrix<DF,dimw,dim>    jac = eg.geometry().jacobianInverseTransposed(it->position());
            std::vector<Dune::FieldVector<RF,dim> > gradphi(lfsu.size());
            for (size_type i=0; i<lfsu.size(); i++)
                jac.mv(js[i][0],gradphi[i]);

            // compute gradient of u
            Dune::FieldVector<RF,dim> gradu(0.0);
            for (size_type i=0; i<lfsu.size(); i++)
                gradu.axpy(x[i],gradphi[i]);

            // evaluate parameters
            Dune::FieldVector<RF,dimw> globalpos = eg.geometry().global(it->position());
            RF f  = globalpos.two_norm2()<0.5 ? -10.0/*+8.*(1+cos(10.*M_PI*globalpos.two_norm2()))*/ : 10.0;
            RF a  = 10.0;

            // integrate grad u * grad phi_i + a*u*phi_i - f phi_i
            RF factor = it->weight()*eg.geometry().integrationElement(it->position());
            for (size_type i=0; i<lfsu.size(); i++) {
                r[i] += ( gradu*gradphi[i] + a*u*phi[i] - f*phi[i] )*factor;
            }
        }
    }

private:
    unsigned int intorder;
};



template< typename SetupTraits >
class FemLocalEvalOperator /*:
    public Dune::PDELab::LocalOperatorDefaultFlags*/
{
public:
    struct Result {
        static constexpr unsigned dim     = SetupTraits::dim;
        static constexpr unsigned dimw    = SetupTraits::dimw;

        typedef typename SetupTraits::Coord BT;

        math::ShortVector<BT, dimw>  x;
        BT                           u;
        math::ShortVector<BT, dimw> du;
    };

protected:
    typedef typename Dune::PDELab::LocalFunctionSpace<typename SetupTraits::GridFunctionSpace > LFSU;
    LFSU lfsu;

public:
    FemLocalEvalOperator ( const typename SetupTraits::GridFunctionSpace& gfs_ ) : lfsu(gfs_) {}

    // pointwise eval
    template< typename IT, typename X, class FieldU >
    const Result eval ( IT& it, const X& x, const FieldU& field )
    {
        // extract some types
        typedef typename SetupTraits::Real Real;
        typedef typename LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::DomainFieldType   DF;
        typedef typename LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType    RF;
        typedef typename LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::JacobianType      JacobianType;
        typedef typename LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType         RangeType;
        typedef typename LFSU::Traits::SizeType                                                             size_type;


        // dimensions
        const int dim   = SetupTraits::dim;
        const int dimw  = SetupTraits::dimw;

        Dune::PDELab::LocalVector<typename SetupTraits::FieldU::ElementType, Dune::PDELab::TrialSpaceTag> ul;
        ul.resize(lfsu.size());

        const typename SetupTraits::GridType::template Codim<0>::Entity& e = *it;
        lfsu.bind(e);
        lfsu.vread(field,ul);

        //compute u at integration point
        Real u = 0.0;
        std::vector<RangeType> phi(lfsu.size());
        for ( size_type i = 0; i < lfsu.size(); i++ )
            u += x[i]*phi[i];

        //evaluate gradient of basis functions on reference element
        std::vector<JacobianType> js(lfsu.size());
        lfsu.finiteElement().localBasis().evaluateJacobian(x,js);

        //transform gradients from reference element to real element
        const Dune::FieldMatrix<DF,dimw,dim>    jac = it->geometry().jacobianInverseTransposed(x);
        std::vector<Dune::FieldVector<RF,dim> > gradphi(lfsu.size());
        for (size_type i=0; i<lfsu.size(); i++)
            jac.mv(js[i][0],gradphi[i]);

        //compute gradient of u
        Dune::FieldVector<RF,dim> gradu(0.0);
        for (size_type i=0; i<lfsu.size(); i++)
            gradu.axpy(ul[i],gradphi[i]);

        Result res;
        res.u = u;
        for ( int k = 0; k < dimw; k++ ) {
            res.x(k)     = x[k];
            res.du(k)    = gradu[k];
        }

        return res;
    }

};



class FemFunctionOperator :
    public Dune::PDELab::NumericalJacobianApplyVolume< FemLocalOperator >,
    public Dune::PDELab::NumericalJacobianVolume< FemLocalOperator >,
    public Dune::PDELab::FullVolumePattern,
    public Dune::PDELab::LocalOperatorDefaultFlags
{
public:
    // pattern assembly flags
    enum { doPatternVolume = true };

    // residual assembly flags
    enum { doAlphaVolume = true };

    FemFunctionOperator(unsigned int intorder_=2) : intorder(intorder_) {}

    // volume integral depending on test and ansatz functions
    template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
    void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
    {
        // extract some types
        typedef typename LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::DomainFieldType   DF;
        typedef typename LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType    RF;
        typedef typename LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::JacobianType      JacobianType;
        typedef typename LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType         RangeType;
        typedef typename LFSU::Traits::SizeType                                                             size_type;


        // dimensions
        const int dim   = EG::Geometry::dimension;
//         const int dimw  = EG::Geometry::dimensionworld;

        // select quadrature rule
        Dune::GeometryType gt = eg.geometry().type();
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);

        // loop over quadrature points
        for (typename Dune::QuadratureRule<DF,dim>::const_iterator
            it=rule.begin(); it!=rule.end(); ++it)
        {
            // evaluate basis functions on reference element
            std::vector<RangeType> phi(lfsu.size());
            lfsu.finiteElement().localBasis().evaluateFunction(it->position(),phi);

            // compute u at integration point
            RF u = 0.0;
            for ( size_type i = 0; i < lfsu.size(); i++ )
                u += x[i]*phi[i];

            // integrate grad u * grad phi_i + a*u*phi_i - f phi_i
            RF factor = it->weight()*eg.geometry().integrationElement(it->position());
            for (size_type i=0; i<lfsu.size(); i++) {
                r[i] += u*phi[i]*factor;
            }
        }
    }

private:
    unsigned int intorder;
};






template< typename T, unsigned dim >
inline void asShortVector( const Dune::FieldVector< T, dim >& f, math::ShortVector< T, dim >& v ) {
    for ( int k = 0; k < dim; k++ )
        v.data[k] = f[k];
}

template< typename T, unsigned dim >
inline void asShortVector( Dune::FieldVector< T, dim >& f, math::ShortVector< T, dim >& v ) {
    for ( int k = 0; k < dim; k++ )
        v.data[k] = f[k];
}

// template< typename BT, unsigned dim, class GridView >
// void saveMesh( GridView& gridView, const std::string& path ) {
//     CubeMesh< BT, dim > cmesh;
//
//     BoundingBox< BT, dim > boundingbox;
//
//     cmesh.choleskyFactor = "";
//
//     for ( auto elm = gridView.template begin<0>(); elm!= gridView.template end<0>(); ++elm ) {
//         ShortVector< BT, dim > corner[NUM_CORNERS(dim)];
//
//         for ( int c = 0; c < NUM_CORNERS(dim); c++ ) {
//             asShortVector<BT, dim>( elm->geometry().corner(c), corner[NUM_CORNERS(dim)-c-1] );
//             boundingbox.include(corner[c]);
//         }
//
//         Cube< BT, dim > cube( corner );
//         cmesh.push_back( cube );
//     }
//
//     cmesh.boundingbox    = boundingbox;
//
//     try {
//         hxgeomatch::save( path, cmesh );
//     } catch ( boost::archive::archive_exception& err ) {
//         std::cout << "Boost Archive Error " << err.what() << "[" << boost::archive::Code2String(err.code) << "]" << " --> " << __FILE__ << ":" << __LINE__ << std::endl;
//     }
// }

template< typename BT, unsigned dim >
struct XT {
    math::ShortVector< BT, dim > x;
    BT                           t;

    XT() : x(0.), t(0.) {}
    XT( const math::ShortVector< BT, dim >& x_, const BT t_ ) : x(x_), t(t_) {}
    XT( const XT& xt ) : x(xt.x), t(xt.t) {}
};

template< typename BT, unsigned dim >
class Trajectory {
private:

protected:
    bool                         adepting;
    std::vector< XT<BT, dim> >   data;
    std::deque< XT<BT, dim> >    dq;

public:
    Trajectory() : adepting(false) {}

    typename std::vector< XT<BT, dim> >::iterator vector_begin() { if (adepting) throw; return data.begin(); }
    typename std::vector< XT<BT, dim> >::iterator vector_end()   { if (adepting) throw; return data.end(); }

    typename std::deque< XT<BT, dim> >::iterator deque_begin() { if (!adepting) throw; return dq.begin(); }
    typename std::deque< XT<BT, dim> >::iterator deque_end()   { if (!adepting) throw; return dq.end(); }

    void preAdept() {
        dq.clear();
        dq.reserve( data.size() );
        std::copy( data.begin(), data.end(), dq );
        adepting = true;
    }

    void insert( const XT<BT, dim>& xt ) {
        if (!adepting) throw;
    }

    void insert( const XT<BT, dim>& xt, const typename std::deque< XT<BT, dim> >::iterator& it ) {
        if (!adepting) throw;
        dq.insert(it, xt);
    }

    void postAdept() {
        data.clear();
        data.reserve( dq.size() );
        std::copy( dq.begin(), dq.end(), data );
        dq.clear();
        adepting = false;
    }

    void push_back( const XT<BT, dim>& xt ) {
        if ( adepting ) {
            dq.push_back( xt );
        } else {
            data.push_back( xt );
        }

    }

    void writeVTK ( const std::string path ) const {
        if (adepting) throw;
        //Create points and add a vertex at each point. Really what you are doing is adding
        //cells to the polydata, and the cells only contain 1 element, so they are, by definition,
        //0-D topology (vertices).
        vtkSmartPointer<vtkPoints>      points   = vtkSmartPointer<vtkPoints>::New();
        vtkSmartPointer<vtkCellArray>   vertices = vtkSmartPointer<vtkCellArray>::New();
        vtkSmartPointer<vtkLine>        line     = vtkSmartPointer<vtkLine>::New();
        vtkSmartPointer<vtkCellArray>   lines    = vtkSmartPointer<vtkCellArray>::New();

        const double T = 1./std::abs( data.front().t - data.back().t );
        for ( auto xt = data.begin(); xt != data.end()-1; ++xt ) {
            //Declare a variable to store the index of the point that gets added. This behaves just like an unsigned int.
            vtkIdType pid[2];

            //Add a point to the polydata and save its index, which we will use to create the vertex on that point.
            if ( dim == 2 ) {
                auto a = *xt;
                auto b = *(xt+1);
                pid[0] = points->InsertNextPoint ( a.x(0), a.x(1), T*a.t );
                pid[1] = points->InsertNextPoint ( b.x(0), b.x(1), T*b.t );
            } else if ( dim == 3 ) {
                auto a = *xt;
                auto b = *(xt+1);
                pid[0] = points->InsertNextPoint ( a.x(0), a.x(1), a.x(2) );
                pid[1] = points->InsertNextPoint ( b.x(0), b.x(1), b.x(2) );
            } else
                if (adepting) throw;

            //create a vertex cell on the point that was just added.
            vertices->InsertNextCell(1,pid);

            line->GetPointIds()->SetId(0,pid[0]);
            line->GetPointIds()->SetId(1,pid[1]);
            lines->InsertNextCell ( line );
        }

        //create a polydata object
        vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();

        //set the points and vertices we created as the geometry and topology of the polydata
        polydata->SetPoints( points   );
        polydata->SetVerts ( vertices );
        polydata->SetPolys ( lines    );

        //write the polydata to a file
        vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
        writer->SetFileName ( path.c_str() );
        writer->SetInput ( polydata );
        writer->Write();
    }
};

template< typename SetupTraits >
class FemTest {
public:
    typedef SetupTraits                                  Traits;

    typedef typename SetupTraits::Real                   Real;
    typedef typename SetupTraits::Coord                  Coord;
    typedef typename SetupTraits::GridType               GridType;
    typedef typename SetupTraits::GridView               GridView;
    typedef typename SetupTraits::FEM                    FEM;
    typedef typename SetupTraits::BCFunc                 BCFunc;
    typedef typename SetupTraits::BCExt                  BCExt;
    typedef typename SetupTraits::Constraints            Constraints;
    typedef typename SetupTraits::ConstraintsContainer   ConstraintsContainer;
    typedef typename SetupTraits::GridFunctionSpace      GridFunctionSpace;
    typedef typename SetupTraits::LocalOperator          LocalOperator;
    typedef typename SetupTraits::FunctionOperator       FunctionOperator;
    typedef typename SetupTraits::GridOperator           GridOperator;
    typedef typename SetupTraits::Solver                 Solver;
    typedef typename SetupTraits::FieldU                 FieldU;
    typedef typename SetupTraits::LinearProblemSolver    LinearProblemSolver;
    typedef typename SetupTraits::DiscreteGridFunction   DiscreteGridFunction;

    typedef typename SetupTraits::Projector              Projector;
    typedef typename SetupTraits::ErrorEstimation        ErrorEstimation;
    typedef typename SetupTraits::EstimationAdaptation   EstimationAdaptation;
    typedef typename SetupTraits::GridAdaptor            GridAdaptor;
    typedef FemLocalEvalOperator< SetupTraits >          FemEvalLOP;

protected:
    GridType&   grid;
    GridView    view;



    FEM                             fem;
    BCFunc                          bf;
    Constraints                     ce;
    ConstraintsContainer            cc;
    GridFunctionSpace               gfs;
    LocalOperator                   lop;
    FunctionOperator                funcLop;
    GridOperator                    gos;
    Solver                          solver;
    FieldU                          fieldL;
    FieldU                          fieldH;
    LinearProblemSolver             lpSolverL;
    LinearProblemSolver             lpSolverH;
    FemEvalLOP                      fleo;
    tree::PointLocator< GridView >  root;


public:
    FemTest( GridType& g, int maxIter, Real tol ) :
        grid     ( g ),
        view     ( g.leafView() ),
        fem      ( ),
        bf       ( ),
        ce       ( g, true, bf ),
        cc       ( ),
        gfs      ( view, fem, ce),
        lop      ( 2 ),
        funcLop  ( 2 ),
        gos      ( gfs, gfs, lop ),
        solver   ( maxIter, true ),
        fieldL   ( gfs, .0 ),
        fieldH   ( gfs, .0 ),
        lpSolverL( gos, fieldL, solver, tol ),
        lpSolverH( gos, fieldH, solver, tol ),
        fleo     ( gfs ),
        root     ( view, false )
    {
    }

    void updateDOF( GridAdaptor& gra, const std::vector< FieldU* > field ) {
        // prepare the grid for refinement
        grid.preAdapt();

        // save u
        std::vector< typename GridAdaptor::MapType > transferMap( field.size() );
        auto tm = transferMap.begin();
        for ( auto f = field.begin(); f != field.end(); ++f, ++tm )
            gra.backupData( **f, *tm );

        // adapt the grid
        grid.adapt();

        // update the function spaces
        gfs.update();

        // reset u
        for ( auto f : field )
            *f = FieldU(gfs,0.0);

        tm = transferMap.begin();
        for ( auto f = field.begin(); f != field.end(); ++f, ++tm )
            gra.replayData( **f , *tm );

        // clean up
        grid.postAdapt();

        std::cout << CE_STATUS <<  "building k-d-Tree ..."<< CE_RESET <<  std::endl;
        ProfilerStart("rebuild.prof");
        root.rebuild();
        ProfilerStop();
        std::cout << CE_STATUS <<  "k-d-Tree statistics"<< CE_RESET <<  std::endl;
        root.printTreeStats( std::cout );
    }

    void globalRefine( GridAdaptor& gra, const std::vector< FieldU* > field ) {
        for ( auto cell = view.template begin<0>(); cell != view.template end<0>(); ++cell )
            grid.mark( 1, *cell );

        updateDOF( gra, field );
    }

    void globalCoarsen( GridAdaptor& gra, const std::vector< FieldU* > field ) {
        for ( auto cell = view.template begin<0>(); cell != view.template end<0>(); ++cell )
            grid.mark(-1, *cell );

        updateDOF( gra, field );
    }

    void localCoarsen( GridAdaptor& gra, FieldU& cf, const std::vector< FieldU* > field ) {
        // mark local
        gra.markGrid( cf );

        // coarsen non-refined
        for ( auto cell = view.template begin<0>(); cell != view.template end<0>(); ++cell ) {
            const int m = grid.getMark(*cell);
            if ( m > 0 )
                grid.mark( 0, *cell );
            else
                grid.mark(-1, *cell );
        }

        updateDOF( gra, field );
    }

    void localRefine( GridAdaptor& gra, FieldU& cf, const std::vector< FieldU* > field ) {
        // mark local
        gra.markGrid( cf );
        updateDOF( gra, field );
    }

    void interpolate( const BCExt& g, const std::vector< FieldU* > field ) {
        Dune::PDELab::constraints(bf,gfs,cc);
        for ( auto f : field )
            Dune::PDELab::interpolate( g, gfs, *f );
    }

    void compute( unsigned maxLevel = 2 ) {
        Projector               proj;
        ErrorEstimation         ree(gfs, fieldH, funcLop );
        EstimationAdaptation    ea(grid,gfs,ree,0.5,0.0,1,maxLevel);
        GridAdaptor             gra(grid,gfs,ea,proj);

        BCExt                   g( view );

        for ( unsigned k = 0; k < maxLevel; k++ ) {
            std::cout << CE_STATUS <<  "Grid information LEVEL "<< k <<  CE_RESET <<  std::endl;
            Dune::gridinfo( grid );

            interpolate( g, {&fieldL, &fieldH} );
            lpSolverL.apply();
            globalRefine( gra, {&fieldL, &fieldH} );

            interpolate( g, {&fieldL, &fieldH} );
            lpSolverH.apply();
            if ( k < maxLevel-1 )
                localCoarsen( gra, fieldL, {&fieldL, &fieldH} );
        }

        ProfilerStart("integrate.prof");

        Timer t;
        t.tic();
//         #pragma omp parallel for
//         for ( unsigned k = 0; k < 8; k++ ) {
            integrate( view, fieldH );

//         }
        t.toc();

        ProfilerStop();

        root.printTreeStats( std::cout );
//         std::cout << CE_STATUS << "time elapsed " << t.toc() << ",     " <<  t.toc()/omp_get_num_procs() <<  CE_RESET << std::endl;


        benchmark();
    }

    void writeVTK( std::string path ) {
        DiscreteGridFunction        udgfL( gfs, fieldH );
        Dune::SubsamplingVTKWriter<GridView>   vtkwriterL( view, 2 );
        vtkwriterL.addVertexData( new Dune::PDELab::VTKGridFunctionAdapter<DiscreteGridFunction>( udgfL, "solution subsampling") );
        vtkwriterL.write( "hi_"+path, Dune::VTKOptions::ascii );

        DiscreteGridFunction        udgfH( gfs, fieldH );
        Dune::VTKWriter<GridView>   vtkwriterH( view, Dune::VTKOptions::conforming );
        vtkwriterH.addVertexData( new Dune::PDELab::VTKGridFunctionAdapter<DiscreteGridFunction>( udgfH, "solution") );
        vtkwriterH.write( "lo_"+path, Dune::VTKOptions::ascii );
    }

    typename FemLocalEvalOperator< SetupTraits >::Result rhs ( math::ShortVector<typename SetupTraits::Coord, SetupTraits::dimw>& x, const FieldU field ) {
        auto e = root.findEntity( x );
        const auto res = fleo.eval( e.pointer, e.xl, field );
        return res;
    }

    typename FemLocalEvalOperator< SetupTraits >::Result rhs ( Dune::FieldVector<typename SetupTraits::Coord, SetupTraits::dimw>& x, const FieldU field ) {
        auto e = root.findEntity( asShortVector( x) );
        const auto res = fleo.eval( e.pointer, e.xl, field );
        return res;
    }

    void integrate ( typename SetupTraits::GridView& gv, const typename SetupTraits::FieldU& v ) {
        typedef Dune::FieldVector<typename SetupTraits::Coord, SetupTraits::dimw>   Coord;
        typedef Dune::FieldVector<typename SetupTraits::Coord, 1>                   Scalar;
        typedef typename SetupTraits::Coord                                         Real;

        Trajectory< Real, Traits::dimw > traj;

        Real    dt = .004;                                          // time step
        Real    fr = .02;                                           // friction
        math::ShortVector<Real, Traits::dimw> xo( .65 );
        math::ShortVector<Real, Traits::dimw> xn( .65 );
        math::ShortVector<Real, Traits::dimw> vo( .07 );
        math::ShortVector<Real, Traits::dimw> vn( .07 );
        vo(0) = .0;
        vn(0) = .0;


        Timer t;
        t.tic();

        try {

            typename FemLocalEvalOperator<SetupTraits>::Result du0;
            typename FemLocalEvalOperator<SetupTraits>::Result du1;
            for ( Real t = 0.; t < 500. + .1*dt; t+=dt ) {
                du0 = rhs( xo, fieldH );

                vo = (1.-fr*dt)*vn - .1*dt*du0.du;
                xo = xn +    dt*vn;

                du1 = rhs( xo, fieldH );

                vn = (1.-.5*fr*dt)*vn - .5*.1*dt*(du0.du+du1.du);
                xn = xn + .5*   dt*(vo+vn);

                traj.push_back( XT<Real, Traits::dimw>(xn, t) );
//                 std::cout << "\t t = " << t << "\tx = ( " << xn << " )\n";
            }

        } catch ( GridError& err ) {
            std::cout << CE_ERROR << err.what() << CE_RESET << std::endl;
        }

        std::cout << CE_STATUS << "time elapsed " << t.toc() <<  CE_RESET << std::endl;

//         root.printTreeStats( std::cout );

        std::cout << CE_STATUS << "Write Trajectory to VTK" << CE_RESET << std::endl;
        traj.writeVTK( "traj.vtp" );
    }

    
    void benchmark() {
        Dune::HierarchicSearch< GridType, typename GridType::LeafIndexSet > _hr_locator( grid, grid.leafIndexSet() );
        
        const unsigned nV = 100;
        const unsigned nL = 2000;
        std::cout << CE_STATUS << "Fill random global coord buffer" << CE_RESET << std::endl;
        std::vector< typename Traits::LinaVector  > lv; lv.reserve(nV);
        std::vector< typename Traits::FieldVector > fv; fv.reserve(nV);
        for ( unsigned k = 0; k < nV; k++ ) {
            typename Traits::LinaVector     l;
            typename Traits::FieldVector    f;
            for ( unsigned d = 0; d < Traits::dim; d++ ) {
                l(d) = 2.*drand48()-1.;
                f[d] = l(d);
            }
            lv.push_back(l);
            fv.push_back(f);
        }
        
        Timer t;
        t.tic();
        std::cout << CE_STATUS << "kd-tree " << CE_RESET;
        for ( unsigned l = 0; l < nL; l++ ) {
            for ( unsigned k = 0; k < nV; k++ ) {
                auto ed = root.findEntity( lv[k] );
            }
        }
        const Real ta = t.toc();
        std::cout << ta << std::endl;
        std::cout << CE_STATUS << "hr-tree " << CE_RESET;
        t.tic();
        for ( unsigned l = 0; l < nL/200; l++ ) {
            for ( unsigned k = 0; k < nV; k++ ) {
                typename Traits::EntityPointer ep( _hr_locator.findEntity( fv[k] ) );
                const auto&     e   = *ep;
                const auto&     geo = e.geometry();
                const auto&     gre = Dune::GenericReferenceElements< Real, Traits::dim >::general(geo.type());
                const auto      xl  = geo.local( fv[k] );
            }
        }
        const Real tb = t.toc();
        std::cout << 200.*tb << std::endl;
        std::cout << CE_STATUS << "SPEED-UP  " << CE_RESET << 200.*tb/ta << "x" << std::endl;
    }

};


template< typename SetupTraits >
inline void compute() {
    Dune::FieldVector<typename SetupTraits::Coord, SetupTraits::dimw> lowerLeft (-1. );
    Dune::FieldVector<typename SetupTraits::Coord, SetupTraits::dimw> upperRight( 1. );
    Dune::array<unsigned int, SetupTraits::dim>                        elements;
    elements.fill(9 - 2*SetupTraits::dim);

    Dune::shared_ptr< typename SetupTraits::GridType >    pgrid    = SetupTraits::createGrid(lowerLeft, upperRight, elements);

    Dune::FieldVector<typename SetupTraits::Coord, SetupTraits::dim> x(.1);

    std::cout << CE_STATUS <<  "Setup PDE\n" <<  CE_RESET;
    FemTest< SetupTraits > femTest( *pgrid, 5000, 1e-9 );

    std::cout << CE_STATUS <<  "Solve PDE\n" <<  CE_RESET;
    femTest.compute(3);
    std::cout << CE_STATUS <<  "Grid information FINAL "<< CE_RESET <<  std::endl;
    Dune::gridinfo( *pgrid );

    std::cout << CE_STATUS <<  "Write solution to VTK\n" <<  CE_RESET;
    femTest.writeVTK( "hang_test" );
}



int main ( int argc, char **argv ) {
    Dune::MPIHelper::instance( argc, argv );

    srand((unsigned)std::time(NULL));

    std::cout.setf( std::ios::scientific );
    std::cout.precision( 4 );

    try {
// #define USE_CMD_PARAM
#ifdef USE_CMD_PARAM
        if ( (argc < 1) || (std::string(argv[1]) == "-p1d2") ) {
#endif
            typedef ALUSimplexP1Traits< double, 3, FemLocalOperator, FemFunctionOperator>   SetupTraits;
            compute< SetupTraits >();
#ifdef USE_CMD_PARAM
        } else if (std::string(argv[1]) == "-p1d3") {
            typedef ALUSimplexP1Traits< double, 3, FemLocalOperator, FemFunctionOperator>   SetupTraits;
            compute< SetupTraits >();
        } else if (std::string(argv[1]) == "-q1d2") {
            typedef ALUCubeQ1Traits< double, 2, FemLocalOperator, FemFunctionOperator>      SetupTraits;
            compute< SetupTraits >();
        } else if (std::string(argv[1]) == "-q1d3") {
            typedef ALUCubeQ1Traits< double, 3, FemLocalOperator, FemFunctionOperator>      SetupTraits;
            compute< SetupTraits >();
        } else if ((std::string(argv[1]) == "-h") || (std::string(argv[1]) == "--help")) {
            std::cout << "Test program using DUNE" << std::endl;
            std::cout << std::endl;
            std::cout << "-p1           Simplex P1-fem (default)" << std::endl;
            std::cout << "-q1           Cube    Q1-fem"           << std::endl;
            std::cout << "-h, --help    This help."               << std::endl << std::endl;
        }
#endif

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
