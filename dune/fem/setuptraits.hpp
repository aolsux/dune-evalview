#pragma once

#include <fem/dune.h>
#include <fem/boundary.hpp>
#include <fem/hierarchical.hpp>
#include <fem/pointlocator.hpp>

template< typename BT, unsigned dim_, class LocalOperator_, class FunctionOperator_ >
struct ALUSimplexP1Traits {
    enum {dim   = dim_, 
          dimw  = dim_};
    typedef typename Dune::ALUSimplexGrid< dim, dimw >  GridType;
    typedef typename GridType::LeafGridView             GridView;
    typedef typename GridType::LevelGridView            LevelView;
    typedef typename Dune::DGFWriter< GridView >        GridWriter;
    typedef typename GridView::ctype                    Coord;
    typedef BT                                          Real;
    typedef LocalOperator_                              LocalOperator;
    typedef FunctionOperator_                           FunctionOperator;
    
    typedef BCTypeParam                                 BCFunc; 
    typedef BCExtension<GridView, Real>                 BCExt; 

    typedef typename Dune::PDELab::P1LocalFiniteElementMap< Coord, Real, dim >                                      FEM;    
    typedef typename Dune::PDELab::HangingNodesConstraintsAssemblers::SimplexGridP1Assembler                        HangingNodeAssembler;   
    typedef typename Dune::PDELab::HangingNodesDirichletConstraints< GridType, HangingNodeAssembler, BCTypeParam >  Constraints;
    
    
    typedef typename Dune::PDELab::ISTLVectorBackend<1>                                                             VectorBackend;
    typedef typename Dune::PDELab::GridFunctionSpace<GridView,FEM,Constraints,VectorBackend >                       GridFunctionSpace;
    typedef typename GridFunctionSpace::template ConstraintsContainer<Real>::Type                                   ConstraintsContainer;
    typedef typename Dune::PDELab::ISTLBCRSMatrixBackend<1, 1>                                                      MatrixBackend;
    
    typedef typename Dune::PDELab::GridOperator< GridFunctionSpace, GridFunctionSpace, 
                                                 LocalOperator, MatrixBackend, Real, Real, Real, 
                                                 ConstraintsContainer, ConstraintsContainer>                        GridOperator;

    typedef typename Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR                                                        Solver;

    typedef typename GridOperator::Traits::Domain                                                                   FieldU;
    
    typedef typename Dune::PDELab::StationaryLinearProblemSolver<GridOperator,Solver,FieldU>                        LinearProblemSolver;
    typedef typename Dune::PDELab::DiscreteGridFunction<GridFunctionSpace,FieldU>                                   DiscreteGridFunction;
  
    typedef typename Dune::PDELab::L2Projection<GridFunctionSpace,FieldU>                                           Projector;
    typedef HierarchicalErrorEstimation<GridFunctionSpace,FieldU,FunctionOperator>                                  ErrorEstimation;
    typedef HierarchicalEstimationAdaptation<GridType,GridFunctionSpace,FieldU,ErrorEstimation>                     EstimationAdaptation;
    typedef typename Dune::PDELab::GridAdaptor<GridType,GridFunctionSpace,FieldU,EstimationAdaptation,Projector>    GridAdaptor;
    
    static typename Dune::shared_ptr<GridType> createGrid( const typename Dune::FieldVector<Coord,dimw>& lowerLeft,
                                                           const typename Dune::FieldVector<Coord,dimw>& upperRight,
                                                           const typename Dune::array<unsigned int,dim>& elements) 
    {
        return Dune::StructuredGridFactory<GridType>::createSimplexGrid(lowerLeft, upperRight, elements);
    }
};

template< typename BT, unsigned dim_, class LocalOperator_, class FunctionOperator_ >
struct ALUCubeQ1Traits {
    enum {dim   = dim_, 
          dimw  = dim_};
    typedef typename Dune::ALUCubeGrid< dim, dimw > GridType;
    typedef typename GridType::LeafGridView         GridView;
    typedef typename GridType::LevelGridView        LevelView;
    typedef typename Dune::DGFWriter< GridView >    GridWriter;
    typedef typename GridView::ctype                Coord;
    typedef BT                                      Real;
    typedef LocalOperator_                          LocalOperator;
    typedef FunctionOperator_                       FunctionOperator;
    
    typedef BCTypeParam                             BCFunc; 
    typedef BCExtension<GridView, Real>             BCExt; 

    typedef typename Dune::PDELab::Q1LocalFiniteElementMap< Coord, Real, dim >                                      FEM;    
    typedef typename Dune::PDELab::HangingNodesConstraintsAssemblers::CubeGridQ1Assembler                           HangingNodeAssembler;   
    typedef typename Dune::PDELab::HangingNodesDirichletConstraints< GridType, HangingNodeAssembler, BCTypeParam >  Constraints;
    
    
    typedef typename Dune::PDELab::ISTLVectorBackend<1>                                                             VectorBackend;
    typedef typename Dune::PDELab::GridFunctionSpace<GridView,FEM,Constraints,VectorBackend >                       GridFunctionSpace;
    typedef typename GridFunctionSpace::template ConstraintsContainer<Real>::Type                                   ConstraintsContainer;
    typedef typename Dune::PDELab::ISTLBCRSMatrixBackend<1, 1>                                                      MatrixBackend;
    
    typedef typename Dune::PDELab::GridOperator< GridFunctionSpace, GridFunctionSpace, 
                                                 LocalOperator, MatrixBackend, Real, Real, Real, 
                                                 ConstraintsContainer, ConstraintsContainer>                        GridOperator;

    typedef typename Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR                                                        Solver;

    typedef typename GridOperator::Traits::Domain                                                                   FieldU;
    
    typedef typename Dune::PDELab::StationaryLinearProblemSolver<GridOperator,Solver,FieldU>                        LinearProblemSolver;
    typedef typename Dune::PDELab::DiscreteGridFunction<GridFunctionSpace,FieldU>                                   DiscreteGridFunction;
  
    typedef typename Dune::PDELab::L2Projection<GridFunctionSpace, FieldU>                                          Projector;
    typedef HierarchicalErrorEstimation<GridFunctionSpace,FieldU, FunctionOperator>                                 ErrorEstimation;
    typedef HierarchicalEstimationAdaptation<GridType,GridFunctionSpace,FieldU,ErrorEstimation>                     EstimationAdaptation;
    typedef typename Dune::PDELab::GridAdaptor<GridType,GridFunctionSpace,FieldU,EstimationAdaptation,Projector>    GridAdaptor;
    
    
    
    static typename Dune::shared_ptr<GridType> createGrid( const typename Dune::FieldVector<Coord,dimw>& lowerLeft,
                                                           const typename Dune::FieldVector<Coord,dimw>& upperRight,
                                                           const typename Dune::array<unsigned int,dim>& elements) 
    {
        return Dune::StructuredGridFactory<GridType>::createCubeGrid(lowerLeft, upperRight, elements);
    }
};
