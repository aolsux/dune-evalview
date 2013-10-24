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
