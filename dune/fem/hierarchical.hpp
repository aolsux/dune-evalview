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






template<class GFSU, class U, class LOP = Dune::PDELab::GradientSmoothnessOperator, bool nonoverlapping_mode = false>
class HierarchicalErrorEstimation
{
    typedef typename GFSU::Traits::GridViewType             GV;
    typedef typename GV::Grid::template Codim<0>::Entity    Element;
    typedef typename GV::template Codim<0>::template Partition<Dune::Interior_Partition>::Iterator LeafIterator;
    typedef typename GV::IntersectionIterator               IntersectionIterator;
    typedef typename IntersectionIterator::Intersection     Intersection;
    typedef typename GV::ctype                              Coord;
    typedef typename GFSU::Traits::FiniteElementMapType     FEM;
    typedef typename FEM::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType RT;
    typedef Dune::PDELab::GridFunctionSpace<GV,FEM>         GFSV;
    typedef typename Dune::PDELab::BackendVectorSelector<GFSV,RT>::Type V;
    typedef typename GV::Grid::LeafIndexSet                 IndexSet;
    typedef typename IndexSet::IndexType                    IndexType;

protected:
    const GFSU& gfsu;
    U&          refU;
    const LOP&  lop;
    Dune::GeometryType gt;

public:

    /*! @brief The constructor.
     *
     * @param[in] gfsu_ The GridFunctionSpace of the problem
     * @param[in] lop_  The LocalOperator of the problem
     * @param[in] gsop_ A LocalOperator to evaluate on the Intersections, defaults to GradientSmoothnessOperator
     */
    HierarchicalErrorEstimation(const GFSU& gfsu_, U& refU_, const LOP& lop_ )
  : gfsu(gfsu_), refU(refU_), lop(lop_), gt(Dune::GeometryType::simplex,GV::dimension) {}

    HierarchicalErrorEstimation(Dune::GeometryType gt_, const GFSU& gfsu_, U& refU_, const LOP& lop_ )
  : gfsu(gfsu_), refU(refU_), lop(lop_), gt(gt_) {}

    /*! @brief Calculate an estimate of the error.
     *
     * @param[in]  u        The current solution, tested against the operators
     * @param[out] estimate The error estimate, one value per Elem (i.e. P0)
     */
    void apply(const U& u, V& estimate)
    {
        using namespace Dune::PDELab;
        FEM fem;
        GFSV gfsv(gfsu.gridView(),fem);

        // make local function spaces
        typedef LocalFunctionSpace<GFSU> LFSU;
        LFSU lfsu(gfsu);
        LFSU lfsu_n(gfsu);
        typedef LocalFunctionSpace<GFSV> LFSV;
        LFSV lfsv(gfsv);
        LFSV lfsv_ref(gfsv);

        LocalVector<typename U::ElementType, TrialSpaceTag> ul;
        LocalVector<typename U::ElementType, TrialSpaceTag> ul_ref;

        typedef LocalVector<typename V::ElementType, TestSpaceTag> LV;
        typedef typename LV::WeightedAccumulationView LVView;
        LV vl;
        LV vl_ref;

        const IndexSet& indexset = gfsv.gridView().indexSet();

        // traverse grid view
        for (LeafIterator it = gfsu.gridView().template begin<0,Dune::Interior_Partition>();
            it!=gfsu.gridView().template end<0,Dune::Interior_Partition>(); ++it)
        {
            const Element& e = *it;
            const IndexType index = indexset.index(e);

            // bind local function spaces to element
            lfsu.bind(e);
            lfsv.bind(e);
            lfsv_ref.bind(e);

            ul.resize(lfsu.size());
            ul_ref.resize(lfsu.size());

            // read coefficents
            lfsu.vread(u,ul);
            lfsu.vread(refU,ul_ref);

            vl.assign(lfsu.size(),0.0);
            LVView vlview = vl.weightedAccumulationView(1.0);
            // volume evaluation
            LocalAssemblerCallSwitch<LOP,LOP::doAlphaVolume>::
            alpha_volume(lop,ElementGeometry<Element>(e),lfsu,ul,lfsv,vlview);

            vl_ref.assign(lfsu.size(),0.0);
            LVView vlview_ref = vl_ref.weightedAccumulationView(1.0);
            LocalAssemblerCallSwitch<LOP,LOP::doAlphaVolume>::
            alpha_volume(lop,ElementGeometry<Element>(e),lfsu,ul_ref,lfsv_ref,vlview_ref);
            for ( unsigned k = 0; k < vl.size(); k++ ) {
                const auto aux  = (vl[k] - vl_ref[k]);
                vl[k] = 1.-std::abs(aux);
            }

            lfsv.vadd(vl,estimate);
        }                                                 // end of element
    }

};







template<class Grid, class GFSU, class U, class Estimation>
class HierarchicalEstimationAdaptation : public Dune::PDELab::AdaptationInterface<Grid,U>
{
    typedef typename Grid::template Codim<0>::Entity                    Element;
    typedef typename Grid::LeafGridView                                 LeafGridView;
    typedef typename LeafGridView::template Codim<0>::template Partition<Dune::Interior_Partition>::Iterator LeafIterator;
    typedef typename Grid::LeafIndexSet                                 IndexSet;
    typedef typename IndexSet::IndexType                                IndexType;
    typedef typename GFSU::Traits::GridViewType                         GV;
    typedef typename GV::ctype                                          Coord;
    typedef typename GFSU::Traits::FiniteElementMapType                 FEM;
    typedef typename FEM::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType RT;
    typedef Dune::PDELab::GridFunctionSpace<GV,FEM>                     GFSV;
    typedef Dune::PDELab::LocalFunctionSpace<GFSV>                      LFSV;
    typedef typename Dune::PDELab::BackendVectorSelector<GFSV,RT>::Type V;
    typedef std::multimap<typename V::ElementType, const IndexType>     MapType;

public:

    HierarchicalEstimationAdaptation(Grid& grid_, const GFSU& gfsu_, Estimation& estimation_,
        RT refine_, RT coarsen_ = 0., int min_ = 0, int max_ = std::numeric_limits<int>::max(), bool doCommunicate_ = true)
  : grid(grid_), gfsu(gfsu_), estimation(estimation_),refine(refine_), coarsen(coarsen_),
    min(min_), max(max_), doCommunicate(doCommunicate_), refinementMap(), localEstimate(0.),
    gt(Dune::GeometryType::simplex,GV::dimension) {}

    HierarchicalEstimationAdaptation(Dune::GeometryType gt_, Grid& grid_, const GFSU& gfsu_, Estimation& estimation_,
        RT refine_, RT coarsen_ = 0., int min_ = 0, int max_ = std::numeric_limits<int>::max(), bool doCommunicate_ = true)
  : grid(grid_), gfsu(gfsu_), estimation(estimation_),refine(refine_), coarsen(coarsen_),
    min(min_), max(max_), doCommunicate(doCommunicate_), refinementMap(), localEstimate(0.), gt(gt_) {}

    /*! @brief Prepare information before marking any of the Elements
     *
     * @param[in] u The current solution
     */
    void prepare (const U& u)
    {
        using namespace Dune::PDELab;
        FEM         fem;
        GFSV        gfsv(gfsu.gridView(),fem);
        V           estimate(gfsv,0.);
        estimation.apply(u,estimate);
        LFSV lfsv(gfsv);
        const LeafGridView leafView = grid.leafView();
        std::vector<typename V::ElementType> vl;
        const IndexSet& indexset = grid.leafIndexSet();
        refinementMap.clear();
        localEstimate = 0.;

        std::multimap<typename V::ElementType, const IndexType> tempMultiMap;
        for (LeafIterator it = leafView.template begin<0,Dune::Interior_Partition>();
            it!=leafView.template end<0,Dune::Interior_Partition>(); ++it)
        {
            const Element& e = *it;
            const IndexType index = indexset.index(e);
            lfsv.bind(e);
            lfsv.vread(estimate,vl);
            tempMultiMap.insert(std::pair<typename V::ElementType,const IndexType>(vl[0],index));
            localEstimate += vl[0];                           // * e.geometry().volume();
//             for (unsigned k = 0; k < vl.size(); k++ )
//                 std::cout <<  vl[0] <<  " ,  ";
//              std::cout << std::endl;
        }

        coarsenNumber = coarsen * leafView.size(0);
        refineNumber = (1. - refine) * leafView.size(0);

        unsigned int count = 0;
        RT accumulated_error = 0.0;
        eta.resize(leafView.size(0));
        for (typename std::multimap<typename V::ElementType, const IndexType>::const_iterator it = tempMultiMap.begin();
            it!=tempMultiMap.end(); ++it)
        {
            refinementMap.insert(std::pair<const IndexType, unsigned int>((*it).second,count++));
            accumulated_error += std::abs((*it).first);
            eta[(*it).second] = accumulated_error;
        }

        globalEstimate = localEstimate;
        if (doCommunicate && grid.comm().size() > 1) communicate();
        std::cout << "prepare: sqrt of sum of local estimators = " << sqrt(globalEstimate) << std::endl;
    }

    /*! @brief Estimate the error and mark Elems for refinement
     *
     * @param[in] e The Elem which is marked (or not)
     * @param[in] u The current solution
     */
    void mark (const Element& e, const U& u)
    {
        using namespace Dune::PDELab;
        const int level = e.level();
        const IndexSet& indexset = grid.leafIndexSet();
        const IndexType index = indexset.index(e);
//         const int i = (refinementMap.find(index))->second;

        // implement true bulk criterion for elliptic problems
        if ((eta[index]>=(1.0-refine)*globalEstimate) && (level < max))
        {
            grid.mark( 1, e);
        }
    }

private:

    /*! @brief Communicate the number of Elems each processor wants to mark,
     * so they can shift the quota of refined Elems to where they
     * are the most useful.
     */
    void communicate ()
    {
        using namespace Dune::PDELab;
        globalEstimate = grid.comm().sum(localEstimate);
        //! @todo get rid of ghosts and overlap
        int    localNumElem = grid.leafView().size(0);
        int    globalNumElem  = grid.comm().sum(localNumElem);

        const RT estimateQuotient = localEstimate / globalEstimate;
        refineNumber = localNumElem - estimateQuotient * refine * globalNumElem;
    }

    Grid& grid;
    const GFSU& gfsu;
    Estimation& estimation;
    const RT refine, coarsen;
    int refineNumber, coarsenNumber;
    const int min, max;
    const bool doCommunicate;
    MapType refinementMap;
    RT localEstimate;
    std::vector<RT> eta;
    RT globalEstimate;
    Dune::GeometryType gt;
};




