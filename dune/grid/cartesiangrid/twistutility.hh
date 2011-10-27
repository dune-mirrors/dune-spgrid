#ifndef DUNE_CARTESIANGRID_TWISTUTILITY_HH
#define DUNE_CARTESIANGRID_TWISTUTILITY_HH

#include <cassert>

#if HAVE_DUNE_FEM
#include <dune/grid/cartesiangrid/hostgridaccess.hh>

#include <dune/fem/quadrature/caching/twistutility.hh>

namespace Dune
{

  // TwistUtility for CartesianGrid 
  // ------------------------------
  
  /** \brief Specialization of TwistUtility for CartesianGrid */
  template< class HostGrid >
  struct TwistUtility< CartesianGrid< HostGrid > >
  {
    typedef CartesianGrid< HostGrid > GridType;

    typedef typename GridType::Traits::LeafIntersectionIterator LeafIntersectionIterator;
    typedef typename LeafIntersectionIterator::Intersection LeafIntersection;
    typedef typename GridType::Traits::LevelIntersectionIterator LevelIntersectionIterator;
    typedef typename LevelIntersectionIterator::Intersection LevelIntersection;

    static const int dimension = GridType::dimension;

  private:
    typedef Dune::HostGridAccess< GridType > HostGridAccess;
    typedef TwistUtility< HostGrid > HostTwistUtility;

  public:
    //! \brief return twist for inner face 
    template< class Intersection >
    static int twistInSelf ( const GridType &grid, const Intersection &intersection )
    {
      return HostTwistUtility::twistInSelf( HostGridAccess::hostGrid( grid ), HostGridAccess::hostIntersection( it ) );
    }
    
    //! \brief return twist for outer face 
    template< class Intersection >
    static int twistInNeighbor ( const GridType &grid, const Intersection &intersection )
    {
      return HostTwistUtility::twistInNeighbor( HostGridAccess::hostGrid( grid ), HostGridAccess::hostIntersection( intersection ) );
    }

    /** \brief return element geometry type of inside or outside entity */
    template< class Intersection >
    static inline GeometryType
    elementGeometry ( const Intersection &intersection, const bool inside )
    {
      return HostTwistUtility::elementGeometry( HostGridAccess::hostIntersection( intersection ), inside );
    }
  };
  
} // namespace Dune 

#endif // #if HAVE_DUNE_FEM

#endif // #ifndef DUNE_TWISTUTILITY_HH
