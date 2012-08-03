#ifndef DUNE_SPGRID_GEOMETRYREFERENCE_HH
#define DUNE_SPGRID_GEOMETRYREFERENCE_HH

#include <dune/common/typetraits.hh>

#include <dune/geometry/type.hh>

#include <dune/grid/spgrid/geometry.hh>

namespace Dune
{

  // Internal Forward Declarations
  // -----------------------------

  template< int mydim, int cdim, class Grid >
  class SPLocalGeometryReference;



  // FacadeOptions
  // -------------

  namespace FacadeOptions
  {

    template< int mydim, int cdim, class Grid >
    struct StoreGeometryReference< mydim, cdim, Grid, SPLocalGeometryReference >
    {
      static const bool v = false;
    };

  } // namespace FacadeOptions



  // SPLocalGeometryReference
  // ------------------------

  template< int mydim, int cdim, class Grid >
  class SPLocalGeometryReference
  {
    typedef SPLocalGeometryReference< mydim, cdim, Grid > This;

    typedef SPLocalGeometry< mydim, cdim, Grid > Implementation;

  public:
    static const int mydimension = Implementation::mydimension;
    static const int coorddimension = Implementation::coorddimension;

    typedef typename Implementation::ctype ctype;

    typedef typename Implementation::LocalCoordinate LocalCoordinate;
    typedef typename Implementation::GlobalCoordinate GlobalCoordinate;
    
    typedef typename Implementation::Jacobian Jacobian;
    typedef typename Implementation::JacobianTransposed JacobianTransposed;

    SPLocalGeometryReference ( const Implementation &impl )
    : impl_( &impl )
    {}

    GeometryType type () const { return impl().type(); }

    bool affine() const { return impl().affine(); }

    int corners () const { return impl().corners(); }
    GlobalCoordinate corner ( int i ) const { return impl().corner( i ); }
    GlobalCoordinate center () const { return impl().center(); }

    GlobalCoordinate global ( const LocalCoordinate &local ) const
    {
      return impl().global( local );
    }

    LocalCoordinate local ( const GlobalCoordinate &global ) const
    {
      return impl().local( global );
    }

    ctype integrationElement ( const LocalCoordinate &local ) const
    {
      return impl().integrationElement( local );
    }
   
    ctype volume () const { return impl().volume(); }

    const JacobianTransposed &jacobianTransposed ( const LocalCoordinate &local ) const
    {
      return impl().jacobianTransposed( local );
    }

    const Jacobian &jacobianInverseTransposed ( const LocalCoordinate &local ) const
    {
      return impl().jacobianInverseTransposed( local );
    }

    const Implementation &impl () const { return *impl_; }

  private:
    const Implementation *impl_;
  };



  // Definitions of GeometryReference
  // --------------------------------

  template< int mydim, int cdim, class Grid >
  const int SPLocalGeometryReference< mydim, cdim, Grid >::mydimension;

  template< int mydim, int cdim, class Grid >
  const int SPLocalGeometryReference< mydim, cdim, Grid >::coorddimension;

} // namespace Dune

#endif // #ifndef DUNE_SPGRID_GEOMETRYREFERENCE_HH
