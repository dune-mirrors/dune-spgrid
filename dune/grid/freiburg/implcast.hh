#ifndef DUNE_FREIBURG_IMPLCAST_HH
#define DUNE_FREIBURG_IMPLCAST_HH

#include <dune/grid/common/grid.hh>

namespace Dune
{

  // Internal Forward Declarations
  // -----------------------------

  template< class Interface >
  struct ImplCastHelper;



  // External Forward Declarations
  // -----------------------------

  template< int mydim, int cdim, class Grid, template< int, int, class > class Impl >
  class Geometry;

  template< int codim, int dim, class Grid, template< int, int, class > class Impl >
  class Entity;

  template< class Grid, template< class > class Impl >
  class Intersection;

  template< class ViewTraits >
  class GridView;



  // ImplCastDefaultHelper
  // ---------------------

  template< class Grid, class Interface, class I >
  struct ImplCastDefaultHelper
  : public GridDefaultImplementation< Grid::dimension, Grid::dimensionworld, typename Grid::ctype, typename Grid::GridFamily >
  {
    typedef GridDefaultImplementation< Grid::dimension, Grid::dimensionworld, typename Grid::ctype, typename Grid::GridFamily > Base;
    typedef I Impl;

    static const Impl &impl ( const Interface &interface )
    {
      return Base::getRealImplementation( interface );
    }

    static Impl &impl ( Interface &interface )
    {
      return Base::getRealImplementation( interface );
    }
  };



  // ImplCastHelper
  // --------------

  template< class Interface >
  struct ImplCastHelper< const Interface >
  : public ImplCastHelper< Interface >
  {
    typedef const typename ImplCastHelper< Interface >::Impl Impl;
  };

  template< int mydim, int cdim, class Grid, template< int, int, class > class Impl >
  struct ImplCastHelper< Geometry< mydim, cdim, Grid, Impl > >
  : public ImplCastDefaultHelper< Grid, Geometry< mydim, cdim, Grid, Impl >, Impl< mydim, cdim, Grid > >
  {};

  template< int codim, int dim, class Grid, template< int, int, class > class Impl >
  struct ImplCastHelper< Entity< codim, dim, Grid, Impl > >
  : public ImplCastDefaultHelper< Grid, Entity< codim, dim, Grid, Impl >, Impl< codim, dim, Grid > >
  {};

  template< class Grid, template< class > class Impl >
  struct ImplCastHelper< Intersection< Grid, Impl > >
  : public ImplCastDefaultHelper< Grid, Intersection< Grid, Impl >, Impl< Grid > >
  {};

  template< class ViewTraits >
  struct ImplCastHelper< GridView< ViewTraits > >
  : public ImplCastDefaultHelper< typename ViewTraits::Grid, GridView< ViewTraits >, typename ViewTraits::GridViewImp >
  {};



  // ImplType
  // --------

  template< class Interface >
  struct ImplType
  {
    typedef typename ImplCastHelper< Interface >::Impl Type;
  };



  // impl_cast
  // ---------

  template< class Interface >
  inline typename ImplType< Interface >::Type &impl_cast ( Interface &interface )
  {
    return ImplCastHelper< Interface >::impl( interface );
  }

}

#endif // #ifndef DUNE_FREIBURG_IMPLCAST_HH
