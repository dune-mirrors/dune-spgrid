#ifndef DUNE_GRID_SPGRID_PICKLE_HH
#define DUNE_GRID_SPGRID_PICKLE_HH

#include <tuple>

#if HAVE_DUNE_COREPY
#include <dune/corepy/common/pickle.hh>
#endif // #if HAVE_DUNE_COREPY

#include <dune/grid/spgrid/gridview.hh>

namespace Dune
{

  namespace CorePy
  {

#if HAVE_DUNE_COREPY
    template< class G >
    struct Pickler< GridView< SPGridViewTraits< G > >, void >
    {
      typedef GridView< SPGridViewTraits< G > > Self;

      typedef typename Self::Grid Grid;
      typedef std::tuple< const Grid &, int > State;

      static State getState ( const Self &self )
      {
        if( &self.indexSet() == &self.grid().leafGridView().indexSet() )
          return State( self.grid(), -1 );
        else
          return State( self.grid(), self.impl().gridLevel().level() );
      }

      static void setState ( Self &self, const State &state )
      {
        if( std::get< 1 >( state ) < 0 )
          new (&self) Self( std::get< 0 >( state ).leafGridView() );
        else
          new (&self) Self( std::get< 0 >( state ).levelGridView( std::get< 1 >( state ) ) );
      }
    };
#endif // #if HAVE_DUNE_COREPY

  } // namespace CorePy

} // namespace Dune

#endif // #ifndef DUNE_GRID_SPGRID_PICKLE_HH
