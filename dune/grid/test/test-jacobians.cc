#include <config.h>

#include <type_traits>
#include <utility>

#include <dune/common/fmatrix.hh>
#include <dune/common/hybridutilities.hh>

#include <dune/grid/spgrid/geometrycache.hh>

int main ( int argc, char *argv[] )
try
{
  Dune::Hybrid::forEach( std::make_integer_sequence< int, 4 >(), [] ( auto &&codim ) {
      Dune::SPJacobianTransposed< double, 3, std::decay_t< decltype( codim ) >::value > jt;
      Dune::FieldMatrix< double, std::decay_t< decltype( codim ) >::value, 3 > fm_jt( jt );
      if( fm_jt.frobenius_norm2() > 0.0 )
        DUNE_THROW( Dune::Exception, "Default initialized JacobianTransposed is not zero." );

      Dune::SPJacobianInverseTransposed< double, 3, std::decay_t< decltype( codim ) >::value > jit;
      Dune::FieldMatrix< double, 3, std::decay_t< decltype( codim ) >::value > fm_jit( jit );
      if( fm_jt.frobenius_norm2() > 0.0 )
        DUNE_THROW( Dune::Exception, "Default initialized JacobianInverseTransposed is not zero." );
    } );

  return 0;
}
catch( const Dune::Exception &e )
{
  std::cerr << e << std::endl;
  return 1;
}
