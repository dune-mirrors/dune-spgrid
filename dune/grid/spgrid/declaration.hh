#ifndef DUNE_SPGRID_DECLARATION_HH
#define DUNE_SPGRID_DECLARATION_HH

namespace Dune
{

  // External Forward Declaration
  // ----------------------------

  template< int >
  class SPIsotropicRefinement;

  template< int >
  class SPAnisotropicRefinement;

  template< int >
  class SPBisectionRefinement;

  template< class, int, template< int > class, class >
  class SPGrid;

} // namespace Dune

#endif // #ifndef DUNE_SPGRID_DECLARATION_HH
