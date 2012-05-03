#ifndef DUNE_SPGRID_DECLARATION_HH
#define DUNE_SPGRID_DECLARATION_HH

namespace Dune
{

  // SPRefinementStrategy
  // --------------------

  /** \enum  SPRefinementStrategy
   *  \brief possible refinement strategies for \ref Dune::SPGrid "SPGrid"
   */
  enum SPRefinementStrategy
  {
    /** Each element is split into 2<sup>dim</sup> children.
     *  This is the default refinement technique.
     */
    SPIsotropicRefinement,
    /** The user may choose freely along which axes the grid should be refined.
     *  By default, this coincides with SPIsotropicRefinement.
     */
    SPAnisotropicRefinement,
    /** Each element is split into 2 children.
     *  The axis along which the elements are split may be chosen by the user.
     *  By default, the axes are cycled periodically.
     */
    SPBisectionRefinement
  };



  // External Forward Declaration
  // ----------------------------

  template< class, int, SPRefinementStrategy, class >
  class SPGrid;

} // namespace Dune

#endif // #ifndef DUNE_SPGRID_DECLARATION_HH
