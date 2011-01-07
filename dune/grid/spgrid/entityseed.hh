#ifndef DUNE_SPGRID_ENTITYSEED_HH
#define DUNE_SPGRID_ENTITYSEED_HH

/** \file
 *  \author Martin Nolte
 *  \brief  entity seed for \ref Dune::SPGrid "SPGrid"
 */

namespace Dune
{

  /** \brief entity seed for \ref Dune::SPGrid "SPGrid"
   *
   *  The entity seed contains the minimal information required to rebuild the
   *  entity, if the grid is still known.
   *  It should be preferred over the entity pointer when storing entities.
   *  
   *  \note While rebuilding an entity pointer from an entity seed is assumed to
   *        be cheap, it introduces non-negligible overhead.
   *        If a large quantity of the grid is to be stored and the order is of
   *        no import, filtering an iterator might be a better option.
   *
   *  \tparam  codim  codimension of the seeded entity
   *  \tparam  Grd    Grid (must be a \ref Dune::SPGrid "SPGrid")
   */
  template< int codim, class Grd >
  struct SPEntitySeed
  {
    /** \brief type of grid this entity seed belongs to */
    typedef typename remove_const< Grd >::type Grid;

    typedef typename Grid::Traits::ReferenceCube ReferenceCube;

    /** \brief dimension of the grid */
    static const int dimension = ReferenceCube::dimension;
    /** \brief codimension of the seeded entity */
    static const int codimension = codim;
    /** \brief dimension of the seeded entity */
    static const int mydimension = dimension - codimension;
    /** \brief world dimension of the grid */
    static const int dimensionworld = dimension;

    /** \brief type of the seeded entity */
    typedef typename Grid::Traits::template Codim< codimension >::Entity Entity;

    typedef typename ReferenceCube::MultiIndex MultiIndex;

    /** \brief constructor
     *
     *  \note This constructor is an implementational detail.
     *
     *  \param[in]  level            level of the seeded entity
     *  \param[in]  id               internal id of the seeded entity
     *  \param[in]  partitionNumber  number of the partition, the seeded entity
     *                               belongs to
     */
    SPEntitySeed ( const int level, const MultiIndex &id, const unsigned int partitionNumber )
    : level_( level ), id_( id ), partitionNumber_( partitionNumber )
    {}

    /** \brief obtain level of the seeded entity
     *
     *  \note This method is an implementational detail.
     *
     *  \returns the level of the seeded entity
     */
    int level () const
    {
      return level_;
    }

    /** \brief obtain the internal id of the seeded entity
     * 
     *  \note This method is an implementational detail.
     *
     *  \returns the internal id of the seeded entity
     */
    MultiIndex id () const
    {
      return id_;
    }

    /** \brief obtain number of the partition, the seeded entity belongs to
     *
     *  \note This method is an implementational detail.
     *
     *  \returns the number of the partition, the seeded entity belongs to
     */
    unsigned int partitionNumber () const
    {
      return partitionNumber_;
    }

  private:
    int level_;
    MultiIndex id_;
    unsigned int partitionNumber_;
  };

} // namespace Dune

#endif // #ifndef DUNE_SPGRID_ENTITYSEED_HH
