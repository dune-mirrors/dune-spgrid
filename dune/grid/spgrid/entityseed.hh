#ifndef DUNE_SPGRID_ENTITYSEED_HH
#define DUNE_SPGRID_ENTITYSEED_HH

#include <cassert>

#include <dune/common/typetraits.hh>

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
    typedef typename std::remove_const< Grd >::type Grid;

    /** \internal \brief type of reference cube */
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

    /** \internal \brief type of multi index */
    typedef typename ReferenceCube::MultiIndex MultiIndex;

    /** \brief default constructor
     *
     *  \note A default constructed entity seed is invalid.
     */
    SPEntitySeed ()
      : level_( -1 )
    {}

    /** \internal
     *  \brief constructor
     *
     *  \param[in]  level            level of the seeded entity
     *  \param[in]  id               multi index of the seeded entity
     *  \param[in]  partitionNumber  number of the partition, the seeded entity
     *                               belongs to
     */
    SPEntitySeed ( const int level, const MultiIndex &id, const unsigned int partitionNumber )
      : level_( level ), id_( id ), partitionNumber_( partitionNumber )
    {}

    /** \brief check whether this seed generates a valid entity */
    bool isValid () const { return (level_ >= 0); }

    /** \internal
     *  \brief obtain level of the seeded entity
     *
     *  \returns the level of the seeded entity
     */
    int level () const
    {
      assert( isValid() );
      return level_;
    }

    /** \internal
     *  \brief obtain the multi index of the seeded entity
     * 
     *  \returns the multi index of the seeded entity
     */
    MultiIndex id () const
    {
      assert( isValid() );
      return id_;
    }

    /** \internal
     *  \brief obtain number of the partition, the seeded entity belongs to
     *
     *  \returns the number of the partition, the seeded entity belongs to
     */
    unsigned int partitionNumber () const
    {
      assert( isValid() );
      return partitionNumber_;
    }

  private:
    int level_;
    MultiIndex id_;
    unsigned int partitionNumber_;
  };

} // namespace Dune

#endif // #ifndef DUNE_SPGRID_ENTITYSEED_HH
