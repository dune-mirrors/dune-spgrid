#ifndef DUNE_SPGRID_DECOMPOSITION_HH
#define DUNE_SPGRID_DECOMPOSITION_HH

#include <algorithm>

#include <dune/grid/spgrid/mesh.hh>
#include <dune/grid/spgrid/multiindex.hh>

namespace Dune
{

  // SPDecomposition
  // ---------------

  template< int dim >
  class SPDecomposition
  {
    typedef SPDecomposition< dim > This;

  public:
    static const int dimension = dim;

    typedef SPMultiIndex< dimension > MultiIndex;
    typedef SPMesh< dimension > Mesh;

  private:
    struct Node
    {
      Node ( const Mesh &mesh, const unsigned int size );
      ~Node ();

      const Mesh &mesh () const;
      const Mesh &subMesh ( const unsigned int rank ) const;
      void subMeshes ( std::vector< Mesh > &meshes ) const;

      unsigned int size () const;

    private:
      Mesh mesh_;
      unsigned int size_;
      Node *left_, *right_;
    };

  public:
    SPDecomposition ( const Mesh &mesh, const unsigned int size );
    SPDecomposition ( const MultiIndex &width, const unsigned int size );

    const Mesh &mesh () const;
    const Mesh &subMesh ( const unsigned int rank ) const;
    std::vector< Mesh > subMeshes () const;

    unsigned int size () const;

  private:
    Node root_;
  };



  // Implementation of SPDecomposition::Node
  // ---------------------------------------

  template< int dim >
  inline SPDecomposition< dim >::Node::Node ( const Mesh &mesh, const unsigned int size )
  : mesh_( mesh ),
    size_( size ),
    left_( 0 ),
    right_( 0 )
  {
    if( size_ > 1 )
    {
      const int leftWeight = size_/2;
      const int rightWeight = size_ - leftWeight;

      const MultiIndex &width = mesh.width();
      const std::pair< Mesh, Mesh > split
        = mesh_.split( std::max_element( width.begin(), width.end() ) - width.begin(), leftWeight, rightWeight );
      left_ = new Node( split.first, leftWeight );
      right_ = new Node( split.second, rightWeight );
    }
  }


  template< int dim >
  inline SPDecomposition< dim >::Node::~Node ()
  {
    delete left_;
    delete right_;
  }


  template< int dim >
  inline const typename SPDecomposition< dim >::Mesh &
  SPDecomposition< dim >::Node::mesh () const
  {
    return mesh_;
  }


  template< int dim >
  inline const typename SPDecomposition< dim >::Mesh &
  SPDecomposition< dim >::Node::subMesh ( const unsigned int rank ) const
  {
    assert( rank < size_ );
    if( size_ > 1 )
    {
      assert( (left_ != 0) && (right_ != 0) );
      if( rank < size_/2 )
        return left_->subMesh( rank );
      else
        return right_->subMesh( rank - size_/2 );
    }
    else
      return mesh();
  }


  template< int dim >
  inline void
  SPDecomposition< dim >::Node::subMeshes ( std::vector< Mesh > &meshes ) const
  {
    if( size_ > 1 )
    {
      assert( (left_ != 0) && (right_ != 0) );
      left_->subMeshes( meshes );
      right_->subMeshes( meshes );
    }
    else
      meshes.push_back( mesh() );
  }


  template< int dim >
  inline unsigned int SPDecomposition< dim >::Node::size () const
  {
    return size_;
  }



  // Implementation of SPDecomposition
  // ---------------------------------

  template< int dim >
  inline SPDecomposition< dim >
    ::SPDecomposition ( const Mesh &mesh, const unsigned int size )
  : root_( mesh, size )
  {}


  template< int dim >
  inline SPDecomposition< dim >
    ::SPDecomposition ( const MultiIndex &width, const unsigned int size )
  : root_( Mesh( width ), size )
  {}


  template< int dim >
  inline const typename SPDecomposition< dim >::Mesh &
  SPDecomposition< dim >::mesh () const
  {
    return root_.mesh();
  }


  template< int dim >
  inline const typename SPDecomposition< dim >::Mesh &
  SPDecomposition< dim >::subMesh ( const unsigned int rank ) const
  {
    return root_.subMesh( rank );
  }


  template< int dim >
  inline std::vector< typename SPDecomposition< dim >::Mesh >
  SPDecomposition< dim >::subMeshes () const
  {
    std::vector< Mesh > meshes;
    meshes.reserve( root_.size() );
    root_.subMeshes( meshes );
    return meshes;
  }


  template< int dim >
  inline unsigned int SPDecomposition< dim >::size () const
  {
    return root_.size();
  }

} // namespace Dune

#endif // #ifndef DUNE_SPGRID_DECOMPOSITION_HH
