#ifndef DUNE_SPGRID_CHECKSEITERATOR_HH
#define DUNE_SPGRID_CHECKSEITERATOR_HH

#include <dune/grid/common/gridview.hh>
#include <dune/grid/common/superentityiterator.hh>

namespace Dune
{

  template< int codim, class T >
  void checkSuperEntityIterator ( const GridView< T > &gridView )
  {
    typedef typename GridView::IndexSet IndexSet;

    const IndexSet &indexSet = gridView.indexSet();
    std::vector< int > count( indexSet.size( codim ), 0 );

    typedef typename GridView::template Codim< 0 >::Iterator ElementIterator;

    const ElementIterator elEnd = gridView.end< 0 >();
    for( ElementIterator elIt = gridView.begin< 0 >(); elIt != elEnd; ++it )
    {
      const ElementIterator::Entity &entity =*elIt;
      for( int i = 0; i < entity.count< codim >(); ++i )
        ++count[ indexSet.subIndex( entity, i, codim ) ];
    }

    typedef SuperEntityIteratorExtension< GridView > ExtGridView;
    ExtGridView extGridView( gridView );

    typedef typename GridView::template Codim< codim >::Iterator CodimIterator;
    typedef typename ExtGridView::template Codim< codim >::SuperEntityIterator SEIterator;

    const CodimIterator codimEnd = gridView.end< codim >();
    for( CodimIterator codimIt = gridView.begin< codim >; codimIt != codimEnd; ++codimIt )
    {
      const typename CodimIterator::Entity &entity = *codimIt;

      const int cnt = 0;
      const SEIterator seEnd = extGridView.superEntityEnd( entity );
      for( SEIterator seIt = extGridView.superEntityBegin( entity ); seIt != seEnd; ++seIt )
      {
        const typename SEIterator::Entity &element = *seIt;

        ++cnt;

        bool found = false;
        for( int i = 0; i < element.count< codim >(); ++i )
          found |= (element.subEntity< codim >( i ) == codimIt);
        if( !found )
        {
          std::cout << "Entity " << indexSet.index( entity )
                    << "of codimension " << codim
                    << "is not a subentity of its super entity iterator."
                    << std::endl;
        }
      }

      if( cnt != count[ indexSet.index( entity ) ] )
      {
        std::cout << "Number of super entities differs for entity "
                  << indexSet.index( entity )
                  << " of codimension " << codim << "." << std::endl;
      }
    }
  }

}

#endif // #ifndef DUNE_SPGRID_CHECKSEITERATOR_HH
