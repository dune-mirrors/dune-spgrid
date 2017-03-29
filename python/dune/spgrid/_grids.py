from __future__ import absolute_import, division, print_function, unicode_literals

import dune.common.checkconfiguration as checkconfiguration


def spBisectionGrid(constructor, dimgrid, ctype="double"):
    from ..grid.grid_generator import module

    checkconfiguration.have("HAVE_DUNE_SPGRID")

    typeName = "Dune::SPGrid< " + ctype + ", " + str(dimgrid) + ", Dune::SPBisectionRefinement >"
    includes = ["dune/grid/spgrid.hh", "dune/grid/spgrid/dgfparser.hh", "dune/grid/spgrid/pickle.hh"]
    gridModule = module(includes, typeName)

    return gridModule.LeafGrid(gridModule.reader(constructor))


def spIsotropicGrid(constructor, dimgrid, ctype="double"):
    from ..grid.grid_generator import module

    checkconfiguration.have("HAVE_DUNE_SPGRID")

    typeName = "Dune::SPGrid< " + ctype + ", " + str(dimgrid) + ", Dune::SPIsotropicRefinement >"
    includes = ["dune/grid/spgrid.hh", "dune/grid/spgrid/dgfparser.hh", "dune/grid/spgrid/pickle.hh"]
    gridModule = module(includes, typeName)

    return gridModule.LeafGrid(gridModule.reader(constructor))


if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS)
