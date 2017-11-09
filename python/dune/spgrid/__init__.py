from __future__ import absolute_import, division, print_function, unicode_literals

from dune.common.checkconfiguration import assertHave, ConfigurationError

try:
    assertHave("HAVE_DUNE_SPGRID")
except ConfigurationError:
    raise ImportError("DUNE module dune-spgrid was not found.")


def spBisectionGrid(domain, dimgrid=None, ctype="double"):
    from ..grid.grid_generator import module, getDimgrid

    if dimgrid is None:
        dimgrid = getDimgrid(domain)

    typeName = "Dune::SPGrid< " + ctype + ", " + str(dimgrid) + ", Dune::SPBisectionRefinement >"
    includes = ["dune/grid/spgrid.hh", "dune/grid/spgrid/dgfparser.hh"]
    gridModule = module(includes, typeName)

    return gridModule.LeafGrid(gridModule.reader(domain))


def spIsotropicGrid(domain, dimgrid=None, ctype="double"):
    from ..grid.grid_generator import module, getDimgrid

    if dimgrid is None:
        dimgrid = getDimgrid(domain)

    typeName = "Dune::SPGrid< " + ctype + ", " + str(dimgrid) + ", Dune::SPIsotropicRefinement >"
    includes = ["dune/grid/spgrid.hh", "dune/grid/spgrid/dgfparser.hh"]
    gridModule = module(includes, typeName)

    return gridModule.LeafGrid(gridModule.reader(domain))


registry = {}
registry["grid"] = {
        "SPBisection": spBisectionGrid,
        "SPIsotropic": spIsotropicGrid
    }


if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS)
