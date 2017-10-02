from __future__ import absolute_import, division, print_function, unicode_literals

import dune.common.checkconfiguration as checkconfiguration

registry = { "grid": dict() }

try:
    checkconfiguration.assertHave("HAVE_DUNE_SPGRID")

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

    registry["grid"]["SPBisection"] = spBisectionGrid
    registry["grid"]["SPIsotropic"] = spIsotropicGrid

except checkconfiguration.ConfigurationError:
    pass

if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS)
