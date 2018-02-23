# Master (will become 2.7)

# Release 2.6

- Optional support for the brand-new `dune-python` module was added.

- MPI communication has been hardened:
  - It detects inconsistencies between `dataHandle.size` and `dataHandle.gather`.
  - It checks whether `dataHandle.scatter` reads bytes the correct amoutn of bytes.
  - It makes the tag number independent of the template arguments `Grid` and `DataHandle`.

- The `SPGrid` class is now move constructible (but not copy constructible).

- The `SPGrid` class is no longer default constructible.

- `Dune::FieldTraits` is now specialized for `SPJacobianTransposed` and
  `SPJacobianInverseTransposed`.

- The `JacobianTransposed` and `JacobianInverseTransposed` implementations are now
  top-level classes called `SPJacobianTransposed` and `SPJacobianInverseTransposed`
  to allow for template specialization.

- The `JacobianTransposed` and `JacobianInverseTransposed` now provide a method
  `determinant` to be compatible with the `DenseMatrix` interface.
