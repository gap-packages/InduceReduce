This file describes changes in the `InduceReduce` package.

## Unreleased

- The test suite no longer assumes an order of the conjugacy classes of a
  character table, or that `Irr` selects the method of this package.
  Both assumptions fail when the `CCTab` package is loaded.
- Fixed `IrrUnger( G, opts )` and `CharacterTableUnger( G, opts )`, which
  returned wrong characters or raised an error when called a second time
  for the same group. Every call with options now runs the algorithm anew.
  `UngerRecord( G )` is now the finished record of the run with the default
  options.

## 1.4 (2026-07-31)

- Dropped the dependency on the `GAPDoc` package, which was only needed
  for building the manual.
- The test suite no longer requires the `SmallGrp` package.

## 1.3 (2025-10-16)

- Set power maps in the character tables whose irreducibles have been
  computed by the `Irr` method that calls `IrrUnger`.

## 1.2 (2025-10-15)

- In the `Irr` method that calls `IrrUnger`,
  delegate to `IrrBaumClausen` if the group is abelian-by-supersolvable;
  this is usually faster.

## 1.1 (2025-07-16)

- Added `IrrUnger` function which is also installed as method for `Irr`,
  so that the code in this package is automatically used when applicable
- Some code cleanup & refactoring, minor speed ups for some functionality
- Fixed the test suite and set up continuous integration tests
- Added Thomas Breuer and Max Horn as co-maintainers
- Moved package to `gap-packages` organization
- Various janitorial updates

## 1.0 (2018-09-04)

- First public release
