# Functional Encryption for [[alloc] init]'s Cryptography Suite

[![Build Status](https://travis-ci.com/NilFoundation/hash.svg?branch=master)](https://travis-ci.com/NilFoundation/hash)

Functional Encryption for [[alloc] init]'s cryptography suite offering different
state-of-the-art implementations of functional encryption schemes, specifically
FE schemes for _linear_ polynomials (e.g. _inner products_).

## Building

This library uses Boost CMake build modules (https://github.com/BoostCMake/cmake_modules.git). To actually include this
library in a project it is required to:

1. Add [CMake Modules](https://github.com/BoostCMake/cmake_modules.git) as submodule to target project repository.
2. Add all the internal dependencies using [CMake Modules](https://github.com/BoostCMake/cmake_modules.git) as
   submodules to target project repository.
3. Initialize parent project with [CMake Modules](https://github.com/BoostCMake/cmake_modules.git) (Look
   at [crypto3](https://github.com/alloc-init/crypto3.git) for the example)

## Dependencies

### Internal

* [Algebra](https://github.com/alloc-init/crypto3-algebra.git)
* [Block Ciphers](https://github.com/alloc-init/crypto3-block.git)

### External
* [Boost](https://boost.org) (>= 1.74)
