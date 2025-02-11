# Zero-Knowledge Cryptography Schemes for [[alloc] init]'s Cryptography Suite

Zero-Knowledge cryptography schemes for [[alloc] init]'s cryptography suite. 
SNARK-alike schemes for now. More trivial Pedersen commitment schemes, STARKs, 
IOP-based SNARKs, Bulletproofs etc in future.

[![Run tests](https://github.com/alloc-init/crypto3-zk/actions/workflows/run_tests.yml/badge.svg)](https://github.com/alloc-init/crypto3-zk/actions/workflows/run_tests.yml)

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

* [Multiprecision](https://github.com/alloc-init/crypto3-multiprecision.git)
* [Algebra](https://github.com/alloc-init/crypto3-algebra.git)
* [FFT](https://github.com/alloc-init/crypto3-fft.git)

### External

* [Boost](https://boost.org) (>= 1.74)
