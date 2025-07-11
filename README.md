# [[alloc] init] C++ Cryptography Suite

[![Twitter](https://img.shields.io/twitter/follow/alloc_init_)](https://twitter.com/alloc_init_)
[![Telegram](https://img.shields.io/badge/Telegram-2CA5E0?style=flat-square&logo=telegram&logoColor=dark)](https://t.me/alloc_init)

Crypto3 cryptography suite's purpose is:

1. To provide a secure, fast and architecturally clean C++ generic cryptography schemes implementation.
2. To provide a developer-friendly, modular suite, usable for novel schemes implementation and further
   extension.
3. To provide a Standard Template Library-alike C++ interface and concept-based architecture implementation.

Libraries are designed to be state of the art, highly performant and providing a one-stop solution for
all cryptographic operations. They are supported on all operating systems (*nix, windows, macOS)
and architectures(x86/ARM).

Initially developed by [=nil; Crypto3](https://crypto3.nil.foundation), part
of [=nil; Foundation](https://nil.foundation) and now supported by [[[alloc] init]](https://allocin.it).

Rationale, tutorials and references are available [here](https://docs.allocin.it/crypto3)

## Contents

1. [Repository Structure](#repository-structure)
2. [Installation](#installation)
3. [Usage](#usage)
3. [Contributing](#contributing)
4. [Community](#community)

## Repository Structure

This repository is an umbrella-repository for the whole suite. Single-purposed libraries repositories (e.g. [block
](https://github.com/alloc-init/block) or [hash](https://github.com/alloc-init/hash)) are not advised to be
used outside this suite or properly constructed CMake project and should be handled with great care.

```
root
├── cmake: cmake sub-module with helper functions/macros to build crypto3 library umbrella-repository
├── docs: documentation , tutorials and guides
├── libs: all directories added as submodules which are independent projects.
│   ├── algebra: algebraic operations and structures being used for elliptic-curve cryptography
│   ├── block: block ciphers
│   ├── blueprint: components and circuits for zk schemes
│   ├── codec: encoding/decoding algorithms
│   ├── containers: containers and generic commitment schemes for accumulating data, includes Merkle Tree
│   ├── hash: hashing algorithms
│   ├── mac: message authentication codes
│   ├── marshalling: marshalling libraries for types in crypto3 library
│   ├── math: set of Fast Fourier Transforms evaluation algorithms and Polynomial Arithmetics
│   ├── modes: cipher modes
│   ├── multiprecision: integer, rational, floating-point, complex and interval number types. 
│   ├── passhash: password hashing operations 
│   ├── pbkdf: password based key derivation functions
│   ├── pkmodes: threshold, aggregation modes for public key schemes
│   ├── pkpad: padding module for public key schemes
│   ├── pubkey: pubkey signing APIs
│   ├── random: randomisation primitives 
│   ├── stream: stream ciphers
│   ├── vdf: verifiable delay functions 
│   ├── zk: zk cryptography schemes
```

## Installation

### Dependencies

- [clang](https://clang.llvm.org/) (>= 17.0)/GCC (>= 10.0)/MSVC (>= 14.20)
- [cmake](https://cmake.org) (>= 3.6)
- [boost](https://boost.org) (>= 1.87)
- [cmake_modules](https://github.com/BoostCMake/cmake_modules) (57639741ecf018835deb97a04db2200241d7fbd3)

### Apple Silicon (Note)

For apple silicone, some of the requirements differ:

Homebrew packages:
* Required: you need to download the latest clang from homebrew `brew install llvm`
* Optional: Doxygen (if you want to use it) `brew install doxygen`
* Optional: GraphViz (needed for doxygens) `brew install graphviz`

#### Tested on the following systems:
* Apple Silicon: it's tested on:Homebrew clang version 20.1.5

### Clone & Build

```
git clone --recurse-submodules https://github.com/alloc-init/crypto3.git 
cd crypto3 && mkdir build && cd build
cmake ..
make tests
```

> Note that you might need to set `-DCMAKE_CXX_COMPILER` to point to an up-to-date clang and `-DBOOST_ROOT` with  `-DBoost_NO_SYSTEM_PATHS=ON` to point to correct version of boost.

## Nix support

This repository provides Nix flake, so once you have installed Nix with flake support, you can use single command to fetch all the dependencies and build:

```bash
nix build
```

To activate Nix development environment:

```bash
nix develop
```

To run all tests:

```bash
nix flake check
```

To run single test:

```bash
nix develop . -c cmake -B build -DCMAKE_CXX_STANDARD=20 -DCMAKE_BUILD_TYPE=Debug -DBUILD_SHARED_LIBS=FALSE -DCMAKE_ENABLE_TESTS=TRUE -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_BUILD_TYPE=Debug  -DCMAKE_CXX_FLAGS=-ggdb
nix develop -c cmake --build build -t <test_target> // for example multiprecision_modular_adaptor_fixed_test
```

## Usage

Cryptography suite can be used as follows:

1. Generic.
2. Selective.

The suite is used as a header-only and is currently statically linked. Future versions will allow dynamic linking.

#### Generic

Generic usage of cryptography suite consists of all modules available at
[GitHub =nil; Crypto3 Team Repositories](https://github.com/orgs/NilFoundation/teams/nil-crypto3/repositories).
This is an umbrella-repository where Modules
are added as submodules emplaced in `libs` directory. A developer can thus add this  
project as a submodule and would not need to resolve dependencies.
See [crypto3-template](https://github.com/alloc-init/crypto3-template) as an example of usage.

The generic module can be added to your c++ project as follows

``` git submodule add https://github.com/alloc-init/crypto3.git <dir>```

### Selective

Developer can select to include a one or more modules to reduce the sources of resulting project and dependencies tree
height. This however
does require the developer to manually resolve all required dependencies and stay upto date regarding
compatibilities across modules.

Example of such embedding is =nil; Foundation's [Actor Library](https://github.com/alloc-init/actor). It uses only
[hashes](https://github.com/alloc-init/hash) so the dependency graph requires
for the project to submodule [block ciphers library](https://github.com/alloc-init/block) and optional
[codec library](https://github.com/alloc-init/codec) for testing purposes. So,
the root Actor repository has only related libraries submoduled:
[block](https://github.com/alloc-init/mtl/libs/block),
[codec](https://github.com/alloc-init/mtl/libs/codec) and
[hash](https://github.com/alloc-init/mtl/hash).

Selective modules can be added to your project as follows:

``` git submodule add https://github.com/alloc-init/crypto3-<lib>.git <dir>```

## Contributing

See [contributing](./docs/manual/contributing.md) for contribution guidelines.

## Support

This cryptography suite is maintained by [[alloc] init], which can be contacted in several ways:

* E-Mail. Just drop a line to [nemo@allocin.it](mailto:nemo@allocin.it).
* Telegram Group. Join our Telegram group [@alloc-init](https://t.me/alloc-init) and ask any question in there.

[//]: # ( * Discord [channel]&#40;https://discord.gg/KmTAEjbmM3&#41; for discussions.)

* Issue. Issue which does not belong to any particular module (or you just don't know where to put it) can be
  created in this repository. The team will answer that.
* Discussion Topic (proposal, tutorial request, suggestion, etc). Would be happy to discuss that in the repository's
  GitHub [Discussions](https://github.com/alloc-init/crypto3/discussions)

## Licence

The software is provided under [MIT](LICENSE) Licence.

