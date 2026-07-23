# [[alloc] init] C++ Cryptography Suite

[![Twitter](https://img.shields.io/twitter/follow/alloc_init_)](https://twitter.com/alloc_init_)
[![Telegram](https://img.shields.io/badge/Telegram-2CA5E0?style=flat-square&logo=telegram&logoColor=dark)](https://t.me/alloc_init)

Crypto3 is a modular, header-only C++ cryptography suite. It provides generic,
STL-style interfaces for cryptographic primitives and the algebraic and
mathematical structures used to build them.

Initially developed by [=nil; Crypto3](https://crypto3.nil.foundation), the
project is now maintained by [[[alloc] init]](https://allocin.it). Guides are
available under [`docs/manual`](docs/manual/).

Continuous integration currently covers Ubuntu x86-64 and macOS. Other
operating systems and architectures are not continuously tested.

## Repository Structure

Crypto3 is maintained as a monorepo. The source trees under `libs` are regular
directories in this repository; `cmake/modules` is the only Git submodule.

```text
root
|-- cmake
|   `-- modules              BoostCMake helper submodule
|-- docs                     Doxygen configuration and guides
`-- libs
    |-- algebra              Finite fields, curves, and pairings
    |-- benchmark_tools      Shared benchmark utilities
    |-- block                Block ciphers
    |-- codec                Encoding and decoding algorithms
    |-- containers           Merkle trees and related containers
    |-- hash                 Hash functions
    |-- kdf                  Key derivation functions, including PBKDF2
    |-- mac                  Message authentication codes
    |-- marshalling          Algebra, container, math, and numeric marshalling
    |-- math                 Polynomial arithmetic and FFT algorithms
    |-- modes                Cipher mode headers and tests
    |-- multiprecision       Extended multiprecision types and utilities
    |-- parallelization-utils
    |-- pkpad                Public-key padding schemes
    |-- pubkey               Public-key schemes and secret sharing
    |-- random               Randomization primitives
    `-- stream               Stream ciphers
```

## Requirements

- CMake 3.22 or newer
- A C++23-capable C++ compiler
- Boost with the `container`, `random`, `filesystem`, `log`, `log_setup`,
  `program_options`, `thread`, `unit_test_framework`, and `timer` components
- Git
- Doxygen and Graphviz when building the API documentation

The build does not enforce a specific Boost version. CI currently builds on
Ubuntu 26.04 with its distribution compiler and on macOS 26 with Homebrew LLVM.

## Clone And Build

The submodule URL uses GitHub SSH. If you do not have GitHub SSH credentials,
configure this checkout to fetch GitHub dependencies over HTTPS:

```sh
git clone https://github.com/alloc-init/crypto3.git
cd crypto3
git config url."https://github.com/".insteadOf git@github.com:
git submodule update --init --recursive

cmake -S . -B build -DBUILD_TESTS=ON
cmake --build build --target tests --parallel
ctest --test-dir build --output-on-failure
```

Useful CMake options include:

- `BUILD_TESTS`: configure tests; defaults to `ON`
- `BUILD_BENCH_TESTS`: configure performance benchmarks; defaults to `OFF`
- `BUILD_DOCS`: configure Doxygen documentation support; defaults to `OFF`
- `BUILD_EXAMPLES`: build component examples; defaults to `OFF`
- `BUILD_WITH_TARGET_ARCHITECTURE`: override automatic target architecture detection

Generate the API documentation under `build/docs` with:

```sh
cmake -S . -B build -DBUILD_DOCS=ON
cmake --build build --target docs
```

On Apple Silicon, install the dependencies with `brew install boost llvm` and
select Homebrew Clang when configuring:

```sh
cmake -S . -B build \
  -DBUILD_TESTS=ON \
  -DCMAKE_CXX_COMPILER=/opt/homebrew/opt/llvm/bin/clang++
```

## Usage

Crypto3 components are header-only CMake `INTERFACE` libraries. In a CMake
project that includes Crypto3, link only the components you use:

```cmake
add_subdirectory(path/to/crypto3)
target_link_libraries(my_target PRIVATE crypto3::hash crypto3::algebra)
```

Component dependencies are expressed by their CMake targets. Buildable examples
are located in component-specific `libs/*/example` directories. See the
[quickstart](docs/manual/quickstart.md) for an in-tree example.

## Contributing

See the [contribution guidelines](docs/manual/contributing.md).

## Support

- Email: [nemo@allocin.it](mailto:nemo@allocin.it)
- Telegram: [@alloc-init](https://t.me/alloc-init)
- [GitHub issues](https://github.com/alloc-init/crypto3/issues)
- [GitHub Discussions](https://github.com/alloc-init/crypto3/discussions)

## License

The software is provided under the [MIT License](LICENSE).
