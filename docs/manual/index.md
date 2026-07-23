# Crypto3 C++ Cryptography Suite {#mainpage}

@tableofcontents

Crypto3 is a modular, header-only C++ cryptography suite with generic,
STL-style interfaces. It includes cryptographic primitives as well as the
algebraic and mathematical structures used to construct them.

Continuous integration currently covers Ubuntu x86-64 and macOS. Other
operating systems and architectures are not continuously tested.

See the [contribution guidelines](contributing.md) to contribute to the project.

## Components

The monorepo contains components for algebra, block and stream ciphers, codecs,
containers, hashes, key derivation, message authentication, marshalling,
polynomial arithmetic, multiprecision arithmetic, public-key padding,
public-key schemes, randomization, and supporting utilities.

Source trees are maintained directly under `libs`; they are not Git submodules.
The repository's only Git submodule is `cmake/modules`.

API groups are listed on the generated [Modules page](modules.html). Buildable
examples are located in component-specific `libs/*/example` directories.

## Usage

Crypto3 components are CMake `INTERFACE` libraries and do not produce static or
shared library artifacts. A project that includes this repository can link the
required component targets:

```cmake
add_subdirectory(path/to/crypto3)
target_link_libraries(my_target PRIVATE crypto3::hash crypto3::algebra)
```

Dependencies between components are carried by their CMake targets.

## Next Steps

Follow the [Quickstart](quickstart.md) to build the project, run the tests, and
compile an example. Continue with the [BLS signing guide](intermediate.md) for a
small public-key example.
