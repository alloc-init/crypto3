# [[alloc] init] C++ Cryptography Suite

[![Twitter](https://img.shields.io/twitter/follow/alloc_init_)](https://twitter.com/alloc_init_)
[![Telegram](https://img.shields.io/badge/Telegram-2CA5E0?style=flat-square&logo=telegram&logoColor=dark)](https://t.me/alloc_init)

Crypto3 cryptography suite's purpose is:
1. To provide a secure, fast and architecturally clean C++ generic cryptography schemes implementation.
2. To provide a developer-friendly, modular suite, usable for novel schemes implementation and further
   extension.
3. To provide a Standard Template Library-alike C++ interface and concept-based architecture implementation.

Libraries are designed to be state of the art, highly performant and providing a one-stop solution for
all cryptographic operations. They are supported on Linux operating system and architectures (x86/ARM).

Developed and supported by [=nil; Foundation](https://nil.foundation).

## Contents
1. [Structure](#structure)
2. [Build & test](#build_&_test)
3. [Usage](#usage)

## Structure
This folder contains the whole suite. Single-purposed libraries (e.g. [algebra
](https://github.com/NilFoundation/placeholder/tree/master/crypto3/libs/algebra) or
[hash](https://github.com/NilFoundation/placeholder/tree/master/crypto3/libs/hash)) are not advised to be
used outside this suite or properly constructed CMake project and should be handled with great care.

```
crypto3
├── benchmarks
├── cmake: cmake sub-module with helper functions/macros to build crypto3 library umbrella-repository
├── docs: documentation, tutorials and guides
├── libs
│   ├── algebra: algebraic operations and structures being used for elliptic-curve cryptography
│   ├── benchmark_tools: utilities to run benchmarks
│   ├── blueprint: components and circuits for zk schemes
│   ├── containers: containers and generic commitment schemes for accumulating data, includes Merkle Tree
│   ├── hash: hashing algorithms
│   ├── marshalling: marshalling libraries for types in crypto3 library
│   ├── math: set of Fast Fourier Transforms evaluation algorithms and Polynomial Arithmetics
│   ├── multiprecision: integer, rational, floating-point, complex and interval number types. 
│   ├── random: randomisation primitives 
│   ├── transpiler
│   └── zk: zk cryptography schemes
```

## Build & test

WRITE ME for cmake & docker

## Usage

The suite is used as header-only libraries.

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

