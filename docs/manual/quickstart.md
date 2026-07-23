\page Quickstart Quickstart

# Quickstart

This guide builds Crypto3, runs its tests, and builds the in-tree BLS example.

## Requirements

- CMake 3.22 or newer
- A C++23-capable C++ compiler
- Boost with the components listed in the [root README](../../README.md#requirements)
- Git

On Ubuntu, the required packages can be installed with:

```sh
sudo apt update
sudo apt install build-essential libboost-all-dev cmake git
```

On macOS, install Boost and LLVM with Homebrew:

```sh
brew install boost llvm
```

## Get The Source

```sh
git clone https://github.com/alloc-init/crypto3.git
cd crypto3
```

`cmake/modules` is a Git submodule whose configured URL uses GitHub SSH. If you
do not have GitHub SSH credentials, use an HTTPS rewrite for this checkout:

```sh
git config url."https://github.com/".insteadOf git@github.com:
git submodule update --init --recursive
```

## Build And Test

```sh
cmake -S . -B build -DBUILD_TESTS=ON
cmake --build build --target tests --parallel
ctest --test-dir build --output-on-failure
```

To select a compiler explicitly, add a compiler path during configuration:

```sh
cmake -S . -B build \
  -DBUILD_TESTS=ON \
  -DCMAKE_CXX_COMPILER=/path/to/clang++
```

## Build The BLS Example

Examples are disabled by default. Enable them and build the public-key BLS
example with:

```sh
cmake -S . -B build -DBUILD_EXAMPLES=ON
cmake --build build --target pubkey_bls_example --parallel
./build/libs/pubkey/example/pubkey_bls_example
```

The example creates a key, signs `hello world`, and asserts that verification
succeeds. It exits silently with status zero on success.

## Next Steps

- Read the [BLS signing guide](intermediate.md).
- Browse generated API groups on the [Modules page](modules.html).
- Explore component-specific examples under `libs/*/example`.
