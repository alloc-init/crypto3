# Contributing

Bug reports, fixes, documentation improvements, and proposals are welcome.

## Branches

The default integration branch is `master`. Create a descriptively named
feature or fix branch from the current `master` branch and open a pull request
targeting `master`.

## Issues And Discussions

Open issues in the [Crypto3 issue tracker](https://github.com/alloc-init/crypto3/issues).
Include the affected path, a minimal reproducer when possible, and details about
your compiler, Boost version, and operating system.

Use [GitHub Discussions](https://github.com/alloc-init/crypto3/discussions) for
design proposals, questions, and tutorial requests.

## Pull Requests

Components in this checkout are maintained directly in the monorepo. Submit
changes to [alloc-init/crypto3](https://github.com/alloc-init/crypto3), targeting
`master`, and mention the affected path, such as `libs/hash`, `libs/pubkey`, or
`libs/marshalling/core`, in the pull request description.

Before opening a pull request, configure, build, and run the relevant tests:

```sh
cmake -S . -B build -DBUILD_TESTS=ON
cmake --build build --target tests --parallel
ctest --test-dir build --output-on-failure
```
