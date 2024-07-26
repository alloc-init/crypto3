# Benchmark framework for Crypto3 suite

The whole thing revolves around `CRYPO3_RUN_BENCHMARK` macro.

It should be invoked like this:

```
int main()
{
    CRYPTO3_RUN_BENCHMARK(
        "Benchmark name",
        "sample-results.json",
        bench_type_A,
        bench_type_B,
        bench_type_C,
        A = B * C
    );
}
```

Parameters:

* Human-readable description.
* File name for saving results.
* Three type names for parameters.
* Expression to benchmark.

`bench_type_A` (and others) should have following interface:

```
class bench_type_X {
public:
    typedef a_type_name_for_A type_name;
    type_name sample_random_element();
    static std::string name();
    std::string to_string(type_name const& A);
};
```

A handy template with specializations is available:

```
bench_type<nil::crypto3::curves::bls12_381::scalar_field>
```

You can specialize your own, just follow the interface.

