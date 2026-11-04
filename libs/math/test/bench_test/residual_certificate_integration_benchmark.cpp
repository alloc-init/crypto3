// Copyright (c) 2026. Distributed under the MIT License.

#include <algorithm>
#include <bit>
#include <chrono>
#include <cstdint>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <sys/resource.h>

#include <nil/crypto3/marshalling/math/types/polynomial.hpp>
#include <nil/crypto3/math/polynomial/backends/bn254_small_prime_fft_backend.hpp>
#include <nil/crypto3/math/polynomial/factorization/kaltofen_shoup_distinct_degree_factorization.hpp>

namespace {
    namespace math = nil::crypto3::math;
    namespace pa = math::polynomial_arithmetic;
    namespace fs = std::filesystem;
    using backend = pa::bn254_small_prime_fft_backend;
    using polynomial = backend::polynomial_type;
    using value = polynomial::value_type;
    using field = value::field_type;
    using endianness = nil::marshalling::option::little_endian;
    using clock_type = std::chrono::steady_clock;
    using discovery_table = math::detail::kaltofen_shoup_frobenius_precomputation<backend>;
    constexpr std::size_t maximum_coefficients = 16385;
    static_assert(field::arity == 12 && sizeof(std::size_t) == 8);

    void require(bool condition, const char *message) {
        if (!condition)
            throw std::runtime_error(message);
    }

    // Use the same Crypto3 marshalling format and canonical-byte checks as the saved-factor harness.
    std::vector<std::uint8_t> encode(const polynomial &p) {
        auto encoded = nil::crypto3::marshalling::types::fill_polynomial<endianness>(p);
        std::vector<std::uint8_t> bytes(encoded.length());
        auto it = bytes.begin();
        require(encoded.write(it, bytes.size()) == nil::marshalling::status_type::success && it == bytes.end(),
                "polynomial encoding failed");
        return bytes;
    }

    void write(const fs::path &path, const polynomial &p) {
        const auto bytes = encode(p);
        std::ofstream stream;
        stream.exceptions(std::ios::failbit | std::ios::badbit);
        stream.open(path, std::ios::binary | std::ios::noreplace);
        stream.write(reinterpret_cast<const char *>(bytes.data()), bytes.size());
        stream.close();
    }

    polynomial read(const fs::path &path) {
        std::ifstream stream(path, std::ios::binary | std::ios::ate);
        stream.exceptions(std::ios::failbit | std::ios::badbit);
        const auto length = stream.tellg();
        require(length >= 8 && length <= std::streamoff(8 + maximum_coefficients * 384), "invalid polynomial length");
        std::vector<std::uint8_t> bytes(static_cast<std::size_t>(length));
        stream.seekg(0);
        stream.read(reinterpret_cast<char *>(bytes.data()), length);
        std::uint64_t count = 0;
        for (std::size_t i = 0; i < 8; ++i)
            count |= std::uint64_t(bytes[i]) << (8 * i);
        require(count >= 1 && count <= maximum_coefficients && length == std::streamoff(8 + count * 384),
                "invalid Fp12 coefficient count");
        using encoded_type =
            nil::crypto3::marshalling::types::polynomial<nil::marshalling::field_type<endianness>, polynomial>::type;
        encoded_type encoded;
        auto it = bytes.begin();
        require(encoded.read(it, bytes.size()) == nil::marshalling::status_type::success && it == bytes.end(),
                "polynomial decoding failed");
        auto p = nil::crypto3::marshalling::types::make_polynomial<endianness, polynomial>(encoded);
        require(p.size() == count && (count == 1 || !p.back().is_zero()) && encode(p) == bytes,
                "noncanonical polynomial");
        return p;
    }

    std::size_t integer(const std::string &text) {
        require(!text.empty() && text.find_first_not_of("0123456789") == std::string::npos, "invalid integer");
        return std::stoull(text);
    }

    struct options {
        fs::path state, output;
        std::string mode = "next";
        std::size_t width = 246, prepared_bytes = std::size_t(8) << 30;
        std::size_t minimum_excluded = 0, maximum_attempts = 1, certificate_bytes = 0;
    };

    options parse(int argc, char **argv) {
        options result;
        for (int i = 1; i < argc; ++i) {
            const std::string key = argv[i];
            require(i + 1 < argc, "missing option value");
            const std::string argument = argv[++i];
            if (key == "--state")
                result.state = argument;
            else if (key == "--output")
                result.output = argument;
            else if (key == "--mode")
                result.mode = argument;
            else if (key == "--width")
                result.width = integer(argument);
            else if (key == "--prepared-bytes")
                result.prepared_bytes = integer(argument);
            else if (key == "--minimum-excluded-degree")
                result.minimum_excluded = integer(argument);
            else if (key == "--maximum-attempts")
                result.maximum_attempts = integer(argument);
            else if (key == "--certificate-prepared-bytes")
                result.certificate_bytes = integer(argument);
            else
                throw std::invalid_argument("unknown option: " + key);
        }
        require(!result.state.empty() && !result.output.empty(), "--state and --output are required");
        require(result.mode == "next" || result.mode == "all" || result.mode == "certificate", "invalid mode");
        require(result.width >= 1 && result.width <= 512 && result.prepared_bytes <= (std::size_t(16) << 30) &&
                    result.certificate_bytes <= (std::size_t(16) << 30),
                "invalid arithmetic bounds");
        return result;
    }

    pa::polynomial_context_options arithmetic_options(const options &o) {
        pa::polynomial_context_options result;
        result.gcd_half_gcd_cutoff = 128;
        result.half_gcd_basecase_cutoff = 60;
        result.distinct_degree_screening_limit = 32;
        result.modular_composition_block_size = o.width;
        result.modular_composition_prepared_byte_limit = o.prepared_bytes;
        result.modular_composition_simd_lanes = 8;
        result.distinct_degree_certificate_minimum_excluded_degree = o.minimum_excluded;
        result.distinct_degree_certificate_maximum_attempts = o.maximum_attempts;
        result.distinct_degree_certificate_prepared_byte_limit = o.certificate_bytes;
        return result;
    }
    std::size_t peak_rss() {
        rusage usage {};
        require(getrusage(RUSAGE_SELF, &usage) == 0, "getrusage failed");
        std::size_t bytes = usage.ru_maxrss;
#ifndef __APPLE__
        bytes *= 1024;
#endif
        return bytes;
    }

    struct timer {
        clock_type::time_point wall = clock_type::now();
        std::clock_t cpu = std::clock();
        double seconds() const {
            return std::chrono::duration<double>(clock_type::now() - wall).count();
        }
        double cpu_seconds() const {
            return double(std::clock() - cpu) / CLOCKS_PER_SEC;
        }
    };

    struct state {
        polynomial modulus, residual, giant;
        std::vector<polynomial> babies, factors;
        std::size_t block = 0, excluded = 0;
    };

    template<typename Context>
    void verify_state_factors(const state &s, Context &context) {
        require(s.modulus.size() >= s.residual.size() && s.residual.size() > 1 && s.modulus.back() == value::one() &&
                    s.residual.back() == value::one(),
                "invalid monic modulus or residual");
        polynomial product {value::one()};
        for (const auto &factor : s.factors) {
            // Irreducibility is supplied by the pinned captured-factor provenance, not inferred from degree.
            require(factor.size() - 1 > s.excluded && factor.back() == value::one(),
                    "captured factor contradicts the excluded-degree bound");
            math::multiplication(product, product, factor, context);
        }
        require(product == s.residual, "captured factors do not reproduce the residual");
        const math::polynomial_divisor_context<backend> divisor(s.residual, s.modulus.size() - s.residual.size() + 1,
                                                                context);
        polynomial remainder;
        math::remainder(remainder, s.modulus, divisor, context);
        require(math::is_zero(remainder), "residual does not divide the original modulus");
    }

    state load(const fs::path &directory) {
        std::ifstream input(directory / "state.txt");
        std::map<std::string, std::size_t> meta;
        std::string key, text;
        while (input >> key) {
            require(bool(input >> text), "missing saved-state value");
            require(meta.emplace(key, integer(text)).second, "duplicate saved-state key");
        }
        require(input.eof() && meta.size() == 5 && meta.contains("block") && meta.contains("excluded") &&
                    meta.contains("factors") && meta.contains("modulus_degree") && meta.contains("residual_degree"),
                "invalid saved-state metadata");
        state s;
        s.block = meta.at("block");
        s.excluded = meta.at("excluded");
        require(s.block >= 1 && s.block <= 128 && s.excluded >= s.block && s.excluded < maximum_coefficients &&
                    s.excluded % s.block == 0 && meta.at("factors") >= 1 && meta.at("factors") <= 128,
                "invalid saved-state bounds");
        s.modulus = read(directory / "modulus.bin");
        s.residual = read(directory / "residual.bin");
        s.giant = read(directory / "giant.bin");
        require(s.modulus.degree() == meta.at("modulus_degree") && s.residual.degree() == meta.at("residual_degree") &&
                    s.giant.size() < s.modulus.size(),
                "saved-state degree mismatch");
        for (std::size_t i = 0; i <= s.block; ++i) {
            s.babies.push_back(read(directory / ("baby-" + std::to_string(i) + ".bin")));
            require(s.babies.back().size() < s.modulus.size(), "unreduced saved baby step");
        }
        require(s.babies.front() == polynomial {value {}, value::one()}, "invalid zeroth baby step");
        for (std::size_t i = 0; i < meta.at("factors"); ++i)
            s.factors.push_back(read(directory / ("factor-" + std::to_string(i) + ".bin")));
        return s;
    }

}    // namespace

int main(int argc, char **argv) {
    try {
        const auto o = parse(argc, argv);
        require(fs::create_directory(o.output), "output already exists or has no parent");
        const timer process;
        auto saved = load(o.state);
        const auto capacity = std::bit_ceil(2 * saved.modulus.degree());
        pa::polynomial_context<backend> context(backend(capacity), arithmetic_options(o));
        verify_state_factors(saved, context);
        const math::polynomial_divisor_context<backend> divisor(saved.modulus, saved.modulus.degree() - 1, context);
        const discovery_table table(std::move(saved.babies), divisor, context);
        const auto schedule = math::detail::kaltofen_shoup_plan_residual_certificate(saved.residual.degree(),
                                                                                     saved.block, saved.excluded);
        {
            std::ofstream plan(o.output / "plan.json");
            plan << "{\"block\":" << saved.block << ",\"excluded\":" << saved.excluded
                 << ",\"modulus_degree\":" << saved.modulus.degree()
                 << ",\"residual_degree\":" << saved.residual.degree() << ",\"jump_blocks\":" << schedule.jump_blocks
                 << ",\"construction_compositions\":" << schedule.construction_compositions
                 << ",\"remainder_compositions\":" << schedule.remainder_compositions
                 << ",\"jump_compositions\":" << schedule.jump_compositions << ",\"targets\":[";
            for (std::size_t i = 0; i < schedule.targets.size(); ++i)
                plan << (i ? "," : "") << schedule.targets[i];
            plan << "]}\n";
            plan.close();
            require(bool(plan), "cannot write schedule");
        }
        const auto setup_seconds = process.seconds();
        std::vector<math::distinct_degree_factor<polynomial>> groups;
        math::detail::kaltofen_shoup_residual_certificate_state attempts;
        bool certified = false;
        std::cout << "MEASUREMENT_START mode=" << o.mode << " minimum_excluded=" << o.minimum_excluded << '\n'
                  << std::flush;
        const timer measured;
        if (o.mode == "certificate") {
            certified = math::detail::kaltofen_shoup_certify_residual(saved.residual, saved.giant, saved.excluded,
                                                                      table, divisor, context);
        } else {
            const auto status = math::detail::kaltofen_shoup_continue_discovery<backend>(
                nullptr, saved.residual, saved.excluded, saved.giant, table, divisor, context,
                [&](const auto &group) {
                    groups.push_back(group);
                    return o.mode == "next" ? math::factorization_control::stop_factorization :
                                              math::factorization_control::continue_factorization;
                },
                attempts);
            require(status == (o.mode == "next" ? math::factorization_control::stop_factorization :
                                                  math::factorization_control::continue_factorization),
                    "unexpected continuation status");
        }
        const auto seconds = measured.seconds(), cpu_seconds = measured.cpu_seconds();
        const auto rss = peak_rss();
        // The captured factors are a validation oracle only. No factor degree or outcome reaches the production
        // trigger.
        std::map<std::size_t, polynomial> expected;
        for (const auto &factor : saved.factors) {
            auto [entry, inserted] = expected.try_emplace(factor.degree(), polynomial {value::one()});
            math::multiplication(entry->second, entry->second, factor, context);
        }
        if (o.mode == "certificate") {
            require(certified == (saved.factors.size() == 1), "certificate disagrees with captured factors");
        } else {
            require(groups.size() == (o.mode == "next" ? 1 : expected.size()), "unexpected degree group count");
            auto entry = expected.begin();
            for (std::size_t i = 0; i < groups.size(); ++i, ++entry) {
                require(groups[i].irreducible_factor_degree == entry->first && groups[i].polynomial == entry->second,
                        "continued discovery differs from captured factor bytes or degree order");
                write(o.output / ("group-" + std::to_string(i) + ".bin"), groups[i].polynomial);
            }
        }
        std::ofstream report(o.output / "measurement.json");
        report << std::setprecision(17) << "{\"mode\":\"" << o.mode << "\",\"seconds\":" << seconds
               << ",\"cpu_seconds\":" << cpu_seconds << ",\"setup_seconds\":" << setup_seconds
               << ",\"process_seconds\":" << process.seconds() << ",\"peak_rss_bytes\":" << rss
               << ",\"modulus_degree\":" << saved.modulus.degree() << ",\"residual_degree\":" << saved.residual.degree()
               << ",\"excluded_degree\":" << saved.excluded << ",\"width\":" << o.width
               << ",\"prepared_byte_limit\":" << o.prepared_bytes
               << ",\"certificate_prepared_byte_limit\":" << o.certificate_bytes
               << ",\"minimum_excluded_degree\":" << o.minimum_excluded
               << ",\"maximum_attempts\":" << o.maximum_attempts << ",\"attempts\":" << attempts.attempts
               << ",\"last_attempted_degree\":" << attempts.last_attempted_degree
               << ",\"certificate_result\":" << (o.mode == "certificate" ? (certified ? "true" : "false") : "null")
               << ",\"rng_draws\":0,\"groups\":[";
        for (std::size_t i = 0; i < groups.size(); ++i)
            report << (i ? "," : "") << "{\"factor_degree\":" << groups[i].irreducible_factor_degree
                   << ",\"group_degree\":" << groups[i].polynomial.degree() << "}";
        report << "]}\n";
        report.close();
        require(bool(report), "cannot write measurement");
        std::cout << "MEASUREMENT_SECONDS=" << seconds << " ATTEMPTS=" << attempts.attempts << " PEAK_RSS_BYTES=" << rss
                  << '\n';
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "ERROR: " << error.what() << '\n';
        return 1;
    }
}
