//---------------------------------------------------------------------------//
// Copyright (c) 2018-2021 Mikhail Komarov <nemo@nil.foundation>
// Copyright (c) 2020-2021 Nikita Kaskov <nbering@nil.foundation>
// Copyright (c) 2026 Alloc Init
//
// MIT License
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//---------------------------------------------------------------------------//

#ifndef CRYPTO3_ZK_SAP_QUOTIENT_HPP
#define CRYPTO3_ZK_SAP_QUOTIENT_HPP

#include <cassert>
#include <cstddef>
#include <vector>

#include <nil/crypto3/algebra/fields/params.hpp>
#include <nil/crypto3/math/coset.hpp>
#include <nil/crypto3/math/domains/evaluation_domain.hpp>

namespace nil {
    namespace crypto3 {
        namespace zk {
            namespace snark {
                namespace reductions {
                    namespace detail {

                        /**
                         * Compute H = (A^2 - C) / Z using the domain's coset operations.
                         * The distinct buffers a and c contain domain.size() coefficients in ascending powers.
                         * On return, a contains H's coefficients, retaining trailing zeros; c is scratch storage.
                         * The caller must ensure exact divisibility and that the multiplicative-generator coset
                         * is disjoint from the domain. This routine does not check divisibility or add blinding.
                         */
                        template<typename FieldType>
                        void compute_sap_quotient(math::evaluation_domain<FieldType> &domain,
                                                  std::vector<typename FieldType::value_type> &a,
                                                  std::vector<typename FieldType::value_type> &c) {
                            assert(a.size() == domain.size());
                            assert(c.size() == domain.size());
                            assert(&a != &c);

                            const typename FieldType::value_type coset(
                                algebra::fields::arithmetic_params<FieldType>::multiplicative_generator);
                            math::multiply_by_coset(a, coset);
                            domain.fft(a);
                            math::multiply_by_coset(c, coset);
                            domain.fft(c);

#ifdef MULTICORE
#pragma omp parallel for
#endif
                            for (std::size_t i = 0; i < domain.size(); ++i) {
                                a[i] = a[i] * a[i] - c[i];
                            }

                            domain.divide_by_z_on_coset(a);
                            domain.inverse_fft(a);
                            math::multiply_by_coset(a, coset.inversed());
                        }

                    }    // namespace detail
                }    // namespace reductions
            }    // namespace snark
        }    // namespace zk
    }    // namespace crypto3
}    // namespace nil

#endif    // CRYPTO3_ZK_SAP_QUOTIENT_HPP
