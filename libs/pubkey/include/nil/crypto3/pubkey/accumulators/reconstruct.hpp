//---------------------------------------------------------------------------//
// Copyright (c) 2021 Mikhail Komarov <nemo@nil.foundation>
// Copyright (c) 2021 Ilias Khairullin <ilias@nil.foundation>
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

#ifndef CRYPTO3_ACCUMULATORS_PUBKEY_SSS_RECONSTRUCT_HPP
#define CRYPTO3_ACCUMULATORS_PUBKEY_SSS_RECONSTRUCT_HPP

#include <set>
#include <utility>
#include <algorithm>
#include <iterator>

#include <boost/concept_check.hpp>

#include <boost/accumulators/framework/accumulator_base.hpp>
#include <boost/accumulators/framework/parameters/sample.hpp>

#include <nil/crypto3/pubkey/accumulators/parameters/threshold_value.hpp>
#include <nil/crypto3/pubkey/accumulators/parameters/iterator_last.hpp>

#include <nil/crypto3/pubkey/keys/share_sss.hpp>
#include <nil/crypto3/pubkey/keys/public_share_sss.hpp>

#include <nil/crypto3/pubkey/modes/isomorphic.hpp>

namespace nil {
    namespace crypto3 {
        namespace pubkey {
            namespace accumulators {
                namespace impl {
                    template<typename ProcessingMode, typename = void>
                    struct reconstruct_impl;

                    template<typename ProcessingMode>
                    struct reconstruct_impl<ProcessingMode> : boost::accumulators::accumulator_base {
                    protected:
                        typedef ProcessingMode processing_mode_type;
                        typedef typename processing_mode_type::scheme_type scheme_type;
                        typedef typename processing_mode_type::op_type op_type;
                        typedef typename processing_mode_type::accumulator_type accumulator_type;

                    public:
                        typedef typename processing_mode_type::result_type result_type;

                        template<typename Args>
                        reconstruct_impl(const Args &args) : seen_shares(0) {
                        }

                        inline result_type result(boost::accumulators::dont_care) const {
                            return processing_mode_type::process(acc);
                        }

                        template<typename Args>
                        inline void operator()(const Args &args) {
                            resolve_type(args[boost::accumulators::sample],
                                         args[crypto3::accumulators::iterator_last | nullptr]);
                        }

                    protected:
                        inline void resolve_type(const share_sss<scheme_type> &share, std::nullptr_t = nullptr) {
                            processing_mode_type::update(acc, share);
                            seen_shares++;
                        }

                        inline void resolve_type(const public_share_sss<scheme_type> &public_share, std::nullptr_t = nullptr) {
                            processing_mode_type::update(acc, public_share);
                            seen_shares++;
                        }

                        template<typename InputRange>
                        inline void resolve_type(const InputRange &range, std::nullptr_t) {
                            for (const auto &s : range) {
                                resolve_type(s);
                            }
                        }

                        template<typename InputIterator>
                        inline void resolve_type(InputIterator first, InputIterator last) {
                            for (auto it = first; it != last; it++) {
                                resolve_type(*it);
                            }
                        }

                        std::size_t seen_shares;
                        mutable accumulator_type acc;
                    };
                }    // namespace impl

                namespace tag {
                    template<typename ProcessingMode>
                    struct reconstruct : boost::accumulators::depends_on<> {
                        typedef ProcessingMode mode_type;

                        /// INTERNAL ONLY
                        ///

                        typedef boost::mpl::always<impl::reconstruct_impl<mode_type>> impl;
                    };
                }    // namespace tag

                namespace extract {
                    template<typename ProcessingMode, typename AccumulatorSet>
                    typename boost::mpl::apply<AccumulatorSet, tag::reconstruct<ProcessingMode>>::type::result_type
                        reconstruct(const AccumulatorSet &acc) {
                        return boost::accumulators::extract_result<tag::reconstruct<ProcessingMode>>(acc);
                    }
                }    // namespace extract
            }        // namespace accumulators
        }            // namespace pubkey
    }                // namespace crypto3
}    // namespace nil

#endif    // CRYPTO3_ACCUMULATORS_PUBKEY_SSS_RECONSTRUCT_HPP
