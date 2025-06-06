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

#ifndef CRYPTO3_ACCUMULATORS_PUBKEY_SIGN_HPP
#define CRYPTO3_ACCUMULATORS_PUBKEY_SIGN_HPP

#include <iterator>
#include <type_traits>

#include <boost/parameter/value_type.hpp>

#include <boost/accumulators/framework/accumulator_base.hpp>
#include <boost/accumulators/framework/extractor.hpp>
#include <boost/accumulators/framework/depends_on.hpp>
#include <boost/accumulators/framework/parameters/sample.hpp>

#include <nil/crypto3/pubkey/accumulators/parameters/iterator_last.hpp>

#include <nil/crypto3/pubkey/modes/isomorphic.hpp>

namespace nil {
    namespace crypto3 {
        namespace pubkey {
            namespace accumulators {
                namespace impl {
                    template<typename ProcessingMode, typename = void>
                    struct sign_impl;

                    template<typename ProcessingMode>
                    struct sign_impl<ProcessingMode> : boost::accumulators::accumulator_base {
                    protected:
                        typedef ProcessingMode processing_mode_type;
                        typedef typename processing_mode_type::key_type key_type;
                        typedef typename processing_mode_type::accumulator_type accumulator_type;

                    public:
                        typedef typename processing_mode_type::result_type result_type;

                        template<typename Args>
                        sign_impl(const Args &args) : key(args[boost::accumulators::sample]) {
                            processing_mode_type::init_accumulator(key, acc);
                        }

                        template<typename Args>
                        inline void operator()(const Args &args) {
                            resolve_type(args[boost::accumulators::sample | nullptr],
                                         args[crypto3::accumulators::iterator_last | nullptr]);
                        }

                        inline result_type result(boost::accumulators::dont_care) const {
                            return processing_mode_type::process(key, acc);
                        }

                    protected:
                        //
                        // pop prove
                        //
                        inline void resolve_type(std::nullptr_t, std::nullptr_t) {
                        }

                        template<typename InputRange>
                        inline void resolve_type(const InputRange &range, std::nullptr_t) {
                            processing_mode_type::update(key, acc, range);
                        }

                        template<typename InputIterator>
                        inline void resolve_type(InputIterator first, InputIterator last) {
                            processing_mode_type::update(key, acc, first, last);
                        }

                        key_type key;
                        mutable accumulator_type acc;
                    };
                } // namespace impl

                namespace tag {
                    template<typename ProcessingMode>
                    struct sign : boost::accumulators::depends_on<> {
                        typedef ProcessingMode processing_mode_type;

                        /// INTERNAL ONLY
                        ///

                        typedef boost::mpl::always<impl::sign_impl<processing_mode_type>> impl;
                    };
                } // namespace tag

                namespace extract {
                    template<typename ProcessingMode, typename AccumulatorSet>
                    typename boost::mpl::apply<AccumulatorSet, tag::sign<ProcessingMode>>::type::result_type
                    sign(const AccumulatorSet &acc) {
                        return boost::accumulators::extract_result<tag::sign<ProcessingMode>>(acc);
                    }
                } // namespace extract
            } // namespace accumulators
        } // namespace pubkey
    } // namespace crypto3
} // namespace nil

#endif    // CRYPTO3_ACCUMULATORS_PUBKEY_SIGN_HPP
