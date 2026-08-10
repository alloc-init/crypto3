//---------------------------------------------------------------------------//
// Copyright (c) 2026 Riccardo Abbate
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

#ifndef CRYPTO3_ALGEBRA_FIELDS_FIELD_ELEMENT_COORDINATE_TRAITS_HPP
#define CRYPTO3_ALGEBRA_FIELDS_FIELD_ELEMENT_COORDINATE_TRAITS_HPP

namespace nil {
    namespace crypto3 {
        namespace algebra {
            namespace fields {

                /**
                 * Describes how a field element is represented by base-field coordinates.
                 *
                 * The default leaves the representation unavailable. Field-element implementations may specialize
                 * this trait so algorithms can operate on their coordinates without depending on their storage layout.
                 */
                template<typename FieldElement>
                struct field_element_coordinate_traits {
                    constexpr static bool is_supported = false;
                };

            }    // namespace fields
        }    // namespace algebra
    }    // namespace crypto3
}    // namespace nil

#endif    // CRYPTO3_ALGEBRA_FIELDS_FIELD_ELEMENT_COORDINATE_TRAITS_HPP
