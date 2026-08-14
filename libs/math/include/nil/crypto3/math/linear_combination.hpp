//---------------------------------------------------------------------------//
// Copyright (c) 2020-2021 Mikhail Komarov <nemo@nil.foundation>
// Copyright (c) 2020-2021 Nikita Kaskov <nbering@nil.foundation>
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
// @file Declaration of interfaces for:
// - a variable (i.e., x_i),
// - a linear term (i.e., a_i * x_i), and
// - a linear combination (i.e., sum_i a_i * x_i).
//---------------------------------------------------------------------------//

#pragma once

#include <algorithm>
#include <cstddef>
#include <ostream>
#include <type_traits>
#include <vector>

namespace nil::crypto3::math {

    // Describes how a variable assignment represents the distinguished
    // constant-one variable x_0.
    enum class assignment_layout {
        // x_0 is omitted. Term index 0 evaluates to one, while variable i is
        // read from assignment[i - 1]. This is the standard Crypto3 R1CS layout.
        implicit_constant,

        // assignment[0] stores x_0 explicitly, normally as one. Every term
        // index directly addresses the assignment. AADP uses this layout so
        // witness vectors align with matrix columns.
        explicit_constant,
    };

    /**
     * Forward declaration.
     */
    template<typename VariableType, assignment_layout AssignmentLayout = assignment_layout::implicit_constant>
    class linear_combination;

    /****************************** Linear term **********************************/

    /**
     * A linear term represents a formal expression of the form "coeff * x_{index}".
     */
    template<typename VariableType>
    class linear_term {
        typedef typename VariableType::value_type field_value_type;

    public:
        typedef typename VariableType::field_type field_type;
        typedef VariableType variable_type;

        std::size_t index;
        field_value_type coeff;

        linear_term() { };
        linear_term(const VariableType &var) : index(var.index), coeff(field_value_type::one()) {
        }

        linear_term operator*(const field_value_type &field_coeff) const {
            linear_term result(this->index);
            result.coeff = field_coeff * this->coeff;
            return result;
        }

        // linear_combination<VariableType> operator+(
        //     const linear_combination<VariableType> &other) const {
        //     return linear_combination<VariableType>(*this) + other;
        // }

        // linear_combination<VariableType> operator-(
        //     const linear_combination<VariableType> &other) const {
        //     return (*this) + (-other);
        // }

        linear_term operator-() const {
            return linear_term(this->index) * (-this->coeff);
        }

        bool operator==(const linear_term &other) const {
            return (this->index == other.index && this->coeff == other.coeff);
        }
    };

    template<typename VariableType>
    linear_term<VariableType> operator*(const typename VariableType::value_type &field_coeff,
                                        const linear_term<VariableType> &lt) {
        return lt * field_coeff;
    }

    // template<typename VariableType>
    // linear_combination<VariableType> operator+(const typename VariableType::value_type
    // &field_coeff,
    //                                         const linear_term<VariableType> &lt) {
    //     return linear_combination<VariableType>(field_coeff) + lt;
    // }

    // template<typename VariableType>
    // linear_combination<VariableType> operator-(const typename VariableType::value_type
    // &field_coeff,
    //                                         const linear_term<VariableType> &lt) {
    //     return linear_combination<VariableType>(field_coeff) - lt;
    // }

    /***************************** Linear combination ****************************/

    /**
     * A linear combination represents a formal expression of the form "sum_i coeff_i * x_{index_i}".
     */
    template<typename VariableType, assignment_layout AssignmentLayout>
    class linear_combination {
        typedef typename VariableType::value_type field_value_type;

        constexpr static const bool RotationSupport = false;

    public:
        typedef typename VariableType::field_type field_type;
        typedef VariableType variable_type;
        std::vector<linear_term<VariableType>> terms;

        linear_combination() { };

        linear_combination(const field_value_type &field_coeff) {
            this->add_term(linear_term<VariableType>(0) * field_coeff);
        }

        template<typename T = VariableType,
                 std::enable_if_t<std::is_same_v<typename T::field_type, typename T::value_type>, int> = 0>
        linear_combination(std::size_t index) {
            this->add_term(VariableType(index));
        }

        linear_combination(const VariableType &var) {
            this->add_term(var);
        }

        linear_combination(const linear_term<VariableType> &lt) {
            this->add_term(lt);
        }

        linear_combination(const std::vector<linear_term<VariableType>> &all_terms) {
            if (all_terms.empty()) {
                return;
            }

            terms = all_terms;
            std::sort(terms.begin(), terms.end(),
                      [](linear_term<VariableType> a, linear_term<VariableType> b) { return a.index < b.index; });

            auto result_it = terms.begin();
            for (auto it = ++terms.begin(); it != terms.end(); ++it) {
                if (it->index == result_it->index) {
                    result_it->coeff += it->coeff;
                } else {
                    *(++result_it) = *it;
                }
            }
            terms.resize((result_it - terms.begin()) + 1);
        }

        /* for supporting range-based for loops over linear_combination */
        typename std::vector<linear_term<VariableType>>::const_iterator begin() const {
            return terms.begin();
        }

        typename std::vector<linear_term<VariableType>>::const_iterator end() const {
            return terms.end();
        }

        void add_term(const VariableType &var) {
            this->terms.emplace_back(linear_term<VariableType>(var));
        }
        void add_term(const VariableType &var, const field_value_type &field_coeff) {
            this->terms.emplace_back(linear_term<VariableType>(var) * field_coeff);
        }
        void add_term(const linear_term<VariableType> &lt) {
            this->terms.emplace_back(lt);
        }

        // AssignmentLayout determines how term indices map to assignment slots;
        // see assignment_layout above. Coefficients and assignment values may
        // belong to different fields if they define a scalar action such as
        // Fp * Fp12.
        template<typename AssignmentType>
        typename AssignmentType::value_type evaluate(const AssignmentType &assignment) const {
            using result_type = typename AssignmentType::value_type;
            result_type acc = result_type::zero();
            for (const auto &lt : terms) {
                if constexpr (AssignmentLayout == assignment_layout::implicit_constant) {
                    acc += lt.coeff * (lt.index == 0 ? result_type::one() : assignment[lt.index - 1]);
                } else {
                    acc += lt.coeff * assignment[lt.index];
                }
            }
            return acc;
        }

        linear_combination operator*(const field_value_type &field_coeff) const {
            linear_combination result;
            result.terms.reserve(this->terms.size());
            for (const linear_term<VariableType> &lt : this->terms) {
                result.terms.emplace_back(lt * field_coeff);
            }
            return result;
        }

        linear_combination operator+(const linear_combination &other) const {
            linear_combination result;

            auto it1 = this->terms.begin();
            auto it2 = other.terms.begin();

            /* invariant: it1 and it2 always point to unprocessed items in the corresponding linear
             * combinations
             */
            while (it1 != this->terms.end() && it2 != other.terms.end()) {
                if (it1->index < it2->index) {
                    result.terms.emplace_back(*it1);
                    ++it1;
                } else if (it1->index > it2->index) {
                    result.terms.emplace_back(*it2);
                    ++it2;
                } else {
                    /* it1->index == it2->index */
                    result.terms.emplace_back(linear_term<VariableType>(VariableType(it1->index)) *
                                              (it1->coeff + it2->coeff));
                    ++it1;
                    ++it2;
                }
            }

            if (it1 != this->terms.end()) {
                result.terms.insert(result.terms.end(), it1, this->terms.end());
            } else {
                result.terms.insert(result.terms.end(), it2, other.terms.end());
            }

            return result;
        }

        linear_combination operator-(const linear_combination &other) const {
            return (*this) + (-other);
        }

        linear_combination operator-() const {
            return (*this) * (-field_value_type::one());
        }

        bool operator==(const linear_combination &other) const {

            std::vector<linear_term<VariableType>> thisterms = this->terms;
            std::sort(thisterms.begin(), thisterms.end(),
                      [](linear_term<VariableType> a, linear_term<VariableType> b) { return a.index < b.index; });

            std::vector<linear_term<VariableType>> otherterms = other.terms;
            std::sort(otherterms.begin(), otherterms.end(),
                      [](linear_term<VariableType> a, linear_term<VariableType> b) { return a.index < b.index; });

            return (thisterms == otherterms);
        }

        bool is_valid(size_t num_variables) const {
            if (terms.empty()) {
                return true;
            }

            /* check that all terms in linear combination are sorted */
            for (std::size_t i = 1; i < terms.size(); ++i) {
                if (terms[i - 1].index >= terms[i].index) {
                    return false;
                }
            }

            /* check that the variables are in proper range. as the variables
               are sorted, it suffices to check the last term */
            if ((--terms.end())->index >= num_variables) {
                return false;
            }

            return true;
        }
    };

    template<typename VariableType, assignment_layout AssignmentLayout>
    linear_combination<VariableType, AssignmentLayout>
        operator*(const typename VariableType::value_type &field_coeff,
                  const linear_combination<VariableType, AssignmentLayout> &lc) {
        return lc * field_coeff;
    }

    template<typename VariableType, assignment_layout AssignmentLayout>
    linear_combination<VariableType, AssignmentLayout>
        operator+(const typename VariableType::value_type &field_coeff,
                  const linear_combination<VariableType, AssignmentLayout> &lc) {
        return linear_combination<VariableType, AssignmentLayout>(field_coeff) + lc;
    }

    template<typename VariableType, assignment_layout AssignmentLayout>
    linear_combination<VariableType, AssignmentLayout>
        operator-(const typename VariableType::value_type &field_coeff,
                  const linear_combination<VariableType, AssignmentLayout> &lc) {
        return linear_combination<VariableType, AssignmentLayout>(field_coeff) - lc;
    }

    template<typename VariableType, assignment_layout AssignmentLayout>
    void add_scaled(linear_combination<VariableType, AssignmentLayout> &out,
                    const typename VariableType::value_type &coefficient,
                    const linear_combination<VariableType, AssignmentLayout> &in) {
        for (const auto &term : in.terms) {
            out.add_term(VariableType(term.index), coefficient * term.coeff);
        }
    }

    template<typename VariableType, assignment_layout AssignmentLayout>
    std::ostream &operator<<(std::ostream &out, const linear_combination<VariableType, AssignmentLayout> &combination) {
        for (std::size_t i = 0; i < combination.terms.size(); ++i) {
            const auto &term = combination.terms[i];
            out << term.coeff << " * v" << term.index;
            if (i + 1 < combination.terms.size()) {
                out << " + ";
            }
        }
        return out;
    }

}    // namespace nil::crypto3::math
