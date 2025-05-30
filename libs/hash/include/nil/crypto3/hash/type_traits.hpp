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

#ifndef CRYPTO3_HASH_TYPE_TRAITS_HPP
#define CRYPTO3_HASH_TYPE_TRAITS_HPP


#include <boost/config.hpp>


#define STRINGIZE(x) #x

#ifdef __has_include
#if __has_include(<version>)

#include <version>

#ifdef __cpp_lib_is_constant_evaluated
#include <type_traits>
#define CRYPTO3_HAS_IS_CONSTANT_EVALUATED
#endif
#endif
#endif

#ifdef __has_builtin
#if __has_builtin(__builtin_is_constant_evaluated) && !defined(BOOST_NO_CXX14_CONSTEXPR) && \
    !defined(BOOST_NO_CXX11_UNIFIED_INITIALIZATION_SYNTAX)
#define CRYPTO3_HAS_BUILTIN_IS_CONSTANT_EVALUATED
#endif
#endif
//
// MSVC also supports __builtin_is_constant_evaluated if it's recent enough:
//
#if defined(_MSC_FULL_VER) && (_MSC_FULL_VER >= 192528326)
#define CRYPTO3_HAS_BUILTIN_IS_CONSTANT_EVALUATED
#endif
//
// As does GCC-9:
//
#if defined(BOOST_GCC) && !defined(BOOST_NO_CXX14_CONSTEXPR) && (__GNUC__ >= 9) && \
    !defined(CRYPTO3_HAS_BUILTIN_IS_CONSTANT_EVALUATED)
#define CRYPTO3_HAS_BUILTIN_IS_CONSTANT_EVALUATED
#endif

#if defined(CRYPTO3_HAS_IS_CONSTANT_EVALUATED) && !defined(BOOST_NO_CXX14_CONSTEXPR)
#define CRYPTO3_IS_CONST_EVALUATED(x) std::is_constant_evaluated()
#elif defined(CRYPTO3_HAS_BUILTIN_IS_CONSTANT_EVALUATED)
#define CRYPTO3_IS_CONST_EVALUATED(x) __builtin_is_constant_evaluated()
#elif !defined(BOOST_NO_CXX14_CONSTEXPR) && defined(BOOST_GCC) && (__GNUC__ >= 6)
#define CRYPTO3_IS_CONST_EVALUATED(x) __builtin_constant_p(x)
#else
#define CRYPTO3_NO_CONSTEXPR_DETECTION
#endif

#define CRYPTO3_CXX14_CONSTEXPR BOOST_CXX14_CONSTEXPR
//
// Early compiler versions trip over the constexpr code:
//
#if defined(__clang__) && (__clang_major__ < 5)
#undef CRYPTO3_CXX14_CONSTEXPR
#define CRYPTO3_CXX14_CONSTEXPR
#endif
#if defined(__apple_build_version__) && (__clang_major__ < 9)
#undef CRYPTO3_CXX14_CONSTEXPR
#define CRYPTO3_CXX14_CONSTEXPR
#endif
#if defined(BOOST_GCC) && (__GNUC__ < 6)
#undef CRYPTO3_CXX14_CONSTEXPR
#define CRYPTO3_CXX14_CONSTEXPR
#endif
#if defined(BOOST_INTEL)
#undef CRYPTO3_CXX14_CONSTEXPR
#define CRYPTO3_CXX14_CONSTEXPR
#define CRYPTO3_NO_CONSTEXPR_DETECTION
#endif

#ifdef CRYPTO3_NO_CONSTEXPR_DETECTION
#define BOOST_CXX14_CONSTEXPR_IF_DETECTION
#else
#define BOOST_CXX14_CONSTEXPR_IF_DETECTION constexpr
#endif

#define GENERATE_HAS_MEMBER_TYPE(Type)                                                                                 \
    template<class T, typename Enable = void>                                                                          \
    class HasMemberType_##Type {                                                                                       \
    public:                                                                                                            \
        static constexpr bool RESULT = false;                                                                          \
    };                                                                                                                 \
                                                                                                                       \
    template<class T>                                                                                                  \
    class HasMemberType_##Type<T, typename std::enable_if<std::is_class<T>::value || std::is_union<T>::value>::type> { \
    private:                                                                                                           \
        using Yes = char[2];                                                                                           \
        using No = char[1];                                                                                            \
                                                                                                                       \
        struct Fallback {                                                                                              \
            struct Type { };                                                                                           \
        };                                                                                                             \
        struct Derived : T, Fallback { };                                                                              \
                                                                                                                       \
        template<class U>                                                                                              \
        static No &test(typename U::Type *);                                                                           \
        template<typename U>                                                                                           \
        static Yes &test(U *);                                                                                         \
                                                                                                                       \
    public:                                                                                                            \
        static constexpr bool RESULT = sizeof(test<Derived>(nullptr)) == sizeof(Yes);                                  \
    };                                                                                                                 \
                                                                                                                       \
    template<class T>                                                                                                  \
    struct has_##Type : public std::integral_constant<bool, HasMemberType_##Type<T>::RESULT> { };

#define GENERATE_HAS_MEMBER(member)                                                                                  \
    template<class T, typename Enable = void>                                                                        \
    class HasMember_##member {                                                                                       \
    public:                                                                                                          \
        static constexpr bool RESULT = false;                                                                        \
    };                                                                                                               \
                                                                                                                     \
    template<class T>                                                                                                \
    class HasMember_##member<T, typename std::enable_if<std::is_class<T>::value || std::is_union<T>::value>::type> { \
    private:                                                                                                         \
        using Yes = char[2];                                                                                         \
        using No = char[1];                                                                                          \
                                                                                                                     \
        struct Fallback {                                                                                            \
            int member;                                                                                              \
        };                                                                                                           \
        struct Derived : T, Fallback { };                                                                            \
                                                                                                                     \
        template<class U>                                                                                            \
        static No &test(decltype(U::member) *);                                                                      \
        template<typename U>                                                                                         \
        static Yes &test(U *);                                                                                       \
                                                                                                                     \
    public:                                                                                                          \
        static constexpr bool RESULT = sizeof(test<Derived>(nullptr)) == sizeof(Yes);                                \
    };                                                                                                               \
                                                                                                                     \
    template<class T>                                                                                                \
    struct has_##member : public std::integral_constant<bool, HasMember_##member<T>::RESULT> { };

#define GENERATE_HAS_MEMBER_FUNCTION(Function, ...)                                  \
                                                                                     \
    template<typename T>                                                             \
    struct has_##Function {                                                          \
        struct Fallback {                                                            \
            void Function(##__VA_ARGS__);                                            \
        };                                                                           \
                                                                                     \
        struct Derived : Fallback { };                                               \
                                                                                     \
        template<typename C, C>                                                      \
        struct ChT;                                                                  \
                                                                                     \
        template<typename C>                                                         \
        static char (&f(ChT<void (Fallback::*)(##__VA_ARGS__), &C::Function> *))[1]; \
                                                                                     \
        template<typename C>                                                         \
        static char (&f(...))[2];                                                    \
                                                                                     \
        static bool const value = sizeof(f<Derived>(0)) == 2;                        \
    };

#define GENERATE_HAS_MEMBER_CONST_FUNCTION(Function, ...)                                  \
                                                                                           \
    template<typename T>                                                                   \
    struct has_##Function {                                                                \
        struct Fallback {                                                                  \
            void Function(##__VA_ARGS__) const;                                            \
        };                                                                                 \
                                                                                           \
        struct Derived : Fallback { };                                                     \
                                                                                           \
        template<typename C, C>                                                            \
        struct ChT;                                                                        \
                                                                                           \
        template<typename C>                                                               \
        static char (&f(ChT<void (Fallback::*)(##__VA_ARGS__) const, &C::Function> *))[1]; \
                                                                                           \
        template<typename C>                                                               \
        static char (&f(...))[2];                                                          \
                                                                                           \
        static bool const value = sizeof(f<Derived>(0)) == 2;                              \
    };

#define GENERATE_HAS_MEMBER_RETURN_FUNCTION(Function, ReturnType, ...)                       \
                                                                                             \
    template<typename T>                                                                     \
    struct has_##Function {                                                                  \
        struct Dummy {                                                                       \
            typedef void ReturnType;                                                         \
        };                                                                                   \
        typedef typename std::conditional<has_##ReturnType<T>::value, T, Dummy>::type TType; \
        typedef typename TType::ReturnType type;                                             \
                                                                                             \
        struct Fallback {                                                                    \
            type Function(##__VA_ARGS__);                                                    \
        };                                                                                   \
                                                                                             \
        struct Derived : TType, Fallback { };                                                \
                                                                                             \
        template<typename C, C>                                                              \
        struct ChT;                                                                          \
                                                                                             \
        template<typename C>                                                                 \
        static char (&f(ChT<type (Fallback::*)(##__VA_ARGS__), &C::Function> *))[1];         \
                                                                                             \
        template<typename C>                                                                 \
        static char (&f(...))[2];                                                            \
                                                                                             \
        static bool const value = sizeof(f<Derived>(0)) == 2;                                \
    };

#define GENERATE_HAS_MEMBER_CONST_RETURN_FUNCTION(Function, ReturnType, ...)                 \
                                                                                             \
    template<typename T>                                                                     \
    struct has_##Function {                                                                  \
        struct Dummy {                                                                       \
            typedef void ReturnType;                                                         \
        };                                                                                   \
        typedef typename std::conditional<has_##ReturnType<T>::value, T, Dummy>::type TType; \
        typedef typename TType::ReturnType type;                                             \
                                                                                             \
        struct Fallback {                                                                    \
            type Function(##__VA_ARGS__) const;                                              \
        };                                                                                   \
                                                                                             \
        struct Derived : TType, Fallback { };                                                \
                                                                                             \
        template<typename C, C>                                                              \
        struct ChT;                                                                          \
                                                                                             \
        template<typename C>                                                                 \
        static char (&f(ChT<type (Fallback::*)(##__VA_ARGS__) const, &C::Function> *))[1];   \
                                                                                             \
        template<typename C>                                                                 \
        static char (&f(...))[2];                                                            \
                                                                                             \
        static bool const value = sizeof(f<Derived>(0)) == 2;                                \
    };

#include <type_traits>

namespace nil {
    namespace crypto3 {
        namespace hashes {
            GENERATE_HAS_MEMBER(word_bits)

            GENERATE_HAS_MEMBER(block_bits)

            GENERATE_HAS_MEMBER(state_bits)

            GENERATE_HAS_MEMBER_TYPE(digest_type)

            GENERATE_HAS_MEMBER_TYPE(word_type)

            GENERATE_HAS_MEMBER_TYPE(block_type)

            GENERATE_HAS_MEMBER_TYPE(state_type)

            GENERATE_HAS_MEMBER_RETURN_FUNCTION(absorb, digest_type);

            GENERATE_HAS_MEMBER_RETURN_FUNCTION(squeeze, block_type)

            template<typename PolicyType>
            class poseidon;

            template<typename ParamsType, typename HashType, typename GroupType>
            struct find_group_hash;

            template<typename ParamsType, typename BasePointGeneratorHash, typename GroupType>
            struct pedersen_to_point;

            template<typename ParamsType, typename BasePointGeneratorHash, typename GroupType>
            struct pedersen;

            template<typename FieldType, typename HashType, typename ParamsType>
            struct h2f;

            template<typename GroupType, typename HashType, typename ParamsType>
            struct h2c;

            template<typename HashType>
            struct is_find_group_hash : std::integral_constant<bool, false> {
            };

            template<typename ParamsType, typename HashType, typename GroupType>
            struct is_find_group_hash<find_group_hash<ParamsType, HashType, GroupType>>
                    : std::integral_constant<bool, true> {
            };

            template<typename HashType>
            struct is_pedersen : std::integral_constant<bool, false> {
            };

            template<typename ParamsType, typename BasePointGeneratorHash, typename GroupType>
            struct is_pedersen<pedersen_to_point<ParamsType, BasePointGeneratorHash, GroupType>>
                    : std::integral_constant<bool, true> {
            };

            template<typename ParamsType, typename BasePointGeneratorHash, typename GroupType>
            struct is_pedersen<pedersen<ParamsType, BasePointGeneratorHash, GroupType>>
                    : std::integral_constant<bool, true> {
            };

            template<typename HashType>
            struct is_h2f : std::integral_constant<bool, false> {
            };

            template<typename FieldType, typename HashType, typename ParamsType>
            struct is_h2f<h2f<FieldType, HashType, ParamsType>> : std::integral_constant<bool, true> {
            };

            template<typename HashType>
            struct is_h2c : std::integral_constant<bool, false> {
            };

            template<typename GroupType, typename HashType, typename ParamsType>
            struct is_h2c<h2c<GroupType, HashType, ParamsType>> : std::integral_constant<bool, true> {
            };

            // TODO: change this to more generic type trait to check for all sponge based hashes.
            template<typename HashType, typename Enable = void>
            struct is_poseidon {
            public:
                constexpr static const bool value = false;
            };

            template<typename HashType>
            struct is_poseidon<HashType, typename std::enable_if_t<std::is_same<nil::crypto3::hashes::poseidon<typename
                        HashType::policy_type>, HashType>::value>> {
            public:
                constexpr static const bool value = true;
                typedef HashType type;
            };

            template<template<typename...> class PrimaryTemplate, typename T>
            struct is_specialization_of : std::false_type {
            };

            template<template<typename...> class PrimaryTemplate, typename... Args>
            struct is_specialization_of<PrimaryTemplate, PrimaryTemplate<Args...>> : std::true_type {
            };

            template<typename HashType>
            struct is_sponge_construction {
                constexpr static const bool value =
                        has_digest_type<HashType>::value && has_block_type<HashType>::value &&
                        has_state_type<HashType>::value && has_word_type<HashType>::value &&
                        has_word_bits<HashType>::value && has_block_bits<HashType>::value &&
                        has_state_bits<HashType>::value && has_absorb<HashType>::value && has_squeeze<HashType>::value;
                typedef HashType type;
            };
        } // namespace hashes
    } // namespace crypto3
} // namespace nil

#endif    // CRYPTO3_HASH_TYPE_TRAITS_HPP
