//---------------------------------------------------------------------------//
// Copyright (c) 2024 Mikhail Komarov <nemo@allocin.it>
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

#ifndef CRYPTO3_FUNCTIONAL_FH_MULTI_IPE_H
#define CRYPTO3_FUNCTIONAL_FH_MULTI_IPE_H

#include <nil/crypto3/algebra/dlog.hpp>
#include <nil/crypto3/algebra/random_element.hpp>
#include <nil/crypto3/algebra/algorithms/pair.hpp>
#include <nil/crypto3/algebra/matrix/matrix.hpp>
#include <nil/crypto3/algebra/matrix/math.hpp>

#include <nil/crypto3/pubkey/keys/public_key.hpp>
#include <nil/crypto3/pubkey/keys/private_key.hpp>

#include <nil/crypto3/functional/keys/functional_key.hpp>

/**
 * \file
 * \ingroup fullysec
 * \brief // FH-Multi-IPE represents a Function Hiding Multi-client Inner
 * Product Encryption scheme based on the paper by P. Datta, T. Okamoto, and
 * J. Tomida:
 * "Full-Hiding (Unbounded) Multi-Input Inner Product Functional Encryption
 * from the 𝒌-Linear Assumption".
 * It allows clients to encrypt vectors {x_1,...,x_m} and derive a secret key
 * based on an inner product vectors {y_1,...,y_m} so that a decryptor can
 * decrypt the sum of inner products <x_1,y_1> + ... + <x_m, y_m> without
 * revealing vectors x_i or y_i. The scheme is slightly modified from the
 * original one to achieve a better performance. The difference is in
 * storing the secret master key as matrices B_hat, B_hat_star, instead of matrices
 * of elliptic curve elements g_1^B_hat, g_2^B_hat_star. This replaces elliptic
 * curves operations with matrix multiplications.
 */

namespace nil {
    namespace crypto3 {
        namespace functional {
            /**
             * fh_multi_ipe contains the shared choice for parameters on which
             * the functionality of the scheme depend.
             * @param SecLevel The parameter defines the security assumption of the scheme,
             * so called K-Lin assumption, where K is the specified sec_level
             * @param Clients Number of clients participating in the scheme
             * @param CiphertextSize Length of the vectors that each client will encrypt
             * @param BoundX Bound on the inputs of the vectors that will be encrypted
             * @param BoundY Bound on the inputs of the inner product vectors for which the functional keys will be
             * generated
             * .
             */
            template<typename CurveType, std::size_t K, std::size_t Clients, std::size_t CiphertextSize,
                     std::size_t BoundX, std::size_t BoundY>
            struct fh_multi_ipe {
                constexpr static const std::size_t sec_level = K;
                constexpr static const std::size_t clients = Clients;

                constexpr static const std::size_t bound_x = BoundX;
                constexpr static const std::size_t bound_y = BoundY;

                typedef CurveType curve_type;

                typedef typename curve_type::base_field_type base_field_type;
                typedef typename base_field_type::value_type base_field_value_type;

                constexpr static const std::size_t schedule_size = CiphertextSize;
                typedef std::array<typename curve_type::template g1_type<>, 2 * schedule_size + 2 * sec_level + 1>
                    schedule_type;

                constexpr static const std::size_t plaintext_size = BoundX * base_field_value_type::modulus_bits;
                typedef std::array<base_field_value_type, BoundX> plaintext_type;

                constexpr static const std::size_t digest_bits = curve_type::gt_type::value_type::modulus_bits;
                typedef typename curve_type::gt_type::value_type digest_type;

                constexpr static const typename algebra::fields::arithmetic_params<base_field_type>::integral_type
                    group_order = algebra::fields::arithmetic_params<base_field_type>::group_order;

                static_assert(BoundX * BoundY * schedule_size * Clients < group_order);
            };

            template<typename CurveType, std::size_t K, std::size_t Clients, std::size_t CiphertextSize,
                     std::size_t BoundX, std::size_t BoundY>
            struct functional_key<fh_multi_ipe<CurveType, K, Clients, CiphertextSize, BoundX, BoundY>> {
                typedef fh_multi_ipe<CurveType, K, Clients, CiphertextSize, BoundX, BoundY> scheme_type;

                typedef CurveType curve_type;
                typedef typename curve_type::base_field_type field_type;

                constexpr static const std::size_t sec_level = scheme_type::sec_level;
                constexpr static const std::size_t clients = scheme_type::clients;
                constexpr static const std::size_t Q = BoundY;

                constexpr static const std::size_t schedule_size = scheme_type::schedule_size;
                typedef typename scheme_type::schedule_type schedule_type;

                constexpr static const std::size_t plaintext_size = scheme_type::plaintext_size;
                typedef typename scheme_type::plaintext_type plaintext_type;

                constexpr static const std::size_t digest_bits = scheme_type::digest_bits;
                typedef typename scheme_type::digest_type digest_type;

                functional_key() {
                    m.fill(typename curve_type::template g2_type<>::value_type::zero());
                }

                algebra::matrix<typename curve_type::template g2_type<>, Clients, 2 * schedule_size + 2 * sec_level + 1>
                    m;
            };
        }    // namespace functional

        namespace pubkey {
            template<typename CurveType, std::size_t K, std::size_t Clients, std::size_t CiphertextSize,
                     std::size_t BoundX, std::size_t BoundY>
            struct public_key<functional::fh_multi_ipe<CurveType, K, Clients, CiphertextSize, BoundX, BoundY>> {
                typedef functional::fh_multi_ipe<CurveType, K, Clients, CiphertextSize, BoundX, BoundY> scheme_type;

                constexpr static const std::size_t sec_level = K;
                constexpr static const std::size_t clients = Clients;
                constexpr static const std::size_t bound_y = scheme_type::bound_y;
                constexpr static const std::size_t bound_x = scheme_type::bound_x;

                typedef CurveType curve_type;
                typedef typename curve_type::base_field_type field_type;
                typedef typename field_type::value_type field_value_type;

                constexpr static const std::size_t schedule_size = scheme_type::schedule_size;
                typedef typename scheme_type::schedule_type ciphertext_type;

                constexpr static const std::size_t plaintext_size = scheme_type::plaintext_size;
                typedef typename scheme_type::plaintext_type plaintext_type;

                constexpr static const std::size_t digest_bits = scheme_type::digest_bits;
                typedef typename scheme_type::digest_type digest_type;

                typedef algebra::matrix<field_type, schedule_size + sec_level + 1,
                                        2 * schedule_size + 2 * sec_level + 1>
                    client_key_type;

                template<typename DistributionType =
                             boost::random::uniform_int_distribution<typename field_type::integral_type>,
                         typename UniformRandomBitGenerator = boost::random::random_device>
                public_key(const field_value_type &mu =
                               algebra::random_element<field_type, DistributionType>(UniformRandomBitGenerator())) :
                    pkey(algebra::final_exponentiation<CurveType>(
                             pair(algebra::fields::detail::fp2_extension_params<field_type>::arithmetic_generator,
                                  algebra::fields::arithmetic_params<field_type>::arithmetic_generator))
                             .pow(mu)) {
                    // TODO: a G2 generator is supposed to be ECP2_CURVE_NAME generator and possibly is not correct in
                    // here
                }

                /**
                 * Accepts the encrypted vectors and functional encryption key. It returns the
                 * inner product of x and y, i.e. <x_1,y_1> + ... + <x_m, y_m> where x_i is i-th
                 * encrypted vector and y_i is i-th inner product vector (i-th row of y).
                 * If decryption failed, an error is returned.
                 *
                 * @param res The result of the decryption (the value will be stored here)
                 * @param ciphers An array of the ciphertexts
                 * @param fe_key A pointer to the functional encryption key
                 * @return Error code
                 */
                digest_type decrypt(const std::array<ciphertext_type, Clients> &ciphers,
                                    const functional::functional_key<scheme_type> &fe_key) {
                    typename curve_type::gt_type::value_type sum = typename curve_type::gt_type::value_type::one(),
                                                             paired;

                    for (size_t i = 0; i < Clients; i++) {
                        for (size_t j = 0; j < 2 * schedule_size + 2 * sec_level + 1; j++) {
                            paired = algebra::pair<curve_type>(fe_key.m[i][j], ciphers[i][j]);
                            paired = algebra::final_exponentiation<CurveType>(paired);
                            sum *= paired;
                        }
                    }

                    return algebra::baby_giant_dlog<field_type>(sum, pkey, BoundY * BoundX * Clients * schedule_size);
                }

                typename curve_type::gt_type::value_type pkey;
            };

            /**
             * fh_multi_ipe_sec_key represents a master secret key in fh_multi_ipe scheme.
             */
            template<typename CurveType, std::size_t K, std::size_t Clients, std::size_t CiphertextSize,
                     std::size_t BoundX, std::size_t BoundY>
            struct private_key<functional::fh_multi_ipe<CurveType, K, Clients, CiphertextSize, BoundX, BoundY>>
                : public public_key<functional::fh_multi_ipe<CurveType, K, Clients, CiphertextSize, BoundX, BoundY>> {
                typedef typename public_key<functional::fh_multi_ipe<CurveType, K, Clients, CiphertextSize, BoundX,
                                                                     BoundY>>::scheme_type scheme_type;

                typedef CurveType curve_type;

                typedef typename curve_type::base_field_type field_type;
                typedef typename curve_type::base_field_type::value_type field_value_type;

                constexpr static const std::size_t sec_level = K;
                constexpr static const std::size_t clients = Clients;
                constexpr static const std::size_t Q = BoundY;

                constexpr static const std::size_t schedule_size = scheme_type::schedule_size;
                typedef typename scheme_type::schedule_type schedule_type;

                constexpr static const std::size_t plaintext_size = scheme_type::plaintext_size;
                typedef typename scheme_type::plaintext_type plaintext_type;

                constexpr static const auto group_order = scheme_type::group_order;

                typedef algebra::matrix<field_type, schedule_size + sec_level + 1,
                                        2 * schedule_size + 2 * sec_level + 1>
                    client_key_type;

                /**
                 * Generates a master secret key and a public key for the scheme. It returns an
                 * error if generating one of the parts of the secret key failed.
                 *
                 * @param sec_key A pointer to a fh_multi_ipe_sec_key struct (the
                 * master secret key will be stored here)
                 * @param pub_key A pointer to a FH12_BN256 struct (the public key will
                 * be stored here)
                 * @param c A pointer to an instance of the scheme (*initialized* fh_multi_ipe struct)
                 * @return Error code
                 */
                template<typename DistributionType =
                             boost::random::uniform_int_distribution<typename field_type::integral_type>,
                         typename UniformRandomBitGenerator = boost::random::random_device>
                private_key(const field_value_type &mu =
                                algebra::random_element<field_type, DistributionType>(UniformRandomBitGenerator())) :
                    public_key<functional::fh_multi_ipe<CurveType, sec_level, Clients, schedule_size, BoundX, BoundY>>(
                        mu) {
                    algebra::matrix<client_key_type, Clients, Clients> B, B_star;

                    for (size_t i = 0; i < Clients; i++) {
                        random_OB(B[i], B_star[i], mu, group_order);

                        for (size_t j = 0; j < schedule_size + sec_level + 1; j++) {
                            if (j < schedule_size) {
                                B_hat[i][j] = B[i][j];
                                B_star_hat[i][j] = B_star[i][j];
                            } else if (j == schedule_size) {
                                B_hat[i][j] = B[i][j + schedule_size];
                                B_star_hat[i][j] = B_star[i][j + schedule_size];
                            } else if (j < schedule_size + sec_level) {
                                B_hat[i][j] = B[i][j - 1 + schedule_size + sec_level];
                                B_star_hat[i][j] = B_star[i][j + schedule_size];
                            } else {
                                B_hat[i][j] = B[i][j - 1 + schedule_size + sec_level];
                            }
                        }
                    }
                }

                /**
                 * The function is called by a client that encrypts input vector x with
                 * the provided part master secret key. It returns a ciphertext struct.
                 * If encryption failed, an error is returned.
                 *
                 * @param cipher A pointer to an initialized vec_G1 struct
                 * (the resulting ciphertext will be stored here)
                 * @param x A pointer to the plaintext vector
                 * @param part_sec_key A pointer to a matrix representing a part
                 * of the master secret key (i-th client gets i-th matrix in array B_hat)
                 * @return Error code
                 */
                schedule_type encrypt(const plaintext_type &x, const client_key_type &part_sec_key) {
                    field_value_type s = field_value_type::zero();
                    std::array<field_value_type, sec_level> phi;
                    std::array<field_value_type, 2 * schedule_size + 2 * sec_level + 1> key_vec;
                    std::array<field_value_type, 2 * schedule_size + 2 * sec_level + 1> tmp_vec;

                    key_vec.fill(field_value_type::zero());
                    tmp_vec.fill(field_value_type::zero());

                    boost::random::mt19937 gen;
                    boost::random::uniform_int_distribution<> dist(0, group_order);

                    for (const auto &v : phi) {
                        v = dist(gen);
                    }

                    for (size_t j = 0; j < schedule_size + sec_level + 1; j++) {
                        if (j < schedule_size) {
                            s = x[j];
                        } else if (j == schedule_size) {
                            s = field_value_type::one();
                        } else {
                            s = phi[j - schedule_size - 1];
                        }

                        tmp_vec = part_sec_key.B_hat[j] * s;
                        key_vec += tmp_vec;
                    }

                    return key_vec * algebra::fields::arithmetic_params<field_type>::arithmetic_generator;
                }

                /**
                 * Takes a master secret key and input matrix y, and derives the functional
                 * encryption key. In case the key could not be derived, it returns an error.
                 *
                 * @param fe_key A pointer to a mat_G2 struct (the functional
                 * encryption key will be stored here)
                 * @param y A pointer to the inner product matrix
                 * @return Error code
                 */
                template<typename MatrixType,
                         typename = typename std::enable_if<std::is_same<
                             typename MatrixType::value_type, typename curve_type::gt_type::value_type>::value>::type>
                functional::functional_key<scheme_type> derive(const MatrixType &y) {
                    field_value_type s = field_value_type::zero();

                    MatrixType gamma(sec_level, Clients), key_mat(Clients, 2 * schedule_size + 2 * sec_level + 1);
                    key_mat.fill(field_value_type::zero());

                    boost::random::mt19937 gen;

                    for (const auto &v : gamma) {
                        v = algebra::random_element<field_type, boost::random::uniform_int_distribution<>>(gen) %
                            group_order;
                    }

                    gamma[0][Clients - 1] = field_value_type::zero();
                    for (size_t i = 0; i < Clients - 1; i++) {
                        gamma[0][Clients - 1] = gamma[0][Clients - 1] + gamma[0][i];
                    }
                    gamma[0][Clients - 1] = -gamma[0][Clients - 1];

                    std::vector<field_value_type> tmp_vec(2 * schedule_size + 2 * sec_level + 1);

                    for (size_t i = 0; i < Clients; i++) {
                        for (size_t j = 0; j < schedule_size + sec_level; j++) {
                            if (j < schedule_size) {
                                s = y[i][j];
                            } else {
                                s = gamma[j - schedule_size][i];
                            }

                            tmp_vec = B_star_hat * s;

                            key_mat[i] = key_mat[i] * tmp_vec;
                        }
                    }

                    return key_mat * algebra::fields::arithmetic_params<field_type>::arithmetic_generator;
                }

            protected:
                // random_OB is a helping function used in fh_multi_ipe_generate_keys.
                template<typename FieldType, template<typename T> typename MatrixType>
                void random_OB(MatrixType<typename FieldType::value_type> &B,
                               MatrixType<typename FieldType::value_type> &B_star,
                               const typename FieldType::value_type &mu,
                               const typename FieldType::value_type &p) {
                    std::fill(B.begin(), B.end(),
                              []() -> typename FieldType::value_type { return algebra::random_element<FieldType>(); });

                    B_star = algebra::transpose(algebra::inverse(B));
                    B_star = B_star * mu;
                }

            public:
                algebra::matrix<client_key_type, Clients, Clients> B_hat;
                algebra::matrix<client_key_type, Clients, Clients> B_star_hat;
            };
        }    // namespace pubkey
    }    // namespace crypto3
}    // namespace nil

#endif
