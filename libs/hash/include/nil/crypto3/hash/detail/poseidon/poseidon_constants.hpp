//---------------------------------------------------------------------------//
// Copyright (c) 2020 Ilias Khairullin <ilias@nil.foundation>
// Copyright (c) 2020 Mikhail Komarov <nemo@nil.foundation>
//
// Distributed under the Boost Software License, Version 1.0
// See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt
//---------------------------------------------------------------------------//

#ifndef CRYPTO3_HASH_POSEIDON_CONSTANTS_HPP
#define CRYPTO3_HASH_POSEIDON_CONSTANTS_HPP

#include <nil/crypto3/algebra/matrix/matrix.hpp>
#include <nil/crypto3/algebra/matrix/math.hpp>
#include <nil/crypto3/algebra/matrix/operators.hpp>
#include <nil/crypto3/algebra/vector/vector.hpp>
#include <nil/crypto3/algebra/vector/math.hpp>
#include <nil/crypto3/algebra/vector/operators.hpp>

#include <nil/crypto3/algebra/random_element.hpp>

#include <nil/crypto3/algebra/fields/bls12/scalar_field.hpp>
#include <nil/crypto3/algebra/fields/bls12/scalar_field.hpp>
#include <nil/crypto3/algebra/fields/alt_bn128/scalar_field.hpp>
#include <nil/crypto3/algebra/fields/pallas/base_field.hpp>
#include <nil/crypto3/algebra/fields/vesta/base_field.hpp>
#include <nil/crypto3/algebra/fields/alt_bn128/scalar_field.hpp>

#include <nil/crypto3/hash/detail/poseidon/poseidon_policy.hpp>

#include <boost/assert.hpp>

#include <type_traits>

namespace nil {
    namespace crypto3 {
        namespace hashes {
            namespace detail {
                using namespace boost::multiprecision;

                // Uses Grain-LFSR stream cipher for constants generation.
                template<typename PolicyType>
                class poseidon_constants {
                public:
                    typedef PolicyType policy_type;
                    typedef typename policy_type::field_type field_type;

                    constexpr static const std::size_t state_words = policy_type::state_words;
                    constexpr static const std::size_t word_bits = policy_type::word_bits;
                    constexpr static const std::size_t full_rounds = policy_type::full_rounds;
                    constexpr static const std::size_t half_full_rounds = policy_type::half_full_rounds;
                    constexpr static const std::size_t part_rounds = policy_type::part_rounds;
                    constexpr static const std::size_t rate = policy_type::block_words;

                    typedef typename field_type::value_type element_type;
                    typedef typename field_type::integral_type integral_type;
                    constexpr static const integral_type modulus = field_type::modulus;

                    constexpr static const std::size_t lfsr_state_bits = 80;
                    typedef number<backends::cpp_int_modular_backend<lfsr_state_bits>> lfsr_state_type;

                    typedef algebra::vector<element_type, state_words> state_vector_type;
                    typedef algebra::matrix<element_type, state_words, state_words> mds_matrix_type;

                    typedef algebra::matrix<element_type, full_rounds + part_rounds, state_words> round_constants_type;
                    typedef std::pair<mds_matrix_type, round_constants_type> constants_type;

                    inline constexpr static const element_type &round_constant(std::size_t round, std::size_t i) {
                        return round_constants[round][i];
                    }

                    inline constexpr static void product_with_mds_matrix(state_vector_type &A_vector) {
                        A_vector = algebra::vectmatmul(A_vector, mds_matrix);
                    }

                private:
                    constexpr static inline mds_matrix_type generate_mds_matrix() {
                        mds_matrix_type mtx, transposed_mtx;

                        state_vector_type x;
                        state_vector_type y;
                        bool secure_mds_found = false;
                        while (!secure_mds_found) {
                            secure_mds_found = true;
                            for (std::size_t i = 0; i < state_words; i++) {
                                x[i] = algebra::random_element<field_type>();
                                y[i] = algebra::random_element<field_type>();
                            }

                            for (std::size_t i = 0; i < state_words; i++) {
                                for (std::size_t j = 0; j < state_words; j++) {
                                    if ((i != j && x[i] == x[j]) ||
                                        (i != j && y[i] == y[j]) ||
                                        (x[i] == y[j])) {
                                        secure_mds_found = false;
                                        break;
                                    }
                                    // We use minus in the next line, as is done in Mina implementation.
                                    // Original implementation uses + instead, but it doesn't matter,
                                    // since X and Y are random elements.
                                    mtx[i][j] = (x[i] - y[i]).inversed();
                                }
                                if (!secure_mds_found) {
                                    break;
                                }
                            }
                            // Determinant of the matrix must not be 0.
                            if (det(mtx) == 0)
                                secure_mds_found = false;

                            // TODO(martun): check that mds has NO eignevalues.
                            // if len(new_mds_matrix.characteristic_polynomial().roots()) == 0:
                            // return new_mds_matrix
                            // The original matrix security check is here: https://extgit.iaik.tugraz.at/krypto/hadeshash/-/blob/master/code/generate_params_poseidon.sage
                        }

                        // Transpose the matrix.
                        for (std::size_t i = 0; i < state_words; i++) {
                            for (std::size_t j = 0; j < state_words; j++) {
                                transposed_mtx[i][j] = mtx[j][i];
                            }
                        }

                        return transposed_mtx;
                    }

                    constexpr static round_constants_type generate_round_constants() {
                        round_constants_type rc;

                        integral_type constant = 0;
                        lfsr_state_type lfsr_state = get_lfsr_init_state();

                        for (std::size_t r = 0; r < full_rounds + part_rounds; r++) {
                            for (std::size_t i = 0; i < state_words; i++) {
                                while (true) {
                                    constant = 0;
                                    for (std::size_t j = 0; j < word_bits; j++) {
                                        lfsr_state = update_lfsr_state(lfsr_state);
                                        constant = set_new_bit<integral_type>(
                                                constant, get_lfsr_state_bit(lfsr_state, lfsr_state_bits - 1));
                                    }
                                    if (constant < modulus) {
                                        rc[r][i] = element_type(constant);
                                        break;
                                    }
                                }
                            }
                        }
                        return rc;
                    }

                    static constexpr lfsr_state_type get_lfsr_init_state() {
                        lfsr_state_type state = 0;
                        int i = 0;
                        for (i = 1; i >= 0; i--)
                            state = set_new_bit(state, (1U >> i) & 1U); // field - as in filecoin
                        for (i = 3; i >= 0; i--)
                            state = set_new_bit(state, (1U >> i) & 1U); // s-box - as in filecoin
                        for (i = 11; i >= 0; i--)
                            state = set_new_bit(state, (word_bits >> i) & 1U);
                        for (i = 11; i >= 0; i--)
                            state = set_new_bit(state, (state_words >> i) & 1U);
                        for (i = 9; i >= 0; i--)
                            state = set_new_bit(state, (full_rounds >> i) & 1U);
                        for (i = 9; i >= 0; i--)
                            state = set_new_bit(state, (part_rounds >> i) & 1U);
                        for (i = 29; i >= 0; i--)
                            state = set_new_bit(state, 1);
                        // idling
                        for (i = 0; i < 160; i++)
                            state = update_lfsr_state_raw(state);
                        return state;
                    }

                    static constexpr lfsr_state_type update_lfsr_state(lfsr_state_type state) {
                        while (true) {
                            state = update_lfsr_state_raw(state);
                            if (get_lfsr_state_bit(state, lfsr_state_bits - 1))
                                break;
                            else
                                state = update_lfsr_state_raw(state);
                        }
                        return update_lfsr_state_raw(state);
                    }

                    static constexpr inline lfsr_state_type update_lfsr_state_raw(lfsr_state_type state) {
                        bool new_bit = get_lfsr_state_bit(state, 0) != get_lfsr_state_bit(state, 13) !=
                                       get_lfsr_state_bit(state, 23) != get_lfsr_state_bit(state, 38) !=
                                       get_lfsr_state_bit(state, 51) != get_lfsr_state_bit(state, 62);
                        return set_new_bit(state, new_bit);
                    }

                    static constexpr inline bool get_lfsr_state_bit(const lfsr_state_type &state, std::size_t pos) {
                        return bit_test(state, lfsr_state_bits - 1 - pos);
                    }

                    template<typename T>
                    static constexpr inline T set_new_bit(T var, bool new_bit) {
                        return (var << 1) | (new_bit ? 1 : 0);
                    }

                public:
                    /*!
                     * @brief Randomly generates all the constants required, using the correct generation rules.
                     * If called multiple times, will return DIFFERENT constants.
                     */

                    constexpr static const round_constants_type round_constants = generate_round_constants();
                    constexpr static const mds_matrix_type mds_matrix = generate_mds_matrix();
                    constexpr static const constants_type constants = {mds_matrix, round_constants};
                };

                template<>
                class poseidon_constants<mina_poseidon_policy<algebra::fields::pallas_base_field>> {
                public:
                    typedef mina_poseidon_policy<algebra::fields::pallas_base_field> policy_type;
                    typedef typename policy_type::field_type field_type;

                    constexpr static const std::size_t state_words = policy_type::state_words;
                    constexpr static const std::size_t word_bits = policy_type::word_bits;
                    constexpr static const std::size_t full_rounds = policy_type::full_rounds;
                    constexpr static const std::size_t half_full_rounds = policy_type::half_full_rounds;
                    constexpr static const std::size_t part_rounds = policy_type::part_rounds;
                    constexpr static const std::size_t rate = policy_type::block_words;

                    typedef typename field_type::value_type element_type;
                    typedef typename field_type::integral_type integral_type;
                    constexpr static const integral_type modulus = field_type::modulus;

                    constexpr static const std::size_t lfsr_state_bits = 80;
                    typedef number<backends::cpp_int_modular_backend<lfsr_state_bits>> lfsr_state_type;

                    typedef algebra::vector<element_type, state_words> state_vector_type;

                    typedef algebra::matrix<element_type, full_rounds + part_rounds, state_words> round_constants_type;
                    constexpr static const round_constants_type round_constants = {
                            0x2ec559cd1a1f2f6889fc8ae5f07757f202b364429677c8ff6603fd6d93659b47_cppui_modular253,
                            0x2553b08c788551bfe064d91c17eb1edb8662283229757711b2b30895f0aa3bad_cppui_modular253,
                            0x25a706fb0f35b260b6f28d61e082d36a8f161be1f4d9416371a7b65f2bfafe4e_cppui_modular253,
                            0x37c0281fda664cc2448d0e7dd77aaa04752250817a945abeea8cfaaf3ee39ba0_cppui_modular253,
                            0x140488321291998b8582eaceeb3fa9ca3980eb64a453573c5aaa2910405936b6_cppui_modular253,
                            0x3a73fe35b1bdd66b809aad5eab47b5c83b0146fd7fc632dfb49cd91ae1169378_cppui_modular253,
                            0x21b7c2b35fd7710b06245711f26c0635d3e21de4db10dd3a7369f59f468d7be6_cppui_modular253,
                            0x1803a068d25fef2ef652c8a4847aa18a29d1885e7bf77fd6a34d66536d09cad7_cppui_modular253,
                            0x291de61c5e6268213772cf7e03c80c2e833eb77c58c46548d158a70fbbd9724b_cppui_modular253,
                            0x230043a0dc2dfab63607cbe1b9c482fdd937fdefecc6905aa5012e89babead13_cppui_modular253,
                            0x218af77a05c502d3fa3144efcf47a0f2a0292498c10c6e2368565674e78764f4_cppui_modular253,
                            0x223e2d94c177d27e071d55729d13a9b216955c7102cc9a95ea40058efb506117_cppui_modular253,
                            0x2a18257c15ad9b6fe8b7c5ad2129394e902c3c3802e738f24ce2f585ae5f6a38_cppui_modular253,
                            0xa6f7ba75f216403d2e4940469d199474a65aa5ef814e36400bddef06158dcf8_cppui_modular253,
                            0x169be41c6227956efef5b4cdde65d00d5e04fe766178bdc731615c6e5b93e31e_cppui_modular253,
                            0x2e28f50a9a55d2e91774083072734544417e290a1cfebc01801b94d0728fe663_cppui_modular253,
                            0xfdedf8da8654a22831040cfc74432464b173ee68628fd90498480b9902f2819_cppui_modular253,
                            0x46a3ed9863d2d739dd8bc9e90a746fda1197162d0a0bec3db1f2f6042cf04e2_cppui_modular253,
                            0x219e08b460c305b428670bacab86ac1e9458075778d35c3619ae7ba1f9b2ed76_cppui_modular253,
                            0x38bb36a12ebcec4d4e8728eb43e3f12a6e33b1ffa1463379018d4e12424e62ca_cppui_modular253,
                            0x1e9aa3fe25d116ccfbd6a8fccdae0aa9bc164a03ab7e951704ee9a715fbedee6_cppui_modular253,
                            0x30f33ed70da4c2bfb844ff1a7558b817d1ec300da86a1694f2db45047d5f18b_cppui_modular253,
                            0x282b04137350495ab417cf2c47389bf681c39f6c22d9e370b7af75cbcbe4bb1_cppui_modular253,
                            0x9b1528dea2eb5bd96905b88ff05fdf3e0f220fe1d93d1b54953ac98fec825f0_cppui_modular253,
                            0x30083dbbb5eab39311c7a8bfd5e55567fa864b3468b5f9200e529cda03d9ef71_cppui_modular253,
                            0x17eace73cf67c6112239cbf51dec0e714ee4e5a91dbc9209dc17bbea5bcd094_cppui_modular253,
                            0x37af1de8f5475ba165b90f8d568683d54e215df97e9287943370cf4118428097_cppui_modular253,
                            0x16ff7592836a45340ec6f2b0f122736d03f0bcb84012f922a4baa73ea0e66f51_cppui_modular253,
                            0x1a5985d4b359d03de60b2edabb1853f476915febc0e40f83a2d1d0084efc3fd9_cppui_modular253,
                            0x255a9d4beb9b5ea18ab9782b1abb267fc5b773b98ab655fd4d469698e1e1f975_cppui_modular253,
                            0x34a8d9f45200a9ac28021712be81e905967bac580a0b9ee57bc4231f5ecb936a_cppui_modular253,
                            0x979556cb3edcbe4f33edd2094f1443b4b4ec6c457b0425b8463e788b9a2dcda_cppui_modular253,
                            0x2a4d028c09ad39c30666b78b45cfadd5279f6239379c689a727f626679272654_cppui_modular253,
                            0xc31b68f6850b3bd71fe4e89984e2c87415523fb54f24ec8ae71430370154b33_cppui_modular253,
                            0x1a27ca0b953d3dba6b8e01cf07d76c611a211d139f2dff5ac023ed2454f2ed90_cppui_modular253,
                            0x109ae97c25d60242b86d7169196d2212f268b952dfd95a3937916b9905303180_cppui_modular253,
                            0x3698c932f2a16f7bb9abac089ec2de79c9965881708878683caf53caa83ad9c4_cppui_modular253,
                            0x3c7e25e0ac8fba3dc1360f8a9a9fa0be0e031c8c76a93497b7cac7ed32ade6c0_cppui_modular253,
                            0x2fc5023c5e4aed5aa7dfca0f5492f1b6efab3099360ec960237512f48c858a79_cppui_modular253,
                            0x2c124735f3f924546fb4fdfa2a018e03f53063d3a2e87fd285ba8d647eda6765_cppui_modular253,
                            0x12c875c9b79591acf9033f8b6c1e357126c44b23f3486fbee0d98340a3382251_cppui_modular253,
                            0x3cda935e895857d39a7db8476aeda5a5131cb165a353073fd3e473fd8855528d_cppui_modular253,
                            0x218eb756fa5f1df9f1eb922ef80b0852588779a7368e3d010def1512815d8759_cppui_modular253,
                            0x23bcf1032957015ef171fbb4329bca0c57d59885522f25f4b082a3cf301cfbc6_cppui_modular253,
                            0x17474c3b6a9bc1057df64b9e4d62badbc7f3867b3dd757c71c1f656205d7bceb_cppui_modular253,
                            0x19826c0ee22972deb41745d3bd412c2ae3d4c18535f4b60c9e870edffa3d550_cppui_modular253,
                            0x30bcb17dfd622c46f3275f698319b68d8816bed0368ded435ed61992bc43efa9_cppui_modular253,
                            0x3bd816c214c66410229cfbd1f4a3a42e6a0f82f3c0d49b09bc7b4c042ff2c94b_cppui_modular253,
                            0x8943ec01d9fb9f43c840757738979b146c3b6d1982280e92a52e8d045633ea1_cppui_modular253,
                            0x2670bf8c01822e31c70976269d89ed58bc79ad2f9d1e3145df890bf898b57e47_cppui_modular253,
                            0xdd53b41599ae78dbd3e689b65ebcca493effa94ed765eeec75a0d3bb20407f9_cppui_modular253,
                            0x68177d293585e0b8c8e76a8a565c8689a1d88e6a9afa79220bb0a2253f203c3_cppui_modular253,
                            0x35216f471043866edc324ad8d8cf0cc792fe7a10bf874b1eeac67b451d6b2cf5_cppui_modular253,
                            0x1fd6efb2536bfe11ec3736e7f7448c01eb2a5a9041bbf84631cc83ee0464f6af_cppui_modular253,
                            0x2c982c7352102289fc1b48dafcd9e3cc364d5a4324575e4721daf0af10033c67_cppui_modular253,
                            0x352f7e8c7662d86db9c722d4d07778858771b832af5bb5dc3b13cf94851c1b45_cppui_modular253,
                            0x18e3c0c1caa5e3ed66ee1ab6f55a5c8063d8c9b034ae47db43435147149e37d5_cppui_modular253,
                            0x3124b12deb37dcbb3d96c1a08d507523e30e03e0919559bf2daaab238422eade_cppui_modular253,
                            0x143bf0def31437eb21095200d2d406e6e5727833683d9740b9bfc1713215dc9a_cppui_modular253,
                            0x1ebee92143f32b4f9d9a90ad62b8483c977480767b53c71f6bde934a8ef38f17_cppui_modular253,
                            0xff6c794ad1afaa494088d5f8ee6c47bf9e83013478628cf9f41f2e81383ebeb_cppui_modular253,
                            0x3d0a10ac3ee707c62e8bdf2cdb49ac2cf4096cf41a7f214fdd1f8f9a24804f17_cppui_modular253,
                            0x1d61014cd3ef0d87d037c56bdfa370a73352b95d472ead1937bed06a31801c91_cppui_modular253,
                            0x123e185b2ec7f072507ac1e4e743589bb25c8fdb468e329e7de169875f90c525_cppui_modular253,
                            0x30b780c0c1cb0609623732824c75017da9799bdc7e08b527bae7f409ebdbecf2_cppui_modular253,
                            0x1dfb3801b7ae4e209f68195612965c6e37a2ed5cf1eeee3d46edf655d6f5afef_cppui_modular253,
                            0x2fdee42805b2774064e963c741552556019a9611928dda728b78311e1f049528_cppui_modular253,
                            0x31b2b65c431212ed36fdda5358d90cd9cb51c9f493bff71cdc75654547e4a22b_cppui_modular253,
                            0x1e3ca033d8413b688db7a543e62ac2e69644c0614801379cfe62fa220319e0ef_cppui_modular253,
                            0xc8ef1168425028c52a32d93f9313153e52e9cf15e5ec2b4ca09d01730dad432_cppui_modular253,
                            0x378c73373a36a5ed94a34f75e5de7a7a6187ea301380ecfb6f1a22cf8552638e_cppui_modular253,
                            0x3218aeec20048a564015e8f221657fbe489ba404d7f5f15b829c7a75a85c2f44_cppui_modular253,
                            0x3312ef7cbbad31430f20f30931b070379c77119c1825c6560cd2c82cf767794e_cppui_modular253,
                            0x356449a71383674c607fa31ded8c0c0d2d20fb45c36698d258cecd982dba478c_cppui_modular253,
                            0xcc88d1c91481d5321174e55b49b2485682c87fac2adb332167a20bcb57db359_cppui_modular253,
                            0x1defccbd33740803ad284bc48ab959f349b94e18d773c6c0c58a4b9390cc300f_cppui_modular253,
                            0x2d263cc2e9af126d768d9e1d2bf2cbf32063be831cb1548ffd716bc3ee7034fe_cppui_modular253,
                            0x111e314db6fb1a28e241028ce3d347c52558a33b6b11285a97fffa1b479e969d_cppui_modular253,
                            0x27409401e92001d434cba2868e9e371703199c2372d23ef329e537b513f453e_cppui_modular253,
                            0x24a852bdf9cb2a8fedd5e85a59867d4916b8a57bdd5f84e1047d410770ffffa0_cppui_modular253,
                            0x205d1b0ee359f621845ac64ff7e383a3eb81e03d2a2966557746d21b47329d6e_cppui_modular253,
                            0x25c327e2cc93ec6f0f23b5e41c931bfbbe4c12da7d55a2b1c91c79db982df903_cppui_modular253,
                            0x39df3e22d22b09b4265da50ef175909ce79e8f0b9599dff01cf80e70884982b9_cppui_modular253,
                            0x9b08d58853d8ac908c5b14e5eb8611b45f40faaa59cb8dff98fb30efcdfaa01_cppui_modular253,
                            0x1ece62374d79e717db4a68f9cddaaf52f8884f397375c0f3c5c1dbaa9c57a0a6_cppui_modular253,
                            0x3bd089b727a0ee08e263fa5e35b618db87d7bcce03441475e3fd49639b9fa1c1_cppui_modular253,
                            0x3fedea75f37ad9cfc94c95141bfb4719ee9b32b874b93dcfc0cc12f51a7b2aff_cppui_modular253,
                            0x36dfa18a9ba1b194228494a8acaf0668cb43aca9d4e0a251b20ec3424d0e65cd_cppui_modular253,
                            0x119e98db3f49cd7fcb3b0632567d9ccaa5498b0d411a1437f57c658f41931d0c_cppui_modular253,
                            0x1100b21c306475d816b3efcd75c3ae135c54ad3cc56ca22abd9b7f45e6d02c19_cppui_modular253,
                            0x15791f9bbea213937208c82794eb667f157f003c65b64aa9800f4bbee4ea5119_cppui_modular253,
                            0x1adbeb5e9c4d515ecfd250ebee56a2a816eb3e3dc8d5d440c1ab4285b350be64_cppui_modular253,
                            0x1fbf4738844a9a249aec253e8e4260e4ab09e26bea29ab0020bf0e813ceecbc3_cppui_modular253,
                            0x3418a929556ec51a086459bb9e63a821d407388cce83949b9af3e3b0434eaf0e_cppui_modular253,
                            0x9406b5c3af0290f997405d0c51be69544afb240d48eeab1736cda0432e8ff9e_cppui_modular253,
                            0x23ece5d70b38ccc9d43cd923e5e3e2f62d1d873c9141ef01f89b6de1336f5bc7_cppui_modular253,
                            0x1852d574e46d370a0b1e64f6c41eeb8d40cf96c524a62965661f2ef87e67234d_cppui_modular253,
                            0xa657027cce8d4f238ea896dde273b7537b508674a366c66b3789d9828b0ce90_cppui_modular253,
                            0x3482f98a46ec358108fbbb68fd94f8f2baa73c723baf21922a850e45511f5a2d_cppui_modular253,
                            0x3f62f164f8c905b335a6cbf76131d2430237e17ad6abc76d2a6329c1ec5463ee_cppui_modular253,
                            0x7e397f503f9c1cea028465b2950ea444b15c5eab567d5a69ea2925685694df0_cppui_modular253,
                            0x405f1fc711872373d6eb50a09fbfb05b2703ae0a0b4edb86aedb216db17a876_cppui_modular253,
                            0xbe0848eb3e09c7027110ad842c502441c97afa14a844406fcfec754a25658c1_cppui_modular253,
                            0x26b78788fd98ac020bac92d0e7792bb5ffed06b697d847f61d984f905d9ba870_cppui_modular253,
                            0x38fd5318d39055c82fef9bdd33315a541c0ec4363e6cc0687005871355dfa573_cppui_modular253,
                            0x380bd03b840c48c8ba3830e7cace72f91a5002218c617294e8c8bc687d5216de_cppui_modular253,
                            0x2c6e57ddc1d7c81a0299ed49c3d74759416bc8426f30e2af5622895c531b4e1c_cppui_modular253,
                            0x11d3a81b262fc76ef506ee6d88e5991d0de8cb9dd162d97c58b175e3bc4584f3_cppui_modular253,
                            0x9b6b283ebaf45fbb1e448969ace9be62adf67ddf58614925741deb6a1ba7def_cppui_modular253,
                            0x15d5095164c885763fa83cdf776d436382821a17bc5563a5b6f6dfcdac504ade_cppui_modular253,
                            0x3427fdbfca3cea23063eb138c5055c6cad9c4252b23d12c12293308eff7d9124_cppui_modular253,
                            0x272f12e731077b74317ef2543c33b86194db1da5f6a7e1eee0656672c81685fe_cppui_modular253,
                            0x5323f85deb8c07c193c37a73d76f6114967913a2bdce11995f183e769f42967_cppui_modular253,
                            0x3d5ce415ecae4ba42b417ea3a501b44694f46efddff2fcca952b097f3852d3d8_cppui_modular253,
                            0xe8ec18c7b52c514d42047f1f0b2a90cb8c0c7391cf9479cd7fd5bfe1d3db8f2_cppui_modular253,
                            0x1591c865ea7065d54304519f8bb268bddbeaf3afae54edcd01a833ed0a9ef1a_cppui_modular253,
                            0x3eddbeeee5eca5deee4bf1789c435e1241e0d71186d8f0f62d74729dfc3119fb_cppui_modular253,
                            0x23691c7009b9283b268766e8d491716d3c1993e6ecf458def8f762af3e355707_cppui_modular253,
                            0x26cdab2c837ebeac5bea4be1d6f0488034907374d81a61a34f1c4db397d4c09b_cppui_modular253,
                            0x2d2206730664d58be0676dad1fee0e990c264a7410a2cdb6b55653c1df72ef56_cppui_modular253,
                            0x2bb74bb185372334a4ef5f6d18e2ece54086e62b04985dd794b7117b0be9217f_cppui_modular253,
                            0x366250fe928c45d8d5aa35f0a142754907ff3c598410199b589b28cd851b2204_cppui_modular253,
                            0x1868f8118482c6b4a5a61a81c8aaca128953179c20f73a44022d9976bdc34af1_cppui_modular253,
                            0xb7901c670e1d75d726eb88d000950b3c963f0f7a6ca24994bdc07ae2f78b4d3_cppui_modular253,
                            0x32c4bd8ab70e1f25af77af57dd340c8e6c8a101dfc5e8dd03314566db90b870_cppui_modular253,
                            0x1ce36db31fe6ea3cd9308db9aa43a8af5c41a8f0a6509bfe00f0e7b486c0ab8a_cppui_modular253,
                            0x26596ea9e1915e53da3479e9d13c3c920505e2449e325810ff6ca855fe4b7c6e_cppui_modular253,
                            0x30f296a269868a7fca8f5b1e269c0116304df31729559a270e713509d3a6d5dc_cppui_modular253,
                            0x2588961eff7897d87eb6ac72350ef9f52640647cbd23136919a994dfd1979d5_cppui_modular253,
                            0x16a49e69721e80690d41e06229e9bc2dbaf9a2abf4b89388db2485595409d62b_cppui_modular253,
                            0x3d7aca02c051fcad8073cfd67210cd423a31888afc4a444d9d3adf3d6c5da7bf_cppui_modular253,
                            0x299bd48a740b7790075268312ab8072c72421de5a6437fa5e25431ef951847b4_cppui_modular253,
                            0x11a69b867d9ea22ec1b2f28e96617129e36eefaea9e8126bdc6a42b99072902b_cppui_modular253,
                            0x25bc1af391f3c1f2284a95da92b5883d1b3a40794b2358b2e7a70fca22da64ce_cppui_modular253,
                            0x361ab3843f4d8ddadede39d82bb1a8109f89b6d9aa117b8f365de43895de0baa_cppui_modular253,
                            0x38ef3ab5b61c117a3465a017a9c8ba4c227659b41fdf145206d5c960f49dd45b_cppui_modular253,
                            0x3992f83f26143dbdbd335604a1a14daf238ae43c249783f694feaf560aaae20f_cppui_modular253,
                            0x350287977eb71c81b10ecd039aad99cfa9ed84a04301cb30869e1dc7fa1dc638_cppui_modular253,
                            0x3afb5bc126020586dcccba32dd054cd9a3f3b834ca9678d6802c48b1da97d6ed_cppui_modular253,
                            0x172b7c2d8e7e4b06d183a2575b790749d0970c54966407fa8f59072c729de671_cppui_modular253,
                            0x2eb53fe3a278688a70494569e54a0f0d269935aec6c897bef4d368c1f67d57e4_cppui_modular253,
                            0x375ae56b8d9310d553ed77d406dedc3f0393e5a321b71caee6a5bb7078b5035_cppui_modular253,
                            0x1d49a0d53bc2993cbf1fb5d1da9bb76fe46a7031d5e5d43fadbf54bc17c1ef38_cppui_modular253,
                            0x132d17b87cab6d707ddfa1f01df1724ad37957e989c44f1ff71426367f953160_cppui_modular253,
                            0x62da5280948d8c6c4acc7e6a1aa421f0f9ec179a44146750060be4be6755f85_cppui_modular253,
                            0xa4b4d5cde54a974ea4e57ee4132d2ab2510c300f21930d6bbbf211d1add80f9_cppui_modular253,
                            0x3356f1fbeac493ccab752b70bbed821ce49965c19284d7aacd78fbf3ff864e91_cppui_modular253,
                            0x42721e8a9cc32557851feb0e0190c5dfbf4cb1b8f47d37e7e653ec6ff8a4059_cppui_modular253,
                            0x53d9b2633fff31ca4fc5724ce6b4422318128cdf01897d321e86f47cdf748b1_cppui_modular253,
                            0x267d96caeafde5dbd3db1f0668b09ccd532a22f0205494716a786219fb4c801c_cppui_modular253,
                            0x39316997737610193c3f9ffcfd4e23d38aac12cd7b95b8d256d774101650a6ca_cppui_modular253,
                            0x191e377462986563fdabf9b23529f7c84c6b200b9101b3a5096bca5f377981fb_cppui_modular253,
                            0x20f89af9722f79c860d2059a0ec209cf3a7925ad0798cab655eca62fe73ff3d9_cppui_modular253,
                            0x1ca568aeddb2ef391a7c78ecf104d32d785b9ca145d97e35879df3534a7d1e0b_cppui_modular253,
                            0x25de9ba0a37472c3b4c0b9c3bc25cbbf78d91881b6f94ee70e4abf090211251c_cppui_modular253,
                            0x3393debd38d311881c7583bee07e605ef0e55c62f0508ccc2d26518cd568e1ef_cppui_modular253,
                            0x38df2fd18a8d7563806aa9d994a611f642d5c397388d1dd3e78bc7a4515c5b1_cppui_modular253,
                            0x5c6503ff1ee548f2435ad9148d7fb94c9222b0908f445537a6667047f6d501c_cppui_modular253,
                            0x104c88d6d0682d82d3d664826dc9565db101a220aa8f90572eb798468a82a2ab_cppui_modular253,
                            0x2caad6108c09ee6aee7851b4a2d2d3b7c3ca3c56a80003c8471f90bfa4ac628b_cppui_modular253,
                            0xa57dbd4c327826c8a97bc7285f94bcddb966177346f1792c4bd7088aa0353f3_cppui_modular253,
                            0x3c15552f9124318b8433d01bb53ba04ba1cc9eb91d83b918e32fea39fbe908fa_cppui_modular253,
                            0xe10c10cbbe1717a9441c6299c4fc087c222208bd4fa8f3be66d2075f623b513_cppui_modular253,
                            0x1e8b254cbff2c92a83dff1728c81dd22a9570f590e497cb2d640042cb879a930_cppui_modular253,
                            0x1812dbcd70c440610057bbfdd0cc4d31d1faf5786419b53841c4adc43f2b2352_cppui_modular253
                    };

                    typedef algebra::matrix<element_type, state_words, state_words> mds_matrix_type;
                    constexpr static const mds_matrix_type mds_matrix = {
                            0x1a9bd250757e29ef4959b9bef59b4e60e20a56307d6491e7b7ea1fac679c7903_cppui_modular253,
                            0x384aa09faf3a48737e2d64f6a030aa242e6d5d455ae4a13696b48a7320c506cd_cppui_modular253,
                            0x3d2b7b0209bc3080064d5ce4a7a03653f8346506bfa6d076061217be9e6cfed5_cppui_modular253,
                            0x9ee57c70bc351220b107983afcfabbea79868a4a8a5913e24b7aaf3b4bf3a42_cppui_modular253,
                            0x20989996bc29a96d17684d3ad4c859813115267f35225d7e1e9a5b5436a2458f_cppui_modular253,
                            0x14e39adb2e171ae232116419ee7f26d9191edde8a5632298347cdb74c3b2e69d_cppui_modular253,
                            0x174544357b687f65a9590c1df621818b5452d5d441597a94357f112316ef67cb_cppui_modular253,
                            0x3ca9263dc1a19d17cfbf15b0166bb25f95dffc53212db207fcee35f02c2c4137_cppui_modular253,
                            0x3cf1fbef75d4ab63b7a812f80b7b0373b2dc21d269ba7c4c4d6581d50aae114c_cppui_modular253
                    };

                    typedef std::pair<mds_matrix_type, round_constants_type> constants_type;
                    constexpr static const constants_type constants = {mds_matrix, round_constants};

                    inline constexpr static const element_type &round_constant(std::size_t round, std::size_t i) {
                        return round_constants[round][i];
                    }

                    inline constexpr static void product_with_mds_matrix(state_vector_type &A_vector) {
                        A_vector = algebra::vectmatmul(A_vector, mds_matrix);
                    }
                };

                template<>
                class poseidon_constants<mina_poseidon_policy<algebra::fields::vesta_base_field>> {
                public:
                    typedef mina_poseidon_policy<algebra::fields::vesta_base_field> policy_type;
                    typedef typename policy_type::field_type field_type;

                    constexpr static const std::size_t state_words = policy_type::state_words;
                    constexpr static const std::size_t word_bits = policy_type::word_bits;
                    constexpr static const std::size_t full_rounds = policy_type::full_rounds;
                    constexpr static const std::size_t half_full_rounds = policy_type::half_full_rounds;
                    constexpr static const std::size_t part_rounds = policy_type::part_rounds;
                    constexpr static const std::size_t rate = policy_type::block_words;

                    typedef typename field_type::value_type element_type;
                    typedef typename field_type::integral_type integral_type;
                    constexpr static const integral_type modulus = field_type::modulus;

                    constexpr static const std::size_t lfsr_state_bits = 80;
                    typedef number<backends::cpp_int_modular_backend<lfsr_state_bits>> lfsr_state_type;

                    typedef algebra::vector<element_type, state_words> state_vector_type;

                    inline constexpr static const element_type &round_constant(std::size_t round, std::size_t i) {
                        return round_constants[round][i];
                    }

                    inline constexpr static void product_with_mds_matrix(state_vector_type &A_vector) {
                        A_vector = algebra::vectmatmul(A_vector, mds_matrix);
                    }

                    typedef algebra::matrix<element_type, full_rounds + part_rounds, state_words> round_constants_type;
                    constexpr static const round_constants_type round_constants = {
                            0x590ef2a14ba3cef7e8f93a6dde4d481057d5d0547f6f09341b6b8be19c00ee6_cppui_modular255,
                            0x77faa77ed78ff8b695859df34db5157f6b491567f5f382a8fce538f0e5ffe6f_cppui_modular255,
                            0x3e54b7c94955c8994ed16ec9950d59aca4c9b6e419ef4935682528c2eba2de50_cppui_modular255,
                            0x37d991dc8d4de3912355745c7d78f8b04516b14d30e29324bb5dd075ca0f0c1d_cppui_modular255,
                            0xc0614dd1cff6c6817aff09d82ef828e80caed4da023823088fd021020f81f0e_cppui_modular255,
                            0x3335e335a3fed44842359528b3e88e1824a173da819d7ee6905e82eed054243_cppui_modular255,
                            0xb2202aa54d42f4f07693766723b9624c9fca4d33a2b9ee40f1c809a15a48a1d_cppui_modular255,
                            0x290253e0e1d2c72b32a5b272137a0892b5934b0b8f26b4fc25ea00d63a70e9df_cppui_modular255,
                            0x3e99873e73025d7c8b71fd209d13dba7a1021013f0815ea33a42ae94b63d00f3_cppui_modular255,
                            0x164682f55ec314f639f5f8062a4ddf11ed80d5822591a22ff54f340d90165d85_cppui_modular255,
                            0x309ba21093c9d04c81bd5273ad1064e1bd9067312d3269dddadf74c2eb1d3e01_cppui_modular255,
                            0x159e72bb030cb8994b2eac1d4ee7d0f06b0b092e7611d460605b3d8c60a274d9_cppui_modular255,
                            0xd743dbfc6f3c833ce2ef4956bead3c118fd3198652038781903ac929218fdd6_cppui_modular255,
                            0x18cb5a9230eb74045ede834ac6dd129bd2a0462dca1d96d167b9be0e1e96a688_cppui_modular255,
                            0x2d82f85fc222b215902d61c85c968b39759d6c2e9aa0e11fd08881bfae311e66_cppui_modular255,
                            0x2920828be5972cb8ff8023386a90a837bbfcca99be240137f7d211ecb72521e6_cppui_modular255,
                            0x3101774e1c3d72d010efb29c16c476e988bdb47321af3f82e05cc9c6b0360853_cppui_modular255,
                            0x327b4e6353c099e41a8ffab9103996b9d29d07da0f1a191aa6fb55c0720c6f54_cppui_modular255,
                            0x71c29018dd48d5c557379ea9d4afd80b92788ed509ced6bac47a65ba8b475c_cppui_modular255,
                            0x25efdeef6c5ad56834b24cfe03d57360b4335ec902c78ee9348ebaceab726038_cppui_modular255,
                            0x109ffe5cd918fcd7da7fdb40d32ac406f453874fda431c35c9e35601bcf708e9_cppui_modular255,
                            0x1f4de5d78b4378e0eca49ed94999d8bc91489fadfd896c8affbaa6e2654d18bf_cppui_modular255,
                            0x173185e1eaad0664ba1c01b8e417a4422c22a43d622c5df98c11481e205e499e_cppui_modular255,
                            0x161a0e8b31a6fd42727dc0a37ae4f715683af35873bd37e78e10abcb0e21fabd_cppui_modular255,
                            0x3decab3f42934acc644cc227315ecd6bcee79e5d92dc686823f60e6a3c40a7cd_cppui_modular255,
                            0x29d7541d2a4fcdf9c7f144ce1e957a5e5c6d5d064618416817d0ad39708b2807_cppui_modular255,
                            0x1d0525558685977d321fe86c05f462ae2e569e6d202bd5c62b0815320454114a_cppui_modular255,
                            0x27d1aec0ccc80f71d09d2a9c0b76ee5fe9a87516f0e691a9f5fba360cb79f32_cppui_modular255,
                            0x1c28ed68159e54df8296e654b0c1b5872de41557b7b02adc256dcc1600229ba8_cppui_modular255,
                            0x15c9cbe29bf4e7d8bae22dd2213c86724e9944ea4b9e34b6681beb1b0972215e_cppui_modular255,
                            0xd479e19db4686f5cb1ef9a8331a1ab680c5d3770e9a9a8a7a6ac58f8006c38a_cppui_modular255,
                            0x3494f6ecf12d5c3d758c5380652154e26f7f3c888d362ea512da8dc265fc32b0_cppui_modular255,
                            0x37ed9343bcc46adb4300f3d8cb88c311383061710836351ded0a146de837966_cppui_modular255,
                            0x35548be14e1cbcbd7d2c0e8c4a95e5fc2893daba34197ef41c350ae7072cde4e_cppui_modular255,
                            0x34e58327efe8d41b81b66b6c3fad424b2ff9008392909bb90eb10f08462b998b_cppui_modular255,
                            0xf55c1223abf50500c4ac4103f679dcfea4eebb368cf64ef3a63ee27146846f_cppui_modular255,
                            0x11dd4ab1734f7069498cc390a41b7de375d8968cec91b5c74cef9812e8ee7ce7_cppui_modular255,
                            0x1e344f255d7c5e537439e75f9c4ea64dd1fda1b0988e5c83626055859369b43c_cppui_modular255,
                            0x147db9afad2d2f7c4249357587faba99a6a38da16fe9ba74ef2f3fc5a0878f44_cppui_modular255,
                            0x31774ce29d00f566bd499f181517df231be7205c05e7527d71a1c89cb0e841a7_cppui_modular255,
                            0x32bdf60a6685665871f654169996f508be8710c99f3fa6f44a7bc4d2c25fbfd8_cppui_modular255,
                            0x2f567f84ec13720611900c4b9e8303f04c8cc5c57daa4d95d9ee009514205e65_cppui_modular255,
                            0x2dbd279621e591da57f54459f4160dde2f5c78e478d20f2f4763832e013bc07f_cppui_modular255,
                            0x1275fb5ba53b7d2b5322e63f09a48026d684369c8e12241a808085a78ab3a369_cppui_modular255,
                            0x1dd0beba925fe1df13f732b03287cad943569d62ec9059afc2c8120655e97d78_cppui_modular255,
                            0xa37d78e392a5c8441f98e9dbd51a9151e78fb877885ecb885b0834c50cfea4d_cppui_modular255,
                            0x1ebb7e2592122cd16d27e13410b2b48d520d8e99d38c1d86af0ac13565dfeb88_cppui_modular255,
                            0x24a6454b0a69c59916d64f532b56226f8d49969432b7d0efc675f599c3bdb64f_cppui_modular255,
                            0x269668b3e7835df2f85b82e9ef8647c43205e799135ce669256bf55f07448209_cppui_modular255,
                            0x15c87375d4514bbdddbfd84e51f246446f1b16bb58bd4bd9fa2ff57e6aa66057_cppui_modular255,
                            0x11ce62bbe1242334c260a67817be908a9422d9b9c6ee96c00772fcc8fc501db6_cppui_modular255,
                            0x20348b7d6b381bfd3ac923d60b965086d281d8a654ad5f3210d277789641fe98_cppui_modular255,
                            0x1398d090fd1144d1e84798e3a0efa942abfe650947e4a3cfa409ff14b541fae9_cppui_modular255,
                            0x2461a1a2d6e3a0b2e5185ae6c844fe2a3b2d85dfb1cf891efc79ae80dd776bed_cppui_modular255,
                            0x3e1f1de94c4af008188ba5eaef1da9ab9792ce54eda56bb5a519a65cd808885b_cppui_modular255,
                            0x1dee6ead07fbc0fe883f4d397994d75ba3c4f90720e74ae2da13066bc3a7dc3b_cppui_modular255,
                            0x287d06396bcb63555cb2ff408ea075cf402b10a3c608043d0cf2e3685ec6e2ad_cppui_modular255,
                            0x36d84c953d584607478da6183dc4da71bdbf737d45fb57d5a53badc123ae071c_cppui_modular255,
                            0x24c8fd13d2687a9f90c61da26823d4934b350cfa488d528482399e106a70ac74_cppui_modular255,
                            0x52e052a6a493457c9476ccc4fd9924e5c7247b98e58a3cfa688c0f8314bea68_cppui_modular255,
                            0x2fd32bae8a40ab498f6ba290733bb82504de1be782c1cdf039e2fbc843a01e52_cppui_modular255,
                            0x4e8e7d3413c8c8ccfe154dc51f31c7682627c71fa4b50daab27f2a4d2623ea6_cppui_modular255,
                            0x20c16d0097cebeb385508b606487baaf3bad515ba8a0b977f15cb50239418e38_cppui_modular255,
                            0x34f1df6035aac75204368125b0c4cec107e2f9eb0005517d26d6113e1f366271_cppui_modular255,
                            0x375973b59ed7b4bdb33642d20e6364f37a942f9018f6bca5abc10705481425e0_cppui_modular255,
                            0x269e8c978803e51d43439b7c18c4260e819e09e7d8c8d38706463bbb811c698c_cppui_modular255,
                            0x21be1913f874f3edb88a1f60cd157fcb76ff20b4eb139aae205b5a2764098782_cppui_modular255,
                            0x37a0a8ba83db884f721c25027d188c7ab7c7840b7860675b33e1c93e4023927f_cppui_modular255,
                            0x56d0e67fde779b7be5f308a3ce119e23e0503e6dabdbbd5189bb44dc6a6f0a4_cppui_modular255,
                            0x144723436a329da5644cce96fee4952b066092c36bd12838b4ffd4283cfe82c4_cppui_modular255,
                            0xec0b5f14ba50aa2b022d06fbb920a2aafb465b8c7f81fc119371a4cbb6acff7_cppui_modular255,
                            0x685de18d9a346a35c44a2a4ac7283d6fe2e4a9dc058bd537700bc2495271721_cppui_modular255,
                            0x178dcb74b546adea41afd5d93ef564cb3adb0ef5200201daea0faa5026bb8cbc_cppui_modular255,
                            0x1c1dcb1ef6cf5f036ae0030bf78f1643c439843959dd74fa28ea3663735cc923_cppui_modular255,
                            0xcfae6c99994c5f702cba3b32a4e38f3764207bfe7cd9bf577633b41843ea138_cppui_modular255,
                            0x2838a02558716d2b49c06fb34c49cd820ec71e861caa935f4a303e42030ae3df_cppui_modular255,
                            0x2c1944f3ec2852ed6b50fbc4abbc8f284797b36a23b321d2763ef48b1a5a0212_cppui_modular255,
                            0x30a218acd109f04657954e82f9faccc477731f4a954cf8ac12d15ebd450e9dcb_cppui_modular255,
                            0x2488defa4553fa5bd5afbb5fd28a1e99c585c5f541c6242e702215b2212c1d23_cppui_modular255,
                            0x3d0c9d7282245c776daa1655697fa879e470a26fcbb3bea62fa8ff32a4f04e50_cppui_modular255,
                            0x33aac46524f32f3556ed16a0912ef27482c2afcacbfb99ced98394b6c0e3765f_cppui_modular255,
                            0x1858a5f543ab0a70cb3957e0884b146b42cc3863fba4e034145ab09cc77d428d_cppui_modular255,
                            0x2d9d6fae68eff2e79396617207e28dba3d793b1e3739d30e9e9b10644e9f99cd_cppui_modular255,
                            0x1747fab074b37cc1ca7dbf7d6dc136740f5d26e30319b3577fc8987f1247caae_cppui_modular255,
                            0x38f905db5128f24e498e36a84df5a58ed3c0b0ed8f39336eb792cb634f86b87_cppui_modular255,
                            0xfffe42ce4a87a0b3a9ebe7eedf16c0cdb29c959b6e594faa69c0727c6e825f_cppui_modular255,
                            0x314c3090cd0a465da95afd515c0771703e4ee2a8eabe8fa405daf8bd49bce458_cppui_modular255,
                            0x3e5fb71d9071c658c39fe64392e90bac65bdaf8f723b6790cce7dd7440ce06aa_cppui_modular255,
                            0x3e9fe7b8fd0aaa379fa7be0dbd64309607cc5b00474ef6670370e631902e98cd_cppui_modular255,
                            0x33ee4f76ff95bd735ec602ee6f4d1664caec27a7c435ead3b4c8df6cb51f010e_cppui_modular255,
                            0x1670c2080f2965bed3f49db0b63aef5f562b347235645b921c0132b01cc82130_cppui_modular255,
                            0x210565224e2ee64dd479be3a969dc30c65933352ba9b2271a0942bf1bf485743_cppui_modular255,
                            0x9a7c6dd48dfbf50b13055b30fe85f934be9518b8af074b88f9de4b1df689616_cppui_modular255,
                            0x1f9116811eaadf677e6cb50fb59ce0fab11fa9f0ddf1432403610e1932a7aa1c_cppui_modular255,
                            0x19b51a48c225daf9b34611ccc5ba077ebbc0a19cfc9bbbd78ade11cfa655075f_cppui_modular255,
                            0x3286d29eb60c3d6204eb534d13f40d1af6364f0fe1622a12ba5fa069886f31fe_cppui_modular255,
                            0x9bd403d05db137ea793f10b6dd087a74a78c9b01bcd6f9daf39af2ef57d346e_cppui_modular255,
                            0x3a71654023e43363e60889eac50eb1f17c044606886771eaaf851bb2d00b3aeb_cppui_modular255,
                            0x3415b94f62c59466f102442b4bae7d6bb348987154cce16bd187525a6fb5b443_cppui_modular255,
                            0x3ca35f0fc660092b81f15dd6f0b3d17a16a053480ef2f935fce806dd0d9a3466_cppui_modular255,
                            0x26e1360af7fdc62e9be08651c2c5900ed5aefcb0d84b3aa88e354c6658a07863_cppui_modular255,
                            0x30d05884174d7a1de9d34c89224d17f3b9dbdfb0793b54c0d2aaaeedcc357bd6_cppui_modular255,
                            0x2c7f66f8b0580236f025dd626520049a09e1bfff0e5fd9f69cbc70daf0ac56c4_cppui_modular255,
                            0xc5cb9a350d2dc463dd05dbd696e122c6917b76654180c323937dee44c6beb93_cppui_modular255,
                            0x14d4d799d43d91b4d09d9c2bfdc13a64b48d18750503324361f9bf7267ec9b92_cppui_modular255,
                            0x60c56a884cd6a1d3514f2895816b84e7160df5106e8d031710769be1ac5c04c_cppui_modular255,
                            0x23e15f37c21266c86ead998a46e42f6e97fbd5d1c384f51d8b54d051a80d753d_cppui_modular255,
                            0x25eb2911034ab6bef4a969653f5cc33e6914b8b6411f064ec01bcf157fea4e55_cppui_modular255,
                            0x1e95c04c5057abd1b43a2fbc942b2391d0e0daef873838b3494e6d5fb067a117_cppui_modular255,
                            0x1547602fc83558aa1327221fd220fa22bcb1f6ec42edb7cc05eff508c65883cb_cppui_modular255,
                            0x16b669eac31e72a9e739fb03fd7ea3882fc5791b157143929ae12fc2fefe8b3d_cppui_modular255,
                            0x7034f4e251a65c4423479dc9d5287a341c108e0b56e29a391f9a07a0ca822f1_cppui_modular255,
                            0x3fdf9d5731ba040dc568e61b8571ea95ead2e89f0a9856b2d12a7e87e43f5683_cppui_modular255,
                            0x33f2cdf6960139a0fb4a3a8127992e2abbd42847728425228a35ee72bd5b01c7_cppui_modular255,
                            0x35616d55033d8fc092398f6c58bfc6eaaf2ec9dd500122516f489dbc631457b_cppui_modular255,
                            0x1eca80189643df1473e98da93fe58a9576def0d192d4153faebcd1b210c1603f_cppui_modular255,
                            0x26223ca4af2d8d878ca5530c3e67ff1c95b50b9c5b8295e19150bc31ef90ba98_cppui_modular255,
                            0x19180fa5facb64ee9b4827ccd766622adf12fe80ab17c7395075368e10a2a361_cppui_modular255,
                            0x169f165855e097501f25d6b3aae815ce6e8a1c289850936d956657f0ed99446_cppui_modular255,
                            0x363a8f891de5974f06bae043bc6a26b4518d217af6590e9318e325fb215cda00_cppui_modular255,
                            0x122aaa7c330ddcb57180749e659600a4dfac5dda7b9b68ab0f8b2ee6de350ced_cppui_modular255,
                            0xed203defca13ebdf6af805a9f5dbdfef90007df2ad32fb1c83165e837ab5e3f_cppui_modular255,
                            0x11cce94bbc7a96e9708e99d3666c0a275329ac4bff42634a5f989ddcfc28fd68_cppui_modular255,
                            0x1705663587a03cb11485ac9d01fd10cb1138be1820d26a14e4ab7b1c0fdec8d2_cppui_modular255,
                            0x12ad28a60485a2d911639051971f43dd15a0dfd2f8a0de756f0c847fed63ed7d_cppui_modular255,
                            0xa9e61cc35eba9374eea117753aaaa93d6b29f550c2c54bce0a6078e05db9475_cppui_modular255,
                            0x72c3d62cf006a95dc8b2a53f878bb26fcaf3c28d709a91634f3a09f525054ad_cppui_modular255,
                            0x1ce8f168b446f7e797b91677fc46a975d2caa63dc359132c7c9729f5be24a7c_cppui_modular255,
                            0xe846a7211efda3d8115b5bf76aab7eac2b6099026fc7504fb81ac4a77c5560d_cppui_modular255,
                            0xabb8fd9d6fa3772022fa88800c12bdcbb1234473022cd141213d452255a0f55_cppui_modular255,
                            0x1c5d9938bc35a4832e8375dc307dba7a116d2a566e406ab31e8b03a36ec807cf_cppui_modular255,
                            0x35bea7ac6f40e0f50f08d325be9f051fd75ada8c03461f4d15b2c5e1a3d72431_cppui_modular255,
                            0x419357c205a7e1e028c0f49cbdeab85b82f4db78f1afb1b5568ec1bd2e48cb0_cppui_modular255,
                            0x1933e424c788e7466a159f1fe015ac7210f47044d9df6872cdfa227ae4d2190a_cppui_modular255,
                            0xde27ccdda95abb3d98db76d6f7f152a08d37ba81758beaf2eddbc58d13e560f_cppui_modular255,
                            0x35a312d5d6cbf00d55f097febaf9bd5eac5f2881ebf0afa377e2ba7cdcf2f51_cppui_modular255,
                            0xce6f415449ca515e4da9177527c9242adcc988de5e1846d07cdd5284f39f9d0_cppui_modular255,
                            0x38fd71543da5c4c0447dc22aa2c1e3744cb84eb1ff17040640b50f5ddf8c8e61_cppui_modular255,
                            0x158de859aad53c6a17de455ab067a09ad6cba22f4101d19e77d8a2975c0dc965_cppui_modular255,
                            0x2c300588eeae8cbc3814bd1d7646f472ef6b44a60c710bf6100937504e532c8b_cppui_modular255,
                            0xb198cf742a029409ac02397b91e2704fa94ecf147909fa8d71ece5087e2cfc3_cppui_modular255,
                            0x100b375c21d357d5679d8e6d9eb7bff8edd4575535bf651ba0b1bd83cfb54598_cppui_modular255,
                            0x15a474d44590e2b23b8bb1e79f5613f1659e7ae2bce10def0ce1a101eb3e3ce5_cppui_modular255,
                            0x2aa20e6642a989e1e6f9814c24f022991c23a7e40af505d4b931079025b7ed4d_cppui_modular255,
                            0x196597f2d65c5692706795bf46eb7be96b31647c23441213642ccceedc01ebc4_cppui_modular255,
                            0x248291aa516daa0a6cd191c1c651a82f7d1b5f087dcb7cee91a27c488483e2bd_cppui_modular255,
                            0x36c02b98ad2722b774aeb131b31bfd087c6a7f2d0a3faa40bd9899e5f270877f_cppui_modular255,
                            0x1240e06949a1ad92bd8ae90772b5d8505174182c87a23227aa74b7630dba4195_cppui_modular255,
                            0x3b83f7e36f30939a78ec63cb2554aa0669a1bfc1b8b8714c6b8a3958beb6a163_cppui_modular255,
                            0x1668b0582ce04f7f5b1e35e1b7cc3e05be23cc2c9e0be9436559193f2a8d102e_cppui_modular255,
                            0x26d6a708e9464c85e9c7605e87fb96036fd1fe87379ac43ad560885582e4026d_cppui_modular255,
                            0x594fccf1863993b43ad0a13c5fc7a53f59f7d622e7b206d425907243a69e62d_cppui_modular255,
                            0x78e4c588b6ddd0fe7ed53a9f25b6ac3c2eac1c63faecc7e916f4d4599051940_cppui_modular255,
                            0xf44ea3e14c3e4849ee7a525fe77170b8658a6753680e269c9fd1d12932af69d_cppui_modular255,
                            0x2e8567bc9e8e369bdf7748d6c7f677837c601455d4651a2f102b94ff1f951379_cppui_modular255,
                            0x37c35b056171982cc7d74e6081fcac2f764f1fe30ee985db306a22b097d51bae_cppui_modular255,
                            0x29dbcffd5b55d671c85ca42037ac5e64d2ef42d2704af47a20877e3a5e5f1d9d_cppui_modular255,
                            0x201098422e054c1ddcc465411d002d2bc5a824e1c7f4f2ded9443c37bd04a520_cppui_modular255,
                            0x7de32ed4c5143430ef43aef100f948ef859ab3793aa52640156f5e7d92cdc84_cppui_modular255,
                            0x34e95adcc0c5c34fd38ab9246a04cc1029f678ba53c0f6fd27f8805094e36199_cppui_modular255,
                            0x1d5faf157126c599232982356ca0ea7b81d875c01d842b5cd1998a5c470fa623_cppui_modular255,
                            0x160a80176bd281e3fa9b82e44063cc7bf86eb81397e51e41fe4745e27c57e1d2_cppui_modular255,
                            0x17ecc7f5deb148c542a22d02b098439724910a3bbd4903428c8fc680f31b2406_cppui_modular255,
                            0x20a6aae17f822bc7035da3b8931896c82152346f2a43ab4e0029dbf0101b3d_cppui_modular255,
                            0x9ea0ec10c0e77b9385a58ccd5ecc3c88b5bed58af72a6d87bb446e14fa7c8d6_cppui_modular255
                    };

                    typedef algebra::matrix<element_type, state_words, state_words> mds_matrix_type;
                    constexpr static const mds_matrix_type mds_matrix = {
                            0x3e28f7dd17f47a7e304a54d377dd7aeead6b92027d60baf300246cf023dd594e_cppui_modular255,
                            0x30db06abb696fccb92b28ac214f4893d3fd84b3d4a9018754975e24477c32600_cppui_modular255,
                            0x174110bc1b058c6016ff5e8152ab3ffb6e2e6c4d01e66aba302659c51b7f563a_cppui_modular255,
                            0x12d36fa83503146980c05a1d48bcd50d2e9d4390e353a158a0fe387e2b4aeb0c_cppui_modular255,
                            0x2ab17c8eb369bea76e9f0c385e8bafc71536bedc8e06d06fd65c1670e94d9c55_cppui_modular255,
                            0xcc915328165c13986af127e108b9e5d9a60c5dc92e3e7636b8c3da5b4a8537_cppui_modular255,
                            0x4d9a6d270696688eb4346153b380c613a3dcaf0fb5a1e8380409ae0a143d31b_cppui_modular255,
                            0x2a805eee3317c8bae1f7d15abe4d27fee5fabcf9a3334d18b1932a33774c324_cppui_modular255,
                            0x19b092e9c6dffd1eb1b6df2dbc00bb2283b9a787273dcbad9b8d89cd502b7bbd_cppui_modular255
                    };

                    typedef std::pair<mds_matrix_type, round_constants_type> constants_type;
                    constexpr static const constants_type constants = {mds_matrix, round_constants};
                };

                template<>
                class poseidon_constants<poseidon_policy<
                        algebra::fields::alt_bn128_scalar_field<254>, 128, 2>> {
                public:
                    typedef poseidon_policy<algebra::fields::alt_bn128_scalar_field<254>, 128, 2> policy_type;
                    typedef typename policy_type::field_type field_type;

                    constexpr static const std::size_t state_words = policy_type::state_words;
                    constexpr static const std::size_t word_bits = policy_type::word_bits;
                    constexpr static const std::size_t full_rounds = policy_type::full_rounds;
                    constexpr static const std::size_t half_full_rounds = policy_type::half_full_rounds;
                    constexpr static const std::size_t part_rounds = policy_type::part_rounds;
                    constexpr static const std::size_t rate = policy_type::block_words;

                    typedef typename field_type::value_type element_type;
                    typedef typename field_type::integral_type integral_type;
                    constexpr static const integral_type modulus = field_type::modulus;

                    constexpr static const std::size_t lfsr_state_bits = 80;
                    typedef number<backends::cpp_int_modular_backend<lfsr_state_bits>> lfsr_state_type;

                    typedef algebra::vector<element_type, state_words> state_vector_type;

                    inline constexpr static const element_type &round_constant(std::size_t round, std::size_t i) {
                        return round_constants[round][i];
                    }

                    inline constexpr static void product_with_mds_matrix(state_vector_type &A_vector) {
                        A_vector = algebra::vectmatmul(A_vector, mds_matrix);
                    }

                    typedef algebra::matrix<element_type, full_rounds + part_rounds, state_words> round_constants_type;
                    constexpr static const round_constants_type round_constants = {
                            0x0ee9a592ba9a9518d05986d656f40c2114c4993c11bb29938d21d47304cd8e6e_cppui_modular254,
                            0x00f1445235f2148c5986587169fc1bcd887b08d4d00868df5696fff40956e864_cppui_modular254,
                            0x08dff3487e8ac99e1f29a058d0fa80b930c728730b7ab36ce879f3890ecf73f5_cppui_modular254,
                            0x2f27be690fdaee46c3ce28f7532b13c856c35342c84bda6e20966310fadc01d0_cppui_modular254,
                            0x2b2ae1acf68b7b8d2416bebf3d4f6234b763fe04b8043ee48b8327bebca16cf2_cppui_modular254,
                            0x0319d062072bef7ecca5eac06f97d4d55952c175ab6b03eae64b44c7dbf11cfa_cppui_modular254,
                            0x28813dcaebaeaa828a376df87af4a63bc8b7bf27ad49c6298ef7b387bf28526d_cppui_modular254,
                            0x2727673b2ccbc903f181bf38e1c1d40d2033865200c352bc150928adddf9cb78_cppui_modular254,
                            0x234ec45ca27727c2e74abd2b2a1494cd6efbd43e340587d6b8fb9e31e65cc632_cppui_modular254,
                            0x15b52534031ae18f7f862cb2cf7cf760ab10a8150a337b1ccd99ff6e8797d428_cppui_modular254,
                            0x0dc8fad6d9e4b35f5ed9a3d186b79ce38e0e8a8d1b58b132d701d4eecf68d1f6_cppui_modular254,
                            0x1bcd95ffc211fbca600f705fad3fb567ea4eb378f62e1fec97805518a47e4d9c_cppui_modular254,
                            0x10520b0ab721cadfe9eff81b016fc34dc76da36c2578937817cb978d069de559_cppui_modular254,
                            0x1f6d48149b8e7f7d9b257d8ed5fbbaf42932498075fed0ace88a9eb81f5627f6_cppui_modular254,
                            0x1d9655f652309014d29e00ef35a2089bfff8dc1c816f0dc9ca34bdb5460c8705_cppui_modular254,
                            0x04df5a56ff95bcafb051f7b1cd43a99ba731ff67e47032058fe3d4185697cc7d_cppui_modular254,
                            0x0672d995f8fff640151b3d290cedaf148690a10a8c8424a7f6ec282b6e4be828_cppui_modular254,
                            0x099952b414884454b21200d7ffafdd5f0c9a9dcc06f2708e9fc1d8209b5c75b9_cppui_modular254,
                            0x052cba2255dfd00c7c483143ba8d469448e43586a9b4cd9183fd0e843a6b9fa6_cppui_modular254,
                            0x0b8badee690adb8eb0bd74712b7999af82de55707251ad7716077cb93c464ddc_cppui_modular254,
                            0x119b1590f13307af5a1ee651020c07c749c15d60683a8050b963d0a8e4b2bdd1_cppui_modular254,
                            0x03150b7cd6d5d17b2529d36be0f67b832c4acfc884ef4ee5ce15be0bfb4a8d09_cppui_modular254,
                            0x2cc6182c5e14546e3cf1951f173912355374efb83d80898abe69cb317c9ea565_cppui_modular254,
                            0x005032551e6378c450cfe129a404b3764218cadedac14e2b92d2cd73111bf0f9_cppui_modular254,
                            0x233237e3289baa34bb147e972ebcb9516469c399fcc069fb88f9da2cc28276b5_cppui_modular254,
                            0x05c8f4f4ebd4a6e3c980d31674bfbe6323037f21b34ae5a4e80c2d4c24d60280_cppui_modular254,
                            0x0a7b1db13042d396ba05d818a319f25252bcf35ef3aeed91ee1f09b2590fc65b_cppui_modular254,
                            0x2a73b71f9b210cf5b14296572c9d32dbf156e2b086ff47dc5df542365a404ec0_cppui_modular254,
                            0x1ac9b0417abcc9a1935107e9ffc91dc3ec18f2c4dbe7f22976a760bb5c50c460_cppui_modular254,
                            0x12c0339ae08374823fabb076707ef479269f3e4d6cb104349015ee046dc93fc0_cppui_modular254,
                            0x0b7475b102a165ad7f5b18db4e1e704f52900aa3253baac68246682e56e9a28e_cppui_modular254,
                            0x037c2849e191ca3edb1c5e49f6e8b8917c843e379366f2ea32ab3aa88d7f8448_cppui_modular254,
                            0x05a6811f8556f014e92674661e217e9bd5206c5c93a07dc145fdb176a716346f_cppui_modular254,
                            0x29a795e7d98028946e947b75d54e9f044076e87a7b2883b47b675ef5f38bd66e_cppui_modular254,
                            0x20439a0c84b322eb45a3857afc18f5826e8c7382c8a1585c507be199981fd22f_cppui_modular254,
                            0x2e0ba8d94d9ecf4a94ec2050c7371ff1bb50f27799a84b6d4a2a6f2a0982c887_cppui_modular254,
                            0x143fd115ce08fb27ca38eb7cce822b4517822cd2109048d2e6d0ddcca17d71c8_cppui_modular254,
                            0x0c64cbecb1c734b857968dbbdcf813cdf8611659323dbcbfc84323623be9caf1_cppui_modular254,
                            0x028a305847c683f646fca925c163ff5ae74f348d62c2b670f1426cef9403da53_cppui_modular254,
                            0x2e4ef510ff0b6fda5fa940ab4c4380f26a6bcb64d89427b824d6755b5db9e30c_cppui_modular254,
                            0x0081c95bc43384e663d79270c956ce3b8925b4f6d033b078b96384f50579400e_cppui_modular254,
                            0x2ed5f0c91cbd9749187e2fade687e05ee2491b349c039a0bba8a9f4023a0bb38_cppui_modular254,
                            0x30509991f88da3504bbf374ed5aae2f03448a22c76234c8c990f01f33a735206_cppui_modular254,
                            0x1c3f20fd55409a53221b7c4d49a356b9f0a1119fb2067b41a7529094424ec6ad_cppui_modular254,
                            0x10b4e7f3ab5df003049514459b6e18eec46bb2213e8e131e170887b47ddcb96c_cppui_modular254,
                            0x2a1982979c3ff7f43ddd543d891c2abddd80f804c077d775039aa3502e43adef_cppui_modular254,
                            0x1c74ee64f15e1db6feddbead56d6d55dba431ebc396c9af95cad0f1315bd5c91_cppui_modular254,
                            0x07533ec850ba7f98eab9303cace01b4b9e4f2e8b82708cfa9c2fe45a0ae146a0_cppui_modular254,
                            0x21576b438e500449a151e4eeaf17b154285c68f42d42c1808a11abf3764c0750_cppui_modular254,
                            0x2f17c0559b8fe79608ad5ca193d62f10bce8384c815f0906743d6930836d4a9e_cppui_modular254,
                            0x2d477e3862d07708a79e8aae946170bc9775a4201318474ae665b0b1b7e2730e_cppui_modular254,
                            0x162f5243967064c390e095577984f291afba2266c38f5abcd89be0f5b2747eab_cppui_modular254,
                            0x2b4cb233ede9ba48264ecd2c8ae50d1ad7a8596a87f29f8a7777a70092393311_cppui_modular254,
                            0x2c8fbcb2dd8573dc1dbaf8f4622854776db2eece6d85c4cf4254e7c35e03b07a_cppui_modular254,
                            0x1d6f347725e4816af2ff453f0cd56b199e1b61e9f601e9ade5e88db870949da9_cppui_modular254,
                            0x204b0c397f4ebe71ebc2d8b3df5b913df9e6ac02b68d31324cd49af5c4565529_cppui_modular254,
                            0x0c4cb9dc3c4fd8174f1149b3c63c3c2f9ecb827cd7dc25534ff8fb75bc79c502_cppui_modular254,
                            0x174ad61a1448c899a25416474f4930301e5c49475279e0639a616ddc45bc7b54_cppui_modular254,
                            0x1a96177bcf4d8d89f759df4ec2f3cde2eaaa28c177cc0fa13a9816d49a38d2ef_cppui_modular254,
                            0x066d04b24331d71cd0ef8054bc60c4ff05202c126a233c1a8242ace360b8a30a_cppui_modular254,
                            0x2a4c4fc6ec0b0cf52195782871c6dd3b381cc65f72e02ad527037a62aa1bd804_cppui_modular254,
                            0x13ab2d136ccf37d447e9f2e14a7cedc95e727f8446f6d9d7e55afc01219fd649_cppui_modular254,
                            0x1121552fca26061619d24d843dc82769c1b04fcec26f55194c2e3e869acc6a9a_cppui_modular254,
                            0x00ef653322b13d6c889bc81715c37d77a6cd267d595c4a8909a5546c7c97cff1_cppui_modular254,
                            0x0e25483e45a665208b261d8ba74051e6400c776d652595d9845aca35d8a397d3_cppui_modular254,
                            0x29f536dcb9dd7682245264659e15d88e395ac3d4dde92d8c46448db979eeba89_cppui_modular254,
                            0x2a56ef9f2c53febadfda33575dbdbd885a124e2780bbea170e456baace0fa5be_cppui_modular254,
                            0x1c8361c78eb5cf5decfb7a2d17b5c409f2ae2999a46762e8ee416240a8cb9af1_cppui_modular254,
                            0x151aff5f38b20a0fc0473089aaf0206b83e8e68a764507bfd3d0ab4be74319c5_cppui_modular254,
                            0x04c6187e41ed881dc1b239c88f7f9d43a9f52fc8c8b6cdd1e76e47615b51f100_cppui_modular254,
                            0x13b37bd80f4d27fb10d84331f6fb6d534b81c61ed15776449e801b7ddc9c2967_cppui_modular254,
                            0x01a5c536273c2d9df578bfbd32c17b7a2ce3664c2a52032c9321ceb1c4e8a8e4_cppui_modular254,
                            0x2ab3561834ca73835ad05f5d7acb950b4a9a2c666b9726da832239065b7c3b02_cppui_modular254,
                            0x1d4d8ec291e720db200fe6d686c0d613acaf6af4e95d3bf69f7ed516a597b646_cppui_modular254,
                            0x041294d2cc484d228f5784fe7919fd2bb925351240a04b711514c9c80b65af1d_cppui_modular254,
                            0x154ac98e01708c611c4fa715991f004898f57939d126e392042971dd90e81fc6_cppui_modular254,
                            0x0b339d8acca7d4f83eedd84093aef51050b3684c88f8b0b04524563bc6ea4da4_cppui_modular254,
                            0x0955e49e6610c94254a4f84cfbab344598f0e71eaff4a7dd81ed95b50839c82e_cppui_modular254,
                            0x06746a6156eba54426b9e22206f15abca9a6f41e6f535c6f3525401ea0654626_cppui_modular254,
                            0x0f18f5a0ecd1423c496f3820c549c27838e5790e2bd0a196ac917c7ff32077fb_cppui_modular254,
                            0x04f6eeca1751f7308ac59eff5beb261e4bb563583ede7bc92a738223d6f76e13_cppui_modular254,
                            0x2b56973364c4c4f5c1a3ec4da3cdce038811eb116fb3e45bc1768d26fc0b3758_cppui_modular254,
                            0x123769dd49d5b054dcd76b89804b1bcb8e1392b385716a5d83feb65d437f29ef_cppui_modular254,
                            0x2147b424fc48c80a88ee52b91169aacea989f6446471150994257b2fb01c63e9_cppui_modular254,
                            0x0fdc1f58548b85701a6c5505ea332a29647e6f34ad4243c2ea54ad897cebe54d_cppui_modular254,
                            0x12373a8251fea004df68abcf0f7786d4bceff28c5dbbe0c3944f685cc0a0b1f2_cppui_modular254,
                            0x21e4f4ea5f35f85bad7ea52ff742c9e8a642756b6af44203dd8a1f35c1a90035_cppui_modular254,
                            0x16243916d69d2ca3dfb4722224d4c462b57366492f45e90d8a81934f1bc3b147_cppui_modular254,
                            0x1efbe46dd7a578b4f66f9adbc88b4378abc21566e1a0453ca13a4159cac04ac2_cppui_modular254,
                            0x07ea5e8537cf5dd08886020e23a7f387d468d5525be66f853b672cc96a88969a_cppui_modular254,
                            0x05a8c4f9968b8aa3b7b478a30f9a5b63650f19a75e7ce11ca9fe16c0b76c00bc_cppui_modular254,
                            0x20f057712cc21654fbfe59bd345e8dac3f7818c701b9c7882d9d57b72a32e83f_cppui_modular254,
                            0x04a12ededa9dfd689672f8c67fee31636dcd8e88d01d49019bd90b33eb33db69_cppui_modular254,
                            0x27e88d8c15f37dcee44f1e5425a51decbd136ce5091a6767e49ec9544ccd101a_cppui_modular254,
                            0x2feed17b84285ed9b8a5c8c5e95a41f66e096619a7703223176c41ee433de4d1_cppui_modular254,
                            0x1ed7cc76edf45c7c404241420f729cf394e5942911312a0d6972b8bd53aff2b8_cppui_modular254,
                            0x15742e99b9bfa323157ff8c586f5660eac6783476144cdcadf2874be45466b1a_cppui_modular254,
                            0x1aac285387f65e82c895fc6887ddf40577107454c6ec0317284f033f27d0c785_cppui_modular254,
                            0x25851c3c845d4790f9ddadbdb6057357832e2e7a49775f71ec75a96554d67c77_cppui_modular254,
                            0x15a5821565cc2ec2ce78457db197edf353b7ebba2c5523370ddccc3d9f146a67_cppui_modular254,
                            0x2411d57a4813b9980efa7e31a1db5966dcf64f36044277502f15485f28c71727_cppui_modular254,
                            0x002e6f8d6520cd4713e335b8c0b6d2e647e9a98e12f4cd2558828b5ef6cb4c9b_cppui_modular254,
                            0x2ff7bc8f4380cde997da00b616b0fcd1af8f0e91e2fe1ed7398834609e0315d2_cppui_modular254,
                            0x00b9831b948525595ee02724471bcd182e9521f6b7bb68f1e93be4febb0d3cbe_cppui_modular254,
                            0x0a2f53768b8ebf6a86913b0e57c04e011ca408648a4743a87d77adbf0c9c3512_cppui_modular254,
                            0x00248156142fd0373a479f91ff239e960f599ff7e94be69b7f2a290305e1198d_cppui_modular254,
                            0x171d5620b87bfb1328cf8c02ab3f0c9a397196aa6a542c2350eb512a2b2bcda9_cppui_modular254,
                            0x170a4f55536f7dc970087c7c10d6fad760c952172dd54dd99d1045e4ec34a808_cppui_modular254,
                            0x29aba33f799fe66c2ef3134aea04336ecc37e38c1cd211ba482eca17e2dbfae1_cppui_modular254,
                            0x1e9bc179a4fdd758fdd1bb1945088d47e70d114a03f6a0e8b5ba650369e64973_cppui_modular254,
                            0x1dd269799b660fad58f7f4892dfb0b5afeaad869a9c4b44f9c9e1c43bdaf8f09_cppui_modular254,
                            0x22cdbc8b70117ad1401181d02e15459e7ccd426fe869c7c95d1dd2cb0f24af38_cppui_modular254,
                            0x0ef042e454771c533a9f57a55c503fcefd3150f52ed94a7cd5ba93b9c7dacefd_cppui_modular254,
                            0x11609e06ad6c8fe2f287f3036037e8851318e8b08a0359a03b304ffca62e8284_cppui_modular254,
                            0x1166d9e554616dba9e753eea427c17b7fecd58c076dfe42708b08f5b783aa9af_cppui_modular254,
                            0x2de52989431a859593413026354413db177fbf4cd2ac0b56f855a888357ee466_cppui_modular254,
                            0x3006eb4ffc7a85819a6da492f3a8ac1df51aee5b17b8e89d74bf01cf5f71e9ad_cppui_modular254,
                            0x2af41fbb61ba8a80fdcf6fff9e3f6f422993fe8f0a4639f962344c8225145086_cppui_modular254,
                            0x119e684de476155fe5a6b41a8ebc85db8718ab27889e85e781b214bace4827c3_cppui_modular254,
                            0x1835b786e2e8925e188bea59ae363537b51248c23828f047cff784b97b3fd800_cppui_modular254,
                            0x28201a34c594dfa34d794996c6433a20d152bac2a7905c926c40e285ab32eeb6_cppui_modular254,
                            0x083efd7a27d1751094e80fefaf78b000864c82eb571187724a761f88c22cc4e7_cppui_modular254,
                            0x0b6f88a3577199526158e61ceea27be811c16df7774dd8519e079564f61fd13b_cppui_modular254,
                            0x0ec868e6d15e51d9644f66e1d6471a94589511ca00d29e1014390e6ee4254f5b_cppui_modular254,
                            0x2af33e3f866771271ac0c9b3ed2e1142ecd3e74b939cd40d00d937ab84c98591_cppui_modular254,
                            0x0b520211f904b5e7d09b5d961c6ace7734568c547dd6858b364ce5e47951f178_cppui_modular254,
                            0x0b2d722d0919a1aad8db58f10062a92ea0c56ac4270e822cca228620188a1d40_cppui_modular254,
                            0x1f790d4d7f8cf094d980ceb37c2453e957b54a9991ca38bbe0061d1ed6e562d4_cppui_modular254,
                            0x0171eb95dfbf7d1eaea97cd385f780150885c16235a2a6a8da92ceb01e504233_cppui_modular254,
                            0x0c2d0e3b5fd57549329bf6885da66b9b790b40defd2c8650762305381b168873_cppui_modular254,
                            0x1162fb28689c27154e5a8228b4e72b377cbcafa589e283c35d3803054407a18d_cppui_modular254,
                            0x2f1459b65dee441b64ad386a91e8310f282c5a92a89e19921623ef8249711bc0_cppui_modular254,
                            0x1e6ff3216b688c3d996d74367d5cd4c1bc489d46754eb712c243f70d1b53cfbb_cppui_modular254,
                            0x01ca8be73832b8d0681487d27d157802d741a6f36cdc2a0576881f9326478875_cppui_modular254,
                            0x1f7735706ffe9fc586f976d5bdf223dc680286080b10cea00b9b5de315f9650e_cppui_modular254,
                            0x2522b60f4ea3307640a0c2dce041fba921ac10a3d5f096ef4745ca838285f019_cppui_modular254,
                            0x23f0bee001b1029d5255075ddc957f833418cad4f52b6c3f8ce16c235572575b_cppui_modular254,
                            0x2bc1ae8b8ddbb81fcaac2d44555ed5685d142633e9df905f66d9401093082d59_cppui_modular254,
                            0x0f9406b8296564a37304507b8dba3ed162371273a07b1fc98011fcd6ad72205f_cppui_modular254,
                            0x2360a8eb0cc7defa67b72998de90714e17e75b174a52ee4acb126c8cd995f0a8_cppui_modular254,
                            0x15871a5cddead976804c803cbaef255eb4815a5e96df8b006dcbbc2767f88948_cppui_modular254,
                            0x193a56766998ee9e0a8652dd2f3b1da0362f4f54f72379544f957ccdeefb420f_cppui_modular254,
                            0x2a394a43934f86982f9be56ff4fab1703b2e63c8ad334834e4309805e777ae0f_cppui_modular254,
                            0x1859954cfeb8695f3e8b635dcb345192892cd11223443ba7b4166e8876c0d142_cppui_modular254,
                            0x04e1181763050e58013444dbcb99f1902b11bc25d90bbdca408d3819f4fed32b_cppui_modular254,
                            0x0fdb253dee83869d40c335ea64de8c5bb10eb82db08b5e8b1f5e5552bfd05f23_cppui_modular254,
                            0x058cbe8a9a5027bdaa4efb623adead6275f08686f1c08984a9d7c5bae9b4f1c0_cppui_modular254,
                            0x1382edce9971e186497eadb1aeb1f52b23b4b83bef023ab0d15228b4cceca59a_cppui_modular254,
                            0x03464990f045c6ee0819ca51fd11b0be7f61b8eb99f14b77e1e6634601d9e8b5_cppui_modular254,
                            0x23f7bfc8720dc296fff33b41f98ff83c6fcab4605db2eb5aaa5bc137aeb70a58_cppui_modular254,
                            0x0a59a158e3eec2117e6e94e7f0e9decf18c3ffd5e1531a9219636158bbaf62f2_cppui_modular254,
                            0x06ec54c80381c052b58bf23b312ffd3ce2c4eba065420af8f4c23ed0075fd07b_cppui_modular254,
                            0x118872dc832e0eb5476b56648e867ec8b09340f7a7bcb1b4962f0ff9ed1f9d01_cppui_modular254,
                            0x13d69fa127d834165ad5c7cba7ad59ed52e0b0f0e42d7fea95e1906b520921b1_cppui_modular254,
                            0x169a177f63ea681270b1c6877a73d21bde143942fb71dc55fd8a49f19f10c77b_cppui_modular254,
                            0x04ef51591c6ead97ef42f287adce40d93abeb032b922f66ffb7e9a5a7450544d_cppui_modular254,
                            0x256e175a1dc079390ecd7ca703fb2e3b19ec61805d4f03ced5f45ee6dd0f69ec_cppui_modular254,
                            0x30102d28636abd5fe5f2af412ff6004f75cc360d3205dd2da002813d3e2ceeb2_cppui_modular254,
                            0x10998e42dfcd3bbf1c0714bc73eb1bf40443a3fa99bef4a31fd31be182fcc792_cppui_modular254,
                            0x193edd8e9fcf3d7625fa7d24b598a1d89f3362eaf4d582efecad76f879e36860_cppui_modular254,
                            0x18168afd34f2d915d0368ce80b7b3347d1c7a561ce611425f2664d7aa51f0b5d_cppui_modular254,
                            0x29383c01ebd3b6ab0c017656ebe658b6a328ec77bc33626e29e2e95b33ea6111_cppui_modular254,
                            0x10646d2f2603de39a1f4ae5e7771a64a702db6e86fb76ab600bf573f9010c711_cppui_modular254,
                            0x0beb5e07d1b27145f575f1395a55bf132f90c25b40da7b3864d0242dcb1117fb_cppui_modular254,
                            0x16d685252078c133dc0d3ecad62b5c8830f95bb2e54b59abdffbf018d96fa336_cppui_modular254,
                            0x0a6abd1d833938f33c74154e0404b4b40a555bbbec21ddfafd672dd62047f01a_cppui_modular254,
                            0x1a679f5d36eb7b5c8ea12a4c2dedc8feb12dffeec450317270a6f19b34cf1860_cppui_modular254,
                            0x0980fb233bd456c23974d50e0ebfde4726a423eada4e8f6ffbc7592e3f1b93d6_cppui_modular254,
                            0x161b42232e61b84cbf1810af93a38fc0cece3d5628c9282003ebacb5c312c72b_cppui_modular254,
                            0x0ada10a90c7f0520950f7d47a60d5e6a493f09787f1564e5d09203db47de1a0b_cppui_modular254,
                            0x1a730d372310ba82320345a29ac4238ed3f07a8a2b4e121bb50ddb9af407f451_cppui_modular254,
                            0x2c8120f268ef054f817064c369dda7ea908377feaba5c4dffbda10ef58e8c556_cppui_modular254,
                            0x1c7c8824f758753fa57c00789c684217b930e95313bcb73e6e7b8649a4968f70_cppui_modular254,
                            0x2cd9ed31f5f8691c8e39e4077a74faa0f400ad8b491eb3f7b47b27fa3fd1cf77_cppui_modular254,
                            0x23ff4f9d46813457cf60d92f57618399a5e022ac321ca550854ae23918a22eea_cppui_modular254,
                            0x09945a5d147a4f66ceece6405dddd9d0af5a2c5103529407dff1ea58f180426d_cppui_modular254,
                            0x188d9c528025d4c2b67660c6b771b90f7c7da6eaa29d3f268a6dd223ec6fc630_cppui_modular254,
                            0x3050e37996596b7f81f68311431d8734dba7d926d3633595e0c0d8ddf4f0f47f_cppui_modular254,
                            0x15af1169396830a91600ca8102c35c426ceae5461e3f95d89d829518d30afd78_cppui_modular254,
                            0x1da6d09885432ea9a06d9f37f873d985dae933e351466b2904284da3320d8acc_cppui_modular254,
                            0x2796ea90d269af29f5f8acf33921124e4e4fad3dbe658945e546ee411ddaa9cb_cppui_modular254,
                            0x202d7dd1da0f6b4b0325c8b3307742f01e15612ec8e9304a7cb0319e01d32d60_cppui_modular254,
                            0x096d6790d05bb759156a952ba263d672a2d7f9c788f4c831a29dace4c0f8be5f_cppui_modular254,
                            0x054efa1f65b0fce283808965275d877b438da23ce5b13e1963798cb1447d25a4_cppui_modular254,
                            0x1b162f83d917e93edb3308c29802deb9d8aa690113b2e14864ccf6e18e4165f1_cppui_modular254,
                            0x21e5241e12564dd6fd9f1cdd2a0de39eedfefc1466cc568ec5ceb745a0506edc_cppui_modular254,
                            0x1cfb5662e8cf5ac9226a80ee17b36abecb73ab5f87e161927b4349e10e4bdf08_cppui_modular254,
                            0x0f21177e302a771bbae6d8d1ecb373b62c99af346220ac0129c53f666eb24100_cppui_modular254,
                            0x1671522374606992affb0dd7f71b12bec4236aede6290546bcef7e1f515c2320_cppui_modular254,
                            0x0fa3ec5b9488259c2eb4cf24501bfad9be2ec9e42c5cc8ccd419d2a692cad870_cppui_modular254,
                            0x193c0e04e0bd298357cb266c1506080ed36edce85c648cc085e8c57b1ab54bba_cppui_modular254,
                            0x102adf8ef74735a27e9128306dcbc3c99f6f7291cd406578ce14ea2adaba68f8_cppui_modular254,
                            0x0fe0af7858e49859e2a54d6f1ad945b1316aa24bfbdd23ae40a6d0cb70c3eab1_cppui_modular254,
                            0x216f6717bbc7dedb08536a2220843f4e2da5f1daa9ebdefde8a5ea7344798d22_cppui_modular254,
                            0x1da55cc900f0d21f4a3e694391918a1b3c23b2ac773c6b3ef88e2e4228325161_cppui_modular254
                    };

                    typedef algebra::matrix<element_type, state_words, state_words> mds_matrix_type;
                    constexpr static const mds_matrix_type mds_matrix = {
                            0x109b7f411ba0e4c9b2b70caf5c36a7b194be7c11ad24378bfedb68592ba8118b_cppui_modular254,
                            0x16ed41e13bb9c0c66ae119424fddbcbc9314dc9fdbdeea55d6c64543dc4903e0_cppui_modular254,
                            0x2b90bba00fca0589f617e7dcbfe82e0df706ab640ceb247b791a93b74e36736d_cppui_modular254,
                            0x2969f27eed31a480b9c36c764379dbca2cc8fdd1415c3dded62940bcde0bd771_cppui_modular254,
                            0x2e2419f9ec02ec394c9871c832963dc1b89d743c8c7b964029b2311687b1fe23_cppui_modular254,
                            0x101071f0032379b697315876690f053d148d4e109f5fb065c8aacc55a0f89bfa_cppui_modular254,
                            0x143021ec686a3f330d5f9e654638065ce6cd79e28c5b3753326244ee65a1b1a7_cppui_modular254,
                            0x176cc029695ad02582a70eff08a6fd99d057e12e58e7d7b6b16cdfabc8ee2911_cppui_modular254,
                            0x19a3fc0a56702bf417ba7fee3802593fa644470307043f7773279cd71d25d5e0_cppui_modular254
                    };

                    typedef std::pair<mds_matrix_type, round_constants_type> constants_type;
                    constexpr static const constants_type constants = {mds_matrix, round_constants};
                };

                template<>
                class poseidon_constants<poseidon_policy<
                        algebra::fields::alt_bn128_scalar_field<254>, 128, 4>> {
                public:
                    typedef poseidon_policy<algebra::fields::alt_bn128_scalar_field<254>, 128, 4> policy_type;
                    typedef typename policy_type::field_type field_type;

                    constexpr static const std::size_t state_words = policy_type::state_words;
                    constexpr static const std::size_t word_bits = policy_type::word_bits;
                    constexpr static const std::size_t full_rounds = policy_type::full_rounds;
                    constexpr static const std::size_t half_full_rounds = policy_type::half_full_rounds;
                    constexpr static const std::size_t part_rounds = policy_type::part_rounds;
                    constexpr static const std::size_t rate = policy_type::block_words;

                    typedef typename field_type::value_type element_type;
                    typedef typename field_type::integral_type integral_type;
                    constexpr static const integral_type modulus = field_type::modulus;

                    constexpr static const std::size_t lfsr_state_bits = 80;
                    typedef number<backends::cpp_int_modular_backend<lfsr_state_bits>> lfsr_state_type;

                    typedef algebra::vector<element_type, state_words> state_vector_type;

                    inline constexpr static const element_type &round_constant(std::size_t round, std::size_t i) {
                        return round_constants[round][i];
                    }

                    inline constexpr static void product_with_mds_matrix(state_vector_type &A_vector) {
                        A_vector = algebra::vectmatmul(A_vector, mds_matrix);
                    }

                    typedef algebra::matrix<element_type, full_rounds + part_rounds, state_words> round_constants_type;
                    constexpr static const round_constants_type round_constants = {
                            0x0eb544fee2815dda7f53e29ccac98ed7d889bb4ebd47c3864f3c2bd81a6da891_cppui_modular254,
                            0x0554d736315b8662f02fdba7dd737fbca197aeb12ea64713ba733f28475128cb_cppui_modular254,
                            0x2f83b9df259b2b68bcd748056307c37754907df0c0fb0035f5087c58d5e8c2d4_cppui_modular254,
                            0x2ca70e2e8d7f39a12447ac83052451b461f15f8b41a75ef31915208f5aba9683_cppui_modular254,
                            0x1cb5f9319be6a45e91b04d7222271c94994196f12ed22c5d4ec719cb83ecfea9_cppui_modular254,
                            0x2eb4f99c69f966ebf8a42192de7ff61621c7bb47b93750c2b9ea08d18446c122_cppui_modular254,
                            0x224a28e5a35385a7c5198169e405d9ea0fc7da8b93ee13b6d5f7d099e299520e_cppui_modular254,
                            0x0f7411b465e600eed8afdd6afca49c3036f33ecbd9a0f97823796b993bbd82f7_cppui_modular254,
                            0x0f9d0d5aad2c9555a2be7150392d8d9819b208ae3370f99a0626f9ff5d90e4e3_cppui_modular254,
                            0x1e9a96dc8292bb596f52a59538d329229732b25259cf744b6a12d30702d6fba0_cppui_modular254,
                            0x08780514ccd90380887d578c45555e593cfe52eab4b945c6c2cd4d528fb3fe3c_cppui_modular254,
                            0x272498fced686c7ac8149fa3f73ef8c2ced64717e3556d5a59f119d629ccb5fc_cppui_modular254,
                            0x01ef8f9dd7c93aac4b7cb80930bd06eb45bd350aff585f10e3d0ef8a782ef7df_cppui_modular254,
                            0x045b9f59b6595e614dc08f222b469b138e886e64bf3c40aa97ea0ae754934d30_cppui_modular254,
                            0x0ac1e91c57d9da919fd6f59d2a40ff8ea3e41e24e247a387adf2584295d61c66_cppui_modular254,
                            0x028a1621a94054b0c7f9a421353cd89d0fd67061aee99979d12e68f04e62d134_cppui_modular254,
                            0x26b41802c071ea4c9632647ed059236e50c19c3fb3c96d09d02aae2a0dcd9dbc_cppui_modular254,
                            0x2fb5dda8072bb72cbaac2f63e468215e05c9de06758db6a94af34384aedb462b_cppui_modular254,
                            0x2212d3a0f5fccaf244ff3547fd823249ad8ab8ba2a18d383dd05c56ee894d850_cppui_modular254,
                            0x1b041ad5b2f0684258e4dfaeea09be56a3276fdb19f44c015cd0c7eed465e2e3_cppui_modular254,
                            0x0a01776bb22f4b6b8eccff33e76fded3144fb7e3ac14e846a91e64afb1500eff_cppui_modular254,
                            0x2b7b5674aaecc3cbf34d3f275066d549a4f33ae8c15cf827f7936440810ace43_cppui_modular254,
                            0x29d299b80cd4489e4cf75779ed54b48c60b042257b78fc004c1b803381a3bdfd_cppui_modular254,
                            0x1c46831d9a74529357641c219d721a74a427110032b5e1dd19dde30424be401e_cppui_modular254,
                            0x06d7626c953ccb72f37141dc34d578e036296c0657674f80739ae1d883e91269_cppui_modular254,
                            0x28ffddc86f18c136c54002748e0c410edc5c440a3022cd960f108c71cda2930c_cppui_modular254,
                            0x2e67f7ee5e4aa295f85deed09e400b17be67f1b7ed2ab6adb8ec0619f6fbc5e9_cppui_modular254,
                            0x26ce38fa636c90630e97f25114a79a2dca56859ef759e53ce7abf22c24e80f27_cppui_modular254,
                            0x2e6e07c3c95bf7c34dd7a01d00a7ffec42cb3d16a1f72721afacb4c4cfd35db1_cppui_modular254,
                            0x2aa74f7597f0c9f45f91d7961c3a54fb8890d276612e1246384b1470da24d8cc_cppui_modular254,
                            0x287d681a46a2faae2c7c090f668ab45b8a71313c1509183e2ec0ca639b7f73fe_cppui_modular254,
                            0x212bd19df812eaaef4a40600528f3d7da5d3106ff565aa3b11e29f3305e73c04_cppui_modular254,
                            0x1154f7cf519186bf1aafb14b350eb860f97fd9740926dab93809c28404713504_cppui_modular254,
                            0x1dff6385cb31f1c24637810a4bd1b16fbf5152905be36583da747e79661fc207_cppui_modular254,
                            0x0e444582d22b4e76c081d34c44c18e424011a34d5476252863ea3c606b551e5c_cppui_modular254,
                            0x0323c9e433ba66c4abab6638328f02f1815773e9c2846323ff72d3aab7e4eff8_cppui_modular254,
                            0x12746bbd71791059193bba79cdec448f25b8cf002740112db70f2c6876a9c29d_cppui_modular254,
                            0x1173b7d112c2a798fd9b9d3751842c75d466c837cf50d73efd049eb4438a2240_cppui_modular254,
                            0x13d51c1090a1ad4876d1e555d7fed13da8e5713b25026ebe5fdb4808703243da_cppui_modular254,
                            0x00874c1344a4ad51ff8dcb7cbd2d9743cb72743f0394efe7f4a58ebeb956baa1_cppui_modular254,
                            0x22df22131aaab85865ce236b07f244fa0eea48d3546e97d6a32a562074fef08f_cppui_modular254,
                            0x0bf964d2dbd25b908708b437a445fc3e984524a59101e6c18bf5eb05a919f155_cppui_modular254,
                            0x09b18d9b917a55bca302be1f7f181e0e640b9d73a9ab298c69b435b5fc502f32_cppui_modular254,
                            0x094f5534444fae36a4bfc1d5bf3dc05bfbbbc70a6365366dd6745a5067289e43_cppui_modular254,
                            0x2999bab1a5f25210519fa6622af53a15a3e240c0da5701cb784fddc0dc23f01f_cppui_modular254,
                            0x2f6898c07581f6371ca94db73710e88084301bce8a93d13669575a11b03a3d23_cppui_modular254,
                            0x07268eaaba08bc19ec16d7e1318a4740565deb1e8e5742f862174b1a6866fccb_cppui_modular254,
                            0x186279b003454db01339ff77113bc9eb62603e078e1c6689a6c9582c41a0529f_cppui_modular254,
                            0x18a3f736509197d6e4915bdd04d3e5ddb67e2cc5de9a22750768e5524737172c_cppui_modular254,
                            0x0a21fa1988cf38d877cc1e2ed24c808c725e2d4bcb2d3a007b5987b87085671d_cppui_modular254,
                            0x15b285cbe26c467f1faf5ef6a64625228328c184a2c43bc00b36a135e785fba2_cppui_modular254,
                            0x164b7062c4671cf08c08b8c3f9806d560b7775b7c902f5788cd28de3e779f161_cppui_modular254,
                            0x0890ba0819ac0a6f86d9865fe7e50ef361c61d3d43b6e65d7a24f651249baa70_cppui_modular254,
                            0x2fbea4d65d7ed425a42712e5a721e4eaa627ac5cb0eb878ccc2ee0aed543e922_cppui_modular254,
                            0x0492bf383c36fa55540303a3b536f85e7b70a58e854ab9b9103d7f5f379abaaa_cppui_modular254,
                            0x05e91fe944e944104e20251c565142d61d6185a9ce85675f6a969d56292dc24e_cppui_modular254,
                            0x12fe5c2029e4b33893d463cb041acad0995b9621e6e49c3b7e380a76e36e6c1c_cppui_modular254,
                            0x024154adf0255d47958f7723921474131f2629fadc89496906cd01dc6fa0784e_cppui_modular254,
                            0x18824a09e6afaf4a36ed2462a86bd0bad798815644f2bbde8813c13457a45550_cppui_modular254,
                            0x0c8b482dba0ad51be9f255de0c3dbddddf84a630af68d50bbb06983e3d5d58a5_cppui_modular254,
                            0x17325fd0ab635871363e0a1667d3b67c5a4fa67fcd6aaf86441392878fdb05e6_cppui_modular254,
                            0x050ae95f6d2f1519122f5af67b690f31e550773fa8d18bf71cc6d0e911fa402e_cppui_modular254,
                            0x0f0d139a0e81e943038cb288d62636764bbb6295f07569885771ec84edc50c40_cppui_modular254,
                            0x1c0f8697795689cdf70fd2f2c0f93d1a79b39ebc7a1b1c549dbbca7b8e747cd6_cppui_modular254,
                            0x2bd0f940ad936b796d2bc2e048bc979e49be23a4b13598f9fe536a16dc1d81e6_cppui_modular254,
                            0x27eb1be27c9c4e934778c09a0053337fa06ebb275e096d167ce54d1e96ee62cb_cppui_modular254,
                            0x2e4889d830a67e5a8f96bdd3155a7ca3284fbd307d1f71b0f151be62548e2aea_cppui_modular254,
                            0x193fe3db0ab47d3c5d2ec5e9c5bd9983c9891f2cadc165db6064bbe6fcc1e305_cppui_modular254,
                            0x2bf3086e96c36c7bce415907ad0c40ed6e9661c009679e4e37cb13027c83e525_cppui_modular254,
                            0x12f16e2de6d4ad46a98cdb697c6cad5dd5e7e413f741ccf29ff2ea486e59bb28_cppui_modular254,
                            0x2a72147d230119f3a0262e3653ddd19f33f3d5d6ec6c4bf0ad919b0343b92d2f_cppui_modular254,
                            0x21be0e2c4bfd64e56dc47f957806dc5f0a2d9bcc26412e2977df79acc10ba974_cppui_modular254,
                            0x0e2d7e1dc946d70b2749a3b54367b25a71b84fb911aa57ae137fd4b6c21b444a_cppui_modular254,
                            0x2667f7fb5a4fa1246170a745d8a4188cc31adb0eae3325dc9f3f07d4b92b3e2e_cppui_modular254,
                            0x2ccc6f431fb7400730a783b66064697a1550c12b08dfeb72830e107da78e3405_cppui_modular254,
                            0x08888a94fc5a2ca34f0201462420001fae6dbee9e8ca0c242ec50621e38e6e5d_cppui_modular254,
                            0x02977b34eeaa3cb6ad40dd42c9b6fdd7a0d2fbe753af88b36acfcd3ccbc53f2a_cppui_modular254,
                            0x120ccce13d28b75cfd6fb6c9ea13a648bfcfe0d7e6ff8e9610b5e9f971e16b9a_cppui_modular254,
                            0x09fad2269c4a8e93c81e1b9770ea098c92787a4575b2bd73a0bf2af32f86ff3c_cppui_modular254,
                            0x026091fd3d4c44d50a4b310e4ac6f0fa0debdb70775eeb8af630cffb60092d6f_cppui_modular254,
                            0x29404aa2ba565b77bb7fba9dfb6fc3212543cc56afad6afcb904fd2bca893994_cppui_modular254,
                            0x2749475c399aaf39d4e87c2548695b4ef1ffd86590e0827de7201351b7c883f9_cppui_modular254,
                            0x098c842322479f7239912b50424685cba2ebe2dc2e4da70ac7557dab65ffa222_cppui_modular254,
                            0x18cef581222b647e31238e57fead7d5c758ace14c93c4da40191d0c053b51936_cppui_modular254,
                            0x13177839c68a5080d4e746745e43711d3cbc0ca4a108f98d63b2aa681698de60_cppui_modular254,
                            0x020ca696f531e43ec088f56f4b74325626cc4df712c0e5f0a907d88e5f0deffd_cppui_modular254,
                            0x27230eede9cccfc9fa805a30fc548db693d13708c646841d16e028387c7ac022_cppui_modular254,
                            0x01645911c1198b01d64fde34a342a1786497c05969a015439057d2fe75bb281c_cppui_modular254,
                            0x2c323fe16481bf496e439c88341ce25f198971e14487056cfdca4a451a5d8643_cppui_modular254,
                            0x0fc082dfe70728e8450bd2074c3e22e1b022c124d3bffe8b5af88ae6db5085c8_cppui_modular254,
                            0x2052c174800db209d8cdca568dcc25b3be9642116ac4c77efe8a488b423521ee_cppui_modular254,
                            0x28e420e10df2fbb5af96d621d55423190be351ce8129065a8dd9fd05b3ece9c0_cppui_modular254,
                            0x25698ca5e24a1b799f783c4462a24db655d6ae1bdacd1cb549d6e0bc3ae5069a_cppui_modular254,
                            0x160a9981a5c89a57cf8ffbfa57d51049a297b61074422ac134d9b857d6984d35_cppui_modular254,
                            0x21c91a39e145c3bc34d9b694b843f3bf8b7cebf59ddbb0a064642b069997f3d4_cppui_modular254,
                            0x1ac8d80dcd5ee876d2b09345ef112345d6eaa029d93f03b6d10975461e41734c_cppui_modular254,
                            0x0ab3e6ad0ecf8b8e7c1662a4174c52225d822895e2755544b8dbcea5657ce02c_cppui_modular254,
                            0x1c675182512620ae27e3b0b917b3a21ca52ef3ef5909b4e1c5b2237cbdab3377_cppui_modular254,
                            0x2cdbc998dfd7affd3d948d0c85bad2e2e37a4a3e07a7d75d0c8a9092ac2bed45_cppui_modular254,
                            0x23b584a56e2117b0774bf67cc0dee33324337350309dff833e491a133bb63b2e_cppui_modular254,
                            0x1e9e2b310f60ba9f8cb73030a3c9d2a10d133bc6ba4ec1152f3d20de1465e9a5_cppui_modular254,
                            0x0e01e365ba5b3031abc3e720140ae746c9ab5dab987520c460bcd4f1fa5b22db_cppui_modular254,
                            0x040884cdcfc64bfc7b7127340498d5c443382011b61c9a4b1387d85bc1264e68_cppui_modular254,
                            0x190b1ee1205eb9500c74a3998f2bea36353f1724d6067ed0a0a17de311ef9668_cppui_modular254,
                            0x1647c72aec6c4388d04f52fc23cd9c08c1dfcf65ce61e165fc28d1f832bd3b2c_cppui_modular254,
                            0x2430006346a0145f799880cc4c8736269f5494d89fb48b02842e595b71e4541d_cppui_modular254,
                            0x177b9a08343917e1365107a3da3ae7f69d853902bb16bacb3221850252b757af_cppui_modular254,
                            0x04a420e642b11ae94e58862a68f5e32609cd53d0ae29423439b11d04666df4f8_cppui_modular254,
                            0x25d0e0f739fb39fc105a88fab0afd810de2461858e956ccccdfabeddb6a25c8f_cppui_modular254,
                            0x04476d91b7eff2fd85905cbf58651edc320cb15610eaed452c4d4ffa0c740a27_cppui_modular254,
                            0x1090c0b68b3d7d7b8bc9ca2419eb8dea1c28f6d5e1250cb5e9780fd9ca286fae_cppui_modular254,
                            0x25393ce3b9256d50448a725c5c7cd5ad376f2d435855c10ebf2899cb5c6617be_cppui_modular254,
                            0x25931c0c7371f4f1fc862f306e6e5830ed824388d6b9342697d144f0fab46630_cppui_modular254,
                            0x2396cb501700bbe6c82aad51b0fb79cf8a4d353185d5808203f73f22afbf62f6_cppui_modular254,
                            0x26a363483348b58954ea748a7129a7b0a3dc9068c3cca7b5b3f0ce03b8724884_cppui_modular254,
                            0x27ca107ca204f2a18d6f1535b92c5478c99b893334215f6ba7a0e5b45fcd6897_cppui_modular254,
                            0x26da28fc097ed77ce4662bde326b2cceac15f7301178581d8d2d02b3b2d91056_cppui_modular254,
                            0x056ab351691d8bb3703e3055070ac9cc655774c1bb35d57572971ba56ee0cb89_cppui_modular254,
                            0x2638b57f23b754aec76d109a2f481aa3c22547a11ffc50152d729af632376a90_cppui_modular254,
                            0x304754bb8c57d60732f492c2605184fdc33e46a532bdec80ea7bc5519ede7cef_cppui_modular254,
                            0x00d1727f8457ee03514f155b5806cbf748ec6857fc554010752ac93a9b7619ac_cppui_modular254,
                            0x00ee1f3c66fbc05c43ba295a303c72fab5bca86805ec9419c588e50947761fa3_cppui_modular254,
                            0x0afafadcf5b4dd4a4a76b5a1d82415fd10a19fbcfc59078c61f9297eb675d972_cppui_modular254,
                            0x0b2449f39746085e86ce45e8eed108ee65a234835a0a6a5ea8996d124dd04d0a_cppui_modular254,
                            0x206b0ce2f1b2c5b7c9f37b0045227095f6c6f071ec3bdda76a7ddf4823dd5dd6_cppui_modular254,
                            0x0feba4fb87834c7cb696e67433628cd6caffc3a4ef20fea852c7e1029459409c_cppui_modular254,
                            0x254dbfac74c49b0b8926752e084e02513b06f1315e6d70e18173e972336e55d3_cppui_modular254,
                            0x0addb1372cee4e164655168c367559e19606c5bd17910aeb37719edfa0ca8762_cppui_modular254,
                            0x26b25b7e257f3e97c799024fb019f65c6ca4d8d81b1ae16221a589d68831d759_cppui_modular254,
                            0x090995b79acec240413b8d4c658787e5a4657b9ab00bdb5b1960b1059e113ba3_cppui_modular254,
                            0x08dbdc2e21ef11f2c57299687843cea3eb0d8e40e99131f42974178d44f73b7b_cppui_modular254,
                            0x09e8aba671481197679faf752a0f78e342fe9c491596ab6758f170939785179f_cppui_modular254,
                            0x1deb05180e833e45659052a7ebaf816c7efd12a7f9eec94b7bc7c683f1363d5c_cppui_modular254,
                            0x19a70ec6bdfc9098a926efbcc04aa9ee248997e8b2c24af335fd6523e5250879_cppui_modular254,
                            0x21d773660adafb8a879986f9aab4890566353a3777d8a3f1eb93abe10bbf1f64_cppui_modular254,
                            0x09f1890f72e9dc713e20ba637b89d5d397a6b01fcd667347f6f46617841c3901_cppui_modular254,
                            0x05af459361eb454d2a300c61e446998d48fa1f897bf219d608c2145c33b111c3_cppui_modular254,
                            0x0fa1a1d6829f0345664a66dc75a657335f336f15f340756cfa12fc850cc8b513_cppui_modular254,
                            0x02e47a35bcc0c3a0bda0b1c0307ad543f4280fcf87f636f853655cf97a628bb0_cppui_modular254,
                            0x14f773e9834c6bdeb8f90e78bf4c24b7203411460112491036621895204d0f12_cppui_modular254,
                            0x102d98cf502ed843255cf19d29bc7d8e642abe7cfd639992ffb091962fc8f7cc_cppui_modular254,
                            0x043dd5f4aa5a76dd4c47f6c65da7ca2320d4c73ad3294738cba686a7e91373c2_cppui_modular254,
                            0x21833819c3337194a6c0d29a48d4f2676f0e7c79743a306f4cfdb2b26bd11efa_cppui_modular254,
                            0x0f281925cf5ee649b474a6819d116ca3eb4eca246c311ecadc53262a3cff2b53_cppui_modular254,
                            0x0d3e2477a7b10beb44709c7746d6824edf625dd60504d5dc93ce662f15c238d6_cppui_modular254,
                            0x2cd7f641bedbf66956ff8a01be9cde35d80f80ab51e73b49acbfc3eff5aefc44_cppui_modular254,
                            0x29e95b492bf2f95f4d09380f98b74e389149d24045811d7a86dd861310463cf8_cppui_modular254,
                            0x22da66bc62e8f011266efca86a6c810f9ae4c51af6ffeb57f8b3c50df83cc13e_cppui_modular254,
                            0x0fe6d30de7a82d163023491794f4aca3220db79e8129df3643072d841925554a_cppui_modular254,
                            0x0050e842a1299909123c46eff185c23ad312d03fef1adfecc7e07ecb298fd67f_cppui_modular254,
                            0x2130a3a7b3221222be34cc53a42d7733666f9ddf714ed7c5885cbbdb63108c21_cppui_modular254,
                            0x2df9ee294edf99e3d8d5883fe0566c24aa66731f34a93280e1d328e67b33c9fa_cppui_modular254,
                            0x1bf7d6e489ad8c0cf26eb68cc21ff54158132396dc250aeba4b6fc5fc3372762_cppui_modular254,
                            0x0c602fa155be958761eaf739617ab136cf7b807728bf7fe35d4778d311780e54_cppui_modular254,
                            0x2e50e2c5b36aa20532407d86b8d22d7d5154080a24972faeb63faf0121ed7f21_cppui_modular254,
                            0x17c2510982a7b5825710d6290ec4f782f674995ee8409b42b459123b180332e1_cppui_modular254,
                            0x0b0d52f03c8af7276803ecf2465b885b21337b538eabd2f6b2ab255f376b42a8_cppui_modular254,
                            0x0f5633df1972b9455953d88a63f80647a9ac77c6c0f85d4561972dd8fab8bd14_cppui_modular254,
                            0x0ebf7ad29ca13804e1422e939681155124780ff43e76e929035498130a7f1572_cppui_modular254,
                            0x1aff13c81bda47e80b02962173bba343e18f94bee27c8a57661b1103a720ffe2_cppui_modular254,
                            0x210449dbf5cf3061da2465be85505862d3f31de1a3b58ff35713be57efac6c07_cppui_modular254,
                            0x088230c2794e50c57d75cd6d3c7b9dbe19d1e2f1d3001044b93ad1c3ee629817_cppui_modular254,
                            0x1c408c256490b0a1da08dc464138dfc78cce9a9e16c7705617a4d6dbb20e7e3a_cppui_modular254,
                            0x074517e081eb4c1f22d1771200fb07658f7c77654d58440490dd6f557e9e3903_cppui_modular254,
                            0x02d04e9c21df1dbd88524bdb203691b4cee5530559d6cf0fa05adf61e12fdcbf_cppui_modular254,
                            0x2eb7a011b8bce91082e13ebd75de3b58eb9b4650dae9f11aa81db32cf1b67b13_cppui_modular254,
                            0x2efda77ed35f4af0299f75d6e8a849b54d2ac6bf95368304e6030c18f0cf17b5_cppui_modular254,
                            0x09199dcafd50ce642eddbeda65206d4f61a73d10852b8114c51b2440192ae064_cppui_modular254,
                            0x268c5cfc446d399c4dd319db666a75b5cb655d8c1797e9fa76181cb4216e1562_cppui_modular254,
                            0x2303a652c949071826b0e9a36c80578697b44e912cce6687012854eda11a18dc_cppui_modular254,
                            0x27c53563b12a6ee2c3f041f31dc45922bc5353eb110868d237073f4efb35fbdf_cppui_modular254,
                            0x1201a87eaf4ae618f02bd82d0a5109049969b5248cfe90f42c278f22615d2b0e_cppui_modular254,
                            0x2c43169439fcd69ead8214997bb069becafcb1ba2c51e5706cb4b43dab2a443d_cppui_modular254,
                            0x0683597315359040ea03c45d6984c6894f46cbb36d702e3c4fb9847e6304d944_cppui_modular254,
                            0x03545706706eab36afb93b128febd16fb0425e158314197b77795ad3a798d183_cppui_modular254,
                            0x1a33c254ec117619d35f1fc051b31728740bed23a6a37870edb393b71a0c0e6b_cppui_modular254,
                            0x1ffe6968a4470cd567b0c002281caf996e88f71e759b87e6f338e517f1690c78_cppui_modular254,
                            0x0fd66e03ba8808ffecb059c899fd80f4140ddd5d2a5c4483107f4e02e355b393_cppui_modular254,
                            0x263ab69f13b966f8197394552906b17e6c8617a7bdd5d74a7be3396b7fe013ab_cppui_modular254,
                            0x16a425e47d1110625054d5a165de413e3bd87d5aa3958fdd6eb7e03e39ba4046_cppui_modular254,
                            0x2dc510a4719ec10cad752f03c673f0e253cc31d13e39e909fcc5f73af9138d9a_cppui_modular254,
                            0x24df8e8d856c5b5e1bd1cad23d07dda3423c5179329b7a82cb4aa709a94576e5_cppui_modular254,
                            0x2bcc94ff4fc3c76f3cd5c68915a042e87628249a01b09561bdf24a6cdce5620f_cppui_modular254,
                            0x076c1e88dc540c8d8de54e343df7c429d3295f52c38cffe6b48be86852da97df_cppui_modular254,
                            0x09b5f209a451ac431c051fb12d9a5e4fe40ee1601120947da990fb8e12cb46e1_cppui_modular254,
                            0x205f17b0d8729e2eaa88d6a44135a6ab64e9424f55b0f1ea0683af75eb677c07_cppui_modular254,
                            0x281c5c688836f6cf912638c38be046cd091681f0a41761720cdd1edf9f237029_cppui_modular254,
                            0x1a053e6878e900f45f4d67448c471cf3009a44e7a02ea50e4afa44f2592621f5_cppui_modular254,
                            0x100dc7d426debe3007fb7ceac84e4f5468efcb897e7bbee981742839d59e064c_cppui_modular254,
                            0x17022672a016a957bb87e2cfadc8b75fb28905bdb62c82c80b1cb31b411e49c8_cppui_modular254,
                            0x1086db7e2760fc8b71053a87ebe151239fb8b547182b170de0c27203f954f4d2_cppui_modular254,
                            0x15384fe39d73b63302460ae4c2942fac2b41fb65a185536fb85dd24fd7584064_cppui_modular254,
                            0x2ebb599fe9136d424bf4abc5342c6c7447b1a853205fcfb5519e551357709008_cppui_modular254,
                            0x1b4b5e87cfb9262cfec3c0f0542e4c5a4cf278292b4ce3eed996fac6f4d37288_cppui_modular254,
                            0x2465053ae50b6885801f3f82e302cafbbb4a7581bb4fba60b637febe659e5057_cppui_modular254,
                            0x114f32edcdea09cd095c5bb5d38f1b97da9f05e18b3708bf6e0ab9d3d54859ef_cppui_modular254,
                            0x2bc70dfeb2baab2f6b387cd77be779ac2e5e5519f3d18123ee28d8c2543c7148_cppui_modular254,
                            0x01c9bf7a203ce22b775e3a61ad7e77b6a78348b9f6ec68a412e49bfe32c05415_cppui_modular254,
                            0x0514b0fe5909ea887bedb0295fbbcec355cfb575ff6a97cd9f4ad00ccb57ee9b_cppui_modular254,
                            0x267c76ec81934cc81a132a8b058910a12092520b12a201af03e3202d7b6c1b7e_cppui_modular254,
                            0x29170e3322b3d8d5c78c84babbb470adf1622493ce83e95cfb151cf757bde5d6_cppui_modular254,
                            0x019f6a8124b19e33af33e5d3873f9c335c6f09a45486cab536dd596ca41d9519_cppui_modular254,
                            0x1904aa4d6908544a8b348e9db1981c27009ed8ea171518ae5405d036242b60e9_cppui_modular254,
                            0x26f17873949bc679f7f043956694e422b3cee1de9dd6f6473b932a476455ff1a_cppui_modular254,
                            0x1ac668f612b8243c193b33720b8aa54040c476031197131ebdcac9b18bc48f75_cppui_modular254,
                            0x0996d961a75c0d07196dae45bf624766ccfbf8555be9796da52f81568ef0663d_cppui_modular254,
                            0x030c97e1b8cad1d4fd50d1b4383fbe6674d171f99c63febb5425b395c24fc819_cppui_modular254,
                            0x06e3ad6a46900e2d3953370255b68f89b3e523f1fe502642ee226f2d8bd0848f_cppui_modular254,
                            0x1d6b3755331cd0216b6880e42f9880f565cb94b0e0455153a329890588cc916e_cppui_modular254,
                            0x28e4dcba4b96f12a59b041535e730ac8c35189dc0b85ac033dd38c08bae531f2_cppui_modular254,
                            0x08b6086046a835508ccf484f2974b6a6b0712a476260376c7a3b3e4bc4a47a14_cppui_modular254,
                            0x162cd2ca7fe3b5f1444bcec97812019bb6fd85fba6a0536a89643e15b9bb3b52_cppui_modular254,
                            0x28f1e03baaea9bbc05af5b11937e4f5cb5c9a9c1192063d1998c01c64d483a76_cppui_modular254,
                            0x1bdb062778d7c15da395af2734c25faa0127d2aab4aa71366031a0bb6791ce10_cppui_modular254,
                            0x2375839502e09890cb2914e829627e0e0fc98870b2324a8b50329ebdd24749cb_cppui_modular254,
                            0x1fa8662fbcb61fb3ad7c55668dc9423a332dc87cfb2df456e92d33611ed7bb50_cppui_modular254,
                            0x1e4fad2dd6b0a6f1f8707f721716c8a446e2fb2c47a5138f3f7f9736079d7694_cppui_modular254,
                            0x211256d16c7269fd6df6f5fcdd1fa788ba3bd050059f53d261b0f5f13731ffe7_cppui_modular254,
                            0x2e49084b336eceaa4f8e2a2e6af08318f42060e574dda341f4a1079b12bcc5a5_cppui_modular254,
                            0x0ce19f54cdc39f7f3bf35192ac6808211aecea08dfe14cab758d25891fb00bb9_cppui_modular254,
                            0x0011c5d56c390e893cc394221261d8748dc60451e4ae4e1c84a8468bab2c14cb_cppui_modular254,
                            0x17d79ff06b63ac2a8a9e05ee6af3dbb7ca60e17bfa39b47514a8cd8051579b4c_cppui_modular254,
                            0x19a7d3a446cb5393dc74560093592b06b1a8b35cd6416a2ecab00173639015fa_cppui_modular254,
                            0x030c00a0933dcdba2a808b2e1b9282f331f04596d8928da7aa6c3c97237037a6_cppui_modular254,
                            0x16bcb447ce2d50f3ae25ad080695382e935d2d00184c4acc9370be8aab64139c_cppui_modular254,
                            0x12341b46b0150aa25ea4ec8715312997e62124f37cab7b6d39255b7cd66feb1d_cppui_modular254,
                            0x0e86d13917f44050b72a97b2bf610c84002fc28e296d1044dc89212db6a49ff4_cppui_modular254,
                            0x08e6eb4089d37d66d357e00b53d7f30d1052a181f8f2eb14d059025b110c7262_cppui_modular254,
                            0x2ea123856245f6c84738d15dd1481a0c0415ccb351a1e0cee10c48ce97ca7b18_cppui_modular254,
                            0x2dca72b2ebcab8c23446e00330b163104195789025413abf664db0f9c84dfa6f_cppui_modular254,
                            0x06ff9ed50d327e8463329f585ec924b3f2f6b4235f036fa4c64a26cbd42b6a6b_cppui_modular254,
                            0x246a10b7e3e0089947f7c9bda3d54df8e2a60e0cca84ea2ac630a4535afbf730_cppui_modular254,
                            0x22a63501c5f04b9018719ed99d700ee52f846a715ae67ad75c96b39d688b6691_cppui_modular254,
                            0x2f4c50477f7fd9c671799ac5d2e224cdb9164f58351d8aa140ec07e514fae937_cppui_modular254,
                            0x10ffb7aad1f51c7d13b17f4d876d9a1e38f0ba8a4a23d4b50cda32cad851567e_cppui_modular254,
                            0x0e9cefddc3c2d3bea4d39722532d5420784027352187e7af1a056935c35803ae_cppui_modular254,
                            0x07af84a4d3141e7ac23352e6dc6ea4afa1656f96a33c8978a3e83bdd4ba62b41_cppui_modular254,
                            0x2d9e31a10aebc761f8de00d14b1e566d1a39323d6e89b638e940f3ec8a22c3c5_cppui_modular254,
                            0x27f19a6532e66b5333db1afd592f66f1d36034b314dad8447656747be27e64c7_cppui_modular254,
                            0x0058fa3c8454d63354b2024c3b4a577a180ed99f8f3155cd7e4d617d47d07ffd_cppui_modular254,
                            0x041627b6715b780967957c080699343eb0414a205d3a175d708964956816a5d5_cppui_modular254,
                            0x006ac49dd9253edc7f632e57b958ccecd98201471cf1f66589888f12b727c52d_cppui_modular254,
                            0x0131adffd8bd7254b1d8c3616bbe3386ec0c9c0d6d25a9a4ec46a6bf18301398_cppui_modular254,
                            0x1c4a6f52c9fccf7a4138e413ef62a28377977ad7e25e49a3cf030e1cd8f9f5b6_cppui_modular254,
                            0x03f2a6be51ec677f946551b3860ea479fee048ae2078aeb7d1f7958d2c2645f6_cppui_modular254,
                            0x2da770aad2c2eb09391a0cb78ef3a9648a1372d8543119564d7376396b8ddc62_cppui_modular254,
                            0x15278463665f74cddc1802febfab02cec9d45fe866c359c738062afb75d64a03_cppui_modular254,
                            0x12fe278aa36544eac9731027090518d434e38ea966a08a6f8d580638ac54c773_cppui_modular254,
                            0x149b9c802182558a4c45d119d3f4cc7fd8587604ca4f0d6e21b06ff30b6a23b6_cppui_modular254,
                            0x0812e7b4d847bc8517d19319772f3c9855e044fd60dbac9a0adc4959b691dfe4_cppui_modular254,
                            0x02ed8d8ddeafe3d9d8df7f28a0bfaa7f555813c7e7503aea2a66973703a0c61b_cppui_modular254,
                            0x0ebd073ba0537b514deb6029f921029e55e5e4d9a03d6b6ba1304038662d4db8_cppui_modular254,
                            0x15c754d5b14b2c4205c6ba8d2ccd028255b3e792c6afa08b44ee75b62eff9f59_cppui_modular254,
                            0x169515c89ac5479db0ed8fa6fa311b391cc1235270f4cbc5c29e7cbc30e8732a_cppui_modular254,
                            0x25479fbfb3a68f982388f2621001101608bdc29f6ff037696d9161f5cd9a4fef_cppui_modular254,
                            0x14475c4bd520451f3c852cb0311a578ca7f8e6e972182196ce09486e94be6071_cppui_modular254,
                            0x045a691066cc66bec9baf2798833a1dfd3a847502aec8d5f5c4e73363d097799_cppui_modular254,
                            0x26029c0c267c799fb833ac8a11e3a3f0147a8ca037221b90013b8bcb37eba683_cppui_modular254,
                            0x163facb34ff572fbf7c946969c1c260873ce12a6a94a3e45b8101d5b948d1641_cppui_modular254,
                            0x2c714e96e1913b351d969320cc69d5ec13e06a6275e58688af8ee00c4240ee28_cppui_modular254,
                            0x1c1661e2a7ce74b75aba84665ecd2bf9ddd6268f06debfe2d52b804eff1d5fa6_cppui_modular254,
                            0x06a69ae795ee9bfe5e5af3e6619a47d26635b34c2a0889fea8c3c068b7dc2c71_cppui_modular254,
                            0x113d58535d892115c5d28b4c19a3609374dbdbadf54195c731416c85d731d46a_cppui_modular254,
                            0x2ab89102e2b8d5e638ff97d761da6042e534f1ff47f7917a2ca1a74063b46101_cppui_modular254,
                            0x03c11ca79e41fdfe962730c45e699546349031893da2b4fd39804fd6a15ad1b3_cppui_modular254,
                            0x27096c672621403888014ddbbbfc9da1f7f67b4d4cfe846c6adf040faaf2669c_cppui_modular254,
                            0x2de32ad15497aef4d504d4deeb53b13c66db790ce486130caa9dc2b57ef5be0d_cppui_modular254,
                            0x0dc108f2b0a280d2fd5d341310722a2d28c738dddaec9f3d255754448eefd001_cppui_modular254,
                            0x1869f3b763fe8164c96858a1bb9efad5bcdc3eebc409be7c7d34ca50365d832f_cppui_modular254,
                            0x022ed3a2d9ff31cbf82559fe6a911843b616945e16a568d48c6d33767129682d_cppui_modular254,
                            0x2155d6005210169e3944ed1365bd0e7292fca1f27c19c26610c6aec077d026bc_cppui_modular254,
                            0x0de1ba7a562a8f7acae93263f5f1b4bbec0c0556c91af3db3ea5928c8caeae85_cppui_modular254,
                            0x05dbb4406024beabcfce5bf46ec7da38126f740bce8d637b6351dfa7da902563_cppui_modular254,
                            0x05d4149baac413bed4d8dc8ad778d32c00e789e3fcd72dccc97e5427a368fd5e_cppui_modular254,
                            0x01cdf8b452d97c2b9be5046e7397e76ff0b6802fa941c7879212e22172c27b2e_cppui_modular254,
                            0x1fc6a71867027f56af8085ff81adce33c4d7c5015eced8c71b0a22279d46c07c_cppui_modular254,
                            0x1040bef4c642d0345d4d59a5a7a3a42ba9e185b75306d9c3568e0fda96aaafc2_cppui_modular254,
                            0x16b79c3a6bf316e0ff2c91b289334a4d2b21e95676431918a8081475ab8fad0d_cppui_modular254,
                            0x20dff1bc30f6db6b434b3a1387e3c8c6a34070e52b601fc13cbe1cdcd59f474e_cppui_modular254,
                            0x0212ac2ab7a6eaaec254955030a970f8062dd4171a726a8bdfb7fd8512ae060d_cppui_modular254,
                            0x2f29377491474442869a109c9215637cb02dc03134f0044213c8119f6996ae09_cppui_modular254,
                            0x0984ca6a5f9185d525ec93c33fea603273be9f3866aa284c5837d9f32d814bfa_cppui_modular254,
                            0x0d080a6b6b3b60700d299bd6fa81220de491361c8a6bd19ceb0ee9294b24f028_cppui_modular254,
                            0x0e65cd99e84b052f6789530638cb0ad821acc85b6400264dce929ed7c85a4544_cppui_modular254,
                            0x2e208875bc7ac1224808f72c716cd05ee30e3d20380ff6a655975da12736920b_cppui_modular254,
                            0x2989f3ae477c2fd376a0b0ff3d7dfac1ae2e3b894afd29f64a60d1aa8592bad5_cppui_modular254,
                            0x11361ce544e941379222d101e6fac0ce918106a463290a3e3a74c3cea7189459_cppui_modular254,
                            0x1e8d014b86cb5a7da539e10c173f6a75d122a822b8fb366c34c8bd05a2061438_cppui_modular254,
                            0x173f65adec8deee27ba812ad29558e23a0c2324167ef6c91212ee2c28ee98733_cppui_modular254,
                            0x01c36daaf9f01f1bafee8bd0c779ac3e5da5df7ad45499d0991bd695310eddd9_cppui_modular254,
                            0x1353acb08c05adb4aa9ab1c485bb85fff277d1a3f2fc89944a6f5741f381e562_cppui_modular254,
                            0x2e5abd2537207cad1860e71ea1188ee4009d33deb4f93aeb20f1c87a3b064d34_cppui_modular254,
                            0x191d5c5edaef42d3d02eedbb7ab8562513deb4eb34913a13421726ba8f69455c_cppui_modular254,
                            0x11d7f8d1f269264282a263fea6d7599d82a04c74c127de9dee7939dd2dcd089e_cppui_modular254,
                            0x04218fde366829ed90f79ad5e67997973445cb4cd6bc6f951bad085286cac971_cppui_modular254,
                            0x0070772f7cf52453048397ca5f47a202027b73b489301c3227b71c730d76d6dd_cppui_modular254,
                            0x038a389baef5d9a7c865b065687a1d9b67681a98cd051634c1dc04dbe3d2b861_cppui_modular254,
                            0x09a5eefab8b36a80cda446b2b4b59ccd0f39d00966a50beaf19860789015a6e5_cppui_modular254,
                            0x01b588848b8b47c8b969c145109b4b583d9ec99edfacb7489d16212c7584cd8c_cppui_modular254,
                            0x0b846e4a390e560f6e1af6dfc3341419545e5abfa323d817fed91e30d42954a6_cppui_modular254,
                            0x23a6679c7d9adb660d43a02ddb900040eb1513bc394fc4f985cabfe85ce72fe3_cppui_modular254,
                            0x2e0374a699197e343e5caa35f1351e9f4c3402fb7c85ecccf72f31d6fe089254_cppui_modular254,
                            0x0752cd899e52dc4d7f7a08af4cde3ff64b8cc0b1176bb9ec37d41913a7a27b48_cppui_modular254,
                            0x068f8813127299dac349a2b6d57397a50275142b664b802c99e2873dd7ae55a7_cppui_modular254,
                            0x2ba70a102355d549677574167434b3f986872d04a295b5b8b374330f2da202b5_cppui_modular254,
                            0x2c467af88748abf6a334d1df03b5521309f9099b825dd289b8609e70a0b50828_cppui_modular254,
                            0x05c5f20bef1bd82701009a2b448ae881e3a52c2d1a31957296d29e5763e8f497_cppui_modular254,
                            0x0dc6385fdc567be5842a381f6006e2c60cd083a2c649d9f23ac8c9fe61b73871_cppui_modular254,
                            0x142d3983f3dc7f7e19d49911b8670fa70378d5b84150d25ed255baa8114b369c_cppui_modular254,
                            0x29a01efb2f6aa894fd7e6d98c96a0fa0f36f86a7a99aa35c00fa18c1b2df67bf_cppui_modular254,
                            0x0525ffee737d605138c4a5066644ec630ab9e8afc64555b7d2a1af04eb613a76_cppui_modular254,
                            0x1e807dca81d79581f076677ca0e822767e164f614910264ef177cf4238301dc8_cppui_modular254,
                            0x0385fb3f89c74dc993510816472474d34c0223e0f733a52fdba56082dbd8757c_cppui_modular254,
                            0x037640dc1afc0143e1a6298e53cae59fcfabd7016fd6ef1af558f337bab0ea01_cppui_modular254,
                            0x1341999a1ed86919f12a6c5260829eee5fd56cf031da8050b7e4c0de896074b4_cppui_modular254,
                            0x069eb075866b0af356906d4bafb10ad773afd642efdcc5657b244f65bed8ece7_cppui_modular254,
                            0x171c0b81e62136e395b38e8e08b3e646d2726101d3afaa02ea1909a619033696_cppui_modular254,
                            0x2c81814c9453f51cb6eb55c311753e84cbbdcb39bfe696f95575107502acced8_cppui_modular254,
                            0x29d843c0415d35d9e3b33fadcf274b2ab04b39032adca92ce39b8a86a7c3a604_cppui_modular254,
                            0x085d6a1070f3513d8436bccdabb78750d8e15ea5947f2cdaa7669cf3fae7728b_cppui_modular254,
                            0x11820363ed541daa10a44ba665bf302cdbf1dd4e6706b02c9e2a5cda412fc394_cppui_modular254,
                            0x201935a58f5c57fc02b60d61a83785bddfd3150e05f1df5d105840b751a16317_cppui_modular254,
                            0x0a8c2820c56971aae27a952abd33a03d46794eedd686cd8ecfed610e87c02e9a_cppui_modular254,
                            0x180638ff301a64ca04abd6d0bd7500b6650b65ff33e6be1fd50dbc163a281877_cppui_modular254,
                            0x095c716266f1de59044f97114a4158a3f85ca8a937cfbec63e9b321a812dd36b_cppui_modular254,
                            0x17c31ea02fbc378320d86ffed6c7ca1583b618c5c1a687818d4087a497d73490_cppui_modular254,
                            0x05b86c4bb8ef318b6a7227e4192d149d3c17a9764ccd660de4d50a77f192a91b_cppui_modular254,
                            0x265bc95df4a4c4876ff70d7ea2fde2c7ab15f4a6ae0d237cd6ce74ba986c7a7b_cppui_modular254,
                            0x24752b47bc6c6bc8d9bbe48f5fef2f6908701739c5f5b4b3d6c886d4715c7929_cppui_modular254,
                            0x14814a1e0f492a4ea0d86e527a96482178d624b98da96ee5e583b9324d974efe_cppui_modular254,
                            0x10def931073b6479bd60577378f29381997c8e041d3cfb3dc7523bca906f00bd_cppui_modular254,
                            0x14f7ae770bf7e95f7f706c0d8ab4ed03fa0b880d28c69d031b4592c98610175f_cppui_modular254,
                            0x1aef50a0cee751b59f926af40e8035d19decc9d428ebe4e775c5cc9dce1ce589_cppui_modular254,
                            0x041935607172f68eba65ca60068dfe3b086c2a2d57d09602951214b57e73cf5a_cppui_modular254,
                            0x26863e9dd24255d1573bd083959b856c0493fbefe83c819837a151d3bf452cb8_cppui_modular254,
                            0x2036efb6f9830965eb3d7a068bd087c9f5adf251ba62052c652738e63ff8b3af_cppui_modular254,
                            0x0c712a975b74dc9d766b639a029969ca30be4f75a753f854b00fa4f1b4f4ee9b_cppui_modular254,
                            0x08014dab3cd1667e27afc99bfac1e6807afdff6456492ca3375731d387539699_cppui_modular254,
                            0x198d07192db4fac2a82a4a79839d6a2b97c4dd4d37b4e8f3b53009f79b34e6a4_cppui_modular254,
                            0x29eb1de42a3ad381b23b4131426897a32709b29d53bb946dfd15784d1f63e572_cppui_modular254
                    };

                    typedef algebra::matrix<element_type, state_words, state_words> mds_matrix_type;
                    constexpr static const mds_matrix_type mds_matrix = {
                            0x251e7fdf99591080080b0af133b9e4369f22e57ace3cd7f64fc6fdbcf38d7da1_cppui_modular254,
                            0x25fb50b65acf4fb047cbd3b1c17d97c7fe26ea9ca238d6e348550486e91c7765_cppui_modular254,
                            0x293d617d7da72102355f39ebf62f91b06deb5325f367a4556ea1e31ed5767833_cppui_modular254,
                            0x104d0295ab00c85e960111ac25da474366599e575a9b7edf6145f14ba6d3c1c4_cppui_modular254,
                            0x0aaa35e2c84baf117dea3e336cd96a39792b3813954fe9bf3ed5b90f2f69c977_cppui_modular254,
                            0x2a70b9f1d4bbccdbc03e17c1d1dcdb02052903dc6609ea6969f661b2eb74c839_cppui_modular254,
                            0x281154651c921e746315a9934f1b8a1bba9f92ad8ef4b979115b8e2e991ccd7a_cppui_modular254,
                            0x28c2be2f8264f95f0b53c732134efa338ccd8fdb9ee2b45fb86a894f7db36c37_cppui_modular254,
                            0x21888041e6febd546d427c890b1883bb9b626d8cb4dc18dcc4ec8fa75e530a13_cppui_modular254,
                            0x14ddb5fada0171db80195b9592d8cf2be810930e3ea4574a350d65e2cbff4941_cppui_modular254,
                            0x2f69a7198e1fbcc7dea43265306a37ed55b91bff652ad69aa4fa8478970d401d_cppui_modular254,
                            0x001c1edd62645b73ad931ab80e37bbb267ba312b34140e716d6a3747594d3052_cppui_modular254,
                            0x15b98ce93e47bc64ce2f2c96c69663c439c40c603049466fa7f9a4b228bfc32b_cppui_modular254,
                            0x12c7e2adfa524e5958f65be2fbac809fcba8458b28e44d9265051de33163cf9c_cppui_modular254,
                            0x2efc2b90d688134849018222e7b8922eaf67ce79816ef468531ec2de53bbd167_cppui_modular254,
                            0x0c3f050a6bf5af151981e55e3e1a29a13c3ffa4550bd2514f1afd6c5f721f830_cppui_modular254,
                            0x0dec54e6dbf75205fa75ba7992bd34f08b2efe2ecd424a73eda7784320a1a36e_cppui_modular254,
                            0x1c482a25a729f5df20225815034b196098364a11f4d988fb7cc75cf32d8136fa_cppui_modular254,
                            0x2625ce48a7b39a4252732624e4ab94360812ac2fc9a14a5fb8b607ae9fd8514a_cppui_modular254,
                            0x07f017a7ebd56dd086f7cd4fd710c509ed7ef8e300b9a8bb9fb9f28af710251f_cppui_modular254,
                            0x2a20e3a4a0e57d92f97c9d6186c6c3ea7c5e55c20146259be2f78c2ccc2e3595_cppui_modular254,
                            0x1049f8210566b51faafb1e9a5d63c0ee701673aed820d9c4403b01feb727a549_cppui_modular254,
                            0x02ecac687ef5b4b568002bd9d1b96b4bef357a69e3e86b5561b9299b82d69c8e_cppui_modular254,
                            0x2d3a1aea2e6d44466808f88c9ba903d3bdcb6b58ba40441ed4ebcf11bbe1e37b_cppui_modular254,
                            0x14074bb14c982c81c9ad171e4f35fe49b39c4a7a72dbb6d9c98d803bfed65e64_cppui_modular254
                    };

                    typedef std::pair<mds_matrix_type, round_constants_type> constants_type;
                    constexpr static const constants_type constants = {mds_matrix, round_constants};
                };

                template<>
                class poseidon_constants<poseidon_policy<
                        algebra::fields::bls12_scalar_field<381>, 128, 2>> {
                public:
                    typedef poseidon_policy<algebra::fields::bls12_scalar_field<381>, 128, 2> policy_type;
                    typedef typename policy_type::field_type field_type;

                    constexpr static const std::size_t state_words = policy_type::state_words;
                    constexpr static const std::size_t word_bits = policy_type::word_bits;
                    constexpr static const std::size_t full_rounds = policy_type::full_rounds;
                    constexpr static const std::size_t half_full_rounds = policy_type::half_full_rounds;
                    constexpr static const std::size_t part_rounds = policy_type::part_rounds;
                    constexpr static const std::size_t rate = policy_type::block_words;

                    typedef typename field_type::value_type element_type;
                    typedef typename field_type::integral_type integral_type;
                    constexpr static const integral_type modulus = field_type::modulus;

                    constexpr static const std::size_t lfsr_state_bits = 80;
                    typedef number<backends::cpp_int_modular_backend<lfsr_state_bits>> lfsr_state_type;

                    typedef algebra::vector<element_type, state_words> state_vector_type;

                    inline constexpr static const element_type &round_constant(std::size_t round, std::size_t i) {
                        return round_constants[round][i];
                    }

                    inline constexpr static void product_with_mds_matrix(state_vector_type &A_vector) {
                        A_vector = algebra::vectmatmul(A_vector, mds_matrix);
                    }

                    typedef algebra::matrix<element_type, full_rounds + part_rounds, state_words> round_constants_type;
                    constexpr static const round_constants_type round_constants = {
                            0x6c4ffa723eaf1a7bf74905cc7dae4ca9ff4a2c3bc81d42e09540d1f250910880_cppui_modular255,
                            0x54dd837eccf180c92c2f53a3476e45a156ab69a403b6b9fdfd8dd970fddcdd9a_cppui_modular255,
                            0x64f56d735286c35f0e7d0a29680d49d54fb924adccf8962eeee225bf9423a85e_cppui_modular255,
                            0x670d5b6efe620f987d967fb13d2045ee3ac8e9cbf7d30e8594e733c7497910dc_cppui_modular255,
                            0x2ef5299e2077b2392ca874b015120d7e7530f277e06f78ee0b28f33550c68937_cppui_modular255,
                            0x0c0981889405b59c384e7dfa49cd4236e2f45ed024488f67c73f51c7c22d8095_cppui_modular255,
                            0x0d88548e6296171b26c61ea458288e5a0d048e2fdf5659de62cfca43f1649c82_cppui_modular255,
                            0x3371c00f3715d44abce4140202abaaa44995f6f1df12384222f61123faa6b638_cppui_modular255,
                            0x4ce428fec6d178d10348f4857f0006a652911085c8d86baa706f6d7975b0fe1b_cppui_modular255,
                            0x1a3c26d755bf65326b03521c94582d91a3ae2c0d8dfb2a345847aece52070ab0_cppui_modular255,
                            0x02dbb4709583838c35a118742bf482d257ed4dfb212014c083a6b059adda82b5_cppui_modular255,
                            0x41f2dd64b9a0dcea721b0035259f45f2a9066690de8f13b9a48ead411d8ff5a7_cppui_modular255,
                            0x5f154892782617b26993eea6431580c0a82c0a4dd0efdb24688726b4108c46a8_cppui_modular255,
                            0x0db98520f9b97cbcdb557872f4b7f81567a1be374f60fc4281a6e04079e00c0c_cppui_modular255,
                            0x71564ed66b41e872ca76aaf9b2fa0ca0695f2162705ca6a1f7ef043fd957f12d_cppui_modular255,
                            0x69191b1fe6acbf888d0c723f754c89e8bd29cb34b1e43ab27be105ea6b38d8b8_cppui_modular255,
                            0x04e9919eb06ff327152cfed30028c5edc667809ce1512e5963329c7040d29350_cppui_modular255,
                            0x573bc78e3ed162e5edd38595feead65481c991b856178f6182a0c7090ff71288_cppui_modular255,
                            0x102800af87fd92eb1dec942469e076602695a1996a4db968bb7f38ddd455db0b_cppui_modular255,
                            0x593d1894c17e5b626f8779acc32d8f188d619c02902ef775ebe81ef1c0fb7a8f_cppui_modular255,
                            0x66850b1b1d5d4e07b03bac49c9feadd051e374908196a806bd296957fa2fe2b7_cppui_modular255,
                            0x46aaa1206232ceb480d6aa16cc03465d8e96a807b28c1e494a81c43e0faffc57_cppui_modular255,
                            0x2102aab97ce5bd94ffd5db908bf28b7f8c36671191d4ee9ac1c5f2fae4780579_cppui_modular255,
                            0x14387b24d1c0c712bbe720164c4093185fcb546a2a7d481abc94e5b8fb5178b7_cppui_modular255,
                            0x5f2179b3a7845836cfced83e64e206f6a6cef2cf737f020b5cfd713c9550fe9f_cppui_modular255,
                            0x1787986ab56e1b56b5443334562b0bc3657d27323b87e3a8485e68ab96d57188_cppui_modular255,
                            0x39ef4b00deefe7e7451adda44428aa22074c496de2c9ed67dcf4861da65f543a_cppui_modular255,
                            0x7271d384cf5c90fd0c48af190c5c765937c7468088b081a99337e6eae53bb20c_cppui_modular255,
                            0x6669e58d04248ca86024fbc196e5f306e522423aa71f84225435328b37a1dd3d_cppui_modular255,
                            0x0c1f1b492b27539d754cba5e46edc1f1ac1c5696da8eb19416b07420bb321c65_cppui_modular255,
                            0x1c4d41a133b97dc467f1f184cf191f331dfc38e79e7e53516c39848c9bd44692_cppui_modular255,
                            0x369ea8e699181b1cf88be9205ab840180c9288e67a359dc0dda4ac74cf9768e2_cppui_modular255,
                            0x4cfa7d72afed332bf0b8a2a719123f7ebfa714b9e3100eaa533dbde6fb985043_cppui_modular255,
                            0x4e592fcde9f3c360e54c6f34d7a8bd41889942e9fe23d9fd4a9e5b3bfbbb3e45_cppui_modular255,
                            0x032b5885586212fb235570996d3a4c40f54ff91598a948ec2722ed865b8438a5_cppui_modular255,
                            0x03f3178956cfd3e2e6614fb134597d3b3cff0d8a33f3523d825982990c068940_cppui_modular255,
                            0x3126e84dfd67a22bf0ce0d9273d8ad40e6109af5bb2bd78d0ac08a16c6248f74_cppui_modular255,
                            0x3527888062f1e2738d7b928e9af244f0a39011390c2dbbcf56d8e087f4087b6f_cppui_modular255,
                            0x64635758efc701dbbe2eb423bf7b5bf6c3d34c6ff92494f3421182a8b187ecf7_cppui_modular255,
                            0x4d7f71960f03db8a2a428cbf77ddc1916a5f4243dbeb2ddaef7b5b5f9d74546e_cppui_modular255,
                            0x37832ba2da93de3643243eba3b9765d75359310617f3fc06d74ac12db57b29c5_cppui_modular255,
                            0x4dce55879ffd9398f96c9e6556a3bb4fc93147965252cb1d6c94b3282ba3fae6_cppui_modular255,
                            0x4ba85e4d2537972c0fd5a4727a58c3d85d98563697a34c0af845bfecd6dc4b40_cppui_modular255,
                            0x582dc453b4cbf6b1d19734b0f337d3423b503703979689f384d0eb96ff5b02ce_cppui_modular255,
                            0x0e6f127f479ee6113540d69b25420a2682f07b23e799566b091a1c891fa224ba_cppui_modular255,
                            0x39c815508d2995bb8ae5035472944706e900b2fb16d5a779fdfff82306f37dbb_cppui_modular255,
                            0x6591aba215bcf96d8aa03220372179a4c5060cfd7f95724ab300d9459f709051_cppui_modular255,
                            0x221807cb4909d549c546a734ad2cd7f60a69e816ace98fad830452a44a343188_cppui_modular255,
                            0x2766a1e33038004da58bce78722380b22b13b0aecb87f38659f3035e1336b53f_cppui_modular255,
                            0x11b5e993e6a9cdc3b5d2f5336dc9bad5074b661537ff890b1babd7f53cada9e3_cppui_modular255,
                            0x29576176f9a5a10e3d0a2c59af26b51f4c5fc86ec59c0f2492deb60ad49eddcd_cppui_modular255,
                            0x51e72c44f9de491c747d8a6d333fb2b3e16ee7571f1340a9a5f6f72363991e98_cppui_modular255,
                            0x2fb360d959be4aa871e071764a5e41eb264d04f0289f098723b69bab09f4d1a6_cppui_modular255,
                            0x03f46b4c3c77957cb595ed61fe13f9e8739a5009311142b69c1e8c07ae250f47_cppui_modular255,
                            0x4683311e382a99927e0ff672cd0543aaebfc0c33ba96ad937818cec979b57b5e_cppui_modular255,
                            0x7117cc69bf566b1b0ba5486b0f1f9bd60f2f945e3cbf33a2ed17076f4caa0dd6_cppui_modular255,
                            0x3bd670c3ce88ea43f254d61c2a9b56d6a4dff19ab5c4d28989d271f3dd6bee25_cppui_modular255,
                            0x2fd2ed0ba1135575995d15061ddb487f2c5c6005feed28d8a01b9d7bee361a1b_cppui_modular255,
                            0x6a66704e22a81e6b7ad8e2f28edd8c9c9a10abf17e053f4d89665810332600ec_cppui_modular255,
                            0x5cbc378be1db3840b32d8d2ebfe2695f810f932a206aacece707ca693f4f933e_cppui_modular255,
                            0x35b716410b3c9374d42e7d39eaca316b6568f0a14cb14d519967aa3ff9970aac_cppui_modular255,
                            0x231c6db056e47a01c192db40e586ededc929b564667377a10bd1465f3852811f_cppui_modular255,
                            0x4904d5de1f512eb14b0f856acb016c7a43079b2f702303752962f336558b0f32_cppui_modular255,
                            0x56d6bc63f429bb7fec7bdd133581f2abc74406a57607c2ba3302481eddba4074_cppui_modular255,
                            0x519d0daccadfbb0167fa79d1afdf36b25f28b9f74f1e65d21d28ce1022579735_cppui_modular255,
                            0x0576cf2418d6bd88f352bb26da1066637575f85688cdb981c7787f8094e5a71a_cppui_modular255,
                            0x16672be70221dfa20aa110bdce12e1e66ab171db4eadd9935baa0e3aa49e437a_cppui_modular255,
                            0x1e51c73bc2aeb9e877d9c2c18f17b03ea3dfcc04adfc649780ce4bcbc43b0b69_cppui_modular255,
                            0x1271c830507a211c8e2ebdfb372f79c8a42a9e84e4fdb0dcb35d55e4d155e169_cppui_modular255,
                            0x67077397c2b01db4de4b78adf97e0ebceb20cb91647db49a7bc06a5ce1b25544_cppui_modular255,
                            0x2e5454b258106b63f0ab01924767b4aecce371202abc28a260adc45f35570b9d_cppui_modular255,
                            0x440f72769f137a8078f05063cfa4e2b73b2381b72b68e97b1c1e9cd18df36f82_cppui_modular255,
                            0x6ae1478fc162c50032fef2ef79c93ca7ee25b16358704f434f6cddcce2fc9c40_cppui_modular255,
                            0x0c0f3630409a2242a39ebb33c5c7cf18965b8932621aab4ca2c315d4441b6987_cppui_modular255,
                            0x0d1bd84a786a990adf88b51f253bd9032cb50ce4682bafe103893af36d5e75dc_cppui_modular255,
                            0x30ce425059810dd94aae2f255666b0fe8bc52ff701c385c43a998926539dd401_cppui_modular255,
                            0x395a1e753153b56d1a9ec2ca73099425e446dfa668dc73da2ea311abe5e3d96d_cppui_modular255,
                            0x57f09d89e827d00392fdc0c3d21b1a5bae2d689894ced82f58e256a03d20ef91_cppui_modular255,
                            0x1065b71b135e4feb8b3cba3c252daa084cb5624b0ba76f48f6a03854bfdbcacc_cppui_modular255,
                            0x3d5f53bd162f053f045547952a06bc83bc413e17957977e359d9bd4c8883203d_cppui_modular255,
                            0x05f467a5081bd3479d6b49f697b0a75d264b42b95b2bed475cd58ffd05322d85_cppui_modular255,
                            0x6f5ad8e3ed272494c36a5a52a7d034e04b633460c16a512d0d8002f8fa0e3484_cppui_modular255,
                            0x23c293275e282bf15cdbffae1f00a2712e76aa6d62820542159e9d6f115df3b8_cppui_modular255,
                            0x3757e7009ca9bec8bba29308b9922354eeeff3beb4113174bf8cde584722d31b_cppui_modular255,
                            0x406f25e72d0264ed50473ec95a7ec53ebe114898f84deb06e53715ae24725342_cppui_modular255,
                            0x046dcfa2d6d655c7c551f7440772b056e7d3f2c65ac52e4496c4fc753130ad45_cppui_modular255,
                            0x49c2e954d649ee1c4e72ce8c1833c33796ab29dbb0486fe53b04687b2063259f_cppui_modular255,
                            0x2caa8aae247ef83e63dbe8e5efc89d7d28ffd8bf7a5331e245af8aebc872a759_cppui_modular255,
                            0x5efa9f8f32d9ec1d3a3d8cea806e068909b3d3562fdc3f91f2d899f8109bc717_cppui_modular255,
                            0x0df424bdf3b0c60395cd7380029a633692b933250b79371e09122c8c39aa1301_cppui_modular255,
                            0x2d012e3e811cf4b88aed6f38d5cc8c3456dbae1741f501574321906efb474930_cppui_modular255,
                            0x709c043fc648c48a5bfb5ea25d5f0557d03aadff9d6ec1afaf2032f3aadb9dba_cppui_modular255,
                            0x1bb9b23d6805ed1179a1dad95740513dcea114185a8ed34e17dc8077dc830916_cppui_modular255,
                            0x0fab922a838c55af1e2349b1e50b56d0690c200d0f2318aad4b7bd8a38a47f61_cppui_modular255,
                            0x4d58799d4501ee8e89c73db7a4ff48d9f5e80fd5984afc67f3054f59d3dc74d1_cppui_modular255,
                            0x4f130b733cb78f3940da337d187934e48765956ad2ca7b75b7bf8e293b46a758_cppui_modular255,
                            0x03e7812afd6c480faef03c3beadfb882923a743a4e60e58a259e7ed4598cca97_cppui_modular255,
                            0x739ea276a5ef7008fffc02a3c853f4d56eaeee7df395cbee8bbe6b502b81ca1a_cppui_modular255,
                            0x0ae97e00a91a4e761815fde0e9506629373ef7ce765ecb1bc7ba0ca2decd7d01_cppui_modular255,
                            0x6d6c41e1315436781a774555668cc3d41c99c78dc107f443ba0ae60cdb287c16_cppui_modular255,
                            0x18d683776871c1918c2b5c632cb1854dff865c4b1b8bd66e46d2fa2a8d515c34_cppui_modular255,
                            0x3597acab641c21dc5475eb8b04b0e2ae91700acad1b543e8c7e69d574eb5a15a_cppui_modular255,
                            0x63df64938297594b4e8bf2ddd6bcaee6f2b9703e5814ddeca44d341b9e7d24a2_cppui_modular255,
                            0x009ab455f6b4c7755da22615073e9839cd12a88d1f9b583d7ad61bde4009b873_cppui_modular255,
                            0x09e21d43c56b0abfc26d0fb7a3ebfd3a7743bbeea99ac2b8f61cc23d1c673a12_cppui_modular255,
                            0x4db404b9eae6a9f39417be43c93a9f6d136a0784b73789d590ada0a60df0d16c_cppui_modular255,
                            0x0c6f0ecaf32a3d60aaebeaf3f8ccb00a10ee19def3836b78fc905bfeaf2b80a9_cppui_modular255,
                            0x3518d688407ca0e548165b9796a4279d038720408a3c822dc44ce8974ea8ad8d_cppui_modular255,
                            0x27ba9d4584a23881e23aa0340dc266b32b56455c30e6da78b37741de7ac5b185_cppui_modular255,
                            0x63d33e44fda7868d50858e482fbff7c29143d60fe00817cf32e0efab4c3ad6eb_cppui_modular255,
                            0x561a72b93fecdbd83d67a5022d9a221cf21b22cff2d79c114bf01c71f2641ae9_cppui_modular255,
                            0x48a1625a9ee1102971aa28bc07a5ba88ac6424801502ff4fcb6994824c2e5e36_cppui_modular255,
                            0x46a003c184ecf0e00fa8ef7dbb356366be4d63a3847634b46a18ecd47667d1bc_cppui_modular255,
                            0x37d6efb2876f3cba63a60821e50853d0997947b96f633607bb36ded243ded838_cppui_modular255,
                            0x14f96acdb291ed2bf98a5bed063f6911598bdff1f6c0219bbefa447ab1918163_cppui_modular255,
                            0x573d156263dc8edf24efced0c465587cbdd1a2c792cbadd58abf95e037d3c668_cppui_modular255,
                            0x46839e7d70370149b35b3a07d8406acbaff07615747d2101bbad18abb9891f95_cppui_modular255,
                            0x3b74a3420d1b988408fe8d8fcb51a81f16f8d17d082da9ba61fbc8031d8ff59b_cppui_modular255,
                            0x059f3301178a22026798b07a8578611d7c56c16bfbbe6a058f4e44016aaa172d_cppui_modular255,
                            0x467d9ff3508feb318b07acf9184537462e987c58b7ef486873e1de428eaa3f32_cppui_modular255,
                            0x716cac6b0fc8f63d406d38d6b82c8ed4e5665e449f07b572b83f43c9f9ba2004_cppui_modular255,
                            0x7121fa9ca506687b3c49dc2060731c85ae48596be138148d8ea365333b8f03a6_cppui_modular255,
                            0x10000c75e6e03366bba4f59c68f312becb7ae0c30d4aa141940a7531105ef7e0_cppui_modular255,
                            0x375487214c07542fa5b6a5736344466a06c2cb4c1838c9966925cd8c5888c3ca_cppui_modular255,
                            0x2361aaf969f732be06b159772a097f3518ed9485449edcfd367e289f0964c486_cppui_modular255,
                            0x2ddba8679308f327c27023a893c0458d1e73dcd64a39b22b130fd9e4f283f906_cppui_modular255,
                            0x6303e21755b1de4d65495bae9685e05162245106f53d7407ec0883e39695b15c_cppui_modular255,
                            0x5aa3dddf8da369722b2e1c8f2aacf0625d08264f8a0ed320df110ab42f5b0c1f_cppui_modular255,
                            0x3525eb41c2db9cf9cd08652d815d7c91f3294defeee702efedb5f777284cd1fd_cppui_modular255,
                            0x0079ae4df49f78b97cb0e3c3f4b225538d4a0c4827e333d27a29398c17c26c9e_cppui_modular255,
                            0x533c8c1b05e2dd7e7e19ea4b027cc8bd559c2e2a622207b0c13bc7afdd7bc3b7_cppui_modular255,
                            0x4989a01e4fe4b1bd544e5cd4288895068897cba899ddb01779f6e2b08024d3ab_cppui_modular255,
                            0x1c7f5858eabb1e2b8c3104808dc68ae3de05381fc74704a2afbd2fcc42cdd3c8_cppui_modular255,
                            0x55faf16bbea2ee0f35413b9808c135fb1e4729c90b4cce4c345238c6dc557639_cppui_modular255,
                            0x156a82f8e5aea455d9c8c436f89c6f9ecbce0ecaafdd13b93f255e075c72ebd0_cppui_modular255,
                            0x37c7047032df0027d7bc128e9a107582f25ba0b7387230a05864aee420724703_cppui_modular255,
                            0x40ab847795176c24af06d5000ceedb82d87492cbde5c1c262a83a9b6b6f4b264_cppui_modular255,
                            0x5a73bece689545bd2de9ef263d5036152f36e2250c76711e8bc9ed9bda7af685_cppui_modular255,
                            0x1c4a903be5dff4440b4f38e56f988cddacc57371aeebb06cb64ab5d21d9562f5_cppui_modular255,
                            0x5bba81a692e87b51c7c176730fd05cfd100b0bd86d69b4b4f367277a2302b2f8_cppui_modular255,
                            0x2f875bdd6669a8ff920c3d7bedd74c101541d4b184b7e1bc0b90ddb26902319d_cppui_modular255,
                            0x5e89035bbe943f9e6024db13c58bbc748d3f1654050c7ffe084b763efceff3bd_cppui_modular255,
                            0x728cff754d7a76a7f8b00656412ad8874e7bab9827706ca6d6d13c72a0c6812e_cppui_modular255,
                            0x6dcfa6338bfe3569524a968abc95c706801fcc695ee3f5854a79e4689625481c_cppui_modular255,
                            0x24ce56469aeaa4243053bb62c07100002b8f74c4ac74c350beff0c0be47e5a51_cppui_modular255,
                            0x6a72f954f591825caa43c3ba7ccfea7aa1a00de5a681e52de6148252062f8363_cppui_modular255,
                            0x59922ae3f06524d2028e9aa00a136613d4306fd5f4247ad0a6a587be0fb0081c_cppui_modular255,
                            0x50d8b98688f4980b1a0c2b5313f8ac9660b1e9199b5f59ed3709e0f1d9185552_cppui_modular255,
                            0x3184262ef10e9b0ab57cfc898fb68342cb86ed6e25e536fa94caa605b4a3caf1_cppui_modular255,
                            0x69980a1f4b883cac1039fc47dba993503d4ae5ad40ed112a5a5070090006f73e_cppui_modular255,
                            0x1d5a91b930b89934745ba00bd9094b67f95e41e3778fe0420880e80bbf8078e1_cppui_modular255,
                            0x0ddebce4b6ca45d69b2f70c8b54e425615c1aadadccda74e0882eb79c445778f_cppui_modular255,
                            0x68c8362e93a371d7c9551edf3e3f3b14c54c729c1fab0fa6eebae7da09855826_cppui_modular255,
                            0x3dcc6a17e074d0350ffc0e5426e1bb6894e6c958f96f3d7d9c4240b948cde438_cppui_modular255,
                            0x03b8aba0ee959a4e51cb5cfc458b0f4ad3a9b59797394c3d3c9eb57adeca2308_cppui_modular255,
                            0x0f24cc57f3b2fbf25375c71d71bbb97b2d193fc1a203ccc514c074d461001ec4_cppui_modular255,
                            0x71e9bfa7f66afbafbf139a70baedfb1b202a2e51e6b6c420e28dd342a5eb0cd6_cppui_modular255,
                            0x3ac9c11890e96a2dcda6405a6c52a47e803d6674e65117f1a8adf701d68cd02a_cppui_modular255,
                            0x45c00146e1b89ad5ccb8a02202482023751b88997d8fba1af5c0e7a68dadb63c_cppui_modular255,
                            0x1f98bdb8dc318e3e2e28cc3d8b85e334f74b57e15b02e1637ae035b04bda3b5c_cppui_modular255,
                            0x2ec077dbbc7bf2affe7ddd8b8a7f900f3019cddc8ce55cf9782004f65f51257b_cppui_modular255,
                            0x32c377fc988f600a2c2ef5d5376e2e31faf1c2d1a618db011fbfec1ff337568d_cppui_modular255,
                            0x0a820d131da844383bdfc1a053d8aceec7f2eb345ab6c21d38e829db8d05861e_cppui_modular255,
                            0x5bd95df8a933f7b7e263e013f45a92c0e786dba563e210b77d5a40f961092e60_cppui_modular255,
                            0x264cf7b75095fb96b420fb3f31c064299e78e796e8b3735bd0a186cd3817708d_cppui_modular255,
                            0x27d3e47b2f11ada6a9a5d329e00a128c9836be92ee92429ab891e71d11dc29f2_cppui_modular255,
                            0x64354b412c8cfa1319e4afd891e619a8fbbde04d85bef4ad0548689295d2bce2_cppui_modular255,
                            0x0db0f967487ee52e0836fb7135bce37fbd32887e911de52d0b855a5afac1f770_cppui_modular255,
                            0x1c9a155911b36c896475995417197faad870737a9ce5d9d3a5000f5396978e9d_cppui_modular255,
                            0x65ae557151ae9ec7f870fa2804bfb88e669dc0f8865b140f964f1f93180ac531_cppui_modular255,
                            0x52c6f6242517362c066020764fef4a5574749106a6dad534d136e7fe885fcb40_cppui_modular255,
                            0x6e44c5bcd5dc6591e2f84290a313b71a04da8da398dd10135d22bb23df41e883_cppui_modular255,
                            0x2146d3e371040feba8595049a285944bd45a458dccb059c785c2adf032c8b710_cppui_modular255,
                            0x16db9ceb3074a795499a37c20ffc9eaca9b07a5a25824aa6adcdb19fabdff0b9_cppui_modular255,
                            0x5903725fd86fec14c9cf2a273017eb01d3a1785039397060650c4e228a6e6571_cppui_modular255,
                            0x54c75952f908e3f99e05718bd1f59bb6c414bc2aebacd81c47189885cbbc566a_cppui_modular255,
                            0x0dba4abc7f188e33e7f309317b7b9f5c22870ca90bcee7b576dd0b52619a39f6_cppui_modular255,
                            0x3950231611808399ad3ba5b78cad4c6bed6f364b9346541dfffa4d16366d257e_cppui_modular255,
                            0x1a6d8230bb9e8d1af552b9bab8babfe505931dd87e200fc7b3c57160a5bc4ae2_cppui_modular255,
                            0x6b3dd35220ecd616eea4309ac9a8118e9dc65a3f7c1ef52dde7a3d33578c43a0_cppui_modular255,
                            0x6da00240c3505b214c8d8ce3f48914247adb9f0ecf239d7baeada5183d31ba54_cppui_modular255,
                            0x37c3720b132d3a719424e29c37acb7dfbd709ec9497a3162175424bf063c6e18_cppui_modular255,
                            0x500f85a3d06a0b5a05c5e93ae70084802fd499c7e6ed1ee6e26b4bf8fd6838fb_cppui_modular255,
                            0x2b37f70d73366d32d575186d0787fc8ce539b73f83c6e7eaab27be85f4faaaf4_cppui_modular255,
                            0x1d8efd6e52d4f936415e5c4814f3366804e2386857a4befa2a53aab21ddb68de_cppui_modular255,
                            0x33303b8a8f2d811be65a977907d17d133f3a64c59fe2a9c5c2d4517e3eb390e3_cppui_modular255,
                            0x2c1ba860f51e0c2eaf4a9a6bf095c65fab3ee15c145f404fbb0272b5ca14a449_cppui_modular255,
                            0x0b0849c7a3adea03a89d101081c9c9f4f66ef917d09c7957584db9a75aec2378_cppui_modular255,
                            0x41e7e30c77579da7809c3e757821c869b53f103fcb752ac82f8a734d4abdc792_cppui_modular255,
                            0x182e66be60686c8c5e6518430845f98924fe8d7d43e628bf75ff52a716371b9c_cppui_modular255,
                            0x373b2508c2fca1a288fa4f54a6edf02f2661e664dcf4ff2a74f3d06b1a00ddc4_cppui_modular255,
                            0x1735b442b3acaad0bbe630f308e03f1aa6f56bdb029e50c1393533cee1a45c30_cppui_modular255,
                            0x22abe8ea470a0372911bcef1367e10aa220491d76caeaa5959feb5d75f4a1f9f_cppui_modular255,
                            0x5caab387eb997f774f64151ed21abfa5364a83c6f065d92bd9c92f2719b8e80b_cppui_modular255,
                            0x57b33094aeff828377897b56e1c432978d07c668ef25a36bc5e2e835aaeff725_cppui_modular255
                    };

                    typedef algebra::matrix<element_type, state_words, state_words> mds_matrix_type;
                    constexpr static const mds_matrix_type mds_matrix = {
                            0x3d955d6c02fe4d7cb500e12f2b55eff668a7b4386bd27413766713c93f2acfcd_cppui_modular255,
                            0x3798866f4e6058035dcf8addb2cf1771fac234bcc8fc05d6676e77e797f224bf_cppui_modular255,
                            0x2c51456a7bf2467eac813649f3f25ea896eac27c5da020dae54a6e640278fda2_cppui_modular255,
                            0x20088ca07bbcd7490a0218ebc0ecb31d0ea34840e2dc2d33a1a5adfecff83b43_cppui_modular255,
                            0x1d04ba0915e7807c968ea4b1cb2d610c7f9a16b4033f02ebacbb948c86a988c3_cppui_modular255,
                            0x5387ccd5729d7acbd09d96714d1d18bbd0eeaefb2ddee3d2ef573c9c7f953307_cppui_modular255,
                            0x1e208f585a72558534281562cad89659b428ec61433293a8d7f0f0e38a6726ac_cppui_modular255,
                            0x0455ebf862f0b60f69698e97d36e8aafd4d107cae2b61be1858b23a3363642e0_cppui_modular255,
                            0x569e2c206119e89455852059f707370e2c1fc9721f6c50991cedbbf782daef54_cppui_modular255
                    };

                    typedef std::pair<mds_matrix_type, round_constants_type> constants_type;
                    constexpr static const constants_type constants = {mds_matrix, round_constants};
                };

                template<>
                class poseidon_constants<poseidon_policy<
                        algebra::fields::bls12_scalar_field<381>, 128, 4>> {
                public:
                    typedef poseidon_policy<algebra::fields::bls12_scalar_field<381>, 128, 4> policy_type;
                    typedef typename policy_type::field_type field_type;

                    constexpr static const std::size_t state_words = policy_type::state_words;
                    constexpr static const std::size_t word_bits = policy_type::word_bits;
                    constexpr static const std::size_t full_rounds = policy_type::full_rounds;
                    constexpr static const std::size_t half_full_rounds = policy_type::half_full_rounds;
                    constexpr static const std::size_t part_rounds = policy_type::part_rounds;
                    constexpr static const std::size_t rate = policy_type::block_words;

                    typedef typename field_type::value_type element_type;
                    typedef typename field_type::integral_type integral_type;
                    constexpr static const integral_type modulus = field_type::modulus;

                    constexpr static const std::size_t lfsr_state_bits = 80;
                    typedef number<backends::cpp_int_modular_backend<lfsr_state_bits>> lfsr_state_type;

                    typedef algebra::vector<element_type, state_words> state_vector_type;

                    inline constexpr static const element_type &round_constant(std::size_t round, std::size_t i) {
                        return round_constants[round][i];
                    }

                    inline constexpr static void product_with_mds_matrix(state_vector_type &A_vector) {
                        A_vector = algebra::vectmatmul(A_vector, mds_matrix);
                    }

                    typedef algebra::matrix<element_type, full_rounds + part_rounds, state_words> round_constants_type;
                    constexpr static const round_constants_type round_constants = {
                            0x5ee52b2f39e240a4006e97a15a7609dce42fa9aa510d11586a56db98fa925158_cppui_modular255,
                            0x3e92829ce321755f769c6fd0d51e98262d7747ad553b028dbbe98b5274b9c8e1_cppui_modular255,
                            0x7067b2b9b65af0519cef530217d4563543852399c2af1557fcd9eb325b5365e4_cppui_modular255,
                            0x725e66aa00e406f247f00002487d092328c526f2f5a3c456004a71cea83845d5_cppui_modular255,
                            0x72bf92303a9d433709d29979a296d98f147e8e7b8ed0cb452bd9f9508f6e4711_cppui_modular255,
                            0x3d7e5deccc6eb706c315ff02070232127dbe99bc6a4d1b23e967d35205b87694_cppui_modular255,
                            0x13558f81fbc15c2793cc349a059d752c712783727e1443c74098cd66fa12b78b_cppui_modular255,
                            0x686f2c6d24dfb9cddbbf717708ca6e04a70f0e077766a39d5bc5de5155e6fcb2_cppui_modular255,
                            0x582bc59317a001ed75ffe1c225901d67d8d3764a70eb254f810afc895cbf231b_cppui_modular255,
                            0x076df166a42eae40f6df9e5908a54f69a77f4c507ea6dd07d671682cbc1a9534_cppui_modular255,
                            0x531f360b9640e565d580688ee5d09e2635997037e87129303bf8297459ab2492_cppui_modular255,
                            0x30be41b5a9d8af19a5f922794008a263a121837bcbe113d59621ea30beefd075_cppui_modular255,
                            0x39f57e4c8a1178d875210f820977f7fcd33812d444f88e471040676e3e591306_cppui_modular255,
                            0x3514084b13bc0be636482204d9cddb072ee674c5cb1238890ee6206a3e7bf035_cppui_modular255,
                            0x6372b6bc660daf6b04361caff785b46bbe59eb6a34ab93e23d6364e655dc3a36_cppui_modular255,
                            0x422af985e648814bec5af62c142828e002d4b014b702760106b0b90c50d11de5_cppui_modular255,
                            0x3296e51f12e0f5c49747c1beb050ff320e2eb7422807eb0c157a372dba2ea013_cppui_modular255,
                            0x3b76246abaf33b03dd5b589b80a7fac0ae7f1ad8a9623bb7cf7432c90e27358d_cppui_modular255,
                            0x0b40e7e02f5cb836c883c7cef72ec48e87c1808f7d829e2ee0bec0ee709f7409_cppui_modular255,
                            0x2ee81b5c29c93b8a6e8871c01d0380a698e547475359b4a4befc22ed2232690f_cppui_modular255,
                            0x341ff90fc4a8afee9b74c464955ba9b357252e915b8d39ea7c1318eda718f54d_cppui_modular255,
                            0x55eddabde058f3b5e9dae90873ec9bd7b05927da36925e7dfb7bc290c1da125e_cppui_modular255,
                            0x6b34ad8cec56aae4595c403377cd2aa990a2f09b931f832781221965bb081b1c_cppui_modular255,
                            0x707de76df294fb845309d2160e1bdffebefd57a80c8658899e2c95e77254c752_cppui_modular255,
                            0x05e9b152bfd4946b9c109f930eb01892f314597507d28c735a266f4277bb2a32_cppui_modular255,
                            0x1589a5cbcee13b696b6f0a1dbbabc08394ab00ed5a6ae6435020e9e3e2fc909a_cppui_modular255,
                            0x7116a5d027fe73fbc45bfc60fd875c3116fe3a567e830d1d2d38655223dbd7ec_cppui_modular255,
                            0x05382ee6ad97381eb3137f5a90ea13298dac6bc7c2204906044fafc01bfe6ae4_cppui_modular255,
                            0x0900bcfe5e7c1b7d0aa80c714b7b2a0c1df7473362138a9dc5c552d11c1d0015_cppui_modular255,
                            0x0513deb89d2e48fc729440dc08d0256a79cda84d511a04e0d92cce3c7e55a7c2_cppui_modular255,
                            0x6bbb5f1736d499fe3fda42ad40a2b124952ac35fe970ebde38c65cc20ad2afc8_cppui_modular255,
                            0x5782ac68a8da0ba09f4d17e7e4b46caa4411a27e60be92168ce75bed95453e05_cppui_modular255,
                            0x2d83f3324639c5d83a1ffcf6ac693eef98d8ea4877d547c62b304b0a9f4a0c28_cppui_modular255,
                            0x16d3a13700ec503e29ca4d0c6342864595134408b6668bbf1766bb48d7f96cba_cppui_modular255,
                            0x318050e971e075931253b00430d35f89f40a88fc73d62150882a8e87149d7244_cppui_modular255,
                            0x7180760dd839d8bffbf9b1e26826cb4f6de65fa868a8143e1dc8c2b6ac6d1ac2_cppui_modular255,
                            0x5cf2aa95907e59c4725cc17c8cf492f9a7eeef2de337ac227a983c444ae0e80e_cppui_modular255,
                            0x2b8345763484d7ec02d6ee267b7c737ca9de41e2186416bf91c65eb0cd11c0a4_cppui_modular255,
                            0x055aa90aa60ef9b7f3c29c7500c64e6b85929220a6418dfad37ead3928059117_cppui_modular255,
                            0x541d5e4be0967bf49a595c1d8290b750305a334f3347c01b57f8ba313170e1ca_cppui_modular255,
                            0x05c0a1f16f97f582caaf4338f018f869e8dd0fa32f007bad1a1a4780053d5817_cppui_modular255,
                            0x01519e13858591aa93b9c1d7f849276ac1d2011b7fd19a475371c7968d9f52cd_cppui_modular255,
                            0x69c30d5a27f4dffa19c956c348287a704676d999f23044036b9e687a45a1a113_cppui_modular255,
                            0x58c93b899aa53e06e82b6346e36338841ba7279d2b7a0ecd3aa20f292852936f_cppui_modular255,
                            0x06b8a12870a15479d41018fed6f1a29102ae23e13d0fbccec93ace48bdb9dc93_cppui_modular255,
                            0x33eda3c347379e61c2297aa1026682d22f95dc3c7e46e68ab3adb4b0939d76e2_cppui_modular255,
                            0x187728045111275b93a1218a148ada85a1f6e2059c443ac7d61fe81e3130b89b_cppui_modular255,
                            0x397ec485c5a8b0c8a03ff543e9a9e5a4dc0dd4849fe955bb77b452e2e22c4f17_cppui_modular255,
                            0x2f33f8de90f81248455d5a6592667092992be0468372addbaff664caa84cd2d5_cppui_modular255,
                            0x061a1a458994ddf9f38c5edfbd737d3ceb05deaee685058b14943e7e9246ebca_cppui_modular255,
                            0x4b73ab5b9d35f47307b731e3cf1a1a22e7068e2744f2af0ef6bd78bf8aae4845_cppui_modular255,
                            0x5578b7ad5f8d4f3b8e618af7d8d5ec8bf837d2d9486527fe2f9bf7464f8516ad_cppui_modular255,
                            0x50b4f055d860f89e12883209f847a4b1a2395fb419eb53c182dbb555c962255c_cppui_modular255,
                            0x0b2da770936d6c778be289557ddd2ca024b93fa38c5d4541344e883a69611813_cppui_modular255,
                            0x47d8441e1ae7cb8ffc52a18c67afff3cf7543cad51605b2d4e2513f1e1868b68_cppui_modular255,
                            0x619da3bf44b42acd949ed572c9f3c195ed20b0b91bcd9e95ee3750d26f3b0ebd_cppui_modular255,
                            0x6c9e249e89b2b4cf9cd7772950e0cc9d06688d4f051095eafd116371ede49ab7_cppui_modular255,
                            0x210bd3217a141c55877d4528a4e80d5d81d78de7addce85994082281a6250d4b_cppui_modular255,
                            0x4e1d8e4079c14c83847af6394d7dc23f33ebf71593379583ec574bf5c86ea9a6_cppui_modular255,
                            0x699187330fc1d606e8b31b677651a2c7d1c87d4d001018031792cad0ad3f2826_cppui_modular255,
                            0x2946bfc0f45c1f1a0dc4c343a85259f6a6237f064481fe66eda76f01998a01ea_cppui_modular255,
                            0x5543e07588375c6d800e5e42d1bfd8b7a92a2a35d65b234ded85f879f82a3d66_cppui_modular255,
                            0x660e9d0f2f866e8d12b40dd9d9c03cc8b9ca78600bd649f0fffb2c388dcc8b43_cppui_modular255,
                            0x38f06c48d4dc53cb1b69619244cc2a610fdc4229ea316980dffe9131a72b4209_cppui_modular255,
                            0x5c9a73a16521ddf463f9de314dd5f7255bc66add48297615b761f34e4636762d_cppui_modular255,
                            0x310931f0204c9936fe659e9ebbda832c930172130b3f5476c6c6ee5e7fef3e45_cppui_modular255,
                            0x72eb1d833664d8989998af11441ac49654c12210b3465e5ac67a99679634a3af_cppui_modular255,
                            0x6981346585a2a466a9255841f710e1d083bdcc21c0aa6721745e158218767a94_cppui_modular255,
                            0x0370a259836b3766d563ed3cdcf55ace52655111a1017d8c76eaf8f97e81d858_cppui_modular255,
                            0x4f63c45a324b8b974c22a20a6c670eb62d47ef900541b63f1d362b8bbe4ec418_cppui_modular255,
                            0x6a4c7347121c2d4745ecffaad22281cc4d58ea74453b7d2b625b890190fdc7ad_cppui_modular255,
                            0x36d8869bb69a51ee99622af09d6878c5b715084b25f6e4560a7498557fe87fb5_cppui_modular255,
                            0x18faa7f51e1b7a442f9123806872094c0de8a46a6d8402f31f0cde3fcb878394_cppui_modular255,
                            0x3610d022aacbe58593e0d6aa7eefdca767f5ddfe7fa1fb9fb4f80225d82b617b_cppui_modular255,
                            0x3b5f13d6a8bbff31569bc6860087b2a4b361146a04ad5fc7396a3d0c59f68c1c_cppui_modular255,
                            0x40e919335051c6aaaee033745c41b6fa36739a097d94ce6eb075ec03da2a978b_cppui_modular255,
                            0x2f54586ab9b7886340f8ed5254f29128a85e2fb1e3725bf3c9cd8bddadc947f1_cppui_modular255,
                            0x00606231b689a040363e5afc050f9fc9296d6c620a885eeaffe91be387cbe96c_cppui_modular255,
                            0x4b55696db6b0fa327527a76e6ab6b688561c879e53d858e4c90a1122210130e1_cppui_modular255,
                            0x569c39bd78356991953aef4b1a01fdf71710bb05eea1f447c3e5efe13bd62894_cppui_modular255,
                            0x537f73fcaa256497a2582e45105f1dc10f39c7fce9b88cab5523af3f5f82dcd9_cppui_modular255,
                            0x2d58d32120c25995cd0754ab9fdf9ad67d67623cfd1fcbf489f51fa6e6eee4a2_cppui_modular255,
                            0x37cb0f655951fca18a4ccdddd4d8466f8839ba8e320a104cb47a59cd387d322f_cppui_modular255,
                            0x4e29d154430c9bced788d2eed8f3e01b5da24c1d3710e490bc40ee6d5903213c_cppui_modular255,
                            0x47597b7a9018192ef22d6dd24555af1c0c51d8a90b54d8a0bdc2df7967d7a28b_cppui_modular255,
                            0x4e01b43205fca0b4a32582abe600f3a326035fe7e028cb0569bac43c997b98ce_cppui_modular255,
                            0x0172ffdfba7e43ca807d5b5de7727b4e41706c1f2858c1e8a46c27ed3eae5ff2_cppui_modular255,
                            0x2216dd907ab98c0d1e720a46ef83334a236d2c134ccf35ef8e889421e70ebe03_cppui_modular255,
                            0x168709f668b635f03607a39390a0de71306d6430ce2babf7292d789d25c0f8d5_cppui_modular255,
                            0x0ff6a3823440877dfd355dea80595e21115d0dfe3472cec4ad1437572cc6151d_cppui_modular255,
                            0x44e37699b3c72f50ec1a754c72e6fa3f5a074181dd63d189ba36447d34e536ff_cppui_modular255,
                            0x267298d2e46227f7f7f422e3059f18d83a8795731b13f6568ce54730cd3fe9ae_cppui_modular255,
                            0x1ecbe7a60848077203373441a5b09b44693a155fe226442259e37ac47209235a_cppui_modular255,
                            0x31cb23e6b5d7393577d5f5c3368c5bdd5b434ee6319f07e502031cc393d4eccb_cppui_modular255,
                            0x5d4c550c4a6eccd74b74d6279b3d9bc755084588156a1bef673657dc2116ecfc_cppui_modular255,
                            0x226056b5dec9afd19190ac48740c3b5ab1bb429b19f56894a3dec3f104d238c0_cppui_modular255,
                            0x09077c021183dd37ad10451ded70d7ae6ec4819ae76ce23fb2a0be63e69907d9_cppui_modular255,
                            0x53545c868ba0fbf0ed1ed7a24ec11b2ecfba5b37fd5cee80774e1ecdea991ed4_cppui_modular255,
                            0x69521c33d148e678ca10b33103812cd27597c4a6cddbe83f4970d4b96e03304d_cppui_modular255,
                            0x01d5779be7477b96aac6532ef919e61c624072be54587e0698999dd5f460e446_cppui_modular255,
                            0x57875a44441d2f191ac7d8de42691ab55fd3401bbaf04b786ef0603b3edf2927_cppui_modular255,
                            0x1d5c957da0832d5b94e76f7abdb190972774b594ed232810bfcafe5441839d37_cppui_modular255,
                            0x1b678335a80fd045fc7ce1897aa129f67bd55ca9ca801bd88eb7cc868538bd7a_cppui_modular255,
                            0x31e69d706a5c1e011c1cb1809e5bf1857c90f9f50b9e1ae5ad36e4d3dcdbb7ed_cppui_modular255,
                            0x485df8462ed7a18de34aa6e99ecc9bbf2db075a096b56bc2943b76a99c4bb1a0_cppui_modular255,
                            0x1e46fdcbb3705f663a350e78f99024912d80c95779195807aae82cbb494ce9e4_cppui_modular255,
                            0x441d0fa0e9cb86c3a2a1f87151681c603c3e028f1a0670be2149eed4f0a24f08_cppui_modular255,
                            0x02a3caff274f40942062340ec1fae17c1b1e97c2f0fc7e847c90e9317fea2c0c_cppui_modular255,
                            0x4caf281080c0b2f2f638bf0f4859442f4c9da94e9994dada34c5c914130c1a9e_cppui_modular255,
                            0x444470c6c49b5b9a38181c3af20bcfea572450946135baea85cfd6b692fa6464_cppui_modular255,
                            0x6d5e07a13376fc883bea2dcdbad7f80b7780f231cdd33f5b98618f42cc49ec2f_cppui_modular255,
                            0x1b9470418a07d8c88c767d1e63e8d5cc7f810cc530db1340181ecbbb212e0f70_cppui_modular255,
                            0x4134c8666c685b712f4aec72077c540ef4a041dcaa123caabd57b83fc6266f14_cppui_modular255,
                            0x3d5d0489e27362db9bf0cc7217477d81d2a73e1a44edc43e32d43bb544287c9d_cppui_modular255,
                            0x71d7d4a91945e796f538f03b9324497489009ec1a0a403de062ed5bb4d7c2400_cppui_modular255,
                            0x646c3d732a94f722384ac266b41e06cf21bf24fb9426c9556d8ac9514f0875f7_cppui_modular255,
                            0x4f860c9e5d9bb73057d93c207902d9e60fd6a7c779fde1ebf16b853dba1ea9ad_cppui_modular255,
                            0x05801566eb9e119e2f9ace565c9488cd999d66a5753eb4b9887363137baa09ab_cppui_modular255,
                            0x0263bdb8654cf1245ae4589370dfd5eeb109a50944eef54308566055b887ee01_cppui_modular255,
                            0x4cc39561e65eb05cb8c83f9854750a9114a996eb23e6a0bb07d2d61f0baf0a62_cppui_modular255,
                            0x36b544778b2fdb94f808ad8d077b7f0b44f3bba515ecdf026919e2fed09a106d_cppui_modular255,
                            0x3fb1f7aec47cbe990151d4bf703c38349b95f409abdf0504e67c1a55ef82294c_cppui_modular255,
                            0x637e7eb19cf539aada7e48bc6b72e5ccb0e3f6913f18a0d55696dddfcb1b587a_cppui_modular255,
                            0x73bc630fcece6947fb81ac8e0f1f1671ed6042c3ef3bbb12ed554f28b48b46ec_cppui_modular255,
                            0x304b46f52d597b964fbec3fc0dceee442febe6131359e156c194ab7be2a11e6d_cppui_modular255,
                            0x067d85956dcfff7fd9f6a0fec505b7f4998e3d85672623677a6d974d6b111de6_cppui_modular255,
                            0x65830d8053bf8afc0ba5274f1a4c4cce617fa624b480f13ed3eb369fbba78e67_cppui_modular255,
                            0x6c32c101e08a962bd996d759a6c012a4d97aedaab9fc99c1fa735a16cd24dd44_cppui_modular255,
                            0x11fb2d160e41a1845fd14578c617285081fb1a16a21b36cfd5065b30fac574e3_cppui_modular255,
                            0x50aada39348c4736f6c59f7f053c488ed999a33ad23501d9c635aa03baf90db5_cppui_modular255,
                            0x5a5f0e3a32b260fbdfdc8c0eaf3a99396992b50b6dbb63a9d1e1ddf9c91d78d4_cppui_modular255,
                            0x62c9f6d9aea355d358f2986ad487c2ae443122e1edfb076930865608d05c3b39_cppui_modular255,
                            0x520cea06cee20150703a1c8000d4a5f22b3efeb9e34eb90bad0b4ff091b33683_cppui_modular255,
                            0x6da4e4682545c1f4c0076f5845fbbcf48632a9c193a92593d12d248031f2c893_cppui_modular255,
                            0x1ba5502cee2ea2d07a64f68f0a7492d2426382a5b9662d0410e086107399989b_cppui_modular255,
                            0x6ab843ca92240f8a82862da071d53f048272d55425907fc8d0e60dcccd5a1ea4_cppui_modular255,
                            0x3f65c2dfa6bb39c1b291c40f810cc912015384a2a24fd322b6375e27bd069322_cppui_modular255,
                            0x6a2df71a64cb0d9a548e3b65ba4e646ff5e519cab564b5f77b3fe08e038b9c3a_cppui_modular255,
                            0x64776bf2b66bcd09c8661ee6ca6b8251bb4aba5a7ba181464d905db561ca45e1_cppui_modular255,
                            0x6d7bed0d258b518eda13368f00be2cc0a94d71cc203d5905c35b10a3ee53eea8_cppui_modular255,
                            0x371b958b5c79c889d1786edfe404119773f728822637fb4890b8847a93f97af1_cppui_modular255,
                            0x56923182c33cb4dbf0988ba2314378dfd7491b3467b6134e6283c87a1478cbb8_cppui_modular255,
                            0x3c4304994ef664d6aa19e3db492c306534281b5b6f857fa6ffae67bdba99c09e_cppui_modular255,
                            0x0d003bd3068fa94c4f7bbe6ba02993acd341a27ed2fd7ecaa4e6b0b9d0abd85a_cppui_modular255,
                            0x1073cb8c08510e7d88ed4cdf78e96b297cabe9d6677db47289b056c2a640da01_cppui_modular255,
                            0x5c57522580fbc75883658d4b7b8ea07e1a4fc75f453c09edd9d249ff1bd31ae0_cppui_modular255,
                            0x2a5bec9b422b4dc64958f4752d0c091ffa7904e0ce4809728d16235bb41d707f_cppui_modular255,
                            0x379c4a9b4174c5878f72b60fa985f7aa86c1fd868683bdbe8fae194cda2e56c7_cppui_modular255,
                            0x3634e042e79d046adb911d57b338e78f51ac7d212c5a5c6dc4fa1a05ddb58c82_cppui_modular255,
                            0x3ace976310c5040e1484d1a6d42993ac5923d474ce5497a3fac468af25843a01_cppui_modular255,
                            0x3f5a856ab863b7584bc2e6e4c610b9df55a9306eb68894d630ff7d04f243e6f5_cppui_modular255,
                            0x0d52822f5581fe9c5dab0b1f8d04eae183deb87c89504544a3d5558594b3149b_cppui_modular255,
                            0x3c119e173586c22059bb09d2af4fc1044c8fc44f709233f7625e5fffa6696596_cppui_modular255,
                            0x3e154fd5a026d7c6584faf8c089d82fd560f138392a8d4a5fe287859994c96b5_cppui_modular255,
                            0x47251339c44d737b21df0ed1e204a28b68c9abb58f1cf2232f8a2da433e24b0b_cppui_modular255,
                            0x73d84625f38db2f3842d7724d8e79d6d0349a93b8d6142603eea382ba6ed8692_cppui_modular255,
                            0x42929bffc19bf9cd1c53d10440b0760a3be6442db20458b692b4ba3901e6003f_cppui_modular255,
                            0x39b16b0fc3700aa93e0cac53fcaf7e84495ac3b49553b2e1a5ff9f73fe74de50_cppui_modular255,
                            0x2b715e21640cfb6f77b91a4f6d3dcaef9b5faa7c0bfe94c8d80b0824292603bc_cppui_modular255,
                            0x306bef0c637b5d7c8d6486915f6623f4e1ed81971f40772ec60feb5e243d32a0_cppui_modular255,
                            0x5287d6ece65ef5df6e1c65dddf1d97cfa019157a5c90c004527c9d7c7496d814_cppui_modular255,
                            0x0d760a2132c9092b0c8c89cbdf4fb1bd282791ef6284b73a44b313e8118e7d0c_cppui_modular255,
                            0x5e830f4484268a349e4d9f6178ef745460f1f8456b04d0dc7814844052d51eb5_cppui_modular255,
                            0x2468669481610965d8439f60a66aa61fbc7b18e82b35aa4755873ec4db82174e_cppui_modular255,
                            0x23b6ea9e4d1fde701c719c2afab1272ea22b172bf7afe0837364ad9a2f698bd4_cppui_modular255,
                            0x412024b2e86e9d5e903a5fbda26200be47003e3b0dcc322480d3079850606cc0_cppui_modular255,
                            0x1f64c17825c1ce9333d211d45a555b5ceaa4608a354ed3237db56225b3a9459b_cppui_modular255,
                            0x0b66fa87587ab95d5d29dde50cd606a1bc2c45fd223c03d0693c88b13ae23039_cppui_modular255,
                            0x3086c386026698e733e54e5e17f65cb26c17fe64e76f85902cc184d5dd8ef0cf_cppui_modular255,
                            0x72036acd9ef575414d5437327d902da6396cc70c0bcffcef2a82b4c296b5ea93_cppui_modular255,
                            0x53d89e4470b3ea1eb861717e47c08fda42f6e61fc08118b16645ae5e8fdd664f_cppui_modular255,
                            0x4ebea65d1fc5c5167b1412ffcbf8900a8df2096c25d7011e6c74f500662465f8_cppui_modular255,
                            0x5ee6e1e0312e78e2e67b246a95afdd79e2e7a5b9b0ef6ee36c3d1639f9736e65_cppui_modular255,
                            0x1d770c0cc2c2231213624d58b7875b1715554f6100784bb2b545e405c7fcb94e_cppui_modular255,
                            0x2ea5c9837af445988c480fc6a55b1e5640dbe38d5e8cf1ddd85bc42c3658d9ca_cppui_modular255,
                            0x6fb78d12c35235f738b1667749064d0066fa7cfe3a9624cb0944f16d37bc485e_cppui_modular255,
                            0x35b75e89e794282cee1e66991ccfb2499dce4366b88d7be5f7b5775c12103a22_cppui_modular255,
                            0x50e83b08162e7ccfe2d0f19aea4753ba83ef5c40572d6e904cfe2419ee9d901d_cppui_modular255,
                            0x3fc5c93031cbcecf12d5831aaa6b2b3071657cd669f7b377b2fef4a7bfc9adf2_cppui_modular255,
                            0x37895bdfe29a174b98cd4b49104e56ea09e41c7b50f9aa95b400b529c545f5b4_cppui_modular255,
                            0x695e405509a0981035ba77e27cdcf53f3bc15d20fe4e43a335aeb6406ae1837d_cppui_modular255,
                            0x104985a48aa7e0a668d8cc7140c255ed1b8482ac5febbd3d7a1cca0e96cf0682_cppui_modular255,
                            0x118220b30330f1954e7d94d40fb1043a1a79ca83e68e9ef590601a86a4a917a4_cppui_modular255,
                            0x098b3be7845a63543c13d211efac076b94a9528d34cb355faf0ff7a0d5ee9991_cppui_modular255,
                            0x69ca1313dcddd8c2f5c5c7ee93a1d2a94726c0c0bc4a303fcf83109b23bf3621_cppui_modular255,
                            0x570c1bd286b258b8bf11e8b85a2eb0c6dbfc2e4cdf01a0cde5464aa009b5bd43_cppui_modular255,
                            0x4f2921de3696018e0d1ca7cdd5a4064ebf51845ab25b2d395b71c341ea8527da_cppui_modular255,
                            0x19035c69cbaf0e0e7e02c5c524a8cc56de0e52d1936a9a10b7580f0c0555878f_cppui_modular255,
                            0x2b8fdad2064a6f58d01e8c48d49bb25730780055829c1faead0430afcfbc5669_cppui_modular255,
                            0x60ef9a74bbf8b98cb8248856492257f30c7520b3353a6fec9d90d48be46070ba_cppui_modular255,
                            0x4c9a6bc8284e783afd6c425f8cbdab82c0db3eac060a2dc00eca48ed6d1d052b_cppui_modular255,
                            0x68e6d3a83ac8e60c92d2860ff7557e1fbe3b91c38fabbde8b28371dccce2a10b_cppui_modular255,
                            0x56e0e39848046f0305d268b28aa753a41d48586e8579d5f95f12dba60e181d4c_cppui_modular255,
                            0x5176824fd8c92fed23df24c382a9fdf86aeeecab0b6716bef53da57bd3f551eb_cppui_modular255,
                            0x3aaf796b71041e8b2b494bca3b030f56a0c5663149093c8a179c0f3e24d0f718_cppui_modular255,
                            0x101cd65865abc573f5382df3636f4d60bc669aaa70f09ba040d61ef8d09c5296_cppui_modular255,
                            0x2581f83d616d932b438bfe0062082d4e1ed7d34b9a1cf63580199731d44a4b25_cppui_modular255,
                            0x65d74f6d1320dd1dc9412547b130bc7ad03c4e80cd8a44c108f24ec7aa35489a_cppui_modular255,
                            0x0d5cb6e19c9aac7d9f51f176ed42d008317a189dc4f6fc5c36fc6d451a035916_cppui_modular255,
                            0x0e367d17423501e62db9fd487f72076f2d1de6dabd3c175341ce35f925c9941e_cppui_modular255,
                            0x3f3f101f7c8abd6bebe6b81dadf0ff5fa31ec7140e317909a8d2f94ce4adc890_cppui_modular255,
                            0x6d5f212b5f4775095ab1d20fffd41dd73ab69b4ac60e9de11693f8e6bab88e67_cppui_modular255,
                            0x6b11154212e86e185a4cb17dc2b9dc061f72bf9cc3df5f95f7b87f1101d09f1c_cppui_modular255,
                            0x43f4cf980ff1a9101ca3c4601814f8de4124d108be2584ee9ffb9505188d35fd_cppui_modular255,
                            0x5d9be9303e3a25e8fa1abb6f2a7e3250231091100f9d7311b050b52666ec8f02_cppui_modular255,
                            0x1eb3b147885e1261d9034ca89a658817caef5ae629e1265cd32c6ef89ce704e9_cppui_modular255,
                            0x1595d95dac2c4653d32b01c3fbc294b2922140e41b93c5e7f5702212226d7140_cppui_modular255,
                            0x578b22f1f6d6eeb61507f0de1c817bb876b9cd079a18be9e99e2faa8e02618e2_cppui_modular255,
                            0x4de38f88c5e8ba1890b3695c912ccacd63721298c9ba3d3668b44f2a13b40abd_cppui_modular255,
                            0x0b9df0b81af072be21be9f08df336d3babe6ed5bfc199c73f2e97ccc73de80ae_cppui_modular255,
                            0x2a1a8c6d54abda22954e90386d40cc7d5c4f54c592ec2c69a9574601e88b6559_cppui_modular255,
                            0x5c5d96136cd1c4ae8fa1db9273083567345b407eb66f73a313ab8ad1a76cb157_cppui_modular255,
                            0x1ade9e2b734e937fc2fa04ca445236faf24e6d47ad1a4baa3408102c0d1e6363_cppui_modular255,
                            0x49354c394824998704e44eeb2ba6cb6fb431c334b648e6c87565e5fe133e8079_cppui_modular255,
                            0x4ea258f019a8055902a696b85547652519b8d8d92de4bf18e2dbfa41264a9a6e_cppui_modular255,
                            0x008a5162adf5ebd8711fd8139418509e472abc02908084f2e494086232336808_cppui_modular255,
                            0x6badee92872dcc00812a1cbc8081dd65ece0c7d3512af5a9df5fed7428557c81_cppui_modular255,
                            0x324c64ef2693e966965246bb7bb8c04b57a21fbee9be8c4a10096222bc83cc51_cppui_modular255,
                            0x3f14138eee87c93b0fbfe7efdcfa906525b0ce0f3b9a7431a492f8cb71514219_cppui_modular255,
                            0x0db99fa5ce25d50f557415ad181f1399840574f678b2534cae8f774bc8703009_cppui_modular255,
                            0x23d984702589f3275211041a4bde9d79329967723ec029da095bdbe733e97381_cppui_modular255,
                            0x6c5144ace155e976e287f1b95951194895bac2e5d54b07b62c3afe0eeafcbe39_cppui_modular255,
                            0x57a3e420fe7e0638bfb4d0b2c6286c2946166a6eb17836571909da153c3204de_cppui_modular255,
                            0x156621c4691a9240863577f10e29dc66a37d1b94e756869984c22d9f9d284726_cppui_modular255,
                            0x1b1e774a7ec903650adffe34f6aa8201d356e41e0951d38fb83a89413d078e4b_cppui_modular255,
                            0x514b940e5717c1ae53ea29b9a5a15998e294f69c1f553fe56124f66a16a78d53_cppui_modular255,
                            0x16350c6898d04d355d966c1d7827eee076a1ebd90781639e120feab665391ea9_cppui_modular255,
                            0x5b8b30d8c5ae46c4171d40478886c71c28fc86a3ae4a52ad1c05d8bcb9991b52_cppui_modular255,
                            0x5226cdc8a40c229ea4fb08f2c10e0e24cd41f24ca5fa5b5ab73e7340f632e727_cppui_modular255,
                            0x64383db664537c84a0a4030c3318f2f19cbeda46c70460035ad9d9240011639d_cppui_modular255,
                            0x61068a086ab73c87701b2642af25f6a430240936ba473a9a258cbf90db275277_cppui_modular255,
                            0x5bf320a3e8a48c6a85e2dffc4740d1b381ec4aa0771d885dc16adee569403ad3_cppui_modular255,
                            0x2603e0fd03264a856c1a7b8f1c5a22c3b98f4858c345e8e0a68e3f6424dd2dfb_cppui_modular255,
                            0x100d221342e64ed7e4f1520be70f5b0134031f8a31b4790ebb8e0a89e50b42e2_cppui_modular255,
                            0x0e61bad85ce909438ecc028b55085ec2cee0dd3ac5a7bcaa79d96186747a4606_cppui_modular255,
                            0x570a2045ca0fa7288d7f372f36bd075c2517a9743c9baa46503c4396e1f316f4_cppui_modular255,
                            0x1a64e108621e134020ea761d8f2c5bb42f24fab7641b095f1d164d1fc7b8be90_cppui_modular255,
                            0x097f0f28fd299e3597ffd761e9ae8b0fa46526c9d78503dc9dd5f61df3a085d7_cppui_modular255,
                            0x1d1063cb1be0f9f96aca5e5e39be9df69c96ff717c7d0c7cfe179cd6dce27505_cppui_modular255,
                            0x3e30f5d48b3c2475b8f3ba08cba27caed32b1cf67f76ba9223803733e13ad863_cppui_modular255,
                            0x2b30db4198cd832506017fa26430d204476113cc791ee110cf5586af5ce3824c_cppui_modular255,
                            0x2b520e374519be203c022ec51dcf8d972dd01abfaea371de9b1532647fca7bfd_cppui_modular255,
                            0x183b9a8e45fd480e822f8a97a8d2f127d0ef561914903229fbc5602bea46cb37_cppui_modular255,
                            0x4e01e6edf11ef4c94fe8589f9622a70709330a12e68591f6ea7dda994117bdc8_cppui_modular255,
                            0x52ee256fb3031d20fc299de7fabd0d4ef2e7f12539760dafb0fbc8560a40ee16_cppui_modular255,
                            0x327f5e141e4758d3d9a94c1628a57c817cf84fc0082b8dc098adbe84c1430979_cppui_modular255,
                            0x3d0e12036899e5be167de13913901831a714ea5617b94de6de070ddc117bac71_cppui_modular255,
                            0x1d9466d50efd1be3080d0aec4b81dd5cdf1ad4681e3ac04d08057f8fe49cdf0b_cppui_modular255,
                            0x2360abd7728da2dcda3f495a9a4f0f2aaff1d2420b8f6a7fed6592e1463f3d00_cppui_modular255,
                            0x23c1df4ddd6da863a1a2837e5222150278adfd4faf2fae7beaf64ed67a30736c_cppui_modular255,
                            0x1e98ec3b325a2a11738273f94516a9d56107f33062661e571342bc043764cf77_cppui_modular255,
                            0x431de5d108f8f7109df3059abcc16ccbd17e18676ef64f8998498e4a3f331fde_cppui_modular255,
                            0x550937f2bf0f1adb53f412d49ffd2886158703c375f87d059461f740d655e3d0_cppui_modular255,
                            0x1341fa99aca4bfc0f511dc9a9bc57c1e7aeb41ebb3a9140f5f93af1b3aeeb582_cppui_modular255,
                            0x706889448219016f970b32463a87e355b55ce0a34401dbfe4dd19fb3f93dec2e_cppui_modular255,
                            0x28d6207e409ab1c6e8e196d9e363040070b6c6fc4685a5482f80ba38cb792dc5_cppui_modular255,
                            0x6827087ecdf4e6bc7c396c59de859cbf08f92c361b5174e7f681ba0e72f83aaa_cppui_modular255,
                            0x553e112dab620286f6cf2d31325b971a6516dc7776a6e5ef37bcb11d1785299d_cppui_modular255,
                            0x40b44f7413d152f0d46460c54e9572fd91174b4b94a3595d709119e49925354c_cppui_modular255,
                            0x4d324dd7dfdf2380ef9f6d3c4f4bc4c5f90dbbbf2f1fd923256913f33a45cc09_cppui_modular255,
                            0x609b3ae79dcdc8a8379a690394c95805d862bc31068b572ac126bbc082ebf8b7_cppui_modular255,
                            0x33973520a1d9fb67048d64a22ad1b75b081d88c135a556dbc1b6a8479f75eaa7_cppui_modular255,
                            0x3bcb7630fc45d34b78fd253d0b5275ecfa85ce48125ef7275c3a9205d01b85d8_cppui_modular255,
                            0x1287f419048e81322d73bb9333e9b854e4ceac4b993b5342547263a486b42e34_cppui_modular255,
                            0x2a2f5a5a689471d5ef46d669e449ccdc1d37256618722f08cc2c7e75d06fc277_cppui_modular255,
                            0x38c913fdc729a28b7e354947f2b6449029976d442e349bc1c2acf3b0fa28bc92_cppui_modular255,
                            0x421826bc690adac2b1f3637bc5e2333cb5e4bce3f9e8eac1a0a76df32a7ebff7_cppui_modular255,
                            0x30ac2452c3a07bb924b6f7ed47cd6581499d532c5f90bf7fbc69556ff3bf6b09_cppui_modular255,
                            0x40ce93f92b281e538efbe7cec9a22a9c005eef428dde3cdd46191630f563ba04_cppui_modular255,
                            0x4fc3dd6720c87f672f7b6ff129e9b2a3236ec760a71f78aee84925d8e7616e97_cppui_modular255,
                            0x3f3ba6f9f12ca6f934f92b17f4f3bd8ec261e5870610557f687bc734eadaa2d4_cppui_modular255,
                            0x11d9eedda8d94fcbed859f5787fe20b7d4483cd319d8215530e2e316c89ee635_cppui_modular255,
                            0x29981cff92be6c882c89feb59849d014fcd163699b5b4fdafca335552c4581d1_cppui_modular255,
                            0x4c4fe2838d175c666c0d3f20d8dfefdcbcdacebca86e013d8ad29b6a0cf6bb79_cppui_modular255,
                            0x630428a99469c03f9027d3c601864185d360d920771ea950732cf000b869a09a_cppui_modular255,
                            0x46a776fbf1f36d7fdfa7a210cbb2ffa533540068c169e12f127cb14d9b587056_cppui_modular255,
                            0x41a775960677e6c5fdf73c2a409b6e5c08e271cbb8c825f598a1801c84fde5ae_cppui_modular255,
                            0x3086af931c41d791deb57f7f82dc511e4d349f42b52c3e0080097c4e44373dc8_cppui_modular255,
                            0x155516da7a229b61392a39cc10a67112f512203cab706428f5fbbb3a9fd89fbd_cppui_modular255,
                            0x41bdb1e32081ac55f42969658f78e308bdf50175b619c3ca8e3bfdf1ca984684_cppui_modular255,
                            0x01344d21e02b9c20d0d886a02167cf8502c3614ab909ae2fa7929b12d3e88519_cppui_modular255,
                            0x733a3e92f74b793915beab78e87bd88a2227aa5406df54dc9a2c5e80a11f71e5_cppui_modular255,
                            0x6a6cc17a31ba2fe1411cdebeb0809bf4ff0069b0d6ac681edf816ef4c59b6f64_cppui_modular255,
                            0x0a77e0a85b06c1b152098066bd36933264641627192e3acdbf611bd002918820_cppui_modular255,
                            0x3efb107ebed9b44672f679bffec0121fb509d19e97ae1bac3a86384e274c8c94_cppui_modular255,
                            0x3c0c4b441b0ea7ffe03c011db9aab4f86ec4849a0c783a3b7af21b05f5654482_cppui_modular255,
                            0x28072c7bfa64f6cb97e4341cd18809ef5cd083374fbec26370c2b0ac02dcdafe_cppui_modular255,
                            0x1962306e92b3c7295b2f7435ed8f67dda3a15ec6d8b0786d7727d071663ab22b_cppui_modular255,
                            0x594dc533611f7f588838f894a26b1cd27432c63f9fbe03ef2d95d9a2d191ae3f_cppui_modular255,
                            0x3e287fec491c686222949bc16c2308ade64e3a0a1dccdb25d64f9d5b94ead6e7_cppui_modular255,
                            0x2a95d47fb725b3978a7f90e601f2a9ab39074b35594e0bd133f9c5f34d765d42_cppui_modular255,
                            0x29c603ecc031a9750a4d826e4abf3874bc76c76cc7ea306b3b9636f9653ff58c_cppui_modular255,
                            0x0bbff6ba283aa42f01172bb82a2838e50941227abc3a2a5b1215b9a6d68de07c_cppui_modular255,
                            0x73c7ee55aaa453d36ed857353bc227375244a7e554ceeea2018eb9cb39a51e74_cppui_modular255,
                            0x3ff41b13d4cb3140ac8426322e88ff6f16895d88e6de3336cc88c693e0d38175_cppui_modular255,
                            0x03043688d4c991763362912a460be95b668fe9b1823fe90febfb3ffc7652ab24_cppui_modular255,
                            0x33a29a0d56a7a64d36a67da2c691ff3eaf8ec7f0d78b357e7d2254c5b0e28f73_cppui_modular255,
                            0x185db562fc75b43ba2710ad5e9114486b3e9712fe4c88f98b333c0c6211ac882_cppui_modular255,
                            0x147b89a0cff9083b8952b3ef292c683f75d523f932711c6e1db3f28f5163b1fb_cppui_modular255,
                            0x58ebc5d6b50bb1e4fdb4dcdfae1b69027978826f757ee4dc10d34f963f98fb59_cppui_modular255,
                            0x1318791367815809badf1f3ed677e50cef92021c65549b2dabaa52c7b424f5a9_cppui_modular255,
                            0x5bce78553694ba32f793c8d7f8d09ac63d0d7ada32b888d61b87849f3eda9557_cppui_modular255,
                            0x026bebcc38f0b2804ed21f2e2b16af2194375ff2559fbc588a8962caf0b684c0_cppui_modular255,
                            0x494bceff689f9885a3998de0eaaa7ac71a04522700f2e067efdbb037c6e53c66_cppui_modular255,
                            0x03ebaf5f0602347c4ed2bdb9a86eb955cb5cd5378f7a6f369dccb69792de8bd2_cppui_modular255,
                            0x3626d91f9f05334cb32d3a42eed03f7a553a0ed4cada2db08b45b548bd3b3655_cppui_modular255,
                            0x63ee9e5c5cd3c83e93757ed93358ff0583d761e595b62f11df27bd4292ffb6e5_cppui_modular255,
                            0x705dd80b2db4492c8b9984439b823681c4d9c8dcddcc04b9786a90051513a0e1_cppui_modular255,
                            0x2636ac2ac559be8fe509641dbc67e55db47bb051e05ef06301020c9501f110f1_cppui_modular255,
                            0x4781b8da302c7764951730e7ac0892de64537d94db2e19b84eec5a2d9539288e_cppui_modular255,
                            0x197852b9a62e16779725f35cd8daf52ffbc8cc9c902c16923f2ff8873795ca86_cppui_modular255,
                            0x1c3e49f33fd73480b280dba7744cf67e244449048f8fc84f7b6e452b4ede9a35_cppui_modular255,
                            0x41d20cdc6a15c07fd9735c89b155412fcbb7bd3cdfc27abaad2a3a8a90e99743_cppui_modular255,
                            0x0c3a7aaeb5f65d907944d7aa48c27648be3d0371bd97a9c060e8ef4f573521b8_cppui_modular255,
                            0x52ea7c3f75cba07991674295c4e1462108401b9a103736623943d42e4fbe334e_cppui_modular255,
                            0x1106537bf3150b442b0992ee517b69707c3042015e938f97a63d5c924e67f677_cppui_modular255,
                            0x71de967042516a5b990ef18ae9956fab89f361b950e0639963017c237ee2a0cf_cppui_modular255,
                            0x664a4487e02f7bfa07a1db6ab94a0d1ed0f9e74002bde9cfcbb65f6f74dbfca0_cppui_modular255,
                            0x1023721fd7285260935b5a347f167ce721dd6ae5004c4debc68066bac8f2c467_cppui_modular255,
                            0x2d52fbc95404515f5456c74b65186c860a89dcda8c84bf68fbf715f3d58fe3f2_cppui_modular255,
                            0x6d987c9de419fb6e075441fd99606303e765d8696bcfe01a0d11aa0bd47c8601_cppui_modular255,
                            0x422016ce4d744029b1440a288d7988e43d0f29d616c47f70322ff87cfbc69301_cppui_modular255,
                            0x1f82afe8eb16611abc6600f7dc2a72c8e1d39643c189f3caa1ead08241a896c4_cppui_modular255,
                            0x3bb8684cf815ae6d8a789e0e488c6fb2ac46883fe1cfeb8cfa6f3dbca0f954bd_cppui_modular255,
                            0x3d5a1a6e571306fac431b098cdb3c4518f5a8fc436535766fe9e1bb8bda95d1d_cppui_modular255,
                            0x5e36e175c5d7df42b86285f43b1e4c6bfbaca19f1019073d38d04de0d0647669_cppui_modular255,
                            0x2c3b1b86ce90cb3fe74c5c99b20c3314e28e2f07ce8d932030caee4dfe5055f1_cppui_modular255,
                            0x0bfba44d41c49044bce730d8af86fe0397fff85ec10288b847868d0e9834f754_cppui_modular255,
                            0x0b79924b9e44662369c615cc8d7f36fe4a4b2a79045cee61c413eaf91d82e0c2_cppui_modular255,
                            0x048a11ec75eb154b70223a40cc0db9104b13f6a4ca24e7b9707963ee6f9f74ef_cppui_modular255,
                            0x6dd58a400d366014e46b0b9785ce9d78516813ed2eb329dc4531bfbd8e80eec0_cppui_modular255,
                            0x112844b7c50e7e676b616e72539d5751dec5a063456921b6b16f9e930cc35ebc_cppui_modular255,
                            0x217b616b50e729547af8ceef5008d1edf8d90bc9a7f3ce7c9bc71867e1c06471_cppui_modular255,
                            0x3f9a0b8402ffa291bccbb46dcd2522dea790b35a8503da46717c63917dcb7b79_cppui_modular255,
                            0x42a44fc114c0cad9badf62b911610bdc4b1a0ba9f656f66173a5476e63dfce86_cppui_modular255,
                            0x294223972f4c7e9c9ebefebf059eb90f44479956f5337b12a2eb803e313e96cc_cppui_modular255,
                            0x448101837874eb1bda92bc8a632cbf8f70a0664bbcf3a196609b14c53ee4dbcb_cppui_modular255,
                            0x53a26c6e2b3df0b17faf6a259bc5531d3ae79da59eb8fc5f594e0b886d8d97be_cppui_modular255,
                            0x207c7c32631a75fe8e0da895367176d24e32c5573ec91acf235f3c6c307807cd_cppui_modular255,
                            0x20f955773b13b160d3575eb2380b466f7d38cb4a0e12a15d43d147645c3944ca_cppui_modular255
                    };

                    typedef algebra::matrix<element_type, state_words, state_words> mds_matrix_type;
                    constexpr static const mds_matrix_type mds_matrix = {
                            0x354423b163d1078b0dd645be56316e34a9b98e52dcf9f469be44b108be46c107_cppui_modular255,
                            0x44778737e8bc1154aca1cd92054a1e5b83808403705f7d54da88bbd1920e1053_cppui_modular255,
                            0x5872eefb5ab6b2946556524168a2aebb69afd513a2fff91e50167b1f6e4055e0_cppui_modular255,
                            0x43dff85b25129835819bc8c95819f1a34136f6114e900cd3656e1b9e0e13f86a_cppui_modular255,
                            0x07803d2ffe72940596803f244ac090a9cf2d3616546520bc360c7eed0b81cbf8_cppui_modular255,
                            0x45d6bc4b818e2b9a53e0e2c0a08f70c34167fd8128e05ac800651ddfee0932d1_cppui_modular255,
                            0x08317abbb9e5046b22dfb79e64c8184855107c1d95dddd2b63ca10dddea9ff1a_cppui_modular255,
                            0x1bb80eba77c5dcffafb55ccba4ae39ac8f94a054f2a0ee3006b362f709d5e470_cppui_modular255,
                            0x038e75bdcf8be7fd3a1e844c4de7333531bbd5a8d2c3779627df88e7480e7c5c_cppui_modular255,
                            0x2dd797a699e620ea6b31b91ba3fad4a82f40cffb3e8a30c0b7a546ff69a9002b_cppui_modular255,
                            0x4b906f9ee339b196e958e3541b555b4b53e540a113b2f1cabba627be16eb5608_cppui_modular255,
                            0x605f0c707b82ef287f46431f9241fe4acf0b7ddb151803cbcf1e7bbd27c3e974_cppui_modular255,
                            0x100c514bf38f6ff10df1c83bb428397789cfff7bb0b1280f52343861e8c8737e_cppui_modular255,
                            0x2d40ce8af8a252f5611701c3d6b1e517161d0549ef27f443570c81fcdfe3706b_cppui_modular255,
                            0x3e6418bdf0313f59afc5f40b4450e56881110ea9a0532e8092efb06a12a8b0f1_cppui_modular255,
                            0x71788bf7f6c0cebae5627c5629d012d5fba52428d1f25cdaa0a7434e70e014d0_cppui_modular255,
                            0x55cc73296f7e7d26d10b9339721d7983ca06145675255025ab00b34342557db7_cppui_modular255,
                            0x0f043b29be2def73a6c6ec92168ea4b47bc9f434a5e6b5d48677670a7ca4d285_cppui_modular255,
                            0x62ccc9cdfed859a610f103d74ea04dec0f6874a9b36f3b4e9b47fd73368d45b4_cppui_modular255,
                            0x55fb349dd6200b34eaba53a67e74f47d08e473da139dc47e44df50a26423d2d1_cppui_modular255,
                            0x45bfbe5ed2f4a01c13b15f20bba00ff577b1154a81b3f318a6aff86369a66735_cppui_modular255,
                            0x6a008906685587af05dce9ad2c65ea1d42b1ec32609597bd00c01f58443329ef_cppui_modular255,
                            0x004feebd0dbdb9b71176a1d43c9eb495e16419382cdf7864e4bce7b37440cd58_cppui_modular255,
                            0x09f080180ce23a5aef3a07e60b28ffeb2cf1771aefbc565c2a3059b39ed82f43_cppui_modular255,
                            0x2f7126ddc54648ab6d02493dbe9907f29f4ef3967ad8cd609f0d9467e1694607_cppui_modular255
                    };

                    typedef std::pair<mds_matrix_type, round_constants_type> constants_type;
                    constexpr static const constants_type constants = {mds_matrix, round_constants};
                };
            } // namespace detail
        } // namespace hashes
    } // namespace crypto3
} // namespace nil

#endif    // CRYPTO3_HASH_POSEIDON_CONSTANTS_HPP
