//---------------------------------------------------------------------------//
// Copyright (c) 2021 Mikhail Komarov <nemo@nil.foundation>
// Copyright (c) 2021-2023 Nikita Kaskov <nbering@nil.foundation>
// Copyright (c) 2022-2023 Polina Chernyshova <pockvokhbtra@nil.foundation>
// Copyright (c) 2022-2023 Elena Tatuzova <e.tatuzova@nil.foundation>
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
// @file Declaration of interfaces for PLONK unified addition component.
//---------------------------------------------------------------------------//

#ifndef CRYPTO3_TRANSPILER_PUBLIC_INPUT_HPP
#define CRYPTO3_TRANSPILER_PUBLIC_INPUT_HPP



namespace nil {
<<<<<<< HEAD:libs/hash/include/nil/crypto3/hash/detail/pedersen/lookup.hpp
    namespace crypto3 {
        namespace hashes {
            namespace detail {
                template<typename ResultT, std::size_t ChunkBits = 3>
                struct lookup;

                template<typename ResultT>
                struct lookup<ResultT, 3> {
                    typedef ResultT result_type;
                    static constexpr std::size_t chunk_bits = 3;

                    template<typename BitRange>
                    static inline result_type process(const BitRange &bits) {
                        int result = (1 - 2 * bits[2]) * (1 + bits[0] + 2 * bits[1]);
                        if (result > 0)
                            return static_cast<result_type>(unsigned(result));
                        return result_type::modulus - unsigned(-result);
                    }
                };
=======
    namespace blueprint {
        template<typename FieldType> 
        std::string convert_numeric_public_input_to_json(std::string path){
            std::ifstream in(path);
            std::stringstream result;

            result << "[" << std::endl;
            std::size_t i = 0;
            if( i != 0 ){
                result << "," << std::endl;
            }
            for (std::string line; std::getline(in, line); i++) 
            {
                if( i != 0 ){
                    result << "," << std::endl;
                }
                result << line;
>>>>>>> nil/master:libs/transpiler/include/nil/blueprint/transpiler/public_input.hpp
            }
            result << std::endl << "]";
            
            return result.str();
        }
    }
}

#endif   // CRYPTO3_MINIMIZED_PROFILING_PLONK_CIRCUIT_HPP