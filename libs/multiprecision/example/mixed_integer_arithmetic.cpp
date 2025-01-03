///////////////////////////////////////////////////////////////
//  Copyright 2012 John Maddock. Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at https://www.boost.org/LICENSE_1_0.txt

//
// Compare arithmetic results using fixed_int to GMP results.
//

#ifdef _MSC_VER
#  define _SCL_SECURE_NO_WARNINGS
#endif

//[mixed_eg
#include <nil/crypto3/multiprecision/cpp_int_modular.hpp>

int main()
{
   using namespace boost::multiprecision;

   boost::uint64_t i = (std::numeric_limits<boost::uint64_t>::max)();
   boost::uint64_t j = 1;

   uint128_modular_t ui128;
   uint256_t ui256;
   //
   // Start by performing arithmetic on 64-bit integers to yield 128-bit results:
   //
   std::cout << std::hex << std::showbase << i << std::endl;
   std::cout << std::hex << std::showbase << add(ui128, i, j) << std::endl;
   std::cout << std::hex << std::showbase << multiply(ui128, i, i) << std::endl;
   //
   // The try squaring a 128-bit integer to yield a 256-bit result:
   //
   ui128 = (std::numeric_limits<uint128_modular_t>::max)();
   std::cout << std::hex << std::showbase << multiply(ui256, ui128, ui128) << std::endl;

   return 0;
}
//]

/* 

Program output:

//[mixed_output

0xffffffffffffffff
0x10000000000000000
0xFFFFFFFFFFFFFFFE0000000000000001
0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFE00000000000000000000000000000001

//]
*/

