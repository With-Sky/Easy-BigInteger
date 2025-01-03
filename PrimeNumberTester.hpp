/*
MIT License

Copyright (c) 2024-2050 Twilight-Dream & With-Sky

https://github.com/Twilight-Dream-Of-Magic/
https://github.com/With-Sky

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#ifndef TWILIGHT_DREAM_PRIME_NUMBER_TESTER_HPP
#define TWILIGHT_DREAM_PRIME_NUMBER_TESTER_HPP

#include "BigInteger.hpp"
#include "BigFraction.hpp"
#include <thread>
#include <future>

namespace TwilightDream
{
	class AKS_Test
	{
	private:
		using BigInteger = TwilightDream::BigInteger::BigInteger;
		using BigFraction = TwilightDream::BigFraction::BigFraction;
		using HashFunction = TwilightDream::BigInteger::HashFunction;

		struct Polynomial
		{
			// Coefficient vector: coefficient[i] denotes the coefficient of x^i
			std::vector<TwilightDream::BigInteger::BigInteger> coefficients;

			Polynomial(size_t degree_count = 1);

			Polynomial(const BigInteger& coefficient);

			Polynomial(const std::vector<BigInteger>& coefficient);
			
			Polynomial(const Polynomial& other) noexcept;

			Polynomial(Polynomial &&other) noexcept;

			Polynomial& operator=(const Polynomial& other) noexcept;

			Polynomial& operator=(Polynomial &&other) noexcept;

			// Getting the subscript of the highest sub-item
			size_t degree() const;

			// 获取和设置系数
			BigInteger GetCoefficient(size_t i) const;

			void SetCoefficient(const BigInteger &c, size_t i);

			Polynomial MultiplyModulo(const Polynomial& other, const BigInteger& modulo, uint64_t poly_modulo );

			Polynomial PowerModulo(const BigInteger& exp, const BigInteger& modulo, uint64_t poly_modulo );

			bool operator==(const Polynomial& other) const;

			bool operator!=(const Polynomial& other) const;

			void Compact();

			void Clear();
		};

		void CheckValueA(BigInteger lower_bound, const BigInteger& r, size_t index, const BigInteger N, std::atomic_bool& result_signal);

		bool IsPerfectPower( const BigInteger& N );
	public:
		/**
		* Performs the Agrawal–Kayal–Saxena primality test on a BigInteger.
		*
		* The AKS primality test is a deterministic algorithm for proving the primality of a number.
		* It is based on the following polynomial congruences:
		*	Input: integer n > 1.
		*		1. If (n == a^{b} for a ∈ N and b > 1)
					output COMPOSITE.
		*		2. Find the smallest r such that order_r(n) > log_{2}(n).
		*		3. If 1 < gcd(a, n) < n for some a ≤ r
					output COMPOSITE.
		*		4. If n ≤ r, output PRIME.
		*		5. For a = 1 to floor(sqrt(φ(r))) * log( n) do
		*			if ((X + a)^{n} != X^{n} + a mod (X^{r} − 1, n)), output COMPOSITE;
		*		6. Output PRIME;
		*
		* @param n The BigInteger to test for primality.
		* @return true if the BigInteger is prime, false otherwise.
		* @see https://www.cse.iitk.ac.in/users/manindra/algebra/primality_v6.pdf
		*/
		bool operator()( const BigInteger& n );
	};

	inline AKS_Test AKS_Test_Instance;

	struct PrimeNumberTester
	{
		using BigInteger = TwilightDream::BigInteger::BigInteger;

		/**
		 * @brief Miller-Rabin primality test for integer `n`.
		 *
		 * This function applies the Miller-Rabin primality test to check if the given
		 * integer `n` is likely prime. The test is repeated `k` times for increased accuracy.
		 *
		 * @param n The integer to test for primality.
		 * @param k The number of iterations for the Miller-Rabin test.
		 * @return True if `n` is likely prime, false otherwise.
		 */
		bool MillerRabin( const BigInteger& n, int k );

		/**
		 * @brief Miller-Rabin primality test with Montgomery multiplication for integer `n`.
		 *
		 * This function applies the Miller-Rabin primality test using Montgomery multiplication
		 * to check if the given integer `n` is likely prime. The test is repeated `k` times for
		 * increased accuracy.
		 *
		 * @param n The integer to test for primality.
		 * @param k The number of iterations for the Miller-Rabin test.
		 * @return True if `n` is likely prime, false otherwise.
		 */
		bool MillerRabinWithMontgomery( const BigInteger& n, int k );

		bool IsPrime_SlowAlgorithm(const BigInteger& Number);

		bool IsPrime_FastAlgorithm(const BigInteger& Number);

		bool IsPrime(const BigInteger& Number);
	};
}

#endif