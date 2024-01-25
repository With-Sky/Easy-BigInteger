#pragma once

#include "BigInteger.hpp"
#include "BigFraction.hpp"

namespace TwilightDream
{
	class AKS_Test
	{
	private:
		using BigInteger = TwilightDream::BigInteger::BigInteger;
		using BigFraction = TwilightDream::BigFraction::BigFraction;
		using HashFunction = TwilightDream::BigInteger::HashFunction;

		bool CheckValueA(const BigInteger& n, BigInteger& bottom, const BigInteger& top, std::unordered_map<BigInteger, bool, HashFunction>& test_result_map);

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