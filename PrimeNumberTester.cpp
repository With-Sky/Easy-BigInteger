#include "PrimeNumberTester.hpp"

namespace TwilightDream
{
	bool AKS_Test::CheckValueA( const BigInteger& n, BigInteger& bottom, const BigInteger& top, std::unordered_map<BigInteger, bool, HashFunction>& test_result_map )
	{
		BigInteger x = bottom / ( top - bottom );
		if ( bottom.IsZero() )
			bottom = 1;

		BigInteger b = 0;
		for ( BigInteger a = bottom; a < top; a++ )
		{
			b = a;
			b.PowerWithModulo( n, n );

			if ( !( ( b - a ).IsZero() ) )
			{
				test_result_map[ x ] = false;
				return false;
			}
		}

		test_result_map[ x ] = true;
		return true;
	}

	bool AKS_Test::operator()( const BigInteger& n )
	{
		// Step 1: Check if n is a power of b integer.
		// If (n == a^{b} for a ∈ N and b > 1), output COMPOSITE.
		if ( BigFraction::IsPerfectPower( n ) )
		{
			return false;
		}

		const BigInteger ZERO( 0 );
		const BigInteger ONE( 1 );
		const BigInteger TWO( 2 );

		const BigInteger n_minus_one = n - ONE;

		// Step 2: Find the smallest r such that order_r(n) > log_{2}(n).
		BigInteger		 log2_n = n.Log2();
		const BigInteger max_k = log2_n.Power( 2 );
		bool			 next_r = true;
		BigInteger		 r = ONE;
		while ( next_r )
		{
			r += ONE;
			next_r = false;
			BigInteger k = ZERO;
			while ( k <= max_k && !next_r )
			{
				k += ONE;
				BigInteger copy_n = n;
				// n ^{k} mod r
				copy_n.PowerWithModulo( k, r );
				if ( copy_n.IsZero() || copy_n == ONE )
				{
					next_r = true;
				}
			}
		}

		// Step 3: a belongs to the range 1~r plus 1 , if the greatest common factor of a and n is greater than 1 and less than n, n is a composite number.
		// If 1 < gcd(a, n) < n for some a ≤ r, a ∈ Range(1, r + 1), output COMPOSITE.
		for ( BigInteger a = ONE; a <= r; a += ONE )
		{
			if ( ONE < BigInteger::GCD( a, n ) < n )
				return false;
		}

		// Step 4: If n is less than or equal to r, then n is a prime number.
		// If n ≤ r, output PRIME.
		if ( n <= r )
			return true;

		// Step 5: check that for every coeficient (a_{i}) in (x-1)^{n} a_{i} mod n == 0
		// https://fishi.devtail.io/weblog/2015/06/25/computing-large-binomial-coefficients-modulo-prime-non-prime/
		// For a = 1, a to floor(sqrt(φ(r))) * log2(n) do
		// {
		//		if ((X + a)^{n} != X^{n} + a (mod X^{r} − 1, n))
		//			output COMPOSITE;
		// }

#if 1

		BigInteger a_i = ONE;
		for ( BigInteger i = ONE; i < n / TWO + ONE; i += ONE )
		{
			a_i = a_i.MultiplyWithModulo( n + ONE - i, n ) / i;

			if ( a_i % n != ZERO )
			{
				return false;
			}
		}

		return true;

#else

		// Since the Euler's totient function is computed here, `r` should preferably be a prime number
		BigInteger phi_r = 0;
		if ( r.IsEven() )
		{
			for ( BigInteger k = ONE; k <= n; k += ONE )
			{
				if ( BigInteger::GCD( k, n ) == ONE )
					phi_r += ONE;
			}
		}
		else
		{
			// φ(r) = r - 1, if r is prime
			// Calculate the Euler's totient function φ(r), assuming r is odd.
			phi_r = r - ONE;
		}

		BigInteger sqrt_phi_r = phi_r;
		sqrt_phi_r.Sqrt();
		BigInteger r_multiply_n = sqrt_phi_r * log2_n;
		if ( r_multiply_n > n )
			r_multiply_n = n;

		std::vector<std::thread>						   threads;
		std::unordered_map<BigInteger, bool, HashFunction> test_result_map;
		BigInteger										   random_value = r_multiply_n / 8;
		if ( random_value.IsZero() )
			random_value = 1;

		for ( BigInteger a = ONE; a < r_multiply_n; a += random_value )
		{
			std::thread thread = std::thread( &AKS_Test::CheckValueA, this, n, a, a + random_value, std::ref( test_result_map ) );
			threads.push_back( std::move( thread ) );
		}

		for ( const auto& [ key, value ] : test_result_map )
		{
			if ( value == false )
			{
				// n is not prime
				return false;
			}
		}

		// n is prime
		return true;

#endif
	}

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
	bool PrimeNumberTester::MillerRabin( const BigInteger& n, int k )
	{
		const BigInteger ZERO( 0 );
		const BigInteger ONE( 1 );
		const BigInteger TWO( 2 );

		if ( n <= ZERO )
			return false;

		if ( n.IsEven() )
			return false;

		const BigInteger n_minus_one = n - ONE;
		const BigInteger n_minus_two = n - TWO;
		BigInteger		 x = ZERO;

		bool TryAgain = false;

		// Decompose (n - 1) to write it as (2 ** s) * d
		// While d is even, divide it by 2 and increase the exponent.
		BigInteger d = n_minus_one;
		size_t	   s = 0;
		while ( d.IsEven() )
		{
			++s;
			d >>= 1;
		}

		const BigInteger RangeInteger = n_minus_two - TWO + ONE;

		// Perform Miller-Rabin test
		// Test k witnesses.
		for ( size_t i = 0; i < k; ++i )
		{
			// Generate random integer a, where 2 <= a <= (n - 2)
			const BigInteger a = BigInteger::RandomGenerateNBit( n.BitSize() );
			x = TWO + ( a % RangeInteger );

			x.PowerWithModulo( d, n );

			if ( x == ONE || x == n_minus_one )
			{
				continue;
			}

			for ( size_t r = 0; r < s; ++r )
			{
				x.PowerWithModulo( TWO, n );

				if ( x == ONE )
				{
					// n is composite.
					return false;
				}
				else if ( x == n_minus_one )
				{
					// Exit inner loop and continue with next witness.
					TryAgain = true;
					break;
				}
			}

			//x != n_minus_one
			if ( !TryAgain )
				return false;
			TryAgain = false;
		}

		//n is *probably* prime
		return true;
	}

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
	bool PrimeNumberTester::MillerRabinWithMontgomery( const BigInteger& n, int k )
	{
		const BigInteger ZERO( 0 );
		const BigInteger ONE( 1 );
		const BigInteger TWO( 2 );

		if ( n <= ZERO )
			return false;

		if ( n.IsEven() )
			return false;

		const BigInteger n_minus_one = n - ONE;
		const BigInteger n_minus_two = n - TWO;
		BigInteger		 x = ZERO;

		bool TryAgain = false;

		// Decompose (n - 1) to write it as (2 ** s) * d
		// While d is even, divide it by 2 and increase the exponent.
		BigInteger d = n_minus_one;
		size_t	   s = 0;
		while ( d.IsEven() )
		{
			++s;
			d >>= 1;
		}

		// Initialize variables for Montgomery multiplication
		size_t	   rsize = n.values.size();
		BigInteger r, rinv = 1, mprime = 0;

		mprime.values.resize( n.values.size() );
		r.values.resize( rsize > 1 ? rsize : 2 );
		r.values[ 1 ] = 1;

		// Set up Montgomery parameters
		for ( size_t i = 0; i < rsize - 1; ++i )
		{
			r <<= TwilightDream::BigInteger::EXPONENT;
		}

		for ( size_t i = 0; i < rsize * TwilightDream::BigInteger::EXPONENT; ++i )
		{
			if ( ( rinv[ 0 ] & 1 ) == 0 )
			{
				rinv >>= 1;
				mprime >>= 1;
			}
			else
			{
				rinv.Add( n );
				rinv >>= 1;
				if ( i != 0 )
					mprime >>= 1;
				mprime.SetBit( rsize * TwilightDream::BigInteger::EXPONENT - 1 );
			}
		}

		const BigInteger RangeInteger = n_minus_two - TWO + ONE;

		// Perform Miller-Rabin test with Montgomery multiplication
		// Test k witnesses.
		for ( size_t i = 0; i < k; ++i )
		{
			// Generate random integer a, where 2 <= a <= (n - 2)
			const BigInteger a = BigInteger::RandomGenerateNBit( n.BitSize() );
			x = TWO + ( a % RangeInteger );

			x.MontgomeryPower( d, n, mprime, r, rsize );

			if ( x == ONE || x == n_minus_one )
			{
				continue;
			}

			for ( size_t r = 0; r < s; ++r )
			{
				x.PowerWithModulo( 2, n );

				if ( x == ONE )
				{
					return false;
				}
				else if ( x == n_minus_one )
				{
					TryAgain = true;
					break;
				}
			}

			//x != n_minus_one
			if ( !TryAgain )
				return false;
			TryAgain = false;
		}

		return true;
	}

	/**********/

	bool PrimeNumberTester::IsPrime_SlowAlgorithm( const BigInteger& Number )
	{
		uint64_t NormalNumber = Number.ToUnsignedInt();

		if ( NormalNumber < 2 )
		{
			return false;
		}

		// Eratosthenes sieve
		std::vector<bool> SieveTable( 10240000, false );

		SieveTable[ 0 ] = SieveTable[ 1 ] = false;

		for ( uint64_t i = 3; i * i <= NormalNumber; i += 2 )
		{
			// If prime[p] is not changed, then it is a prime
			if ( SieveTable[ i / 2 ] == false )
			{
				// Update all multiples of p greater than or equal to the square of it numbers
				// which are multiple of p and are less than p^2 are already been marked.
				for ( uint64_t j = i * 3; j <= NormalNumber; j += 2 * i )
				{
					SieveTable[ j / 2 ] = true;
				}
			}
		}

		return SieveTable[ Number.ToUnsignedInt() ];
	}

	bool PrimeNumberTester::IsPrime_FastAlgorithm( const BigInteger& Number )
	{
		std::vector<BigInteger> SmallPrimes { BigInteger( 2 ), BigInteger( 3 ), BigInteger( 5 ), BigInteger( 7 ), BigInteger( 11 ) };

		// Check for small numbers.
		if ( Number < 13 )
		{
			for ( const auto& SmallPrime : SmallPrimes )
			{
				if ( Number == SmallPrime )
				{
					return true;
				}
			}
			return false;
		}

		size_t bit_size = Number.BitSize();
		/*
			Returns minimum number of rounds for Miller-Rabing primality testing, on number bitsize.

			According to NIST FIPS 186-4, Appendix C, Table C.3, minimum number of rounds of M-R testing
			using an error probability of 2 ** (-100), for different p, q bitsizes are:
			* p, q bitsize: 512; rounds: 7
			* p, q bitsize: 1024; rounds: 4
			* p, q bitsize: 1536; rounds: 3
			See: http://nvlpubs.nist.gov/nistpubs/FIPS/NIST.FIPS.186-4.pdf
		*/
		auto primality_testing_rounds = [ &bit_size ]( const BigInteger& number ) -> size_t {
			// Set number of rounds.
			if ( bit_size >= 1536 )
				return 3;
			if ( bit_size >= 1024 )
				return 4;
			if ( bit_size >= 512 )
				return 7;
			// For smaller bitsizes, set arbitrary number of rounds.
			return 10;
		};

		if ( bit_size > 4096 )
		{
			return AKS_Test_Instance( Number );
		}
		else if ( bit_size > 2048 )
			return MillerRabinWithMontgomery( Number, primality_testing_rounds( Number ) + 1 );
		else
			return MillerRabin( Number, primality_testing_rounds( Number ) + 1 );
	}

	bool PrimeNumberTester::IsPrime( const BigInteger& Number )
	{
		if ( Number < 10240000 )
		{
			return IsPrime_SlowAlgorithm( Number );
		}
		else
		{
			return IsPrime_FastAlgorithm( Number );
		}
	}
}  // namespace TwilightDream