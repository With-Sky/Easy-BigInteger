#include "PrimeNumberTester.hpp"
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

#include "PrimeNumberTester.hpp"

namespace TwilightDream
{
	AKS_Test::Polynomial::Polynomial(size_t degree_count)
	{
		if(degree_count < 1)
		{
			degree_count = 1;
		}

		coefficients.resize(degree_count, BigInteger(0));
	}

	AKS_Test::Polynomial::Polynomial(const BigInteger& coefficient) 
		: coefficients(1, coefficient)
	{
				
	}

	AKS_Test::Polynomial::Polynomial(const std::vector<BigInteger>& coefficient)
		: coefficients(coefficient) 
	{}

	AKS_Test::Polynomial::Polynomial(const Polynomial& other) noexcept
		: coefficients(other.coefficients)
	{}

	AKS_Test::Polynomial::Polynomial(Polynomial &&other) noexcept
		: coefficients(std::move(other.coefficients)) {}

	AKS_Test::Polynomial& AKS_Test::Polynomial::operator=( const Polynomial& other ) noexcept
	{
		if(this == &other)
			return *this;

		coefficients = other.coefficients;

		return *this;
	}

	AKS_Test::Polynomial& AKS_Test::Polynomial::operator=( Polynomial&& other ) noexcept
	{
		if(this == &other)
			return *this;

		coefficients = std::move(other.coefficients);

		return *this;
	}

	// Getting the subscript of the highest sub-item
	size_t AKS_Test::Polynomial::degree() const
	{
		size_t degree = coefficients.size() - 1;
		while (degree > 0 && coefficients[degree].IsZero())
		{
			degree--;
		}
		return (size_t)degree;
	}

	BigInteger::BigInteger AKS_Test::Polynomial::GetCoefficient(size_t i) const
	{
		if (i < coefficients.size())
			return coefficients[i];
				
		return 0;
	}

	void AKS_Test::Polynomial::SetCoefficient(const BigInteger &c, size_t i)
	{
		if (i >= coefficients.size())
			coefficients.resize(i + 1, BigInteger(0ULL));

		coefficients[i] = c;
	}

	AKS_Test::Polynomial AKS_Test::Polynomial::MultiplyModulo(const Polynomial& other, const BigInteger& modulo, uint64_t poly_modulo )
	{
		Polynomial a = *this;
		Polynomial b = other;
		Polynomial result = *this;
		result.Clear();

		size_t a_degree = a.degree();
		size_t b_degree = b.degree();
		size_t degree = a_degree > b_degree ? a_degree : b_degree;

		uint64_t j = 0;
		for (; j <= poly_modulo; j++)
		{
			BigInteger sum = 0;
			uint64_t i = 0;
			for(; i <= j; i++)
			{
				sum = sum + (a.GetCoefficient(i) * b.GetCoefficient(j - i) + b.GetCoefficient(j + poly_modulo - i));
			}

			for(; i <= j + poly_modulo; i++)
			{
				sum = sum + (a.GetCoefficient(i) * b.GetCoefficient(j + poly_modulo - i));
			}

			result.SetCoefficient(sum % modulo, j);
			
			if(j > degree && sum.IsZero())
				break;
		}
		result.Compact();
		return result;
	}

	AKS_Test::Polynomial AKS_Test::Polynomial::PowerModulo(const BigInteger& e, const BigInteger& modulo, uint64_t poly_modulo )
	{
		Polynomial result = *this;
		result.Clear();
		result.SetCoefficient(BigInteger(1), 0);

		uint64_t i = e.BitLength();
		for (; i>= 0; i--)
		{
			result = result.MultiplyModulo(result, modulo, poly_modulo);

			if (e.GetBit(i))
			{
				result = result.MultiplyModulo(*this, modulo, poly_modulo);
			}

			if(i == 0)
				break;
		}
		result.Compact();
		return result;
	}

	bool AKS_Test::Polynomial::operator==(const Polynomial& other) const 
	{
		if(coefficients.size() != other.coefficients.size())
			return false;
		if(degree() != other.degree())
			return false;

		return std::equal(coefficients.begin(), coefficients.end(), other.coefficients.begin(), other.coefficients.end());
	}

	bool AKS_Test::Polynomial::operator!=(const Polynomial& other) const 
	{
		return !(*this == other);
	}

	void AKS_Test::Polynomial::Compact()
	{
		size_t degree_count = this->degree();
		for(; degree_count > 0; degree_count--)
			if(!coefficients[degree_count].IsZero())
				break;
			else
				coefficients.pop_back();
	}

	void AKS_Test::Polynomial::Clear()
	{
		coefficients.clear();
		coefficients.shrink_to_fit();
		coefficients.emplace_back(0);
	}

	void AKS_Test::CheckValueA(BigInteger lower_bound, const BigInteger& r, size_t index, const BigInteger N, std::atomic_bool& result_signal)
	{
		if(!result_signal.load(std::memory_order_relaxed))
			return;

		// 5.2.1 lhs = (X + a)^n mod (X^r-1, n)
		Polynomial base(r);
		base.SetCoefficient(lower_bound % N, 0);
		base.SetCoefficient(1, 1);
		Polynomial lhs = base.PowerModulo(N, N, r.ToUnsignedInt());

		// 5.2.2 rhs = X^n + a mod (X^r-1, n)
		// n mod r => such that m = n mod r
		Polynomial rhs = r;

		// Coefficients of rhs: x^m = 1
		rhs.SetCoefficient(lower_bound % N, 0);
		rhs.SetCoefficient(1, index);

		// 5.2.3 Comparing lhs == rhs
		if( lhs != rhs )
		{
			result_signal.store(false, std::memory_order_relaxed);
		}
	}

	bool AKS_Test::IsPerfectPower( const BigInteger& N )
	{
		if ( N <= 1 )
			return false;

		#if true

		size_t max_power_count = N.BitLength() + 2;
		BigInteger phi = N;
		BigInteger low, high;
		for (size_t exponent = 2; exponent <= max_power_count; exponent++) {
			low = 2, high = phi;
			while (low <= high) {
				BigInteger middle = low + ((high - low) >> 1);
				BigInteger powered_value = middle.Power(exponent);
				if (powered_value == N)
					return true;
				if (powered_value >= N)
					high = middle - 1;
				else
					low = middle + 1;
			}
			phi = high;
		}

		#else

		BigInteger maxExponent = N.Log2() + 1;

		for (BigInteger b = 2; b < maxExponent; ++b)
		{
			BigFraction exponent(1, b, 1); // 表示 1/b
			BigFraction base(N, 1, 1);      // 表示 N
			BigFraction lower_bound = base.Power(exponent); // 计算 N^(1/b)

			if (!lower_bound.IsNaN() && lower_bound.IsInteger() == 1)
			{
				return true;
			}
		}

		#endif

		return false;
	}

	bool AKS_Test::operator()( const BigInteger& n )
	{
		// Step 1: Check if n is a power of b integer.
		// If (n == a^{b} for a ∈ N and b > 1), output COMPOSITE.
		if ( IsPerfectPower( n ) )
		{
			return false;
		}

		const BigInteger ZERO( 0 );
		const BigInteger ONE( 1 );
		const BigInteger TWO( 2 );

		const BigInteger n_minus_one = n - ONE;

		// Step 2: Find the smallest r such that order_{r}(n) > log_{2}(n).
		BigInteger log2_n = n.Log2();
		BigInteger r(1U);
		while (true)
		{
			r += ONE;
			bool found_equal_one = false;
			for (BigInteger k = ONE; k <= log2_n; k++)
			{
				// compute n^k mod r
				BigInteger copy_n = n;
				copy_n.PowerWithModulo(k, r);
				if (copy_n == ONE)
				{
					// Real: order_r(n) <= k
					// Expected: order_{r}(n) > log2(n)
					found_equal_one = true;
					break;
				}
			}
			if (!found_equal_one)
			{
				// Show that for k=1 ... log2(n), none of the n^k mod r == 1
				// Real: order_{r}(n) > log2(n)
				break;
			}
		}
		//Now Finded r

		// Step 3: a belongs to the range 1~r plus 1 , if the greatest common factor of a and n is greater than 1 and less than n, n is a composite number.
		// If 1 < gcd(a, n) < n for some a ≤ r, a ∈ Range(1, r + 1), output COMPOSITE.
		for ( BigInteger a = ONE; a <= r; a += ONE )
		{
			BigInteger value = BigInteger::GCD( a, n );
			if ( value > ONE && value < n )
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

#if false

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
		// Step 5.1 Since the Euler's totient function is computed here, `r` should preferably be a prime number
		BigInteger phi_r = 0;
		if ( r.IsEven() )
		{
			for ( BigInteger k = ONE; k <= r; k += ONE )
			{
				if ( BigInteger::GCD( k, r ) == ONE )
					phi_r += ONE;
			}
		}
		else
		{
			// φ(r) = r - 1, if r is prime
			// Calculate the Euler's totient function φ(r), assuming r is odd.
			phi_r = r - ONE;
		}

		// Step 5.2 compute upper limit = floor( sqrt(phi_r) ) * floor( log(n) )
		BigInteger sqrt_phi_r = phi_r;
		sqrt_phi_r.Sqrt();
		BigInteger upper_bound = sqrt_phi_r * log2_n;
		
		if ( upper_bound > n )
			upper_bound = n;
		
		if(r.Size() > 1)
		{
			std::cout << "Impossible r convert to subscripts!" << std::endl;
			throw std::out_of_range("Impossible r convert to subscripts!");
		}

		BigInteger m = n % r;
		if(m.Size() > 1)
		{
			std::cout << "Impossible m convert to subscripts!" << std::endl;
			throw std::out_of_range("Impossible m convert to subscripts!");
		}
		size_t index = m.ToUnsignedInt();

		// Step 5.3 Compute the upper bound for the range of a(lower_bound)
		
		// 用一个原子布尔来标记是否仍然“全部通过” => true=暂未发现问题
		std::atomic_bool result_signal(true);

		// 先准备存放 future 的容器
		std::vector<std::future<void>> futures;
		futures.reserve(upper_bound.ToUnsignedInt()); // 仅当 upperBound也不太大

		// 遍历 lower_bound = 1 ... upper_bound，给每个 lower_bound 提交一个异步任务
		// 注意：若 upper_bound 非常大，也可能需要分批提交
		for(BigInteger lower_bound = ONE; lower_bound <= upper_bound; ++lower_bound)
		{
			// 提交任务
			auto future = std::async(std::launch::async,
									&AKS_Test::CheckValueA,
									this,             // 成员函数，第一参数是this
									lower_bound, 
									std::cref(r),     // 传引用要用std::ref or std::cref
									(size_t)index, 
									std::cref(n),
									std::ref(result_signal));

			// 收集 future
			futures.push_back(std::move(future));

			if(futures.size() >= std::thread::hardware_concurrency())
			{
				for(auto &f : futures)
				{
					f.get(); // 如果函数里抛异常，这里也会 rethrow
				}

				futures.clear(); // 清空已完成的 futures
			}
		}

		for(auto &f : futures)
		{
			f.get();
		}

		//all test passed, that n is prime
		// 根据 result_signal 判断是否lower_bound全部通过
		return result_signal.load(std::memory_order_relaxed);

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
			const BigInteger a = BigInteger::RandomGenerateNBit( n.BitLength() + 1 );
			x = a % RangeInteger + TWO;

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
		BigInteger		 x = ZERO;

		auto generate_random_number_with_range = []( const BigInteger& start, const BigInteger& end )
		{
			BigInteger range = end - start;
			BigInteger result = BigInteger::RandomGenerateNBit( range.BitLength() + 1 );
			return result % range + start;
		};

		// Decompose (n - 1) to write it as (2 ** s) * d
		// While d is even, divide it by 2 and increase the exponent.
		BigInteger s = n_minus_one;
		size_t	   t = s.CountTrailingZeros();
		s >>= t;

		// Perform Miller-Rabin test with Montgomery multiplication
		// Test k witnesses.

		BigInteger range_integer = n_minus_one;
		TwilightDream::BigInteger::Montgomery montgomery( n );
		for ( size_t count = 0; count < k; count++ )
		{
			// Generate random integer a, where 2 <= a <= (n - 2)
			BigInteger x = generate_random_number_with_range( 2, range_integer );
			x = montgomery.Power( x, s );
			if ( x != 1 )
			{
				size_t i = 0;
				while ( x != range_integer )
				{
					if ( i == t - 1 )
					{
						return false;
					}
					else
					{
						i++;
						x = x * x % n;
					}
				}
			}
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

		size_t bit_size = Number.BitLength();
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