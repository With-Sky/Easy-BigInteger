//Reference code: https://github.com/NoahBz/Easy-BigInt/blob/main/BigInt.h

#pragma once

#include <iostream>
#include <cmath>
#include <string>
#include <algorithm>
#include <vector>
#include <thread>
#include <cstdint>
#include <intrin.h>
#include <immintrin.h>
#include <random>
#include <bitset>

namespace TwilightDream::BigInteger
{
	using digit_type = uint64_t;

	/*
		Caution!
		Here the BASE constant can only be set to a maximum of 2 to the 32nd power, and the EXPONENT constant is set to 32
		Because we have to reserve half of the range of values as overflow space for each element (digit) of the values array, even though the "digit_type" type we use is a 64-bit number.
		Although this guarantees computational accuracy and is very unlikely to overflow the value range, we have BASE set to 2**32 here, which is actually very extreme.
		So actually the optimal should be set to 2**28, which is 2**28 = 268435456
		Then the same thing.
		The BASE constant can be set to the 64th power of 2 if you want, but the "digit_type" type would have to use a 128 bit number.
		And change the EXPONENT constant to 64, in which case it would actually be optimal to set it to 60, i.e. 64 - 4 = 60
	*/

	// Set max binary bit count the values[index]
	const uint16_t EXPONENT_VALUE = 32;
	const uint16_t EXPONENT_VALUE2 = std::numeric_limits<digit_type>::digits / 2;
	const uint16_t EXPONENT = EXPONENT_VALUE < EXPONENT_VALUE2 ? EXPONENT_VALUE : EXPONENT_VALUE2;

	// Set max value of the values[index]
	// 65536, 536870912, 2147483648, 4294967296, 18446744073709551616
	const uint64_t BASE = digit_type(pow(2, EXPONENT));
	
	constexpr int	   KARATSUBA_LIMIT = 70;
	constexpr int	   BINARY_SEARCH_DIVISION_LIMIT = 4;
	constexpr int	   DONALD_KNUTH_LONG_DIVISION_LIMIT = 20;
	constexpr int	   MULTI_THREAD_LIMIT = 10000;
	constexpr bool	   MULTI_THREAD = true;
	constexpr uint64_t STACK_LIMIT = 4096;

	inline std::random_device						TrueRandomDevice;
	inline std::default_random_engine				RandomEngine( TrueRandomDevice() );
	inline std::uniform_int_distribution<uint64_t> Dist( 0, BASE - 1 );

	static void* MarkAlloct( void* p, uint8_t mark )
	{
		static_cast<uint8_t*>( p )[ 0 ] = mark;
		return static_cast<uint8_t*>( p ) + 32;
	}

	static void _free_alloct( void* p )
	{
		uint8_t* bytes = static_cast<uint8_t*>( p ) - 32;

		if ( bytes[ 0 ] == 1 )
		{
			_aligned_free( bytes );
		}
	}

#define _alloct(size)															\
		__pragma(warning(suppress: 6255))										\
		(size <= STACK_LIMIT													\
		? MarkAlloct((void*)(((uint64_t)_alloca(size + 63) + 31) & ~31), 0)		\
		: MarkAlloct(_aligned_malloc(size + 32, 32), 1))

	class BigInteger
	{
	private:
		void Clean()
		{
			if ( IsZero() )
			{
				values.resize( 1 );
				values[ 0 ] = 0;
				sign = 1;
			}

			size_t count = 0;
			while ( ( values.size() - count ) > 1 && values[ values.size() - ( 1 + count ) ] == 0 )
			{
				count++;
			}
			if ( !count )
			{
				return;
			}

			values.erase( values.end() - count, values.end() );

			if ( IsZero() )
				sign = 1;
		}

		void BitLeftShiftCore( const uint32_t& shift )
		{
			uint64_t k, t;
			k = 0;
			for ( size_t i = 0; i < values.size(); ++i )
			{
				// Shift the current value to the right and store the carry
				t = ( uint64_t )values[ i ] >> ( EXPONENT - shift );
				// Shift the current value to the left and combine with the carry
				// BASE is the maximum value that each element of the values array can represent, and is used here as a modulus to prevent overflow.
				values[ i ] = ( ( ( uint64_t )values[ i ] << shift ) | k ) & ( BASE - 1 );
				// Update the carry for the next iteration
				k = t;
				// If this is the last value and there is a carry, add a new value to the BigInteger
				if ( i == values.size() - 1 && k != 0 )
					values.push_back( 0 );
			}
		}

		void BitRightShiftCore( const uint32_t& shift )
		{
			uint64_t k, t;
			k = 0;
			for ( int64_t i = values.size() - 1; i >= 0; --i )
			{
				// Shift the current value to the left and store the carry
				t = ( uint64_t )values[ i ] << ( EXPONENT - shift );
				// Shift the current value to the right and combine with the carry
				// BASE is the maximum value that each element of the values array can represent, and is used here as a modulus to prevent overflow.
				values[ i ] = ( ( values[ i ] >> shift ) | k ) & ( BASE - 1 );
				// Update the carry for the next iteration
				k = t;
			}
		}

		std::vector<digit_type> values;

	public:
		int8_t				  sign = 1;

		enum class ArithmeticMode : uint8_t
		{
			Addition = 0,
			Subtraction = 1,
			Multiplication = 2,
			Division = 4
		};

		BigInteger() {}

		template <size_t N>
		BigInteger(const std::bitset<N>& bits)
		{
			size_t digitCount = (N + EXPONENT - 1) / EXPONENT;

			values.resize(digitCount, 0);

			for (size_t i = 0; i < N; ++i)
			{
				if(i >= EXPONENT * values.size())
				{
					break;
				}

				if (bits.test(i))
				{
					SetBit(i);
				}
			}
		}

		BigInteger( int64_t num )
		{
			if ( num == 0 )
			{
				values.resize( 1 );
				values[ 0 ] = 0;
				return;
			}
			else if ( num < 0 )
			{
				num *= -1;
				sign = -1;
			}
			if ( num < BASE )
			{
				values.push_back( num );
				return;
			}

			size_t count = size_t( log10( num ) / log10( BASE ) + 1 );
			values.reserve( count );
			while ( num > 0 )
			{
				values.push_back( num % BASE );
				num /= BASE;
			}
		}

		explicit BigInteger( const std::string& number_string )
		{
			FromString( number_string );
		}

		BigInteger( const std::string& number_string, uint32_t base )
		{
			FromString( number_string, base );
		}

		BigInteger( const BigInteger& num )
		{
			values = num.values;
			sign = num.sign;
		}

		BigInteger( BigInteger&& num ) noexcept
		{
			values = std::move( num.values );
			sign = num.sign;
		}

		BigInteger& operator+=( const BigInteger& other )
		{
			if(this->IsZero() && other.sign == 1)
			{
				*this = other;
				return *this;
			}
			else if(other.IsZero())
			{
				return *this;
			}

			if ( sign != other.sign )
			{
				sign = ( sign > other.sign || *this < other ? 1 : -1 );
				if ( *this >= other )
				{
					return Subtract( other );
				}
				else
				{
					BigInteger res = other;
					res.Subtract( *this );
					values = std::move( res.values );
					return *this;
				}
			}

			return Add( other );
		}


		/**
		 * @brief Adds another BigInteger to the current instance.
		 *
		 * This function performs addition of two BigIntegers and updates the current
		 * instance with the result.
		 *
		 * @param other The BigInteger to be added.
		 * @return Reference to the modified current instance after addition.
		 */
		BigInteger& Add( const BigInteger& other )
		{
			// Determine the maximum size for the loop and result storage
			const size_t thisSize = values.size();
			const size_t otherSize = other.values.size();
			const size_t max = ( thisSize > otherSize ? thisSize : otherSize ) + 1;

			// Resize the storage for the result
			values.resize( max );

			// Initialize carry and perform addition
			uint8_t carry = 0;
			for ( size_t i = 0; i < max; ++i )
			{
				uint64_t sum = values[ i ] + ( ( i < otherSize ) ? other.values[ i ] : 0 ) + carry;
				values[ i ] = sum & ( BASE - 1 );  // Update current digit with the sum
				carry = sum >> EXPONENT;		   // Update carry for the next iteration
			}

			// Remove leading zeros and resize the storage
			Clean();

			// Return a reference to the modified current instance
			return *this;
		}

		friend BigInteger operator+( const BigInteger& lhs, const BigInteger& rhs )
		{
			BigInteger result( lhs );
			result += rhs;
			return result;
		}

		BigInteger& operator-=( const BigInteger& other )
		{
			if(this->IsZero() && other.sign == -1)
			{
				*this = other;
				this->sign = 1;
				return *this;
			}
			else if(other.IsZero())
			{
				return *this;
			}

			if ( sign == other.sign )
			{
				sign = ( *this > other ? 1 : -1 );
				sign *= other.sign;

				Difference( other );
			}
			else
			{
				Add( other );
			}

			return *this;
		}

		BigInteger& Difference( const BigInteger& other )
		{
			size_t	   max = values.size();
			BigInteger a, b;
			a.values = values;
			b.values = other.values;

			if ( b > a )
			{
				max = other.values.size();
				std::swap( a, b );
				values.resize( max );
			}

			a.Subtract( b );

			values = std::move( a.values );

			return *this;
		}

		/**
		 * @brief Subtracts another BigInteger from the current instance.
		 *
		 * This function performs subtraction of two BigIntegers and updates the current
		 * instance with the result.
		 *
		 * @param other The BigInteger to be subtracted.
		 * @return Reference to the modified current instance after subtraction.
		 */
		BigInteger& Subtract( const BigInteger& other )
		{
			// Determine the size of current and other BigIntegers
			const size_t thisSize = values.size();
			const size_t otherSize = other.values.size();

			// Initialize variables for borrow and temporary subtraction
			int64_t borrow = 0;
			int64_t temp;

			// Perform subtraction for each digit
			for ( size_t i = 0; i < thisSize; ++i )
			{
				temp = values[ i ] - ( ( i >= otherSize ) ? 0 : other.values[ i ] ) + borrow;
				values[ i ] = temp & ( BASE - 1 );	// Update current digit with the result
				borrow = temp >> EXPONENT;			// Update borrow for the next iteration
			}

			// Remove leading zeros and resize the storage
			Clean();

			// Return a reference to the modified current instance
			return *this;
		}

		friend BigInteger operator-( const BigInteger& lhs, const BigInteger& rhs )
		{
			BigInteger result( lhs );
			result -= rhs;
			return result;
		}

		BigInteger operator++()
		{
			*this += BigInteger(1);
			return *this;
		}
		BigInteger operator--()
		{
			*this -= BigInteger(1);
			return *this;
		}

		BigInteger operator++( int )
		{
			BigInteger copy = *this;
			copy += BigInteger(1);
			return copy;
		}
		BigInteger operator--( int )
		{
			BigInteger copy = *this;
			copy -= BigInteger(1);
			return copy;
		}

		BigInteger& operator*=( const BigInteger& other )
		{
			if ( IsZero() || other.IsZero() )
			{
				sign = 1;
				values.resize( 1 );
				values[ 0 ] = 0;
				return *this;
			}
			sign *= other.sign;

#ifdef __AVX2__
			if ( ( values.size() > 254 && other.values.size() > 254 ) )
				return KaratsubaMultiplication( other );
			else
				return *this == other ? Square_I() : BaseMul_I( other );
#else
			if ( ( values.size() > KARATSUBA_LIMIT && other.values.size() > KARATSUBA_LIMIT ) )
				return KaratsubaMultiplication( other );
			else
				return BaseMultiplication( other );
#endif	// __AVX2__
		}

		friend BigInteger operator*( const BigInteger& lhs, const BigInteger& rhs )
		{
			BigInteger result( lhs );
			result *= rhs;
			return result;
		}


		/**
		 * @brief Performs base multiplication of two BigIntegers.
		 *
		 * This function multiplies the current instance by another BigInteger using
		 * base multiplication and updates the current instance with the result.
		 *
		 * @param other The BigInteger to be multiplied.
		 * @return Reference to the modified current instance after multiplication.
		 */
		BigInteger& BaseMultiplication( const BigInteger& other )
		{
			const size_t thisSize = values.size();
			const size_t otherSize = other.values.size();
			uint64_t	 carry, sum;
			size_t		 i, j;

			// Make a copy of the current instance's values
			std::vector<digit_type> copy = std::move( values );

			// Resize the storage for the result
			values.resize( thisSize + otherSize );

			// Perform base multiplication
			for ( i = 0; i < otherSize; ++i )
			{
				carry = 0;
				for ( j = 0; j < thisSize; ++j )
				{
					sum = values[ i + j ];
					sum += carry;
					sum += copy[ j ] * other.values[ i ];
					values[ i + j ] = sum & ( BASE - 1 );  // Update current digit with the result
					carry = sum >> EXPONENT;			   // Update carry for the next iteration
				}
				values[ i + thisSize ] = carry;	 // Store the final carry in the next position
			}

			// Remove leading zeros and resize the storage
			Clean();

			// Return a reference to the modified current instance
			return *this;
		}

		/**
		 * @brief Performs base multiplication of two BigIntegers using AVX2 instructions.
		 *
		 * This function multiplies the current instance by another BigInteger using
		 * AVX2 vectorized instructions for improved performance and updates the current
		 * instance with the result.
		 *
		 * @param other The BigInteger to be multiplied.
		 * @return Reference to the modified current instance after multiplication.
		 */
		BigInteger& BaseMul_I( const BigInteger& other )
		{
			const size_t thisSize = values.size();
			const size_t otherSize = other.values.size();
			size_t		 i, j;
			__m256i		 a, b, c;

			// Allocate memory for result and copy of the current instance's values
			uint64_t* result = ( uint64_t* )_alloct( ( thisSize + otherSize ) * sizeof( uint64_t ) );
			uint64_t* copy = ( uint64_t* )_alloct( thisSize * sizeof( uint64_t ) );

			// Initialize memory and copy values
			memset( &result[ 0 ], 0, ( thisSize + otherSize ) * sizeof( uint64_t ) );
			memcpy( &copy[ 0 ], &values[ 0 ], thisSize * sizeof( uint64_t ) );

			// Prefetch memory for better performance
			_mm_prefetch( ( const char* )&result, _MM_HINT_T0 );
			_mm_prefetch( ( const char* )&copy, _MM_HINT_T0 );

			// Main multiplication loop using AVX2 instructions
			for ( i = 0; i < otherSize; ++i )
			{
				b = _mm256_set1_epi32( other.values[ i ] );
				for ( j = 0; j + 3 < thisSize; j += 4 )
				{
					a = _mm256_load_si256( reinterpret_cast<__m256i*>( &copy[ j ] ) );
					c = _mm256_mul_epu32( a, b );

					a = _mm256_load_si256( reinterpret_cast<__m256i*>( &result[ i + j ] ) );
					c = _mm256_add_epi64( a, c );

					_mm256_store_si256( reinterpret_cast<__m256i*>( &result[ i + j ] ), c );
				}
				for ( ; j < thisSize; ++j )
				{
					result[ i + j ] += values[ j ] * other.values[ i ];
				}
			}

			// Finalize and handle carry
			uint64_t carry = 0;
			for ( i = 0; i < thisSize + otherSize; ++i )
			{
				result[ i ] += carry;
				carry = result[ i ];
				result[ i ] &= ( BASE - 1 );
				carry >>= EXPONENT;
			}

			// Resize the storage for the result
			size_t size = thisSize + otherSize - ( result[ thisSize + otherSize - 1 ] == 0 ? 1 : 0 );
			values.resize( size );

			// Copy the result back to the current instance's values
			memcpy( &values[ 0 ], &result[ 0 ], size * sizeof( uint64_t ) );

			// Free allocated memory
			_free_alloct( result );

			// Return a reference to the modified current instance
			return *this;
		}

		/**
		 * @brief Squares the current BigInteger.
		 *
		 * This function squares the current instance by performing a specialized
		 * multiplication and updates the current instance with the result.
		 *
		 * @return Reference to the modified current instance after squaring.
		 */
		BigInteger& Square()
		{
			std::vector<digit_type> result;

			/*
				uv: Represents the sum of a product in the algorithm, essentially the temporary value before splitting into two parts (v and u).
				c: Represents the carry from the product calculations.
				v: Represents the lower part of the temporary value (uv % BASE).
				u: Represents the upper part of the temporary value ((uv - v) / BASE).
			*/
			uint64_t uv, c, v, u;
			size_t	 t = values.size();
			result.resize( t * 2 );

			// Perform specialized multiplication for squaring
			for ( size_t i = 0; i < t; ++i )
			{
				// Compute square and initialize variables
				uv = result[ 2 * i ] + values[ i ] * values[ i ];
				v = uv % BASE;
				u = ( uv - v ) / BASE;
				result[ 2 * i ] = v;
				c = u;

				// Compute cross products and update result
				for ( size_t j = i + 1; j < t; ++j )
				{
					uv = result[ i + j ] + 2 * values[ j ] * values[ i ] + c;
					v = uv % BASE;
					u = ( uv - v ) / BASE;
					result[ i + j ] = v;
					c = u;
				}
				result[ i + t ] = u;
			}

			values = std::move( result );
			Clean();
			return *this;
		}

		/**
		 * @brief Squares the current BigInteger using AVX2 instructions.
		 *
		 * This function squares the current instance by performing a specialized
		 * multiplication using AVX2 vectorized instructions for improved performance.
		 * It updates the current instance with the squared result.
		 *
		 * @return Reference to the modified current instance after squaring.
		 */
		BigInteger& Square_I()
		{
			// Sizes of the current instance's values and squared result
			const size_t n = values.size(), n2 = n << 1;
			size_t		 i, j;
			__m256i		 a, b, c, r;

			// Aligned memory allocation for squared result and a copy of values
			uint64_t* result = ( uint64_t* )_alloca( n2 * sizeof( uint64_t ) + 31 );
			result = ( uint64_t* )( ( ( ( uint64_t )result + 31 ) & ~31 ) );
			memset( &result[ 0 ], 0, n2 * sizeof( uint64_t ) );

			uint64_t* copy = ( uint64_t* )_alloca( n * sizeof( uint64_t ) + 31 );
			copy = ( uint64_t* )( ( ( ( uint64_t )copy + 31 ) & ~31 ) );
			memcpy( &copy[ 0 ], &values[ 0 ], n * sizeof( uint64_t ) );

			// Perform specialized multiplication for squaring using AVX2 instructions
			for ( i = 0; i < n; i += 1 )
			{
				result[ i << 1 ] += copy[ i ] * copy[ i ];
				a = _mm256_set1_epi32( values[ i ] );

				for ( j = i + 1; j + 3 < n; j += 4 )
				{
					b = _mm256_load_si256( reinterpret_cast<__m256i*>( &copy[ j ] ) );

					c = _mm256_mul_epu32( a, b );
					c = _mm256_slli_epi64( c, 1 );

					r = _mm256_load_si256( reinterpret_cast<__m256i*>( &result[ i + j ] ) );
					c = _mm256_add_epi64( r, c );

					_mm256_store_si256( reinterpret_cast<__m256i*>( &result[ i + j ] ), c );
				}
				for ( j; j < n; ++j )
				{
					result[ i + j ] += ( copy[ i ] * copy[ j ] ) << 1;
				}
			}

			// Finalize and handle carry
			uint64_t carry = 0;
			for ( i = 0; i < n2; ++i )
			{
				result[ i ] += carry;
				carry = result[ i ];
				result[ i ] &= BASE - 1;
				carry >>= EXPONENT;
			}

			// Resize the storage for the result
			size_t size = n2 - ( result[ n2 - 1 ] == 0 );
			values.resize( size );

			// Copy the result back to the current instance's values
			memcpy( &values[ 0 ], &result[ 0 ], size * sizeof( uint64_t ) );

			// Return a reference to the modified current instance
			return *this;
		}

		BigInteger& MultiplyBase( size_t times )
		{
			values.insert( values.begin(), times, 0 );
			return *this;
		}

		/**
		 * @brief Performs Karatsuba multiplication of two BigIntegers.
		 *
		 * This function multiplies the current instance by another BigInteger using
		 * the Karatsuba algorithm, a divide-and-conquer approach for efficient
		 * multiplication of large numbers. It updates the current instance with the result.
		 *
		 * @param other The BigInteger to be multiplied.
		 * @return Reference to the modified current instance after multiplication.
		 */
		BigInteger& KaratsubaMultiplication( const BigInteger& other )
		{
			// Determine the size for splitting and create temporary BigIntegers
			const size_t m = ( values.size() < other.values.size() ? values.size() : other.values.size() );
			const size_t m2 = m >> 1;
			BigInteger	 high1, low1, high2, low2;
			BigInteger	 z0, z1, z2;

			// Split the BigIntegers into high and low parts
			SplitAt( *this, m2, high1, low1 );
			SplitAt( other, m2, high2, low2 );

			static bool workerId0, workerId1, workerId2;

			// Perform Karatsuba multiplication
			if ( !MULTI_THREAD || m < MULTI_THREAD_LIMIT || ( workerId0 || workerId1 || workerId2 ) )
			{
				z0 = low1 * low2;
				z1 = ( low1 + high1 ) * ( low2 + high2 );
				z2 = high1 * high2;
			}
			else
			{
				// Use multi-threading for parallel computation
				workerId0 = true;
				workerId1 = true;
				workerId2 = true;

				auto lambda0 = [ & ] {
					z0 = low1 * low2;
					workerId0 = false;
				};
				auto lambda1 = [ & ] {
					z1 = ( low1 + high1 ) * ( low2 + high2 );
					workerId1 = false;
				};
				auto lambda2 = [ & ] {
					z2 = high1 * high2;
					workerId2 = false;
				};

				std::thread worker0( lambda0 );
				std::thread worker1( lambda1 );
				std::thread worker2( lambda2 );

				worker1.join();
				worker0.join();
				worker2.join();
			}

			// Combine the partial results
			z1.Subtract( z2 );
			z1.Subtract( z0 );
			z2.MultiplyBase( m2 << 1 );
			z1.MultiplyBase( m2 );
			z2.Add( z1 );
			z2.Add( z0 );

			// Update the current instance with the result
			values = std::move( z2.values );

			return *this;
		}

		void SplitAt( const BigInteger& num, const size_t n, BigInteger& high, BigInteger& low )
		{
			std::vector<digit_type> lowValue( num.values.begin(), num.values.begin() + n );
			std::vector<digit_type> highValue( num.values.end() - ( num.values.size() - n ), num.values.end() );

			low.values = std::move( lowValue );
			high.values = std::move( highValue );
		}

		BigInteger& Power( const size_t p )
		{
			if ( p == 0 )
			{
				*this = 1;
				return *this;
			}

			BigInteger result( *this );
			for ( size_t i = 0; i < p - 1; ++i )
			{
				result *= *this;
			}
			*this = result;
			return *this;
		}

		BigInteger& Sqrt()
		{
			BigInteger s = 0, t, u = *this;
			do
			{
				s = u;
				t = s + ( *this ) / s;
				u = t >> ( uint8_t )1;
			} while ( u < s );

			values = std::move( s.values );
			return *this;
		}

		BigInteger& PowerWithModulo( const BigInteger& e, const BigInteger& m )
		{
			size_t	   bitSize = e.BitSize();
			BigInteger result = *this % m;

			for ( size_t i = 1; i < bitSize; ++i )
			{
				result *= result;
				result %= m;
				if ( e.GetBit( bitSize - ( i + 1 ) ) )
				{
					result *= *this;
					result %= m;
				}
			}

			values = std::move( result.values );
			return *this;
		}

		//Use by MontgomeryReduce
		BigInteger& MontDivR( size_t rsize )
		{
			if ( values.size() > rsize )
			{
				values.erase( values.begin(), values.begin() + rsize );
			}
			else
			{
				*this = 0;
			}
			return *this;
		}

		//Use by MontgomeryReduce
		BigInteger& MontModR( size_t rsize )
		{
			if ( values.size() > rsize )
			{
				values.erase( values.end() - ( values.size() - rsize ), values.end() );
			}
			return *this;
		}

		/**
		 * @brief Montgomery reduction operation for modular arithmetic.
		 *
		 * This function performs Montgomery reduction for modular arithmetic.
		 * It reduces the current instance modulo `m` using a precomputed value `mprime`
		 * for optimization. The size of the Montgomery representation is specified by `rsize`.
		 *
		 * @param rsize The size of the Montgomery representation.
		 * @param m The modulus for the reduction.
		 * @param mprime The precomputed Montgomery constant.
		 * @return Reference to the modified current instance after reduction.
		 */
		BigInteger& MontgomeryReduce( size_t rsize, const BigInteger& m, const BigInteger& mprime )
		{
			// Perform Montgomery reduction steps
			BigInteger n = *this;
			n.MontModR( rsize );

			n *= mprime;
			n.MontModR( rsize );

			n *= m;

			Add( n );
			MontDivR( rsize );

			// Subtract m if the result is greater than or equal to m
			if ( *this >= m )
			{
				*this -= m;
			}

			return *this;
		}

		/**
		 * @brief Montgomery exponentiation for modular arithmetic.
		 *
		 * This function performs Montgomery exponentiation, raising the current instance
		 * to the power of `e` modulo `m`. It computes the necessary Montgomery constants
		 * `r`, `rinv`, and `mprime` before applying the exponentiation.
		 *
		 * @param e The exponent.
		 * @param m The modulus.
		 * @return Reference to the modified current instance after exponentiation.
		 */
		BigInteger& MontgomeryPower( const BigInteger& e, const BigInteger& m )
		{
			// Determine the size of the Montgomery representation
			size_t rsize = m.values.size();

			// Initialize Montgomery constants
			BigInteger r, rinv = 1, mprime = 0;
			mprime.values.resize( m.values.size() );
			r.values.resize( rsize > 1 ? rsize : 2 );
			r.values[ 1 ] = 1;

			// Set up Montgomery constants
			for ( size_t i = 0; i < rsize - 1; ++i )
			{
				r <<= EXPONENT;
			}

			for ( size_t i = 0; i < rsize * EXPONENT; ++i )
			{
				if ( ( rinv[ 0 ] & 1 ) == 0 )
				{
					rinv >>= 1;
					mprime >>= 1;
				}
				else
				{
					rinv.Add( m );
					rinv >>= 1;
					if ( i != 0 )
						mprime >>= 1;
					mprime.SetBit( rsize * EXPONENT - 1 );
				}
			}

			// Perform Montgomery exponentiation using precomputed constants
			MontgomeryPower( e, m, mprime, r, rsize );

			// Return a reference to the modified current instance
			return *this;
		}

		/**
		 * @brief Montgomery exponentiation for modular arithmetic.
		 *
		 * This function performs Montgomery exponentiation, raising the current instance
		 * to the power of `e` modulo `m`. It utilizes a precomputed Montgomery constant `mprime`
		 * and a Montgomery representation of `r` for optimization. The size of the Montgomery
		 * representation is specified by `rsize`.
		 *
		 * @param e The exponent.
		 * @param m The modulus.
		 * @param mprime The precomputed Montgomery constant.
		 * @param r The Montgomery representation.
		 * @param rsize The size of the Montgomery representation.
		 * @return Reference to the modified current instance after exponentiation.
		 */
		BigInteger& MontgomeryPower( const BigInteger& e, const BigInteger& m, const BigInteger& mprime, const BigInteger r, const size_t rsize )
		{
			// Initialize precomputed powers of the base
			const uint8_t			k = ( e.BitSize() < 512 ? 4 : 5 );
			std::vector<BigInteger> g( ( uint64_t )1 << k );

			g[ 0 ] = *this * ( ( r * r ) % m );
			g[ 0 ].MontgomeryReduce( rsize, m, mprime );
			BigInteger g2 = g[ 0 ];

			g2 = g2 * g2;
			g2.MontgomeryReduce( rsize, m, mprime );

			for ( size_t i = 1; i < g.size(); ++i )
			{
				g[ i ] = g[ i - 1 ];
				g[ i ] *= g2;
				g[ i ].MontgomeryReduce( rsize, m, mprime );
			}

			// Perform Montgomery exponentiation
			size_t	   bitSize = e.BitSize();
			BigInteger result = r % m;
			int64_t	   i = bitSize - 1;

			while ( i >= 0 )
			{
				if ( e.GetBit( i ) )
				{
					uint64_t l = ( i - k + 1 > 0 ? i - k + 1 : 0 );
					while ( !e.GetBit( l ) )
					{
						l++;
					}

					for ( size_t j = 0; j < i - l + 1; ++j )
					{
						result *= result;
						result.MontgomeryReduce( rsize, m, mprime );
					}

					uint64_t ndx = 0;
					for ( int64_t j = i; j >= l; --j )
					{
						if ( j < 0 )
						{
							break;
						}
						ndx <<= 1;
						if ( e.GetBit( j ) )
						{
							ndx++;
						}
					}
					ndx >>= 1;
					result *= g[ ndx ];
					result.MontgomeryReduce( rsize, m, mprime );
					i = l - 1;
				}
				else
				{
					result *= result;
					result.MontgomeryReduce( rsize, m, mprime );
					i--;
				}
			}

			// Finalize and update the current instance with the result
			result.MontgomeryReduce( rsize, m, mprime );
			values = std::move( result.values );

			return *this;
		}

		BigInteger& operator/=( const BigInteger& other )
		{
			sign *= other.sign;
			Divide( other );
			return *this;
		}

		BigInteger& Divide( const BigInteger& other, BigInteger* r = nullptr )
		{
			// Handle special case when the numerator is zero
			if (values.size() == 1 && values[0] == 0)
			{
				if (r != nullptr)
				{
					// Return 0 with positive sign if the numerator is zero
					r->sign = 1;
					r->values.resize(1);
					r->values[0] = 0;
				}
				return *this;
			}

			// Select division algorithm based on the size of the divisor
			if (other.values.size() > BINARY_SEARCH_DIVISION_LIMIT)
			{
				// Use BinarySearchDivision for large divisors
				auto pair = BinarySearchDivision(other);
				if (r != nullptr)
				{
					// Store the remainder in the provided pointer
					*r = pair.second;
				}
				*this = pair.first;  // Update the current instance with the quotient
				return *this;
			}
			else if (other.values.size() > DONALD_KNUTH_LONG_DIVISION_LIMIT)
			{
				// Use DonaldKnuthLongDivision for intermediate-sized divisors
				return DonaldKnuthLongDivision(other, r);
			}
			else
			{
				// Use ShortDivision for small divisors
				return ShortDivision(other, r);
			}
		}


		/**
		 * @brief Performs division using binary search based on bit chunks.
		 *
		 * This function implements a binary search division algorithm, which is an efficient
		 * method for performing division on large integers. It iteratively calculates quotient
		 * and remainder by considering bit chunks.
		 *
		 * @param other The divisor.
		 * @return A pair containing the quotient and remainder after division.
		 */
		std::pair<BigInteger, BigInteger> BinarySearchDivision(const BigInteger& other)
		{
			// Make copies of the numerator and denominator
			BigInteger A = *this;
			BigInteger B = other;
			BigInteger Quotient = 0;
			BigInteger Remainder = 0;
			const BigInteger Zero = 0;

			// Lambda function for the main division logic
			auto RunFunction = [&Zero](BigInteger& Remainder, const BigInteger& Divisor)
			{
				digit_type QuotientDigit = 0;

				//Multiplier for Reference Remainder
				BigInteger From = 0;
				//Reference Remainder
				BigInteger&& ToSubtract = std::move(Remainder);

				//Copy Divisor
				BigInteger ToAddition = Divisor;
				// Left shift ToAddition by (EXPONENT - 1) to prepare for binary search
				ToAddition <<= EXPONENT - 1;

				// Binary search loop
				for (digit_type index = digit_type(1) << (EXPONENT - 1); index > 0; index >>= 1)
				{
					if (From + ToAddition <= Remainder)
					{
						From += ToAddition;
						QuotientDigit += index;
						
						// Check for potential overflow
						if (QuotientDigit >= BASE)
							throw std::runtime_error("Overflow detected in division!");

						// Note: The original code throws an exception if QuotientDigit >= BASE,
						// which is reasonable to catch unexpected behavior during debugging.
					}
					ToAddition >>= 1;
				}

				// Update Remainder and return QuotientDigit
				ToSubtract -= From;
				return QuotientDigit;
			};

			// Main loop for calculating quotient and remainder (Binary search division based on bit chunk)
			for (auto iter = A.values.rbegin(); iter != A.values.rend(); ++iter)
			{
				digit_type RemainderDigit = *iter;

				// Insert RemainderDigit to the front of Remainder's values
				if (RemainderDigit != 0)
				{
					if(Remainder.IsZero())
					{
						Remainder.values[0] = RemainderDigit;
					}
					else
					{
						Remainder.values.insert(Remainder.values.begin(), RemainderDigit);
					}
				}

				// Calculate QuotientDigit using the RunFunction based on binary search
				digit_type QuotientDigit = (B <= Remainder) ? RunFunction(Remainder, B) : 0;

				if(Quotient.IsZero())
				{
					Quotient.values[0] = QuotientDigit;
				}
				else
				{
					// Insert QuotientDigit to the front of Quotient's values
					Quotient.values.insert(Quotient.values.begin(), QuotientDigit);
				}
			}

			// Return the result as a pair of Quotient and Remainder
			return std::pair<BigInteger, BigInteger>(Quotient, Remainder);
		}

		BigInteger& ShortDivision( const BigInteger& other, BigInteger* r = nullptr )
		{
			if ( *this < other )
			{
				if ( r != nullptr )
				{
					r->values = std::move( values );
				}
				else
				{
					*this = 0;
				}

				return *this;
			}
			
			size_t	 n = values.size();
			uint64_t d = other.values[ 0 ];

			std::vector<digit_type> q( n );
			uint64_t			  t, k = 0;

			for ( int64_t i = n - 1; i >= 0; --i )
			{
				if ( values[ i ] + k >= d )
				{
					t = values[ i ] + k;
					q[ i ] = t / d;
					k = uint64_t( t % d ) << EXPONENT;
				}
				else
				{
					k = uint64_t( values[ i ] % d ) << EXPONENT;
				}
			}
			if ( r != nullptr )
			{
				*r = k >> EXPONENT;
			}

			values = std::move( q );

			Clean();
			return *this;
		}

		/**
		 * @brief Performs long division using Knuth's algorithm.
		 *
		 * This function divides the current instance by another BigInteger using
		 * Donald Knuth's long division algorithm. The quotient is stored in the current
		 * instance, and the remainder can be optionally stored in the provided pointer `r`.
		 *
		 * @param other The divisor.
		 * @param r Pointer to store the remainder (can be nullptr if not needed).
		 * @return Reference to the modified current instance after division.
		 */
		BigInteger& DonaldKnuthLongDivision( const BigInteger& other, BigInteger* r = nullptr )
		{
			// Check if the divisor is greater than the dividend
			if ( *this < other )
			{
				if ( r != nullptr )
				{
					r->values = std::move( values );
				}
				else
				{
					*this = 0;
				}
				return *this;
			}

			uint64_t   n = other.values.size();
			uint64_t   m = values.size() - n;
			uint64_t   qhat, rhat;
			BigInteger un, vn;

			// Calculate the shift for normalization
			uint32_t shift = LeadingZeros( other.values[ n - 1 ] );

			vn = other;
			un = *this;
			un.values.push_back( 0 );
			vn <<= shift;
			un <<= shift;

			std::vector<digit_type> q;
			q.resize( m + 1 );

			// Perform Knuth's long division algorithm
			for ( int64_t j = m; j >= 0; --j )
			{
				qhat = ( un[ j + n ] * BASE + un[ j + n - 1 ] ) / vn[ n - 1 ];
				rhat = ( un[ j + n ] * BASE + un[ j + n - 1 ] ) % vn[ n - 1 ];

				do
				{
					if ( qhat >= BASE || qhat * vn[ n - 2 ] > BASE * rhat + un[ j + n - 2 ] )
					{
						qhat--;
						rhat += vn[ n - 1 ];
					}
					else
					{
						break;
					}
				} while ( rhat < BASE );

				int64_t borrow = 0;
				int64_t t = 0;

				// Update the dividend with the quotient
				for ( size_t i = 0; i < n; ++i )
				{
					uint64_t p = qhat * vn[ i ];
					t = un[ i + j ] - borrow - ( p & ( BASE - 1 ) );
					un[ i + j ] = t & ( BASE - 1 );
					borrow = ( p >> EXPONENT ) - ( t >> EXPONENT );
				}

				t = un[ j + n ] - borrow;
				un[ j + n ] = t;

				q[ j ] = qhat;

				if ( t < 0 )
				{
					// Adjust the dividend if the quotient resulted in a negative value
					q[ j ]--;
					borrow = 0;
					for ( size_t i = 0; i < n; ++i )
					{
						t = un[ i + j ] + vn[ i ] + borrow;
						un[ i + j ] = t & ( BASE - 1 );
						borrow = t >> EXPONENT;
					}
					un[ j + n ] += borrow;
				}
			}

			// Update the current instance with the quotient
			values = std::move( q );
			Clean();

			if ( r != nullptr )
			{
				// Normalize and store the remainder
				un >>= shift;
				r->values = std::move( un.values );
				r->Clean();
			}

			return *this;
		}

		bool IsPrime() const
		{
			const BigInteger zero( 0 );
			const BigInteger one( 1 );
			const BigInteger two( 2 );
			const BigInteger three( 3 );

			if ( *this <= one || *this % two == zero || *this % three == zero )
				return false;

			BigInteger i = two;
			BigInteger i2 = zero;
			while ( i * i <= *this )
			{
				i2 = i + two;

				if ( *this % i == zero || *this % i2 == zero )
					return false;

				i += 6;
			}

			return true;
		}

		bool IsEven() const
		{
			return (values[ 0 ] & 1) == 0;
		}

		//bool IsPowerOfTwo() const
		//{
		//	//Is the current sign bit state negative?
		//	if(this->sign == -1)
		//		return false;

		//	//If the last binary bit is 0, it must be an even number, otherwise it is an odd number.
		//	if((values[ 0 ] & 1) == 0)
		//	{
		//		const BigInteger zero( 0 );
		//		const BigInteger one( 1 );

		//		//Performs, on the copied entire array data, bitwise operations on whether the binary is a power of 2 or not.
		//		return (!this->IsZero()) && ((*this) & (*this - one) == zero); 
		//	}

		//	return false;
		//}

		bool IsZero() const
		{
			if(values.size() == 1)
			{
				return values.size() <= 1 && values[ 0 ] == 0;
			}

			return this->values == std::vector<digit_type>(this->values.size(), 0);
		}

		bool IsNegative() const
		{
			return ( sign == 1 ? -1 : 1 );
		}

		size_t Size() const
		{
			return values.size();
		}

		void CountLeadingZeros(size_t& result) const
		{
			if(values.empty())
				return;
			
			digit_type value = 0;
			for(auto iter = values.rbegin(); iter != values.rend(); iter++)
			{
				value = *iter;
				if(value == 0)
				{
					result += EXPONENT;
				}
				else
				{
					result += LeadingZeros(value);
					break;
				}
			}
		}

		static uint16_t LeadingZeros( digit_type x )
		{
			if(x == 0)
			{
				return EXPONENT;
			}

			uint16_t n = 0;
			while ( x <= ( BASE - 1 ) / 2 )
			{
				x <<= 1;
				n++;
			}

			return n;
		}

		size_t BitSize() const
		{
			uint64_t count = 0;
			digit_type high = values[ values.size() - 1 ];
			while ( high != 0 )
			{
				high >>= 1;
				count += 1;
			}

			if(count)
				return ( values.size() - 1 ) * EXPONENT + count;
			return ( values.size() ) * EXPONENT;
		}

		bool GetBit( size_t bit_position ) const
		{
			size_t wrapper_index = bit_position / EXPONENT;
			size_t bit_index = bit_position - ( wrapper_index * EXPONENT );
			return values[ wrapper_index ] & ( digit_type( 1 ) << bit_index );
		}

		void SetBit( size_t bit_position )
		{
			size_t wrapper_index = bit_position / EXPONENT;
			size_t bit_index = bit_position - ( wrapper_index * EXPONENT );
			values[ wrapper_index ] = values[ wrapper_index ] | ( digit_type( 1 ) << bit_index );
		}

		void SetBit( bool value, size_t bit_position )
		{
			if ( bit_position >= this->values.size() * EXPONENT )
			{
				throw std::out_of_range( "Index out of range from set bit" );
			}
			size_t wrapper_index = bit_position / EXPONENT;
			size_t bit_index = bit_position - ( wrapper_index * EXPONENT );

			if ( value )
			{
				this->values[ wrapper_index ] = values[ wrapper_index ] | ( digit_type( 1 ) << bit_index );
			}
			else
			{
				this->values[ wrapper_index ] = values[ wrapper_index ] | ( digit_type( 0 ) << bit_index );
			}
		}

		// Add 'count' number of bits with the specified value to the most significant end of the integer.
		void BigInteger::ExtendedLeadingBits(size_t count, bool value)
		{
			// Calculate the number of digits needed to accommodate the new bits.
			size_t newDigits = (count + EXPONENT - 1) / EXPONENT;

			const digit_type BIT_ZEROS = digit_type(0) & BASE;
			const digit_type BIT_ONES = ~digit_type(0) & BASE;

			// If needed, add new digits to the beginning of the values vector.
			if (newDigits > 0)
			{
				values.insert(values.begin(), newDigits, value ? BIT_ONES : BIT_ZEROS);
			}

			// Adjust the bits within the existing most significant digit.
			size_t remainingBits = count % EXPONENT;
			if (remainingBits > 0)
			{
				digit_type& msb = values.back();
				if (value)
				{
					msb |= (BIT_ONES >> (EXPONENT - remainingBits));
				}
				else
				{
					msb &= (BIT_ZEROS << remainingBits) | (BIT_ONES >> (EXPONENT - remainingBits));
				}
			}
		}

		// Remove 'count' number of bits from the most significant end of the integer.
		void BigInteger::SqueezeLeadingBits(size_t count)
		{
			const digit_type BIT_ZEROS = digit_type(0) & BASE;
			const digit_type BIT_ONES = ~digit_type(0) & BASE;

			// Calculate the number of digits needed to remove the specified bits.
			size_t removedDigits = count / EXPONENT;

			// If needed, remove digits from the beginning of the values vector.
			if (removedDigits > 0)
			{
				if (removedDigits >= values.size())
				{
					// Remove all digits if necessary.
					values.clear();

					// Keep the meaning of this large integer as the number 0
					values.push_back(0);
				}
				else
				{
					values.erase(values.begin(), values.begin() + removedDigits);
				}
			}

			// Adjust the remaining bits within the new most significant digit.
			size_t remainingBits = count % EXPONENT;
			if (remainingBits > 0 && !values.empty())
			{
				digit_type& msb = values.front();
				msb &= (BIT_ONES >> remainingBits);
			}
		}

		uint64_t ToInt() const
		{
			int64_t result = 0;
			int64_t power = 1;
			int64_t base = BASE;
			for ( size_t i = 0; i < values.size(); ++i )
			{
				result += values[ i ] * power;
				power *= base;
			}

			return result * sign;
		}

		BigInteger operator&( const BigInteger& other ) const
		{
			BigInteger result = std::max( this->values.size(), other.values.size() ) ? (*this) : other;
			size_t	digit_size = std::min( this->values.size(), other.values.size() );
			for ( size_t i = 0; i < digit_size; i++ )
			{
				result.values[i] = this->values[i] & other.values[i];
			}

			return result;
		}

		BigInteger& operator&=( const BigInteger& other )
		{
			*this = *this & other;
			return *this;
		}

		BigInteger operator|( const BigInteger& other ) const
		{
			BigInteger result = std::max( this->values.size(), other.values.size() ) ? (*this) : other;
			size_t	digit_size = std::min( this->values.size(), other.values.size() );
			for ( size_t i = 0; i < digit_size; i++ )
			{
				result.values[i] = this->values[i] | other.values[i];
			}

			return result;
		}

		BigInteger& operator|=( const BigInteger& other )
		{
			*this = *this | other;
			return *this;
		}

		BigInteger operator~() const
		{
			BigInteger result = *this;
			for ( size_t i = 0; i < this->values.size(); i++ )
			{
				result.values[i] = ~(this->values[i]);
				result.values[i] <<= EXPONENT;
				result.values[i] >>= EXPONENT;
			}

			return result;
		}

		BigInteger operator^( const BigInteger& other ) const
		{
			BigInteger result = std::max( this->values.size(), other.values.size() ) ? (*this) : other;
			size_t	digit_size = std::min( this->values.size(), other.values.size() );
			for ( size_t i = 0; i < digit_size; i++ )
			{
				result.values[i] = this->values[i] ^ other.values[i];
			}

			return result;
		}

		BigInteger& operator^=( const BigInteger& other )
		{
			*this = *this ^ other;
			return *this;
		}

		BigInteger& operator<<=( const uint32_t shift )
		{
			// If the shift is greater than the maximum exponent (EXPONENT)
			if ( shift > EXPONENT )
			{
				/* Shift the binary representation */
				//// Convert the BigInteger to a binary string
				////binary_string[0] is MSB
				////binary_string[binary_string.size() - 1] is LSB
				//std::string binary_string = ToString( 2 );
				//size_t		binary_string_size = binary_string.size();

				//// If the binary string is shorter than the shift amount, pad it with zeros
				//if ( binary_string_size < shift )
				//{
				//	binary_string.append( shift - binary_string_size, '0' );
				//	FromString( binary_string, 2 );
				//	return *this;
				//}

				//// Calculate the amount to shift within the binary string size
				//uint32_t shift_amount = shift % binary_string_size;
				//// Shift the binary string by removing the leading bits and appending zeros
				//binary_string.erase( binary_string.begin(), binary_string.begin() + shift_amount );
				//binary_string.append( shift_amount, '0' );
				//// Convert the modified binary string back to a BigInteger
				//FromString( binary_string, 2 );
				//return *this;
				
				auto shift_amount = shift;
				while ( shift_amount >= EXPONENT )
				{
					BitLeftShiftCore( (EXPONENT / 2) );

					shift_amount -= (EXPONENT / 2);
				}

				*this >>= shift_amount;

				return *this;
			}

			// If the shift is zero or negative, return the current BigInteger unchanged
			if ( shift <= 0 )
			{
				return *this;
			}

			// If the BigInteger is zero, return it unchanged
			if ( IsZero() )
			{
				return *this;
			}

			// If the shift is less and equal the maximum exponent, perform the shift directly on the BigInteger's internal representation
			BitLeftShiftCore( shift );

			// Return the modified BigInteger
			return *this;
		}

		friend BigInteger operator<<( const BigInteger& lhs, const uint32_t shift )
		{
			BigInteger result( lhs );
			result <<= shift;
			return result;
		}

		BigInteger& operator>>=( const uint32_t shift )
		{
			// If the shift is greater than the maximum exponent (EXPONENT)
			if ( shift > EXPONENT )
			{
				/*Shift the binary representation*/
				//// Convert the BigInteger to a binary string
				////binary_string[0] is MSB
				////binary_string[binary_string.size() - 1] is LSB
				//std::string binary_string = ToString( 2 );
				//size_t		binary_string_size = binary_string.size();
				
				//// If the shift is greater than or equal to the size of the binary string, set the BigInteger to zero
				//if ( shift >= binary_string_size )
				//{
				//	*this = BigInteger(0);
				//	return *this;
				//}
				
				//// Calculate the amount to shift within the binary string size
				//uint32_t	shift_amount = shift % binary_string_size;
				//// Shift the binary string by removing the trailing bits and inserting zeros at the beginning
				//binary_string.erase( binary_string.end() - shift_amount, binary_string.end());
				//binary_string.insert( binary_string.begin(), shift_amount, '0' );
				//// Convert the modified binary string back to a BigInteger
				//FromString( binary_string, 2 );
				//return *this;

				auto shift_amount = shift;

				while ( shift_amount >= EXPONENT )
				{
					BitRightShiftCore( (EXPONENT / 2) );

					shift_amount -= (EXPONENT / 2);
				}

				*this >>= shift_amount;

				return *this;
			}

			// If the shift is zero or negative, return the current BigInteger unchanged
			if ( shift <= 0 )
			{
				return *this;
			}

			// If the BigInteger is zero, return it unchanged
			if ( IsZero() )
			{
				return *this;
			}

			// If the shift is less and equal the maximum exponent, perform the shift directly on the BigInteger's internal representation
			BitRightShiftCore( shift );

			// Clean up any leading zeros that may have been introduced by the shift operation
			Clean();

			// Return the modified BigInteger
			return *this;
		}

		friend BigInteger operator>>( const BigInteger& lhs, const uint32_t shift )
		{
			BigInteger result( lhs );
			result >>= shift;
			return result;
		}

		/**
		 * @brief Performs a left rotation on the binary representation of the BigInteger.
		 *
		 * @param shift The number of positions to rotate the bits to the left.
		 * @param reference_bit_capacity The desired bit capacity for the resulting binary string.
		 * @return copy to the modified BigInteger.
		 */
		static BigInteger BitRotateLeft(const BigInteger& bits, uint32_t shift, const uint32_t reference_bit_capacity)
		{
			// Check if the BigInteger is zero, if so, no rotation is needed
			if (bits.IsZero() || reference_bit_capacity == 0)
				return bits;

			shift %= reference_bit_capacity;

			#if 1

			// Convert the BigInteger to a binary string representation
			// The binary string length based on the reference_bit_capacity
			std::string binary_string = bits.ToBinaryString(reference_bit_capacity);

			// Left rotate the binary string
			std::rotate(binary_string.begin(), binary_string.begin() + shift, binary_string.end());

			BigInteger result;

			// Update the BigInteger with the rotated binary string
			result.FromString(binary_string, 2);

			// Ensure the BigInteger has the desired bit capacity
			if (reference_bit_capacity / EXPONENT != result.values.size())
			{
				size_t bit_size = result.BitSize();
				if(bit_size < reference_bit_capacity)
					result.ExtendedLeadingBits(reference_bit_capacity - bit_size, false);
			}
			
			#else

			BigInteger left = *this << shift;
			BigInteger right = *this >> shift;
			
			size_t bit_size = left.BitSize();
			if (bit_size < reference_bit_capacity)
			{
				left.ExtendedLeadingBits(reference_bit_capacity - bit_size, false);
			}
			bit_size = right.BitSize();
			if (bit_size < reference_bit_capacity)
			{
				right.ExtendedLeadingBits(reference_bit_capacity - bit_size, false);
			}

			BigInteger result = left | right;
			bit_size = result.BitSize();
			if (bit_size < reference_bit_capacity)
			{
				this->ExtendedLeadingBits(reference_bit_capacity - bit_size, false);
			}

			#endif

			return result;
		}

		/**
		 * @brief Performs a right rotation on the binary representation of the BigInteger.
		 *
		 * @param shift The number of positions to rotate the bits to the right.
		 * @param reference_bit_capacity The desired bit capacity for the resulting binary string.
		 * @return copy to the modified BigInteger.
		 */
		static BigInteger BitRotateRight(const BigInteger& bits, uint32_t shift, const uint32_t reference_bit_capacity)
		{
			// Check if the BigInteger is zero, if so, no rotation is needed
			if (bits.IsZero() || reference_bit_capacity == 0)
				return bits;

			shift %= reference_bit_capacity;

			#if 1
			
			// Convert the BigInteger to a binary string representation
			// The binary string length based on the reference_bit_capacity
			std::string binary_string = bits.ToBinaryString(reference_bit_capacity);

			// Right rotate the binary string
			std::rotate(binary_string.rbegin(), binary_string.rbegin() + shift, binary_string.rend());

			BigInteger result;
			
			// Update the BigInteger with the rotated binary string
			result.FromString(binary_string, 2);

			// Ensure the BigInteger has the desired bit capacity
			if (reference_bit_capacity / EXPONENT != result.values.size())
			{
				size_t bit_size = result.BitSize();
				if(bit_size < reference_bit_capacity)
					result.ExtendedLeadingBits(reference_bit_capacity - bit_size, false);
			}
			
			#else

			BigInteger right = *this >> shift;
			BigInteger left = *this << shift;
			
			size_t bit_size = right.BitSize();
			if (bit_size < reference_bit_capacity)
			{
				right.ExtendedLeadingBits(reference_bit_capacity - bit_size, false);
			}
			bit_size = left.BitSize();
			if (bit_size < reference_bit_capacity)
			{
				left.ExtendedLeadingBits(reference_bit_capacity - bit_size, false);
			}
			
			BigInteger result = left | right;
			bit_size = result.BitSize();
			if (bit_size < reference_bit_capacity)
			{
				this->ExtendedLeadingBits(reference_bit_capacity - bit_size, false);
			}

			#endif

			return result;
		}

		digit_type& operator[]( const size_t index )
		{
			return values[ index ];
		}

		friend BigInteger operator/( const BigInteger& lhs, const BigInteger& rhs )
		{
			BigInteger result( lhs );
			result /= rhs;
			return result;
		}

		BigInteger& operator%=( const BigInteger& other )
		{
			BigInteger r;
			Divide( other, &r );
			values = std::move( r.values );
			return *this;
		}

		friend BigInteger operator%( const BigInteger& lhs, const BigInteger& rhs )
		{
			BigInteger result( lhs );
			result %= rhs;
			return result;
		}

		BigInteger& operator=( const BigInteger& other )
		{
			values = other.values;
			sign = other.sign;
			return *this;
		}

		explicit operator bool()
		{
			return this->IsZero() ? false : true;
		}

		bool operator&&( const BigInteger& other ) const
		{
			bool left = this->IsZero();
			bool right = other.IsZero();

			if ( ( left == true ) && ( right == true ) )
				return false;
			else if ( ( left == false ) && ( right == false ) )
				return true;
			else
				return false;
		}

		bool operator||( const BigInteger& other ) const
		{
			bool left = this->IsZero();
			bool right = other.IsZero();

			if ( ( left == true ) && ( right == true ) )
				return false;
			else if ( ( left == false ) || ( right == false ) )
				return true;
		}

		bool operator!() const
		{
			return !this->IsZero();
		}

		friend bool operator==( const BigInteger& lhs, const BigInteger& rhs )
		{
			const size_t a = lhs.values.size();
			const size_t b = rhs.values.size();

			if ( a != b )
				return false;
			for ( size_t i = 0; i < a; ++i )
			{
				if ( lhs.values[ i ] != rhs.values[ i ] )
					return false;
			}

			return true;
		}

		friend bool operator!=( const BigInteger& lhs, const BigInteger& rhs )
		{
			return !( lhs == rhs );
		}

		friend bool operator>( const BigInteger& lhs, const BigInteger& rhs )
		{
			size_t a = lhs.values.size();
			size_t b = rhs.values.size();

			if ( a > b )
				return true;
			if ( b > a )
				return false;
			for ( size_t i = 0; i < a; ++i )
			{
				digit_type v1 = lhs.values[ a - ( i + 1 ) ];
				digit_type v2 = rhs.values[ a - ( i + 1 ) ];
				if ( v1 > v2 )
					return true;
				else if ( v1 < v2 )
					return false;
			}

			return false;
		}

		friend bool operator>=( const BigInteger& lhs, const BigInteger& rhs )
		{
			if ( lhs > rhs || lhs == rhs )
				return true;

			return false;
		}

		friend bool operator<( const BigInteger& lhs, const BigInteger& rhs )
		{
			size_t a = lhs.values.size();
			size_t b = rhs.values.size();

			if ( a < b )
				return true;
			if ( b < a )
				return false;
			for ( size_t i = 0; i < a; ++i )
			{
				digit_type v1 = lhs.values[ a - ( i + 1 ) ];
				digit_type v2 = rhs.values[ a - ( i + 1 ) ];
				if ( v1 < v2 )
					return true;
				else if ( v1 > v2 )
					return false;
			}

			return false;
		}

		friend bool operator<=( const BigInteger& lhs, const BigInteger& rhs )
		{
			if ( lhs < rhs || lhs == rhs )
				return true;

			return false;
		}

		static BigInteger GCD( BigInteger a, BigInteger b )
		{
			BigInteger temp, t = 1;

			while ( a.IsEven() && b.IsEven() )
			{
				t <<= 1;
				a >>= 1;
				b >>= 1;
			}
			while ( a.IsEven() )
			{
				a >>= 1;
			}
			while ( b.IsEven() )
			{
				b >>= 1;
			}
			while ( a != b )
			{
				temp = ( a < b ? a : b );
				a.Difference( b );
				b = temp;
				while ( a.IsEven() )
				{
					a >>= 1;
				}
			}

			return t * a;
		}

		static void EGCD( const BigInteger& a, const BigInteger& b, BigInteger* gcd, BigInteger* co1, BigInteger* co2 )
		{
			BigInteger old_r, r;
			BigInteger old_s, s;
			BigInteger old_t, t;
			BigInteger q, temp;

			old_r = a;
			r = b;
			old_s = 1;
			s = 0;
			old_t = 0;
			t = 1;

			while ( !r.IsZero() )
			{
				q = old_r / r;

				temp = r;
				BigInteger testk = q * temp;
				r = old_r - testk;
				old_r = temp;

				temp = s;
				s = old_s - q * temp;
				old_s = temp;

				temp = t;
				t = old_t - q * temp;
				old_t = temp;
			}

			if ( co1 != nullptr )
				( *co1 ) = old_s;
			if ( co2 != nullptr )
				( *co2 ) = old_t;
			if ( gcd != nullptr )
				( *gcd ) = old_r;
		}

		static BigInteger LCM( const BigInteger& a, const BigInteger& b )
		{
			BigInteger result = a * b;
			result.sign = 1;
			result /= GCD( a, b );
			return result;
		}

		static BigInteger ModuloInverse( const BigInteger& a, const BigInteger& b )
		{
			BigInteger gcd, x;

			EGCD( a, b, &gcd, &x, nullptr );
			return ( x % b + b ) % b;
		}

		/**
		 * @brief Pollard's Rho algorithm for integer factorization.
		 *
		 * This function applies Pollard's Rho algorithm to find a non-trivial factor
		 * of the given integer `n`. It returns a factor if found, or -1 if no non-trivial
		 * factor is discovered.
		 *
		 * @param n The integer to factorize.
		 * @return A non-trivial factor of `n` or -1 if none is found.
		 */
		static BigInteger PollardRho( const BigInteger& n )
		{
			const BigInteger zero( 0 );
			const BigInteger one( 1 );
			const BigInteger two( 2 );

			// Handle special cases
			if ( n == 1 )
				return n;
			if ( n % two == zero )
				return 2;

			// Initialize variables for the algorithm
			BigInteger x = 2, y = 2, d = 1, diff;
			uint64_t   i = 0;

			// Perform Pollard's Rho algorithm
			while ( d == one )
			{
				x = x * x;
				x.Add( one );
				x %= n;

				y = y * y;
				y.Add( one );
				y %= n;
				y = y * y;
				y.Add( one );
				y %= n;

				diff = ( x > y ? x : y );
				diff = diff.Subtract( ( x > y ? y : x ) );
				d = GCD( diff, n );
			}

			// Check if a non-trivial factor is found
			if ( d == n )
				return -1;
			else
				return d;
		}

		static BigInteger RandomGenerateNBit( size_t n )
		{
			std::uniform_int_distribution<uint16_t> dist01( 0, 1 );
			BigInteger								ud_prng;
			ud_prng.values.resize( n / EXPONENT + 1 );
			ud_prng.SetBit( n - 1 );
			for ( size_t i = 0; i < n; ++i )
			{
				if ( dist01( RandomEngine ) )
				{
					ud_prng.SetBit( i );
				}
			}

			return ud_prng;
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
		static bool MillerRabin( const BigInteger& n, int k )
		{
			const BigInteger zero( 0 );
			const BigInteger one( 1 );
			const BigInteger two( 2 );

			// Initialize random number generator
			std::uniform_int_distribution<uint64_t> ud_prng( 2, ( ( n - 2 ) >= BigInteger( INT64_MAX ) ? INT64_MAX : ( n.ToInt() - 2 ) ) );

			uint32_t   s = 0;
			BigInteger d = n - 1;
			BigInteger a, x;

			bool tryAgain = false;

			// Factorize n-1 to d * 2^s
			while ( d.IsEven() )
			{
				s++;
				d >>= 1;
			}

			// Perform Miller-Rabin test
			for ( size_t i = 0; i < k; ++i )
			{
				a = ud_prng( RandomEngine );
				x = a.PowerWithModulo( d, n );

				if ( x == one || x == n - 1 )
				{
					continue;
				}

				for ( size_t r = 0; r < s; ++r )
				{
					x.PowerWithModulo( 2, n );

					if ( x == one )
					{
						return false;
					}
					else if ( x == n - 1 )
					{
						tryAgain = true;
						break;
					}
				}

				if ( !tryAgain )
					return false;
				tryAgain = false;
			}

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
		static bool MillerRabinWithMontgomery( const BigInteger& n, int k )
		{
			const BigInteger one( 1 );

			// Initialize random number generator
			std::uniform_int_distribution<uint64_t> ud_prng( 2, ( ( n - 2 ) >= BigInteger( INT64_MAX ) ? INT64_MAX : ( n.ToInt() - 2 ) ) );

			uint32_t   s = 0;
			BigInteger d = n - 1;
			BigInteger a, x;

			bool tryAgain = false;

			// Factorize n-1 to d * 2^s
			while ( d.IsEven() )
			{
				s++;
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
				r <<= EXPONENT;
			}

			for ( size_t i = 0; i < rsize * EXPONENT; ++i )
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
					mprime.SetBit( rsize * EXPONENT - 1 );
				}
			}

			// Perform Miller-Rabin test with Montgomery multiplication
			for ( size_t i = 0; i < k; ++i )
			{
				a = ud_prng( RandomEngine );
				x = a.MontgomeryPower( d, n, mprime, r, rsize );

				if ( x == one || x == n - 1 )
				{
					continue;
				}

				for ( size_t r = 0; r < s; ++r )
				{
					x.PowerWithModulo( 2, n );

					if ( x == one )
					{
						return false;
					}
					else if ( x == n - 1 )
					{
						tryAgain = true;
						break;
					}
				}

				if ( !tryAgain )
					return false;
				tryAgain = false;
			}

			return true;
		}

		/**
		 * @brief Performs modular arithmetic operations on BigIntegers.
		 *
		 * This class provides static methods for performing modular arithmetic operations
		 * on BigIntegers. The supported operations include addition, subtraction, multiplication,
		 * and division. The modular operations are performed with respect to a given modulus.
		 * The modulus can be any non-zero BigInteger, and specific optimizations are applied
		 * when the modulus is a power of two.
		 *
		 * @note The division operation requires the modulus to be non-zero and, in the case of
		 *       division, either the divisor must not be zero or the modulus must be a prime number.
		 *
		 * @param number The BigInteger on which the modulo operation is to be performed.
		 * @param modulus The modulus with respect to which the modulo operation is performed.
		 * @return The result of the modulo operation (number % modulus).
		 */
		static BigInteger Modulo(const BigInteger& number, const BigInteger& modulus)
		{
			if(number.IsZero())
			{
				return BigInteger(0);
			}

			BigInteger result;

			/*if ( modulus.IsPowerOfTwo() )
			{
				result = number & ( modulus - 1 );
				result.Clean();
				return result;
			}*/

			result = number % modulus;
			result.Clean();
			return result;
		}

		/**
		 * @brief Performs modular arithmetic operations on two BigIntegers.
		 *
		 * This method supports modular addition, subtraction, multiplication, and division of
		 * two BigIntegers. The specific operation is determined by the `mode` parameter.
		 *
		 * @param mode The arithmetic mode specifying the operation to be performed
		 *             (Addition, Subtraction, Multiplication, or Division).
		 * @param a The first BigInteger operand.
		 * @param b The second BigInteger operand.
		 * @param modulus The modulus with respect to which the modular operation is performed.
		 * @return The result of the specified modular arithmetic operation (a op b) % modulus.
		 *
		 * @throws std::invalid_argument If the modulus is zero or, in the case of division,
		 *         the divisor is zero or the modulus is not prime.
		 * @throws std::invalid_argument If the modular inverse for division cannot be found,
		 *         ensuring (a * b) mod modulus != (a * inverse(b)) mod modulus.
		 */
		static BigInteger ModuloArithmetic(ArithmeticMode mode, const BigInteger& a , const BigInteger& b, const BigInteger& modulus)
		{
			if (modulus.IsZero())
			{
				throw std::invalid_argument("The modulus cannot be zero!");
			}

			BigInteger temporary;

			switch (mode)
			{
				case ArithmeticMode::Addition:
					temporary = a + b;
					break;
				case ArithmeticMode::Subtraction:
					temporary = a - b;
					break;
				case ArithmeticMode::Multiplication:
					temporary = a * b;
					break;
				case ArithmeticMode::Division:
					if (!b.IsZero())
					{
						BigInteger inverse = BigInteger::ModuloInverse(b, modulus);
						if (!inverse.IsZero())
						{
							temporary = a * inverse;
						}
						else
						{
							throw std::invalid_argument("The b is not have multiplication inverse with modulo the modulus number, This mean (a * b) mod modulus_number != (a * inverse(b)) mod modulus_number");
						}
					}
					else
					{
						throw std::invalid_argument("Division by zero!");
					}
					break;
				default:
					throw std::invalid_argument("Invalid modulo arithmetic mode!");
					break;
			}

			return Modulo(temporary, modulus);
		}

		std::string ToString() const
		{
			if ( *this == 0 )
			{
				return "0";
			}

			uint64_t   m = 0;
			BigInteger a = *this;

			const BigInteger TEN( 10 );

			std::string number_string;
			while ( !a.IsZero() )
			{
				number_string += ( a % TEN ).values[ 0 ] + '0';

				a /= TEN;
				m++;
			}
			if ( sign == -1 )
				number_string += '-';

			std::reverse( number_string.begin(), number_string.end() );
			return number_string;
		}

		/**
		* Convert the BigInteger to a binary string representation.
		*
		* @param reference_bit_capacity The desired bit capacity of the resulting binary string.
		* @param have_leading_zeros     Flag indicating whether to include leading zeros in the binary string.
		* @return                       Binary string representation of the BigInteger.
		*/
		std::string ToBinaryString( const uint32_t reference_bit_capacity, bool have_leading_zeros = true ) const
		{
			if ( reference_bit_capacity == 0 )
				return "";

			if ( have_leading_zeros )
			{
				// Initialize the binary string with leading zeros
				std::string number_string( reference_bit_capacity, '0' );

				// Copy the current BigInteger for modification
				BigInteger result = *this;
				size_t	   bit_size = result.BitSize();

				// Set binary digits based on the BigInteger bits
				for ( size_t i = 0; i < bit_size; i++ )
				{
					if ( result.GetBit( i ) )
						number_string[ i ] = '1';
				}

				// Reverse the string to get correct binary representation
				std::reverse( number_string.begin(), number_string.end() );
				return number_string;
			}
			else
			{
				// Initialize the binary string without leading zeros
				std::string number_string;

				// Copy the current BigInteger for processing
				BigInteger result = *this;
				size_t	   bit_size = result.BitSize();

				// Set binary digits based on the BigInteger bits
				for ( size_t i = 0; i < bit_size; i++ )
				{
					number_string.push_back( result.GetBit( i ) ? '1' : '0' );
				}

				// Reverse the string to get correct binary representation
				std::reverse( number_string.begin(), number_string.end() );
				return number_string;
			}
		}

		/**
		* Convert the BigInteger to a string representation with the specified base.
		*
		* @param base_value The desired base for the string representation.
		* @warning Binary base is unsupported in ToString function.
		* @return String representation of the BigInteger in the specified base.
		*/
		std::string ToString( uint32_t base_value ) const
		{
			if ( *this == 0 )
			{
				return "0";
			}

			uint64_t   m = 0;
			BigInteger a = *this;

			std::string number_string;
			
			// Calculate the number of the current bit
			if ( base_value == 2 )
			{
				std::cout << "Warning: Empty string has been returned! \nThe \"ToString\" function does not support binary generation, please use the \"ToBinaryString\" function." << std::endl;
				return number_string;
			}

			BigInteger STRING_BASE = BigInteger( base_value );
			while ( !a.IsZero() )
			{
				// Calculate the current digit
				uint64_t digit = ( a % STRING_BASE ).values[ 0 ];
				
				// Convert digits(number) to characters
				number_string += ( digit < 10 ? '0' + digit : 'A' + digit - 10 );

				// Divide the number in the current bit by 2 from a.
				a /= STRING_BASE;

				// Update m for possible subsequent negative symbols
				m++;
			}
			if ( sign == -1 )
				number_string += '-';

			std::reverse( number_string.begin(), number_string.end() );
			return number_string;
		}

		void FromString( const std::string& number_string )
		{
			size_t	   i = 0, count = number_string.length();
			BigInteger n, result = 0;

			if ( number_string[ 0 ] == '-' )
			{
				sign = -1;
				i++;
			}

			const BigInteger TEN( 10 );

			for ( i; i < count; ++i )
			{
				result *= TEN;
				n = number_string[ i ] - '0';
				result += n;
			}
			values = std::move( result.values );
		}

		void FromString( const std::string& number_string, uint32_t base_value )
		{
			if ( base_value == 0 )
			{
				throw std::invalid_argument( "Base must be greater than 0." );
			}

			if ( base_value > 36 )
			{
				throw std::invalid_argument( "Base must be less than or equal to 36." );
			}

			size_t	   i = 0, count = number_string.length();
			BigInteger result = 0;
			sign = 1;

			if ( number_string[ 0 ] == '-' )
			{
				sign = -1;
				i++;
			}

			BigInteger STRING_BASE = BigInteger( base_value );
			// If the format conversion is binary-based
			if ( base_value == 2 )
			{
				result.ExtendedLeadingBits(count - result.BitSize(), false);

				for ( i; i < count; ++i )
				{
					result.SetBit(number_string[i] == '1' ? true : false, (count - 1) - i);
				}

				values = std::move( result.values );

				return;
			}

			for ( i; i < count; ++i )
			{
				// Get the current character in the string
				char c = number_string[ i ];
				if ( c >= '0' && c <= '9' )
				{
					// The position of the progression is moved and the current number needs to be multiplied by the current progression.
					result *= STRING_BASE;
					
					// Convert characters to digits(numbers)
					int digit = c - '0';
					if ( digit >= base_value )
					{
						throw std::invalid_argument( "Invalid digit for specified base." );
					}
					
					// Calculate the current digit to this result
					result += BigInteger( digit );
				}
				else if ( c >= 'A' && c <= 'Z' )
				{
					// The position of the progression is moved and the current number needs to be multiplied by the current progression.
					result *= STRING_BASE;
					
					// Convert characters to digits(numbers)
					int digit = 10 + c - 'A';
					if ( digit >= base_value )
					{
						throw std::invalid_argument( "Invalid digit for specified base." );
					}
					
					// Calculate the current digit to this result
					result += BigInteger( digit );
				}
				else
				{
					throw std::invalid_argument( "Invalid character in input string." );
				}
			}

			values = std::move( result.values );
		}

		template <size_t N>
		void BitsetData(std::bitset<N>& bits) const
		{
			bool bit = false;
			for ( size_t i = 0; i < N; i++ )
			{
				if(i >= EXPONENT * values.size())
				{
					break;
				}

				bool bit = GetBit(i);

				if(bit)
				{
					bits[i] = true;
				}
			}
		}

		void Print( bool base10 ) const
		{
			if ( base10 )
			{
				std::cout << ToString( 10 ) << "\n\n";
			}
			else
			{
				if ( sign == -1 )
					std::cout << "-";

				for ( const auto& i : values )
				{
					std::cout << i;
				}
				std::cout << "\n\n" << values.size() << " digits" << std::endl;
			}
		}
	};
}  // namespace TwilightDream::BigInteger