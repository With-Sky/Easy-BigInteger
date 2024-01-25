#include "BigInteger.hpp"

namespace TwilightDream::BigInteger
{
	void BigInteger::Clean()
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

	void BigInteger::BitLeftShiftCore( const uint32_t& shift )
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

	void BigInteger::BitRightShiftCore( const uint32_t& shift )
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

	BigInteger BigInteger::RightShiftBlock( size_t d ) const
	{
		if ( d >= Size() )
			return 0;
		BigInteger tmp;
		tmp.values.assign( values.begin() + d, values.end() );	//shift d element, not bit
		return tmp;
	}

	BigInteger BigInteger::LeftShiftBlock( size_t d ) const
	{
		if ( *this == 0 )
			return 0;
		BigInteger tmp;
		tmp.values.resize( d );
		tmp.values.insert( tmp.values.end(), values.begin(), values.end() );
		return tmp;
	}

	BigInteger BigInteger::DivisionInvertReciprocal( size_t n ) const
	{
		// If the exponent is small enough, use the long division algorithm directly.
		if ( std::min( Size(), n - Size() ) <= 64 )
		{
			BigInteger a;
			a.values.resize( n + 1 );
			a.values[ n ] = 1;
			return a.DonaldKnuthLongDivision( *this );	// a / this
		}

		// Calculate the size of the dividend and the number of bits to shift.
		size_t dividend_size = Size();
		size_t k = ( n - dividend_size + 5 ) >> 1;				// k is the number of bits to shift for the divisor.
		size_t k2 = k > dividend_size ? 0 : dividend_size - k;	// k2 is the number of bits to shift in the dividend.

		// Shift the dividend to the right by k2 bits to prepare for division.
		BigInteger shifted_dividend = RightShiftBlock( k2 );

		// Calculate the size of the intermediate result.
		size_t intermediate_size = k + shifted_dividend.Size();

		// Recursively compute the reciprocal of the shifted dividend.
		BigInteger reciprocal = shifted_dividend.DivisionInvertReciprocal( intermediate_size );

		// Calculate the intermediate results for the division.
		BigInteger doubled_reciprocal = reciprocal + reciprocal;
		BigInteger numerator = ( *this ) * reciprocal * reciprocal;

		// Compute the final reciprocal by adjusting for the shifts and intermediate calculations.
		return doubled_reciprocal.LeftShiftBlock( n - intermediate_size - k2 ) - numerator.RightShiftBlock( 2 * ( intermediate_size + k2 ) - n ) - 1;
	}

	BigInteger::BigInteger( int64_t value )
	{
		FromInt( value );
	}

	BigInteger::BigInteger( const std::string& number_string )
	{
		FromString( number_string );
	}

	BigInteger::BigInteger( const std::string& number_string, uint32_t base )
	{
		FromString( number_string, base );
	}

	BigInteger::BigInteger( const BigInteger& num )
	{
		values = num.values;
		sign = num.sign;
	}

	BigInteger::BigInteger( BigInteger&& num ) noexcept
	{
		values = std::move( num.values );
		sign = num.sign;
	}

	BigInteger& BigInteger::operator+=( const BigInteger& other )
	{
		if ( this == &other )
		{
			BigInteger self( other );
			*this += self;
			return *this;
		}

		if ( this->IsZero() && other.sign == 1 )
		{
			*this = other;
			return *this;
		}
		else if ( other.IsZero() )
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

	BigInteger& BigInteger::Add( const BigInteger& other )
	{
		// Determine the maximum size for the loop and result storage
		const size_t this_size = values.size();
		const size_t other_size = other.values.size();
		const size_t max = ( this_size > other_size ? this_size : other_size ) + 1;

		// Resize the storage for the result
		values.resize( max );

		// Initialize carry and perform addition
		uint8_t carry = 0;
		for ( size_t i = 0; i < max; ++i )
		{
			uint64_t sum = values[ i ] + ( ( i < other_size ) ? other.values[ i ] : 0 ) + carry;
			values[ i ] = sum & ( BASE - 1 );  // Update current digit with the sum
			carry = sum >> EXPONENT;		   // Update carry for the next iteration
		}

		// Remove leading zeros and resize the storage
		Clean();

		// Return a reference to the modified current instance
		return *this;
	}

	BigInteger& BigInteger::operator-=( const BigInteger& other )
	{
		if ( this == &other )
		{
			BigInteger self( other );
			*this -= self;
			return *this;
		}

		if ( this->IsZero() && other.sign == -1 )
		{
			*this = other;
			this->sign = 1;
			return *this;
		}
		else if ( other.IsZero() )
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

	BigInteger& BigInteger::Difference( const BigInteger& other )
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

	BigInteger& BigInteger::Subtract( const BigInteger& other )
	{
		// Determine the size of current and other BigIntegers
		const size_t this_size = values.size();
		const size_t other_size = other.values.size();

		// Initialize variables for borrow and temporary subtraction
		int64_t borrow = 0;
		int64_t temp;

		// Perform subtraction for each digit
		for ( size_t i = 0; i < this_size; ++i )
		{
			temp = values[ i ] - ( ( i >= other_size ) ? 0 : other.values[ i ] ) + borrow;
			values[ i ] = temp & ( BASE - 1 );	// Update current digit with the result
			borrow = temp >> EXPONENT;			// Update borrow for the next iteration
		}

		// Remove leading zeros and resize the storage
		Clean();

		// Return a reference to the modified current instance
		return *this;
	}

	BigInteger BigInteger::operator++()
	{
		*this += BigInteger( 1 );
		return *this;
	}

	BigInteger BigInteger::operator--()
	{
		*this -= BigInteger( 1 );
		return *this;
	}

	BigInteger BigInteger::operator++( int )
	{
		BigInteger copy = *this;
		copy += BigInteger( 1 );
		return copy;
	}

	BigInteger BigInteger::operator--( int )
	{
		BigInteger copy = *this;
		copy -= BigInteger( 1 );
		return copy;
	}

	BigInteger BigInteger::operator-() const
	{
		BigInteger copy = *this;
		copy.sign = ( copy.sign == 1 ) ? -1 : 1;
		return copy;
	}

	BigInteger BigInteger::operator+() const
	{
		BigInteger copy = *this;
		copy.sign = ( copy.sign == 1 ) ? 1 : -1;
		return copy;
	}

	BigInteger& BigInteger::BaseMultiplication( const BigInteger& other )
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

	BigInteger& BigInteger::FHTMultiplication( const BigInteger& other )
	{
		const size_t this_size = values.size();
		const size_t other_size = other.values.size();

		// Make a copy of the current instance's values
		std::vector<FHT_BITS_TYPE> copy1( values.begin(), values.end() );
		std::vector<FHT_BITS_TYPE> copy2( other.values.begin(), other.values.end() );
		std::vector<FHT_BITS_TYPE> result( this_size + other_size );

		// Use FHT to accelerate multiplication
		HyperIntegerFunctions::FHTMul( result.data(), copy1.data(), copy1.size(), copy2.data(), copy2.size() );

		values = std::vector<digit_type>( result.begin(), result.end() );
		// Remove leading zeros and resize the storage
		Clean();

		// Return a reference to the modified current instance
		return *this;
	}

	BigInteger& BigInteger::FHTSquare()
	{
		size_t					   t = values.size();
		std::vector<FHT_BITS_TYPE> copy( values.begin(), values.end() );
		std::vector<FHT_BITS_TYPE> result( t * 2 );

		// Perform FHT algorithm for squaring
		HyperIntegerFunctions::FHTSquare( result.data(), copy.data(), t );

		values = std::vector<digit_type>( result.begin(), result.end() );
		Clean();
		return *this;
	}

	BigInteger& TwilightDream::BigInteger::BigInteger::Square()
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

	BigInteger& BigInteger::MultiplyBase( size_t times )
	{
		values.insert( values.begin(), times, 0 );
		return *this;
	}

	void BigInteger::SplitAt( const BigInteger& num, const size_t n, BigInteger& high, BigInteger& low )
	{
		std::vector<digit_type> lowValue( num.values.begin(), num.values.begin() + n );
		std::vector<digit_type> highValue( num.values.end() - ( num.values.size() - n ), num.values.end() );

		low.values = std::move( lowValue );
		high.values = std::move( highValue );
	}

	BigInteger& TwilightDream::BigInteger::BigInteger::KaratsubaMultiplication( const BigInteger& other )
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

	BigInteger& BigInteger::Power( const size_t exponent )
	{
		if ( exponent == 0 )
		{
			*this = 1;
			return *this;
		}

		BigInteger result( *this );
		for ( size_t i = 0; i < exponent - 1; ++i )
		{
			result *= *this;
		}
		*this = result;
		return *this;
	}

	BigInteger& BigInteger::BigPower( const BigInteger& exponent )
	{
		const BigInteger ZERO = 0;
		const BigInteger ONE = 1;
		BigInteger		 result = ONE;
		BigInteger		 base = *this;

		if ( exponent == ZERO )
		{
			*this = 1;
			return *this;
		}

		size_t index_bits = exponent.BitSize();
		for ( size_t i = 0; i < index_bits; i++ )
		{
			if ( exponent.GetBit( i ) )
			{
				result *= base;
			}

			base.Power( 2 );
		}

		values = std::move( result.values );
		return *this;
	}

	BigInteger& BigInteger::MontDivR( size_t rsize )
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

	BigInteger& BigInteger::MontModR( size_t rsize )
	{
		if ( values.size() > rsize )
		{
			values.erase( values.end() - ( values.size() - rsize ), values.end() );
		}
		return *this;
	}

	BigInteger& BigInteger::MontgomeryReduce( size_t rsize, const BigInteger& m, const BigInteger& mprime )
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

	BigInteger& BigInteger::MontgomeryPower( const BigInteger& e, const BigInteger& m, const BigInteger& mprime, const BigInteger r, const size_t rsize )
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

	BigInteger& BigInteger::MontgomeryPower( const BigInteger& e, const BigInteger& m )
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

	BigInteger& BigInteger::Divide( const BigInteger& other, BigInteger* r )
	{
		// Handle special case when the numerator is zero
		if ( values.size() == 1 && values[ 0 ] == 0 )
		{
			if ( r != nullptr )
			{
				// Return 0 with positive sign if the numerator is zero
				r->sign = 1;
				r->values.resize( 1 );
				r->values[ 0 ] = 0;
			}
			return *this;
		}

		/*
				BinarySearch <= Short <= DonaldKnuthLong <= NewtonIteration
			*/

		// Select division algorithm based on the size of the divisor
		if ( other.values.size() < BINARY_SEARCH_DIVISION_LIMIT )
		{
			// Use BinarySearchDivision for large divisors
			auto pair = BinarySearchDivision( other );
			if ( r != nullptr )
			{
				// Store the remainder in the provided pointer
				*r = pair.second;
			}
			*this = pair.first;	 // Update the current instance with the quotient
			return *this;
		}
		else if ( other.values.size() < SHORT_DIVISION_LIMIT )
		{
			//There are special cases where the computation fails using short division.
			if ( other.values[ 0 ] == 0 )
			{
				DonaldKnuthLongDivision( other, r );
				return *this;
			}
			// Use ShortDivision for small divisors
			return ShortDivision( other, r );
		}
		else if ( other.values.size() < DONALD_KNUTH_LONG_DIVISION_LIMIT )
		{
			// Use DonaldKnuthLongDivision for intermediate-sized divisors
			DonaldKnuthLongDivision( other, r );
			return *this;
		}
		else
		{
			// Use Newton iteration for long-sized divisors
			BigInteger copy( *this );
			auto	   pair = copy.NewtonIterationDivision( other );
			if ( r != nullptr )
			{
				*r = pair.second;
			}
			return *this = pair.first;
		}
	}

	std::pair<BigInteger, BigInteger> BigInteger::BinarySearchDivision( const BigInteger& other )
	{
		// Make copies of the numerator and denominator
		BigInteger		 A = *this;
		BigInteger		 B = other;
		BigInteger		 Quotient = 0;
		BigInteger		 Remainder = 0;
		const BigInteger Zero = 0;

		// Lambda function for the main division logic
		auto RunFunction = [ &Zero ]( BigInteger& Remainder, const BigInteger& Divisor ) {
			digit_type QuotientDigit = 0;

			//Multiplier for Reference Remainder
			BigInteger From = 0;
			//Reference Remainder
			BigInteger&& ToSubtract = std::move( Remainder );

			//Copy Divisor
			BigInteger ToAddition = Divisor;
			// Left shift ToAddition by (EXPONENT - 1) to prepare for binary search
			ToAddition <<= EXPONENT - 1;

			// Binary search loop
			for ( digit_type index = digit_type( 1 ) << ( EXPONENT - 1 ); index > 0; index >>= 1 )
			{
				if ( From + ToAddition <= Remainder )
				{
					From += ToAddition;
					QuotientDigit += index;

					// Check for potential overflow
					if ( QuotientDigit >= BASE )
						throw std::runtime_error( "Overflow detected in division!" );

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
		for ( auto iter = A.values.rbegin(); iter != A.values.rend(); ++iter )
		{
			digit_type RemainderDigit = *iter;

			// Insert RemainderDigit to the front of Remainder's values
			if ( RemainderDigit != 0 )
			{
				if ( Remainder.IsZero() )
				{
					Remainder.values[ 0 ] = RemainderDigit;
				}
				else
				{
					Remainder.values.insert( Remainder.values.begin(), RemainderDigit );
				}
			}

			// Calculate QuotientDigit using the RunFunction based on binary search
			digit_type QuotientDigit = ( B <= Remainder ) ? RunFunction( Remainder, B ) : 0;

			if ( Quotient.IsZero() )
			{
				Quotient.values[ 0 ] = QuotientDigit;
			}
			else
			{
				// Insert QuotientDigit to the front of Quotient's values
				Quotient.values.insert( Quotient.values.begin(), QuotientDigit );
			}
		}

		// Return the result as a pair of Quotient and Remainder
		return std::pair<BigInteger, BigInteger>( Quotient, Remainder );
	}

	BigInteger& BigInteger::ShortDivision( const BigInteger& other, BigInteger* r )
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
		uint64_t				t, k = 0;

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

	BigInteger& BigInteger::DonaldKnuthLongDivision( const BigInteger& other, BigInteger* r )
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

	std::pair<BigInteger, BigInteger> BigInteger::NewtonIterationDivision( const BigInteger& divisor ) const
	{
		// Check if the dividend is smaller than the divisor, return 0 as quotient and the dividend as remainder.
		if ( *this < divisor )
		{
			return std::make_pair( BigInteger( 0 ), *this );
		}

		// Calculate the size of the dividend and divisor.
		size_t dividend_size = this->Size();
		size_t divisor_size = divisor.Size();

		// Determine the initial approximation of the quotient size.
		size_t k = dividend_size - divisor_size + 5;		  // k is the number of digits to shift.
		size_t k2 = k > divisor_size ? 0 : divisor_size - k;  // k2 is the number of digits to shift in the divisor.

		// Shift the divisor to the right by k2 bits to approximate the reciprocal.
		BigInteger approximate_divisor = divisor.RightShiftBlock( k2 );
		if ( k2 != 0 )
		{
			approximate_divisor += 1;  // Adjust for the shift.
		}

		// Calculate the size of the intermediate result.
		size_t intermediate_size = k + approximate_divisor.Size();

		// Multiply the dividend by the approximate reciprocal of the divisor.
		BigInteger product = ( *this ) * approximate_divisor.DivisionInvertReciprocal( intermediate_size );

		// Shift the product to the right by the required number of bits.
		BigInteger quotient = product.RightShiftBlock( intermediate_size + k2 );

		// Calculate the initial remainder.
		BigInteger remainder = ( *this ) - quotient * divisor;

		// Iterate to refine the quotient until the remainder is less than the divisor.
		while ( remainder >= divisor )
		{
			quotient += 1;
			remainder -= divisor;
		}

		// Return the final quotient and remainder.
		return { quotient, remainder };
	}

	BigInteger& BigInteger::operator*=( const BigInteger& other )
	{
		if ( this == &other )
		{
			BigInteger self( other );
			*this *= self;
			return *this;
		}

		if ( IsZero() || other.IsZero() )
		{
			sign = 1;
			values.resize( 1 );
			values[ 0 ] = 0;
			return *this;
		}
		sign *= other.sign;

		/*
			BASE <= Karatsuba <= FHT
		*/

		if ( std::min( values.size(), other.values.size() ) < BASE_MULTIPLY_LIMIT )
		{
			return BaseMultiplication( other );
		}
		else if ( ( values.size() + other.values.size() ) < FHT_MULTIPLY_LIMIT )
		{
			if ( *this == other )
				return FHTSquare();
			else
				return FHTMultiplication( other );
		}
		else
		{
			return KaratsubaMultiplication( other );
		}
	}

	BigInteger& BigInteger::operator/=( const BigInteger& other )
	{
		if ( this == &other )
		{
			BigInteger self( other );
			*this /= self;
			return *this;
		}

		sign *= other.sign;
		Divide( other );
		return *this;
	}

	BigInteger& BigInteger::operator%=( const BigInteger& other )
	{
		BigInteger r;
		Divide( other, &r );
		values = std::move( r.values );
		return *this;
	}

	BigInteger& BigInteger::Sqrt()
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

	BigInteger& BigInteger::PowerWithModulo( const BigInteger& exponent, const BigInteger& modulo )
	{
		const BigInteger ZERO = 0;
		const BigInteger ONE = 1;
		BigInteger		 result = ONE;
		BigInteger		 base = *this;
		BigInteger		 exponent_value = exponent;

		if ( modulo == ONE )
		{
			*this = ZERO;
			return *this;
		}

		size_t index_bits = exponent_value.BitSize();
		for ( size_t i = 0; i < index_bits; i++ )
		{
			if ( exponent_value.GetBit( i ) )
			{
				result *= base;
				result %= modulo;
			}

			base.Power( 2 ) %= modulo;
		}

		values = std::move( result.values );
		return *this;
	}

	BigInteger& BigInteger::MultiplyWithModulo( const BigInteger& other, const BigInteger modulo )
	{
		const BigInteger ZERO = 0;
		const BigInteger ONE = 1;
		const BigInteger TWO = 2;

		if ( modulo == ONE )
		{
			*this = ZERO;
			return *this;
		}

		BigInteger result = ZERO;
		BigInteger self = *this;
		BigInteger other_value = other;

		while ( !other_value.IsZero() )
		{
			if ( !other_value.IsEven() )
			{
				result += self;
				result %= modulo;
			}

			self *= TWO;
			self %= modulo;

			if ( other_value.IsPowerOfTwo() )
				other_value >>= 1;
			else
				other_value /= TWO;
		}

		values = std::move( result.values );
		return *this;
	}

	bool BigInteger::IsEven() const
	{
		return ( values[ 0 ] & 1 ) == 0;
	}

	bool BigInteger::IsPowerOfTwo() const
	{
		if ( this->IsZero() )
			return false;

		//Is the current sign bit state negative?
		if ( this->sign == -1 )
			return false;

		//If the last binary bit is 0, it must be an even number, otherwise it is an odd number.
		if ( ( values[ 0 ] & 1 ) == 0 )
		{
			const BigInteger ONE( 1 );
			const BigInteger TEMPORARY = *this & ( *this - ONE );

			//Performs, on the copied entire array data, bitwise operations on whether the binary is a power of 2 or not.
			return TEMPORARY.values == std::vector<digit_type>( TEMPORARY.values.size(), 0 );
		}

		return false;
	}

	bool BigInteger::IsZero() const
	{
		if ( values.size() == 1 )
		{
			return values.size() <= 1 && values[ 0 ] == 0;
		}

		return this->values == std::vector<digit_type>( this->values.size(), 0 );
	}

	bool BigInteger::IsNegative() const
	{
		if ( sign == 1 )
			return false;
		else if ( sign == -1 )
			return true;

		return false;
	}

	size_t BigInteger::Size() const
	{
		return values.size();
	}

	void BigInteger::CountLeadingZeros( size_t& result ) const
	{
		if ( values.empty() )
			return;

		digit_type value = 0;
		for ( auto iter = values.rbegin(); iter != values.rend(); iter++ )
		{
			value = *iter;
			if ( value == 0 )
			{
				result += EXPONENT;
			}
			else
			{
				result += LeadingZeros( value );
				break;
			}
		}
	}

	uint16_t BigInteger::LeadingZeros( digit_type x )
	{
		if ( x == 0 )
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

	size_t BigInteger::BitSize() const
	{
		uint64_t   count = 0;
		digit_type high = values[ values.size() - 1 ];
		while ( high != 0 )
		{
			high >>= 1;
			count += 1;
		}

		if ( count )
			return ( values.size() - 1 ) * EXPONENT + count;
		return ( values.size() ) * EXPONENT;
	}

	bool BigInteger::GetBit( size_t bit_position ) const
	{
		size_t wrapper_index = bit_position / EXPONENT;
		size_t bit_index = bit_position - ( wrapper_index * EXPONENT );
		return values[ wrapper_index ] & ( digit_type( 1 ) << bit_index );
	}

	void BigInteger::SetBit( size_t bit_position )
	{
		size_t wrapper_index = bit_position / EXPONENT;
		size_t bit_index = bit_position - ( wrapper_index * EXPONENT );
		values[ wrapper_index ] = values[ wrapper_index ] | ( digit_type( 1 ) << bit_index );
	}

	void BigInteger::SetBit( bool value, size_t bit_position )
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

	void BigInteger::ExtendedLeadingBits( size_t count, bool value )
	{
		// Calculate the number of digits needed to accommodate the new bits.
		size_t newDigits = ( count + EXPONENT - 1 ) / EXPONENT;

		const digit_type BIT_ZEROS = digit_type( 0 ) & BASE;
		const digit_type BIT_ONES = ~digit_type( 0 ) & BASE;

		// If needed, add new digits to the beginning of the values vector.
		if ( newDigits > 0 )
		{
			values.insert( values.begin(), newDigits, value ? BIT_ONES : BIT_ZEROS );
		}

		// Adjust the bits within the existing most significant digit.
		size_t remainingBits = count % EXPONENT;
		if ( remainingBits > 0 )
		{
			digit_type& msb = values.back();
			if ( value )
			{
				msb |= ( BIT_ONES >> ( EXPONENT - remainingBits ) );
			}
			else
			{
				msb &= ( BIT_ZEROS << remainingBits ) | ( BIT_ONES >> ( EXPONENT - remainingBits ) );
			}
		}
	}

	void BigInteger::SqueezeLeadingBits( size_t count )
	{
		const digit_type BIT_ZEROS = digit_type( 0 ) & BASE;
		const digit_type BIT_ONES = ~digit_type( 0 ) & BASE;

		// Calculate the number of digits needed to remove the specified bits.
		size_t removedDigits = count / EXPONENT;

		// If needed, remove digits from the beginning of the values vector.
		if ( removedDigits > 0 )
		{
			if ( removedDigits >= values.size() )
			{
				// Remove all digits if necessary.
				values.clear();

				// Keep the meaning of this large integer as the number 0
				values.push_back( 0 );
			}
			else
			{
				values.erase( values.begin(), values.begin() + removedDigits );
			}
		}

		// Adjust the remaining bits within the new most significant digit.
		size_t remainingBits = count % EXPONENT;
		if ( remainingBits > 0 && !values.empty() )
		{
			digit_type& msb = values.front();
			msb &= ( BIT_ONES >> remainingBits );
		}
	}

	BigInteger& BigInteger::FromInt( uint64_t value )
	{
		if ( value == 0 )
		{
			values.resize( 1 );
			values[ 0 ] = 0;
			return *this;
		}
		else if ( value < 0 )
		{
			value *= -1;
			sign = -1;
		}
		if ( value < BASE )
		{
			values.push_back( value );
			return *this;
		}

		size_t count = size_t( log10( value ) / log10( BASE ) + 1 );
		values.reserve( count );
		while ( value > 0 )
		{
			values.push_back( value % BASE );
			value /= BASE;
		}
	}

	int64_t BigInteger::ToInt() const
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

	uint64_t BigInteger::ToUnsignedInt() const
	{
		uint64_t result = 0;
		uint64_t power = 1;
		uint64_t base = BASE;
		size_t	 size = values.size();

		for ( size_t i = 0; i < size; ++i )
		{
			result += values[ i ] * power;
			power *= base;
		}

		return result;
	}

	BigInteger& BigInteger::FromUnsignedInt( uint64_t value )
	{
		const uint64_t base = BASE;
		const uint64_t MaxUint64 = std::numeric_limits<uint64_t>::max();

		if ( value <= MaxUint64 / base )
		{
			return BigInteger( static_cast<uint64_t>( value ) );
		}
		else
		{
			std::vector<digit_type> values;
			while ( value > 0 )
			{
				values.push_back( static_cast<digit_type>( value % base ) );
				value /= base;
			}

			this->values = values;
			return *this;
		}
	}
	BigInteger BigInteger::operator&( const BigInteger& other ) const
	{
		BigInteger result = std::max( this->values.size(), other.values.size() ) ? ( *this ) : other;
		size_t	   digit_size = std::min( this->values.size(), other.values.size() );
		for ( size_t i = 0; i < digit_size; i++ )
		{
			result.values[ i ] = this->values[ i ] & other.values[ i ];
		}

		return result;
	}

	BigInteger& BigInteger::operator&=( const BigInteger& other )
	{
		*this = *this & other;
		return *this;
	}

	BigInteger BigInteger::operator|( const BigInteger& other ) const
	{
		BigInteger result = std::max( this->values.size(), other.values.size() ) ? ( *this ) : other;
		size_t	   digit_size = std::min( this->values.size(), other.values.size() );
		for ( size_t i = 0; i < digit_size; i++ )
		{
			result.values[ i ] = this->values[ i ] | other.values[ i ];
		}

		return result;
	}

	BigInteger& BigInteger::operator|=( const BigInteger& other )
	{
		*this = *this | other;
		return *this;
	}

	BigInteger BigInteger::operator~() const
	{
		BigInteger result = *this;
		for ( size_t i = 0; i < this->values.size(); i++ )
		{
			result.values[ i ] = ~( this->values[ i ] );
			result.values[ i ] <<= EXPONENT;
			result.values[ i ] >>= EXPONENT;
		}

		return result;
	}

	BigInteger BigInteger::operator^( const BigInteger& other ) const
	{
		BigInteger result = std::max( this->values.size(), other.values.size() ) ? ( *this ) : other;
		size_t	   digit_size = std::min( this->values.size(), other.values.size() );
		for ( size_t i = 0; i < digit_size; i++ )
		{
			result.values[ i ] = this->values[ i ] ^ other.values[ i ];
		}

		return result;
	}

	BigInteger& BigInteger::operator^=( const BigInteger& other )
	{
		*this = *this ^ other;
		return *this;
	}

	BigInteger& BigInteger::operator<<=( const uint32_t shift )
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
				BitLeftShiftCore( ( EXPONENT / 2 ) );

				shift_amount -= ( EXPONENT / 2 );
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

	BigInteger& BigInteger::operator>>=( const uint32_t shift )
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
				BitRightShiftCore( ( EXPONENT / 2 ) );

				shift_amount -= ( EXPONENT / 2 );
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

	BigInteger BigInteger::BitRotateLeft( const BigInteger& bits, uint32_t shift, const uint32_t reference_bit_capacity )
	{
		// Check if the BigInteger is zero, if so, no rotation is needed
		if ( bits.IsZero() || reference_bit_capacity == 0 )
			return bits;

		shift %= reference_bit_capacity;

#if 1

		// Convert the BigInteger to a binary string representation
		// The binary string length based on the reference_bit_capacity
		std::string binary_string = bits.ToBinaryString( reference_bit_capacity );

		// Left rotate the binary string
		std::rotate( binary_string.begin(), binary_string.begin() + shift, binary_string.end() );

		BigInteger result;

		// Update the BigInteger with the rotated binary string
		result.FromString( binary_string, 2 );

		// Ensure the BigInteger has the desired bit capacity
		if ( reference_bit_capacity / EXPONENT != result.values.size() )
		{
			size_t bit_size = result.BitSize();
			if ( bit_size < reference_bit_capacity )
				result.ExtendedLeadingBits( reference_bit_capacity - bit_size, false );
		}

#else

		BigInteger left = *this << shift;
		BigInteger right = *this >> shift;

		size_t bit_size = left.BitSize();
		if ( bit_size < reference_bit_capacity )
		{
			left.ExtendedLeadingBits( reference_bit_capacity - bit_size, false );
		}
		bit_size = right.BitSize();
		if ( bit_size < reference_bit_capacity )
		{
			right.ExtendedLeadingBits( reference_bit_capacity - bit_size, false );
		}

		BigInteger result = left | right;
		bit_size = result.BitSize();
		if ( bit_size < reference_bit_capacity )
		{
			this->ExtendedLeadingBits( reference_bit_capacity - bit_size, false );
		}

#endif

		return result;
	}

	BigInteger BigInteger::BitRotateRight( const BigInteger& bits, uint32_t shift, const uint32_t reference_bit_capacity )
	{
		// Check if the BigInteger is zero, if so, no rotation is needed
		if ( bits.IsZero() || reference_bit_capacity == 0 )
			return bits;

		shift %= reference_bit_capacity;

#if 1

		// Convert the BigInteger to a binary string representation
		// The binary string length based on the reference_bit_capacity
		std::string binary_string = bits.ToBinaryString( reference_bit_capacity );

		// Right rotate the binary string
		std::rotate( binary_string.rbegin(), binary_string.rbegin() + shift, binary_string.rend() );

		BigInteger result;

		// Update the BigInteger with the rotated binary string
		result.FromString( binary_string, 2 );

		// Ensure the BigInteger has the desired bit capacity
		if ( reference_bit_capacity / EXPONENT != result.values.size() )
		{
			size_t bit_size = result.BitSize();
			if ( bit_size < reference_bit_capacity )
				result.ExtendedLeadingBits( reference_bit_capacity - bit_size, false );
		}

#else

		BigInteger right = *this >> shift;
		BigInteger left = *this << shift;

		size_t bit_size = right.BitSize();
		if ( bit_size < reference_bit_capacity )
		{
			right.ExtendedLeadingBits( reference_bit_capacity - bit_size, false );
		}
		bit_size = left.BitSize();
		if ( bit_size < reference_bit_capacity )
		{
			left.ExtendedLeadingBits( reference_bit_capacity - bit_size, false );
		}

		BigInteger result = left | right;
		bit_size = result.BitSize();
		if ( bit_size < reference_bit_capacity )
		{
			this->ExtendedLeadingBits( reference_bit_capacity - bit_size, false );
		}

#endif

		return result;
	}

	digit_type& BigInteger::operator[]( const size_t index )
	{
		return values[ index ];
	}

	BigInteger& BigInteger::operator=( const BigInteger& other )
	{
		if ( this == &other )
			return *this;

		values = other.values;
		sign = other.sign;
		return *this;
	}

	bool BigInteger::operator&&( const BigInteger& other ) const
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

	bool BigInteger::operator||( const BigInteger& other ) const
	{
		bool left = this->IsZero();
		bool right = other.IsZero();

		if ( ( left == true ) && ( right == true ) )
			return false;
		else if ( ( left == false ) || ( right == false ) )
			return true;
	}

	bool BigInteger::operator!() const
	{
		return !this->IsZero();
	}

	BigInteger BigInteger::Abs() const
	{
		BigInteger copy = *this;

		if ( this->IsNegative() )
		{
			return -copy;
		}

		return copy;
	}

	BigInteger BigInteger::ModuloArithmetic( ArithmeticMode mode, const BigInteger& a, const BigInteger& b, const BigInteger& modulus )
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

	std::string BigInteger::ToString() const
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

	std::string BigInteger::ToBinaryString( const uint32_t reference_bit_capacity, bool have_leading_zeros ) const
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
				{
					if ( i < reference_bit_capacity )
						number_string[ i ] = '1';
					else
						break;
				}
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

	std::string BigInteger::ToString( uint32_t base_value ) const
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
		BigInteger rem;
		while ( !a.IsZero() )
		{
			// Calculate the current digit, divide the number in the current bit by 2 from a.
			a.ShortDivision( STRING_BASE, &rem );
			uint64_t digit = rem.values[ 0 ];

			// Convert digits(number) to characters
			number_string += ( digit < 10 ? '0' + digit : 'A' + digit - 10 );

			// Update m for possible subsequent negative symbols
			m++;
		}
		if ( sign == -1 )
			number_string += '-';

		std::reverse( number_string.begin(), number_string.end() );
		return number_string;
	}

	void BigInteger::FromString( const std::string& number_string )
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

	void BigInteger::FromString( const std::string& number_string, uint32_t base_value )
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
			result.ExtendedLeadingBits( count - result.BitSize(), false );

			for ( i; i < count; ++i )
			{
				result.SetBit( number_string[ i ] == '1' ? true : false, ( count - 1 ) - i );
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

	void BigInteger::Print( bool base10 ) const
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

	size_t BigInteger::SipHash( const BigInteger& Integer, std::vector<uint8_t>* keys ) const
	{
		if ( Integer.values.empty() )
		{
			return 0;
		}

		std::vector<uint8_t> datas( Integer.values.size() * sizeof( digit_type ), 0 );
		::memcpy( &datas, Integer.values.data(), Integer.values.size() * sizeof( digit_type ) );

		uint64_t state = 0x736f6d6570736575UL;
		uint64_t state2 = 0x646f72616e646f6dUL;
		uint64_t key = 0;
		uint64_t key2 = 0;

		if ( keys != nullptr && !( keys->empty() ) )
		{

			if ( keys->size() < sizeof( uint64_t ) )
			{
				::memcpy( &key, keys->data(), keys->size() * sizeof( uint8_t ) );
			}
			else if ( keys->size() >= sizeof( uint64_t ) && keys->size() < sizeof( uint64_t ) * 2 )
			{
				::memcpy( &key, keys->data(), sizeof( uint64_t ) );
				::memcpy( &key2, keys->data() + sizeof( uint64_t ), keys->size() - sizeof( uint64_t ) );
			}
			else if ( keys->size() > sizeof( uint64_t ) * 2 )
			{
				::memcpy( &key, keys->data(), sizeof( uint64_t ) );
				::memcpy( &key2, keys->data() + sizeof( uint64_t ), sizeof( uint64_t ) );
			}

			state ^= key;
			state2 ^= key2;
		}

		uint64_t v0 = state;
		uint64_t v1 = state2;
		uint64_t v2 = 0x6c7967656e657261 ^ key;
		uint64_t v3 = 0x7465646279746573 ^ key2;

		uint64_t		block = 0;
		const uint64_t* start_block = reinterpret_cast<uint64_t*>( datas[ 0 ] );
		const uint64_t* end_block = start_block + ( datas.size() & ~7 );

		//Update blocks with One-way compress function
		uint64_t* block_pointer = ( uint64_t* )start_block;
		while ( block_pointer < end_block )
		{
			block = *block_pointer;

			v3 ^= block;
			for ( size_t round = 0; round < 2; round++ )
			{
				v0 += v1;
				v2 += v3;
				v1 = v1 << 13 | v1 >> 51;
				v3 = v3 << 16 | v3 >> 48;
				v1 ^= v0;
				v3 ^= v2;
				v0 = v0 << 32 | v0 >> 32;
				v2 += v1;
				v0 += v3;
				v1 = v1 << 17 | v1 >> 47;
				v3 = v3 << 21 | v3 >> 43;
				v1 ^= v2;
				v3 ^= v0;
				v2 = v2 << 32 | v2 >> 32;
			}
			v0 ^= block;
		}

		// Update end blocks with One-way compress function
		std::vector<uint8_t> end_block_bytes = std::vector<uint8_t>( reinterpret_cast<uint8_t*>( &end_block ), reinterpret_cast<uint8_t*>( &end_block + 1 ) );
		uint64_t			 last_block = 0;
		switch ( end_block_bytes.size() )
		{
		case 7:
			last_block |= static_cast<uint64_t>( end_block_bytes[ 6 ] ) << 48;
			[[fallthrough]];
		case 6:
			last_block |= static_cast<uint64_t>( end_block_bytes[ 5 ] ) << 40;
			[[fallthrough]];
		case 5:
			last_block |= static_cast<uint64_t>( end_block_bytes[ 4 ] ) << 32;
			[[fallthrough]];
		case 4:
			last_block |= static_cast<uint64_t>( end_block_bytes[ 3 ] ) << 24;
			[[fallthrough]];
		case 3:
			last_block |= static_cast<uint64_t>( end_block_bytes[ 2 ] ) << 16;
			[[fallthrough]];
		case 2:
			last_block |= static_cast<uint64_t>( end_block_bytes[ 1 ] ) << 8;
			[[fallthrough]];
		case 1:
			last_block |= static_cast<uint64_t>( end_block_bytes[ 0 ] );
		}

		v3 ^= last_block;
		for ( size_t round = 0; round < 2; round++ )
		{
			v0 += v1;
			v2 += v3;
			v1 = v1 << 13 | v1 >> 51;
			v3 = v3 << 16 | v3 >> 48;
			v1 ^= v0;
			v3 ^= v2;
			v0 = v0 << 32 | v0 >> 32;
			v2 += v1;
			v0 += v3;
			v1 = v1 << 17 | v1 >> 47;
			v3 = v3 << 21 | v3 >> 43;
			v1 ^= v2;
			v3 ^= v0;
			v2 = v2 << 32 | v2 >> 32;
		}
		v0 ^= last_block;

		// Finalization.
		v2 ^= 0xff;
		for ( size_t round = 0; round < 4; round++ )
		{
			v0 += v1;
			v2 += v3;
			v1 = v1 << 13 | v1 >> 51;
			v3 = v3 << 16 | v3 >> 48;
			v1 ^= v0;
			v3 ^= v2;
			v0 = v0 << 32 | v0 >> 32;
			v2 += v1;
			v0 += v3;
			v1 = v1 << 17 | v1 >> 47;
			v3 = v3 << 21 | v3 >> 43;
			v1 ^= v2;
			v3 ^= v0;
			v2 = v2 << 32 | v2 >> 32;
		}

		return v0 ^ v1 ^ v2 ^ v3;
	}

	BigInteger BigInteger::GCD( BigInteger a, BigInteger b )
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

	void BigInteger::EGCD( const BigInteger& a, const BigInteger& b, BigInteger* gcd, BigInteger* co1, BigInteger* co2 )
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

	BigInteger BigInteger::LCM( const BigInteger& a, const BigInteger& b )
	{
		BigInteger result = a * b;
		result.sign = 1;
		result /= GCD( a, b );
		return result;
	}

	BigInteger BigInteger::ModuloInverse( const BigInteger& a, const BigInteger& b )
	{
		BigInteger gcd, x;

		EGCD( a, b, &gcd, &x, nullptr );
		return ( x % b + b ) % b;
	}

	BigInteger BigInteger::PollardRho( const BigInteger& n )
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

	BigInteger BigInteger::RandomGenerateNBit( size_t n )
	{
		std::uniform_int_distribution<uint16_t> dist01( 0, 1 );
		BigInteger								ud_prng( 0 );
		size_t									digit_type_count = n % EXPONENT == 0 ? n / EXPONENT : n / EXPONENT + 1;
		ud_prng.values = std::vector<digit_type>( digit_type_count, 0 );
		ud_prng.SetBit( n - 1 );
		for ( size_t i = 0; i < n; ++i )
		{
			if ( dist01( RandomEngine ) == 1 )
			{
				ud_prng.SetBit( true, i );
			}
		}

		return ud_prng;
	}

	BigInteger BigInteger::Log2() const
	{
		if (IsZero())
		{
			// Logarithm of 0 is undefined, return an error code.
			return -1;
		}
		if (*this == 1)
		{
			// Logarithm of 1 is 0.
			return 0;
		}

		BigInteger log2 = 0;
		BigInteger current = *this;
		while (current > 1)
		{
			// Shift current right by one bit.
			current >>= 1;
			log2++;
		}

		return log2;
	}

	BigInteger BigInteger::Modulo( const BigInteger& number, const BigInteger& modulus )
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
}  // namespace TwilightDream::BigInteger