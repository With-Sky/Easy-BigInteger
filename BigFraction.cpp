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

#include "BigFraction.hpp"

namespace TwilightDream::BigFraction
{
	BigFraction::BigFraction() : numerator( 0 ), denominator( 1 ), sign( 1 ) {}

	BigFraction::BigFraction( const BigFraction& other ) noexcept
		: numerator( other.numerator ), denominator( other.denominator ), sign( other.sign ),
		simplify_reduced(other.simplify_reduced ),
		PrecisionMode( other.PrecisionMode ), FixedPrecisionCount( other.FixedPrecisionCount )
	{}

	BigFraction::BigFraction( BigFraction&& other ) noexcept
		: numerator( other.numerator ), denominator( other.denominator ), sign( other.sign ),
		simplify_reduced(other.simplify_reduced ),
		PrecisionMode( other.PrecisionMode ), FixedPrecisionCount( other.FixedPrecisionCount )
	{}

	BigFraction::BigFraction( const BigInteger& numerator )
		: numerator( numerator ), denominator( 1 ), sign( 1 ) {}

	BigFraction::BigFraction( const BigInteger& numerator, const BigInteger& denominator )
		: numerator( numerator ), denominator( denominator ), sign( 1 )
	{
		if ( simplify_reduced )
		{
			this->ReduceSimplify();
		}
	}

	BigFraction::BigFraction( const BigInteger& numerator, const BigInteger& denominator, int32_t sign )
		: numerator( numerator ), denominator( denominator ), sign( sign )
	{
		if ( simplify_reduced )
		{
			this->ReduceSimplify();
		}

		this->sign = (this->sign >= 0 ? 1 : -1);
	}

	BigFraction::BigFraction( const std::string& complexString )
	{
		ComputeAndFromDecimalString( complexString );
	}

	void BigFraction::SetSimplifyReduced( bool value )
	{
		simplify_reduced = value;
	}

	bool BigFraction::IsNaN() const
	{
		if ( numerator.IsZero() && denominator.IsZero() )
		{
			return true;
		}

		return false;
	}

	bool BigFraction::IsInfinity() const
	{
		if ( !numerator.IsZero() && denominator.IsZero() )
		{
			return true;
		}

		return false;
	}

	bool BigFraction::IsZero() const
	{
		return numerator.IsZero() && !denominator.IsZero();
	}

	bool BigFraction::IsNegative() const
	{
		return sign == -1;
	}

	BigFraction::BigInteger BigFraction::GetNumerator() const
	{
		return numerator;
	}

	BigFraction::BigInteger BigFraction::GetDenominator() const
	{
		return denominator;
	}

	void BigFraction::SetNumerator( const BigInteger& number )
	{
		numerator = number;
	}

	void BigFraction::SetDenominator( const BigInteger& number )
	{
		denominator = number;
	}

	BigFraction BigFraction::GetFullPrecision() const
	{
		if(BigFractionFullPrecision.denominator.IsZero())
		{
			throw std::invalid_argument( "BigFractionFullPrecision Undefined: a denominator that is 0 is incorrect." );
		}

		return BigFractionFullPrecision;
	}

	void BigFraction::SetFullPrecision(BigInteger number)
	{
		if (number.IsZero())
		{
			return;
		}

		if (!((number % RADIX).IsZero()))
		{
			throw std::invalid_argument("BigFractionFullPrecision Undefined: denominators that are not multiples of 10 are incorrect.");
		}

		// Ensure that number is of the form 10^n
		BigInteger temp = number;
		while ((temp % RADIX).IsZero())
		{
			temp /= RADIX;
		}

		// After division, if temp is not 1, then number was not 10^n
		if (temp != BigInteger(1))
		{
			throw std::invalid_argument("BigFractionFullPrecision Undefined: denominator must be a power of 10.");
		}

		this->PrecisionMode = DecimalPrecisionMode::Full;

		// Set the denominator of BigFractionFullPrecision to number
		BigFractionFullPrecision.SetDenominator(number);

		// Check if the denominator is not zero before setting
		if (BigFractionFullPrecision.GetDenominator().IsZero())
		{
			throw std::invalid_argument("BigFractionFullPrecision Undefined: a denominator that is 0 is incorrect.");
		}
	}

	void BigFraction::ComputeAndFromDecimalString( const std::string& complexString )
	{
		// Parsing Complex String
		bool   isNegative = false;
		bool   hasFractionalPart = false;
		size_t integerPartEnd = 0;
		size_t fractionalPartStart = 0;

		// Determine if the number is negative
		if ( !complexString.empty() && complexString[ 0 ] == '-' )
		{
			isNegative = true;
			integerPartEnd = 1;
		}

		// Find the end of integer part
		while ( integerPartEnd < complexString.size() && std::isdigit( complexString[ integerPartEnd ] ) )
		{
			integerPartEnd++;
		}

		// Check if there is a fractional part
		if ( integerPartEnd < complexString.size() && complexString[ integerPartEnd ] == '.' )
		{
			hasFractionalPart = true;
			fractionalPartStart = integerPartEnd + 1;
		}

		// Parse integer part
		BigInteger integerPart( complexString.substr( isNegative, integerPartEnd - isNegative ) );

		// Parse fractional part if exists
		BigInteger fractionalPart = 0;
		BigInteger fractionalMultiplier = 1;
		if ( hasFractionalPart )
		{
			uint64_t LoopCount = 0;

			switch ( this->PrecisionMode )
			{
				case DecimalPrecisionMode::Fixed:
					LoopCount = this->FixedPrecisionCount;
					break;
				case DecimalPrecisionMode::Full:
					LoopCount = complexString.size();
					break;
				default:
					break;
			}

			for ( size_t i = fractionalPartStart; i < LoopCount; ++i )
			{
				if ( !std::isdigit( complexString[ i ] ) )
				{
					// Invalid character found in fractional part
					// You may want to handle this error case appropriately
					return;
				}
				fractionalPart = fractionalPart * RADIX + ( complexString[ i ] - '0' );
				fractionalMultiplier *= RADIX;
			}
		}

		// Combine integer and fractional parts into the numerator
		numerator = integerPart * fractionalMultiplier + fractionalPart;
		denominator = fractionalMultiplier;

		// Apply the sign
		if ( isNegative )
		{
			numerator *= -1;
		}

		// Simplify the fraction to reduce it to its simplest form
		if ( simplify_reduced )
		{
			ReduceSimplify();
		}
	}

	std::string BigFraction::ComputeAndToDecimalString() const
	{
		if ( IsNaN() )
		{
			return "NaN";
		}

		std::string signString = ( sign == -1 ) ? "-" : "";

		BigFraction abs = this->Abs();
		BigInteger absNumerator = abs.GetNumerator();
		BigInteger absDenominator = abs.GetDenominator();

		// Reduce the fraction to its simplest form if necessary
		if ( absNumerator > 1 && absDenominator > 1 )
		{
			BigInteger gcd = BigInteger::GCD( absNumerator, absDenominator );
			absNumerator /= gcd;
			absDenominator /= gcd;
		}

		// Compute the integer part
		BigInteger	integerPart = absNumerator / absDenominator;
		std::string result = signString + integerPart.ToString( 10 ) + ".";

		switch ( this->PrecisionMode )
		{
			case DecimalPrecisionMode::Fixed:
			{
				if ( this->FixedPrecisionCount == 0 )
				{
					// Return only the integer part
					result += "0";
					return result;
				}

				BigInteger fractionalPart = absNumerator % absDenominator;
				for ( size_t i = 0; i < this->FixedPrecisionCount; ++i )
				{
					fractionalPart *= RADIX;
					BigInteger nextDigit = fractionalPart / absDenominator;
					fractionalPart %= absDenominator;
					result += nextDigit.ToString( 10 );
					if ( fractionalPart == 0 )
					{
						// Stop if the fractional part becomes zero before reaching the desired precision
						break;
					}
				}

				break;
			}
			
			case DecimalPrecisionMode::Full:
			{
				if (absDenominator.IsZero())
				{
					throw std::invalid_argument("Denominator cannot be zero.");
				}

				BigInteger fractionalPart = absNumerator % absDenominator;
				BigFraction precision = GetFullPrecision(); // Assume GetFullPrecision() returns the minimal precision as a fraction

				// Determine the required decimal precision based on the denominator of the precision fraction.
				// We will count how many times the denominator can be divided by 10 (RADIX).
				BigInteger precisionDenominator = precision.GetDenominator();
				BigInteger requiredPrecisionCount = 0;

				// Count the number of decimal places needed by repeatedly dividing the precision denominator by RADIX (10).
				while ((precisionDenominator % RADIX).IsZero())
				{
					precisionDenominator /= RADIX;
					requiredPrecisionCount++;
				}

				bool hasFractionalPart = false;  // Flag to indicate if there is a fractional part

				// Generate the decimal representation by calculating each digit up to the required precision.
				for (BigInteger i = 0; i < requiredPrecisionCount; ++i)
				{
					if (fractionalPart == 0)
					{
						break; // The fractional part has been completely consumed, no more digits needed
					}

					hasFractionalPart = true;
					fractionalPart *= RADIX; // Shift the fractional part to the next decimal place
					BigInteger nextDigit = fractionalPart / absDenominator;
					fractionalPart %= absDenominator; // Update the fractional part for the next iteration
					result += nextDigit.ToString(10); // Append the next digit to the result string
				}

				// If no fractional digits were added, append a '0' to ensure proper decimal formatting (e.g., "1.0").
				if (!hasFractionalPart)
				{
					result += "0";
				}

				break;
			}
			default:
				break;
		}

		return result;
	}

	BigFraction& BigFraction::operator=( const BigFraction& other )
	{
		if ( this != &other )
		{
			this->numerator = other.numerator;
			this->denominator = other.denominator;
			this->sign = other.sign;
			this->simplify_reduced = this->simplify_reduced;
			this->PrecisionMode = other.PrecisionMode;
			this->FixedPrecisionCount = other.FixedPrecisionCount;
		}
		return *this;
	}

	BigFraction& BigFraction::operator=( BigFraction&& other )
	{
		if ( this != &other )
		{
			this->numerator = std::move(other.numerator);
			this->denominator = std::move(other.denominator);
			this->sign = std::move(other.sign);
			this->simplify_reduced = std::move(this->simplify_reduced);
			this->PrecisionMode = std::move(other.PrecisionMode);
			this->FixedPrecisionCount = std::move(other.FixedPrecisionCount);
		}
		return *this;
	}

	BigFraction& BigFraction::operator+=( const BigFraction& other )
	{
		if ( this->IsNaN() || other.IsNaN() )
		{
			this->sign = 1;
			this->numerator = 0;
			this->denominator = 1;
			return *this;
		}

		if ( IsZero() )
		{
			if ( other.IsZero() )
			{
				sign = 1;
				numerator = 0;
				denominator = 1;
				return *this;
			}
			else if ( other.IsNegative() )
			{
				sign = -1;
				numerator = other.numerator.Abs();
				denominator = other.denominator.Abs();
				return *this;
			}
			else
			{
				sign = 1;
				numerator = other.numerator.Abs();
				denominator = other.denominator.Abs();
				return *this;
			}
		}
		else if ( other.IsZero() )
		{
			return *this;
		}

		// Use BigSignedInteger for intermediate calculations
		BigSignedInteger a(this->numerator * other.denominator, this->sign < 0);
		BigSignedInteger b(other.numerator * this->denominator, other.sign < 0);
		BigSignedInteger resultNumeratorSigned = a + b;
		BigInteger resultDenominatorSigned = this->denominator * other.denominator;

		// Determine the sign and absolute values
		sign = resultNumeratorSigned.IsNegative() ? -1 : 1;
		numerator = static_cast<BigInteger>(resultNumeratorSigned.Abs());
		denominator = resultDenominatorSigned;

		if (simplify_reduced)
		{
			ReduceSimplify();
		}

		return *this;
	}

	BigFraction& BigFraction::operator+=( const BigInteger& other )
	{
		if ( this->IsNaN() )
		{
			this->sign = 1;
			this->numerator = 0;
			this->denominator = 1;
			return *this;
		}

		if ( this->IsZero() )
		{
			this->sign = 1;
			this->numerator = other;
			this->denominator = 1;
			return *this;
		}
		else if ( other.IsZero() )
		{
			return *this;
		}

		// Use BigSignedInteger for intermediate calculations
		BigSignedInteger a(this->numerator, this->sign < 0);
		BigSignedInteger b(other, false);  // other is always positive

		BigSignedInteger resultNumeratorSigned;
		if (this->sign == 1)
		{
			resultNumeratorSigned = a + b;
		}
		else
		{
			resultNumeratorSigned = a - b;
		}

		// Determine the sign and absolute values
		sign = resultNumeratorSigned.IsNegative() ? -1 : 1;
		numerator = static_cast<BigInteger>(resultNumeratorSigned.Abs());
		// Denominator remains unchanged

		if ( simplify_reduced )
		{
			ReduceSimplify();
		}

		return *this;
	}

	BigFraction& BigFraction::operator-=( const BigFraction& other )
	{
		if ( this->IsNaN() || other.IsNaN() )
		{
			this->sign = 1;
			this->numerator = 0;
			this->denominator = 1;
			return *this;
		}

		if ( this->IsZero() )
		{
			if ( other.IsZero() )
			{
				this->sign = 1;
				this->numerator = 0;
				this->denominator = 1;
				return *this;
			}
			else if ( other.IsNegative() )
			{
				this->sign = 1;
				this->numerator = other.numerator.Abs();
				this->denominator = other.denominator.Abs();
				return *this;
			}
			else
			{
				this->sign = -1;
				this->numerator = other.numerator.Abs();
				this->denominator = other.denominator.Abs();
				return *this;
			}
		}
		else if ( other.IsZero() )
		{
			return *this;
		}

		// Use BigSignedInteger for intermediate calculations
		BigSignedInteger a(this->numerator * other.denominator, this->sign < 0);
		BigSignedInteger b(other.numerator * this->denominator, other.sign < 0);
		BigSignedInteger resultNumeratorSigned = a - b;
		BigInteger resultDenominator = this->denominator * other.denominator;

		// Determine the sign and absolute values
		sign = resultNumeratorSigned.IsNegative() ? -1 : 1;
		numerator = static_cast<BigInteger>(resultNumeratorSigned.Abs());
		denominator = resultDenominator.Abs();

		if ( simplify_reduced )
		{
			ReduceSimplify();
		}

		return *this;
	}

	BigFraction& BigFraction::operator-=( const BigInteger& other )
	{
		if ( this->IsNaN() )
		{
			this->sign = 1;
			this->numerator = 0;
			this->denominator = 1;
			return *this;
		}

		if ( this->IsZero() )
		{
			if ( other.IsZero() )
			{
				this->sign = 1;
				this->numerator = 0;
				this->denominator = 1;
				return *this;
			}
			else if ( other.IsNegative() )
			{
				this->sign = 1;
				this->numerator = other.Abs();
				this->denominator = 1;
				return *this;
			}
			else
			{
				this->sign = -1;
				this->numerator = other;
				this->denominator = 1;
				return *this;
			}
		}
		else if ( other.IsZero() )
		{
			return *this;
		}

		// Use BigSignedInteger for intermediate calculations
		BigSignedInteger a(this->numerator, this->sign < 0);
		BigSignedInteger b(other, false);

		BigSignedInteger resultNumeratorSigned;
		if ((this->sign > 0) || (this->sign < 0))
		{
			resultNumeratorSigned = a + b.Abs();
		}
		else
		{
			resultNumeratorSigned = a - b;
		}

		// Determine the sign and absolute values
		sign = resultNumeratorSigned.IsNegative() ? -1 : 1;
		numerator = static_cast<BigInteger>(resultNumeratorSigned.Abs());
		// Denominator remains unchanged

		if ( simplify_reduced )
		{
			ReduceSimplify();
		}

		return *this;
	}

	BigFraction& BigFraction::operator*=( const BigFraction& other )
	{
		if ( this->IsNaN() || other.IsNaN() )
		{
			this->sign = 1;
			this->numerator = 0;
			this->denominator = 1;
			return *this;
		}

		if ( other.IsZero() )
		{
			this->sign = 1;
			this->numerator = 0;
			this->denominator = 1;
			return *this;
		}
		if ( ( other.numerator > 0 && other.denominator > 0 ) && ( other.numerator == other.denominator ) )
		{
			return *this;
		}

		// Use BigSignedInteger for intermediate calculations
		BigSignedInteger signed_numerator(this->numerator, this->sign < 0);
		BigSignedInteger signed_denominator(this->denominator, this->sign < 0);
		BigSignedInteger other_signed_numerator(other.numerator, other.IsNegative());
		BigSignedInteger other_signed_denominator(other.denominator, other.IsNegative());

		signed_numerator *= other_signed_numerator;
		signed_denominator *= other_signed_denominator;

		this->sign = (signed_numerator.IsNegative() != other_signed_denominator.IsNegative()) ? -1 : 1;
		this->numerator = static_cast<BigInteger>(signed_numerator.Abs());
		this->denominator = static_cast<BigInteger>(signed_denominator.Abs());

		if ( simplify_reduced )
		{
			ReduceSimplify();
		}

		return *this;
	}

	BigFraction& BigFraction::operator*=( const BigInteger& other )
	{
		if ( other.IsZero() )
		{
			this->sign = 1;
			this->numerator = 0;
			this->denominator = 1;
			return *this;
		}
		if ( other == 1 )
		{
			return *this;
		}

		// Use BigSignedInteger for intermediate calculations
		BigSignedInteger signed_numerator(this->numerator, this->sign < 0);
		BigSignedInteger signed_other(other, false);

		signed_numerator *= signed_other;

		this->numerator = static_cast<BigInteger>(signed_numerator.Abs());
		this->sign = signed_numerator.IsNegative() ? -1 : 1;

		if ( simplify_reduced )
		{
			ReduceSimplify();
		}

		return *this;
	}

	BigFraction& BigFraction::operator/=( const BigFraction& other )
	{
		if ( this->IsNaN() || other.IsNaN() )
		{
			this->sign = 1;
			this->numerator = 0;
			this->denominator = 1;
			return *this;
		}

		if ( other.denominator.IsZero() )
		{
			throw std::runtime_error( "Division by zero error." );
		}

		// Use BigSignedInteger for intermediate calculations
		BigSignedInteger signed_numerator(this->numerator, this->sign < 0);
		BigSignedInteger signed_denominator(this->denominator, this->sign < 0);
		BigSignedInteger other_signed_numerator(other.numerator, other.IsNegative());
		BigSignedInteger other_signed_denominator(other.denominator, other.IsNegative());

		signed_numerator *= other_signed_denominator; // Multiply the numerator by the denominator of another fraction
		signed_denominator *= other_signed_numerator; // Multiply the denominator by the numerator of another fraction

		this->sign = (signed_numerator.IsNegative() != signed_denominator.IsNegative()) ? -1 : 1;
		this->numerator = static_cast<BigInteger>(signed_numerator.Abs());
		this->denominator = static_cast<BigInteger>(signed_denominator.Abs());

		if ( simplify_reduced )
		{
			ReduceSimplify();
		}

		return *this;
	}

	BigFraction& BigFraction::operator/=( const BigInteger& other )
	{
		if ( this->IsNaN() )
		{
			this->sign = 1;
			this->numerator = 0;
			this->denominator = 1;
			return *this;
		}

		if ( other.IsZero() )
		{
			throw std::runtime_error( "Division by zero error." );
		}
		
		// 使用 BigSignedInteger 进行中间计算
		BigSignedInteger signed_numerator(this->numerator, this->sign < 0);
		BigSignedInteger signed_denominator(this->denominator, this->sign < 0);
		BigSignedInteger other_signed(other, false);

		signed_denominator *= other_signed;

		this->sign = (signed_numerator.IsNegative() != signed_denominator.IsNegative()) ? -1 : 1;
		this->numerator = static_cast<BigInteger>(signed_numerator.Abs());
		this->denominator = static_cast<BigInteger>(signed_denominator.Abs());

		if ( simplify_reduced )
		{
			ReduceSimplify();
		}

		return *this;
	}

	BigFraction BigFraction::operator+( const BigFraction& other ) const
	{
		BigFraction result( *this );
		result += other;

		result.PrecisionMode = this->PrecisionMode;
		result.FixedPrecisionCount = this->FixedPrecisionCount;

		return result;
	}

	BigFraction BigFraction::operator+( const BigInteger& other ) const
	{
		BigFraction result( *this );
		result += other;

		result.PrecisionMode = this->PrecisionMode;
		result.FixedPrecisionCount = this->FixedPrecisionCount;

		return result;
	}

	BigFraction BigFraction::operator-( const BigFraction& other ) const
	{
		BigFraction result( *this );
		result -= other;

		result.PrecisionMode = this->PrecisionMode;
		result.FixedPrecisionCount = this->FixedPrecisionCount;

		return result;
	}

	BigFraction BigFraction::operator-( const BigInteger& other ) const
	{
		BigFraction result( *this );
		result -= other;

		result.PrecisionMode = this->PrecisionMode;
		result.FixedPrecisionCount = this->FixedPrecisionCount;

		return result;
	}

	BigFraction BigFraction::operator*( const BigFraction& other ) const
	{
		BigFraction result( *this );
		result *= other;

		result.PrecisionMode = this->PrecisionMode;
		result.FixedPrecisionCount = this->FixedPrecisionCount;

		return result;
	}

	BigFraction BigFraction::operator*( const BigInteger& other ) const
	{
		BigFraction result( *this );
		result *= other;

		result.PrecisionMode = this->PrecisionMode;
		result.FixedPrecisionCount = this->FixedPrecisionCount;

		return result;
	}

	BigFraction BigFraction::operator/( const BigFraction& other ) const
	{
		BigFraction result( *this );
		result /= other;

		result.PrecisionMode = this->PrecisionMode;
		result.FixedPrecisionCount = this->FixedPrecisionCount;

		return result;
	}

	BigFraction BigFraction::operator/( const BigInteger& other ) const
	{
		BigFraction result( *this );
		result /= other;

		result.PrecisionMode = this->PrecisionMode;
		result.FixedPrecisionCount = this->FixedPrecisionCount;

		return result;
	}

	bool BigFraction::operator<( const BigFraction& other ) const
	{
		BigInteger lhs = numerator * other.denominator;
		BigInteger rhs = other.numerator * denominator;
		return lhs < rhs;
	}

	bool BigFraction::operator<=( const BigFraction& other ) const
	{
		return !( *this > other );
	}

	bool BigFraction::operator>( const BigFraction& other ) const
	{
		BigInteger lhs = numerator * other.denominator;
		BigInteger rhs = other.numerator * denominator;
		return lhs > rhs;
	}

	bool BigFraction::operator>=( const BigFraction& other ) const
	{
		return !( *this < other );
	}

	bool BigFraction::operator==( const BigFraction& other ) const
	{
		return numerator == other.numerator && denominator == other.denominator && sign == other.sign;
	}

	bool BigFraction::operator!=( const BigFraction& other ) const
	{
		return !( *this == other );
	}

	BigFraction BigFraction::Abs() const
	{
		BigFraction result = BigFraction( numerator, denominator, 1 );
		
		result.PrecisionMode = this->PrecisionMode;
		result.FixedPrecisionCount = this->FixedPrecisionCount;
		return result;
	}

	BigFraction BigFraction::Reciprocal() const
	{
		BigFraction result;

		if ( IsNaN() )
		{
			result = BigFraction( 0, 0 );
		}
		else if ( IsZero() )
		{
			result = BigFraction( 0, ONE );
		}
		else
		{
			result = BigFraction( denominator, numerator, sign );
		}

		result.PrecisionMode = this->PrecisionMode;
		result.FixedPrecisionCount = this->FixedPrecisionCount;
		result.sign = this->sign;
		return result;
	}

	BigFraction BigFraction::Sqrt()
	{
		BigFraction result;

		if ( IsNaN() )
		{
			result = BigFraction( 0, ONE );
		}
		else if ( IsZero() )
		{
			result = BigFraction( 0, ONE );
		}
		else if ( IsNegative() )
		{
			result = BigFraction( 0, ONE );
		}
		else
		{
			result = BigFraction(BigInteger(numerator).Sqrt(), BigInteger(denominator).Sqrt());
		}

		result.PrecisionMode = this->PrecisionMode;
		result.FixedPrecisionCount = this->FixedPrecisionCount;
		return result;
	}

	BigFraction BigFraction::Cbrt()
	{
		BigFraction result;

		if ( IsNaN() )
		{
			result = BigFraction( 0, ONE );
		}
		else if ( IsZero() )
		{
			result = BigFraction( 0, ONE );
		}
		else if ( IsNegative() )
		{
			result = BigFraction( 0, ONE );
		}
		else
		{
			result = BigFraction(BigInteger(numerator).Cbrt(), BigInteger(denominator).Cbrt());
		}

		result.PrecisionMode = this->PrecisionMode;
		result.FixedPrecisionCount = this->FixedPrecisionCount;
		return result;
	}

	// Logarithm base e
	BigFraction BigFraction::Log( const BigInteger& value ) const
	{
		if ( value.IsZero() )
		{
			throw std::invalid_argument( "Logarithm of non-positive value is undefined." );
		}

		BigFraction result = LogCF( value );

		result.PrecisionMode = this->PrecisionMode;
		result.FixedPrecisionCount = this->FixedPrecisionCount;
		return result;
	}

	BigFraction BigFraction::Log() const
	{
		if ( IsZero() || IsNegative() )
		{
			throw std::invalid_argument( "Logarithm of non-positive value is undefined." );
		}
		BigFraction numeratorLog = LogCF( this->numerator );
		BigFraction denominatorLog = LogCF( this->denominator );
		BigFraction result = numeratorLog - denominatorLog;

		result.PrecisionMode = this->PrecisionMode;
		result.FixedPrecisionCount = this->FixedPrecisionCount;
		return result;
	}

	// Logarithm base 10
	BigFraction BigFraction::Log10( const BigInteger& fraction ) const
	{
		static const BigFraction LOG_10_E = Log( TEN );
		BigFraction result = LogCF( fraction ) / LOG_10_E;

		result.PrecisionMode = this->PrecisionMode;
		result.FixedPrecisionCount = this->FixedPrecisionCount;
		return result;
	}

	BigFraction BigFraction::Log10() const
	{
		static const BigFraction LOG_10_E = Log( TEN );
		BigFraction result = Log() / LOG_10_E;

		result.PrecisionMode = this->PrecisionMode;
		result.FixedPrecisionCount = this->FixedPrecisionCount;
		return result;
	}

	// Exponential function using Taylor series
	BigFraction BigFraction::Exp( const BigInteger& value, DecimalPrecisionMode precision_mode, uint64_t fixed_precision_count )
	{
		BigFraction result(ONE, ONE);
		BigFraction term(ONE, ONE);
		BigInteger i = ONE;

		switch ( precision_mode )
		{
			case DecimalPrecisionMode::Fixed:
			{
				for(uint64_t round = fixed_precision_count + 1; round > 0; round--)
				{
					term = term * BigFraction(value, i);
					result += term;
					i += ONE;
				}

				break;
			}
			
			case DecimalPrecisionMode::Full:
			{
				while (term > BigFractionFullPrecision)
				{
					term = term * BigFraction(value, i);
					result += term;
					i += ONE;
				}

				break;
			}
			default:
				break;
		}

		return result;
	}

	BigFraction BigFraction::Exp( const BigFraction& value )
	{
		BigFraction result(ONE, ONE);
		BigFraction term(ONE, ONE);
		BigInteger i = ONE;
		
		switch ( value.PrecisionMode )
		{
			case DecimalPrecisionMode::Fixed:
			{
				for(uint64_t round = value.FixedPrecisionCount + 1; round > 0; round--)
				{
					term = term * (value / i);
					result += term;
					i += ONE;
				}

				break;
			}
			
			case DecimalPrecisionMode::Full:
			{
				BigFraction precision = value.GetFullPrecision();
				while (term > precision)
				{
					term = term * (value / i);
					result += term;
					i += ONE;
				}

				break;
			}
			default:
				break;
		}

		result.PrecisionMode = value.PrecisionMode;
		result.FixedPrecisionCount = value.FixedPrecisionCount;

		return result;
	}

	BigFraction BigFraction::Power( const BigInteger& exponent ) const
	{
		BigFraction result;

		if ( exponent.IsZero() )
		{
			result = BigFraction( ONE, ONE );
		}
		else if ( exponent > 0 )
		{
			BigInteger numeratorPower = numerator;
			numeratorPower.BigPower( exponent );
			BigInteger denominatorPower = denominator;
			denominatorPower.BigPower( exponent );
			result = BigFraction( numeratorPower, denominatorPower );
		}
		else
		{
			result = Reciprocal().Power( exponent );
		}

		result.PrecisionMode = this->PrecisionMode;
		result.FixedPrecisionCount = this->FixedPrecisionCount;
		return result;
	}

	BigFraction BigFraction::Power( const BigFraction& exponent ) const
	{
		BigFraction result;

		if ( IsNaN() )
		{
			result = BigFraction( 0, ONE );
		}
		else if ( IsZero() )
		{
			result = BigFraction( ONE, ONE );
		}
		else if ( exponent.IsZero() )
		{
			result = BigFraction( ONE, ONE );;
		}
		else if ( numerator.IsZero() )
		{
			result = BigFraction( 0, ONE );
		}
		else
		{
			// If the numerator is greater than or equal to the denominator, the fraction is greater than or equal to 1.
			if ( this->numerator >= this->denominator || (this->numerator > ONE && this->denominator == ONE) )
			{
				BigFraction newNumerator = Log( numerator ) * exponent;
				BigFraction newDenominator = Log( denominator ) * exponent;

				BigFraction resultNumerator = Exp( newNumerator );
				BigFraction resultDenominator = Exp( newDenominator );

				result = BigFraction( resultNumerator.GetNumerator(), resultDenominator.GetNumerator() );
			}
			else
			{
				long double value = static_cast<long double>(*this);
				long double value2 = static_cast<long double>(exponent);
				long double floating_result = std::pow(value, value2);
				BigFraction result = this->FromFloatingNumber<long double>(floating_result);
			}
		}

		result.PrecisionMode = this->PrecisionMode;
		result.FixedPrecisionCount = this->FixedPrecisionCount;
		return result;
	}

	/**
	 * @brief Computes the nth root of a BigFraction.
	 * 
	 * This function calculates the nth root of a BigFraction. The approach depends
	 * on whether the fraction's numerator is greater than or equal to the denominator 
	 * (i.e., the fraction is greater than or equal to 1) or less than the denominator 
	 * (i.e., the fraction is less than 1).
	 * 
	 * @param n The root to compute (e.g., n = 2 for square root).
	 * @return BigFraction The nth root of the fraction.
	 * @throws std::invalid_argument if the fraction is zero or negative, as the nth root 
	 *         is undefined for these cases.
	 */
	BigFraction BigFraction::NthRoot(const BigInteger& n) const
	{
		// Check for invalid input: zero or negative values are not supported.
		if (IsZero() || IsNegative())
		{
			throw std::invalid_argument("Nth root of zero or negative value is undefined.");
		}

		// Nth Root for BigFraction
		// If the numerator is greater than or equal to the denominator, the fraction is greater than or equal to 1.
		if ( this->numerator >= this->denominator || (this->numerator > ONE && this->denominator == ONE) )
		{
			// Use the ShiftingKthRoot algorithm to compute the nth root 
			// for both the numerator and denominator.
			TwilightDream::BigInteger::ShiftingKthRoot kthRootCalculator(n.ToUnsignedInt());
			TwilightDream::BigInteger::BigInteger kthRootedNumerator = kthRootCalculator(this->numerator);
			TwilightDream::BigInteger::BigInteger kthRootedDenominator = kthRootCalculator(this->denominator);

			// Construct the result as a BigFraction.
			BigFraction result(kthRootedNumerator, kthRootedDenominator);
			result.PrecisionMode = this->PrecisionMode;
			result.FixedPrecisionCount = this->FixedPrecisionCount;
			return result;
		}
		else
		{
			// For small fractions (numerator < denominator), using integer root algorithms can be problematic due to potential precision issues. 
			// Therefore, we convert the fraction to a floating-point number, compute the root, and then convert back to BigFraction.
			
			// Convert the fraction to a long double for root calculation.
			long double value = static_cast<long double>(*this);
			// Compute the nth root using standard library pow function.
			long double floating_result = std::pow(value, long double(1.0) / n.ToUnsignedInt());
			// Convert the floating-point result back to BigFraction.
			BigFraction result = this->FromFloatingNumber<long double>(floating_result);
			result.PrecisionMode = this->PrecisionMode;
			result.FixedPrecisionCount = this->FixedPrecisionCount;
			return result;
		}

		// This point should not be reached, as the function should return from one of the branches above.
		return BigFraction(ONE, ONE);
	}


	// Nth Root for BigFraction with fraction parameter
	BigFraction BigFraction::NthRoot( const BigFraction& fraction, const BigInteger& nth ) const
	{
		if ( fraction.IsZero() || fraction.IsNegative() )
		{
			throw std::invalid_argument( "Nth root of zero or negative value is undefined." );
		}

		auto result = fraction.NthRoot(nth);
		// Return the refined result
		result.PrecisionMode = this->PrecisionMode;
		result.FixedPrecisionCount = this->FixedPrecisionCount;
		return result;
	}

	BigFraction BigFraction::Sine( const BigFraction& x ) const
	{
		BigFraction result( 0, ONE );
		BigFraction term = x;  // x^1 / 1!
		BigInteger	i = ONE;   // Factorial index starts at 1

		bool isNegative = false;

		switch ( this->PrecisionMode )
		{
			case DecimalPrecisionMode::Fixed:
			{
				for ( uint64_t round = this->FixedPrecisionCount + 1; round > 0; --round )
				{
					if ( isNegative )
					{
						result -= term;
					}
					else
					{
						result += term;
					}

					// Prepare for the next term: x^n / n!
					i += 2;
					term *= x * x;
					term = BigFraction( term.GetNumerator(), i * ( i - 1 ) );

					isNegative = !isNegative;  // Alternate the sign for each term
				}
				break;
			}

			case DecimalPrecisionMode::Full:
			{
				BigFraction precision = GetFullPrecision();
				while ( term.Abs() > precision )
				{
					if ( isNegative )
					{
						result -= term;
					}
					else
					{
						result += term;
					}

					// Prepare for the next term: x^n / n!
					i += 2;
					term *= x * x;
					term = BigFraction( term.GetNumerator(), i * ( i - 1 ) );

					isNegative = !isNegative;  // Alternate the sign for each term
				}
				break;
			}

			default:
				throw std::invalid_argument( "Unknown PrecisionMode" );
		}

		result.PrecisionMode = this->PrecisionMode;
		result.FixedPrecisionCount = this->FixedPrecisionCount;

		return result;
	}

	BigFraction BigFraction::Cosine( const BigFraction& x ) const
	{
		BigFraction result( ONE, ONE );
		BigFraction term( ONE, ONE );
		BigInteger	i = 2;	// Factorial index starts at 2

		bool isNegative = true;

		switch ( this->PrecisionMode )
		{
			case DecimalPrecisionMode::Fixed:
			{
				for ( uint64_t round = this->FixedPrecisionCount + 1; round > 0; --round )
				{
					term *= x * x;
					term = BigFraction( term.GetNumerator(), i * ( i - 1 ) );

					if ( isNegative )
					{
						result -= term;
					}
					else
					{
						result += term;
					}
					isNegative = !isNegative;  // Alternate the sign for each term

					i += 2;
				}
				break;
			}

			case DecimalPrecisionMode::Full:
			{
				BigFraction precision = GetFullPrecision();
				while ( term.Abs() > precision )
				{
					term *= x * x;
					term = BigFraction( term.GetNumerator(), i * ( i - 1 ) );

					if ( isNegative )
					{
						result -= term;
					}
					else
					{
						result += term;
					}
					isNegative = !isNegative;  // Alternate the sign for each term

					i += 2;
				}
				break;
			}

			default:
				throw std::invalid_argument( "Unknown PrecisionMode" );
		}

		result.PrecisionMode = this->PrecisionMode;
		result.FixedPrecisionCount = this->FixedPrecisionCount;

		return result;
	}

	BigFraction BigFraction::Tangent(const BigFraction& x) const
	{
		BigFraction sinResult = Sine(x);
		BigFraction cosResult = Cosine(x);

		if (cosResult.IsZero())
		{
			throw std::overflow_error("Tangent undefined for this input (cos(x) = 0).");
		}

		auto result = sinResult / cosResult;
		result.PrecisionMode = this->PrecisionMode;
		result.FixedPrecisionCount = this->FixedPrecisionCount;

		return result;
	}

	BigFraction BigFraction::Arctangent( const BigFraction& x ) const
	{
		BigFraction result = x;
		BigFraction term = x;
		BigInteger	i = 3;

		switch ( this->PrecisionMode )
		{
			case DecimalPrecisionMode::Fixed:
			{
				for ( uint64_t round = this->FixedPrecisionCount; round > 0; --round )
				{
					term *= x * x;
					result -= term / i;
					i += 2;
				}
				break;
			}

			case DecimalPrecisionMode::Full:
			{
				BigFraction precision = GetFullPrecision();
				while ( term.Abs() > precision )
				{
					term *= x * x;
					result -= term / i;
					i += 2;
				}
				break;
			}

			default:
				throw std::invalid_argument( "Unknown PrecisionMode" );
		}

		result.PrecisionMode = this->PrecisionMode;
		result.FixedPrecisionCount = this->FixedPrecisionCount;

		return result;
	}

	BigFraction BigFraction::Arcsine( const BigFraction& x ) const
	{
		if ( x > BigFraction( ONE, ONE ) || x < BigFraction( -1, ONE ) )
		{
			throw std::invalid_argument( "asin(x) is undefined for |x| > 1" );
		}

		BigFraction result = x;
		BigFraction term = x;
		BigInteger	i = 1;

		switch ( this->PrecisionMode )
		{
			case DecimalPrecisionMode::Fixed:
			{
				for ( uint64_t round = 1; round <= this->FixedPrecisionCount; ++round )
				{
					term *= x * x * BigFraction( i, i + 1 );
					i += 2;
					term = BigFraction( term.GetNumerator(), i * ( i - 1 ) );
					result += term;
				}
				break;
			}

			case DecimalPrecisionMode::Full:
			{
				BigFraction precision = GetFullPrecision();
				while ( term.Abs() > precision )
				{
					term *= x * x * BigFraction( i, i + 1 );
					i += 2;
					term = BigFraction( term.GetNumerator(), i * ( i - 1 ) );
					result += term;
				}
				break;
			}

			default:
				throw std::invalid_argument( "Unknown PrecisionMode" );
		}

		result.PrecisionMode = this->PrecisionMode;
		result.FixedPrecisionCount = this->FixedPrecisionCount;

		return result;
	}

	BigFraction BigFraction::Arccosine(const BigFraction& x) const
	{
		auto result = (GenerateSrinivasaRamanujanPI() / BigFraction(TWO, ONE)) - Arcsine(x);

		result.PrecisionMode = this->PrecisionMode;
		result.FixedPrecisionCount = this->FixedPrecisionCount;
		
		return result;
	}

	BigFraction::BigInteger BigFraction::Floor() const
	{
		if ( IsNaN() )
		{
			return BigInteger( 0 );
		}
		else if ( IsZero() )
		{
			return BigInteger( 0 );
		}
		else
		{
			BigInteger result = ( numerator * sign ) / denominator;
			if ( sign == -1 )
			{
				result -= 1;
			}
			return result;
		}
	}

	BigFraction::BigInteger BigFraction::Ceil() const
	{
		if ( IsNaN() )
		{
			return BigInteger( 0 );
		}
		else if ( IsZero() )
		{
			return BigInteger( 0 );
		}
		else
		{
			BigInteger result = ( numerator * sign ) / denominator;
			if ( sign == 1 )
			{
				result += 1;
			}
			return result;
		}
	}

	BigFraction::BigInteger BigFraction::Round() const
	{
		if ( IsNaN() )
		{
			return BigInteger( 0 );
		}
		else if ( IsZero() )
		{
			return BigInteger( 0 );
		}
		else
		{
			BigInteger numeratorTwice = ( numerator << 1 );
			BigInteger denominatorTwice = ( denominator << 1 );
			BigInteger result = ( numeratorTwice + denominator ) / ( denominatorTwice * sign );
			return result;
		}
	}

	BigFraction::operator BigInteger() const
	{
		if ( IsNaN() || IsZero() )
		{
			return 0;
		}

		if ( denominator == 1 )
		{
			return numerator;
		}

		return numerator / denominator;
	}

	long double BigFraction::BigIntegerToLongDouble(const BigInteger& big_integer) const
	{
		long double result = 0.0;
		long double factor = 1.0;

		BigInteger temp = big_integer;
		while (!temp.IsZero()) 
		{
			// Extract the last digit
			BigInteger digit = temp % TEN;

			// Convert this digit to long double and add to the result
			result += static_cast<long double>(digit.ToUnsignedInt()) * factor;

			// Prepare to process the next digit
			temp /= TEN;
			factor *= 10.0;
		}

		return result;
	}

	BigFraction::operator long double() const
	{
		return sign * ( BigIntegerToLongDouble(this->numerator) / BigIntegerToLongDouble(this->denominator) );
	}

	BigFraction::operator double() const
	{
		return static_cast<double>( static_cast<long double>( *this ) );
	}

	BigFraction::operator float() const
	{
		return static_cast<float>( static_cast<long double>( *this ) );
	}

	bool BigFraction::IsPerfectPower( const BigInteger& N )
	{
		if ( N <= 1 )
			return false;

		for ( BigInteger b = TWO; b < N.Log2() + 1; ++b )
		{
			BigFraction exponent( 1, b, 1 );
			BigFraction base( N, 1, 1 );
			BigFraction a = base.Power( exponent );

			if ( !a.IsNaN() )
			{
				return true;
			}
		}

		return false;
	}

	std::istream& operator>>( std::istream& is, BigFraction& fraction )
	{
		std::string input;
		is >> input;
		fraction.ComputeAndFromDecimalString( input );
		return is;
	}

	std::ostream& operator<<( std::ostream& os, const BigFraction& fraction )
	{
		std::string output;
		output = fraction.ComputeAndToDecimalString();
		os << output;
		return os;
	}

	void BigFraction::ReduceSimplify()
	{
		// Handle case where numerator is zero
		if (this->numerator.IsZero())
		{
			this->sign = 1;  // Zero has no sign
			this->denominator = BigInteger(1);  // Denominator becomes 1 (0/1)
			return;
		}

		// Handle case where denominator is zero
		if (this->denominator.IsZero())
		{
			if (this->numerator == ONE)
			{
				//inline BigFraction BigFractionFullPrecision(1, 0);
				return;
			}
			throw std::invalid_argument("Denominator cannot be zero.");
		}

		// Calculate the greatest common divisor (GCD) to simplify the fraction
		BigInteger gcd = BigInteger::GCD(this->numerator, this->denominator);
    
		// Simplify the numerator and denominator using the GCD
		this->numerator /= gcd;
		this->denominator /= gcd;
	}

	/*
		Logarithm of a BigInteger using continued fractions
		The function uses the series expansion of the natural logarithm:
		\[
		\ln(1+x) = 2 \left( \frac{x}{1} - \frac{x^3}{3} + \frac{x^5}{5} - \cdots \right)
		\]
		The function uses \( x = \frac{value + 1}{1} \), ensuring \( |x| > 0 \).
	*/
	BigFraction BigFraction::LogCF( const BigInteger& value ) const
	{
		if (value <= 0)
		{
			throw std::invalid_argument("Logarithm is undefined for non-positive values.");
		}
		
		// Initialize x as (value - 1) / (value + 1)
		BigFraction x(value - ONE, value + ONE);
		
		// Variables to store the current and previous terms
		BigFraction result(0, ONE);
		BigFraction newTerm = x;
		BigFraction oldTerm(0, ONE);

		// Initialize variables for the loop
		BigFraction xPower = x;
		BigInteger n = 1;
		bool round_flag = true;
		
		BigFraction term(0, ONE);

		switch (this->PrecisionMode)
		{
			case DecimalPrecisionMode::Fixed:
			{
				// Main loop to calculate
				for (uint64_t round = this->FixedPrecisionCount + 1; round > 0; round--)
				{
					oldTerm = newTerm;
					n += 2;
					xPower *= x.Power(2);  // xPower *= x^2
					term = BigFraction(xPower, n);

					if (round_flag)
					{
						newTerm -= term;
					}
					else
					{
						newTerm += term;
					}

					round_flag = !round_flag;  // Toggle round_flag for the next term
				}

				break;
			}
			
			case DecimalPrecisionMode::Full:
			{
				// Set precision for the calculation
				BigFraction precision = GetFullPrecision(); // 10^-n precision (1 / n)

				BigFraction diff(x, ONE);

				// Main loop to calculate
				while (diff > precision)
				{
					oldTerm = newTerm;
					n += 2;
					xPower *= x.Power(2);  // xPower *= x^2
					term = BigFraction(xPower, n);

					if (round_flag)
					{
						newTerm -= term;
					}
					else
					{
						newTerm += term;
					}

					round_flag = !round_flag;  // Toggle round_flag for the next term

					diff = (newTerm - oldTerm).Abs();
				}

				break;
			}
			
			default:
				throw std::invalid_argument( "Unknown PrecisionMode" );
		}

		// Scale the result by 2
		result = newTerm * BigFraction(2, 1);

		result.PrecisionMode = this->PrecisionMode;
		result.FixedPrecisionCount = this->FixedPrecisionCount;

		return result;
	}

	/*
		@brief Generate a BigFraction representing the value of π (Pi) using the Srinivasa Ramanujan series.

		@details
		This function precomputes the value of π based on a high-precision calculation using the Ramanujan series. The series is implemented in Python code, which is used to generate the large numerator and denominator that represent π as a fraction. The numerator and denominator are then used to create a BigFraction object, which is simplified before returning. This approach ensures that the value of π is available as a precomputed fraction, avoiding the need for real-time calculation.

		Here's the Python code used for the precomputation:

		```python
		import decimal
		import math
		import sys

		# Increase the maximum allowed string length limit 1000000
		sys.set_int_max_str_digits(1000000)

		def ramanujan_pi(iterations):
			decimal.getcontext().prec = 4000  # Setting High Precision (4000) Floating Number
			total_sum = decimal.Decimal(0)
			factor = decimal.Decimal(2 * math.sqrt(2)) / decimal.Decimal(9801)
    
			for k in range(iterations):
				numerator = decimal.Decimal(math.factorial(4 * k)) * decimal.Decimal(1103 + 26390 * k)
				denominator = (decimal.Decimal(math.factorial(k)) ** 4) * (decimal.Decimal(396) ** (4 * k))
				total_sum += numerator / denominator
    
			pi_inv = factor * total_sum
			pi_value = 1 / pi_inv
    
			return pi_value

		# Calculate 4000 iterations
		iterations = 1000
		pi_value = ramanujan_pi(iterations)

		# Extract the numerator and denominator
		numerator = pi_value.as_integer_ratio()[0]
		denominator = pi_value.as_integer_ratio()[1]

		print(f"numerator: {numerator}")
		print("\n")
		print(f"denominator: {denominator}")
		```

		@returns A BigFraction object representing the value of π, which is precomputed and stored as a fraction.
	*/
	BigFraction BigFraction::GenerateSrinivasaRamanujanPI()
	{
		BigInteger numerator
		(
			"3141592653589793023709380780717486763469833607826720486919923336894437464007835588165042081964422742587908889560773603151602852252160721362604174811891856615759117890464640476568650178481249677573589790684201235888650294353755998754685334341161547519996739479088998029573036973027413119828620943663789725484131700408126609909622220518805590414926449571875506022908156575184488006635859659116582117046254701384522374966358155301272693553353003579630093035050347617539808771920501672319758479656692221558575312251138967002902152245874255769149955239237448955176839985431847020706153743905502399179299646534937041554002789010027961962404319344317286479510693244584354856919020779537843922119464912373073596572509774070303324515084520691906992844853606319349101668803386721721598040581777774166381273204366786202324958675881342695628474223358811039066984925125526569697961907489718469177713381893701508628667491043927677900003085691238567118912952987422332684754462530254905466964591042619559310506524647225038356053503517038314302522459085771358839391487226506563236692702240829730941639779948228827744194136658917165343179231213257234946968543185983639441766802188606222112905767965979508122349492249684602497365713780782464471588731382170390365980830890278333528243982312919799161908263914936503227238359902976561935855017607416271780661568874032546521026334890218008723205438436573308335005064367492649130650170928018997041850225811532040195032354048598916679277743469374504223598001834255935428585386570575108999193655530015169673080359135096125370607287492684886298690718642253474552538575943448760746261986701522809808106907856472453745548018843666778246347149199421242796654197581930022582079252315784443840828099414976684367185149593612697279014481924167803593224961844806925090015708913389587171094241095715819246403640345543172057480914110903244374400952316134238229364711669612939589272446948764986560662622803554056559382613847375895090272960402640799417413219275841599247692537975220369891493369865122789481022660204430836546740679836647812978454578829074438622280447133004537289961496898199570487997874049473858235886677924248439420091039036310792704707026395624936939772365850056540274223425967005733999430968934506361300522095662329222033832903147074107889445193213728511664854793343852929172217730802019849568383847258473682687297207396191198944675921694663985208768130199129729807581992912097840078948767054754182676043961089110174429550908179188926388525337029878502685792221563622587586281868524789668842682653706593305182025066942224885298542175285628276455007702606472975814022406582691276168482734482757548369089586081975059391171129937574473723464153781854038582808728063370266549256305374990709200838108880168637423218889972603741587261798443770711093880615118922854878352301092878360016581833083063796618609368560241339159973874058629627138496986957631689871966254367254692956016379279448874683514281697004585555838564094418793338671368601173283704279073264965621690114064163235967729047883678416938917360002662535574528001196028215302453291901414717614600696089088354944244872431473782146646549712341516997796893078368238399904995381372061217567030618141243103498561744102841802985521647089277639437300673578556292373521599048508968947191382785903832432096668474344277696879602308885004950094917338132091884244421556370764488766399194839163626522234902751110737457211049911174568301928147415503014754416774043957433847893694785409916513148219571185372567798354645049837538654975228391463026784522026597955601234576727383064854311734139929145523402730741308615148210849533160273611408093565534511040013619706768849711043039694441521436006230947022332470959566068470691867582191855083701316378489719328591195689069354659130487014161868594172673350216663277301691150694975210755044585491762698066726264933918570727604759370840023937155558976535205544633307037314341960814595142969678651164854900062052682375651304373110563708342214734229470223681470230598001528372515949904254107" , 
			10
		);

		BigInteger denominator
		(
			"1000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000",
			10
		);

		//π
		auto pi = BigFraction(numerator, denominator);
		
		pi.ReduceSimplify();

		return pi;
	}

	BigFraction BigFraction::GenerateNilakanthaArrayPI( uint64_t iteration )
	{
		using TwilightDream::BigInteger::BigInteger;
		using TwilightDream::BigFraction::BigFraction;

		// Initial value π = 3
		BigFraction pi( 3, 1 );

		// Initial numerator = 4
		BigInteger numerator( 4 );

		// Initial denominator terms 2, 3, 4
		BigInteger denominator1( 2 );
		BigInteger denominator2( 3 );
		BigInteger denominator3( 4 );

		for ( uint64_t i = 0; i < iteration; ++i )
		{
			// constructed fraction term：4 / (denominator1 * denominator2 * denominator3)
			BigFraction term( numerator, denominator1 * denominator2 * denominator3 );

			// alternate addition and subtraction
			if ( i % 2 == 0 )
			{
				pi += term;
			}
			else
			{
				pi -= term;
			}

			// Updating denominator terms: (2n+2), (2n+3), (2n+4)
			denominator1 += 2;
			denominator2 += 2;
			denominator3 += 2;
		}

		// Set a very high precision denominator for the final value
		TwilightDream::BigInteger::BigInteger FullPrecision(
			"1000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000", 
			10
		);
		pi.SetFullPrecision(FullPrecision);

		// You can now output the result if needed (though you asked not to)
		//std::cout << "Pi approximated using " << iteration << " iterations: " << pi << std::endl;
	
		//std::cout << "Pi squared: " << pi.Power(2) << std::endl;
		//std::cout << "Square root of Pi: " << pi.Sqrt() << std::endl;
		//std::cout << "Cube root of Pi: " << pi.Cbrt() << std::endl;
		//std::cout << "Logarithm of Pi: " << pi.Log() << std::endl;
		//std::cout << "Log10 of Pi: " << pi.Log10() << std::endl;

		return pi;
	}

}  // namespace TwilightDream::BigFraction