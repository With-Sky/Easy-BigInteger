#pragma once

#include "BigInteger.hpp"

namespace TwilightDream::BigFraction
{
	class BigFraction
	{
	public:
		using BigInteger = TwilightDream::BigInteger::BigInteger;

		// Define radix constant
		const BigInteger TWO = 2;
		const BigInteger FIVE = 5;
		const BigInteger TEN = 10;
		const BigInteger RADIX = TEN;

		bool SimplifyReduced = true;

		BigFraction() : numerator( 0 ), denominator( 1 ), sign( 1 ) {}

		BigFraction( const BigInteger& numerator, const BigInteger& denominator )
			: numerator( numerator ), denominator( denominator ), sign( 1 )
		{
			// If the numerator or denominator is negative, adjust the sign and values
			if ( numerator.IsNegative() )
			{
				this->numerator = numerator.Abs();
				sign *= -1;
			}
			if ( denominator.IsNegative() )
			{
				this->denominator = denominator.Abs();
				sign *= -1;
			}

			if(SimplifyReduced)
			{
				// Reduce the fraction to its simplest form
				ReduceSimplify();
			}
		}

		BigFraction( const BigInteger& numerator, const BigInteger& denominator, int sign )
			: numerator( numerator ), denominator( denominator ), sign( sign >= 0 ? 1 : -1 )
		{
			// If the numerator or denominator is negative, adjust the sign and values
			if ( numerator.IsNegative() )
			{
				this->numerator = numerator.Abs();
			}
			if ( denominator.IsNegative() )
			{
				this->denominator = denominator.Abs();
			}

			if(SimplifyReduced)
			{
				// Reduce the fraction to its simplest form
				ReduceSimplify();
			}
		}

		BigFraction( const BigFraction& other )
			: numerator( other.numerator ), denominator( other.denominator ), sign( other.sign )
		{
		
		}

		BigFraction( const std::string& complexString, size_t precision )
		{
			ComputeAndFromDecimalString(complexString, precision);
		}

		BigFraction& operator=( const BigFraction& other )
		{
			if ( this != &other )
			{
				numerator = other.numerator;
				denominator = other.denominator;
				sign = other.sign;
			}
			return *this;
		}

		BigFraction& operator=( BigFraction& other )
		{
			if ( this != &other )
			{
				numerator = other.numerator;
				denominator = other.denominator;
				sign = other.sign;
			}
			return *this;
		}

		//This is not a number
		bool IsNaN() const
		{
			// 0÷0
			if(numerator.IsZero() && denominator.IsZero())
			{
				return true;
			}

			// N÷0
			if(numerator.IsZero() && denominator.IsZero())
			{
				return true;
			}

			return false;
		}

		bool IsZero() const
		{
			// 0÷N
			return numerator.IsZero() && !denominator.IsZero();
		}

		bool IsNegative() const
		{
			return sign == -1;
		}

		BigInteger GetNumerator() const
		{
			return numerator;
		}

		void SetNumerator( const BigInteger& number )
		{
			numerator = number;
		}

		BigInteger SetDenominator() const
		{
			return denominator;
		}

		void SetDenominator( const BigInteger& number )
		{
			denominator = number;
		}

		// High-precision fractions are converted to high-precision decimals and output by approximation
		std::string ComputeAndToDecimalString( size_t precision ) const
		{
			if (IsNaN())
			{
				return "NaN";
			}

			if (precision == 0)
			{
				// Return integer part only
				return numerator.ToString(10);
			}

			BigInteger absNumerator = numerator.Abs();
			BigInteger absDenominator = denominator.Abs();

			// Reduce the fraction to its simplest form
			BigInteger gcd = BigInteger::GCD(absNumerator, absDenominator);
			absNumerator /= gcd;
			absDenominator /= gcd;

			// Compute the integer part
			BigInteger integerPart = absNumerator / absDenominator;
			std::string integerPartStr = integerPart.ToString(10);

			// Compute the fractional part
			BigInteger fractionalPart = absNumerator % absDenominator;
			std::string fractionalPartStr = "";

			// Extend the fractional part to the desired precision
			for (size_t i = 0; i < precision; ++i)
			{
				// Multiply remainder by RADIX (10) and append next digit to fractional part
				fractionalPart *= RADIX;
				BigInteger nextDigit = fractionalPart / absDenominator;
				fractionalPart %= absDenominator;
				fractionalPartStr += nextDigit.ToString(10);

				// If fractional part becomes zero, no need to continue
				if (fractionalPart == 0)
					break;
			}

			// Combine the integer and fractional parts
			std::string result = integerPartStr + ".";
			if (!fractionalPartStr.empty())
			{
				result += fractionalPartStr;
			}

			// Apply the sign
			if (sign == -1)
			{
				result = "-" + result;
			}

			return result;
		}

		// High-precision fractions are converted from high-precision decimals and input by approximation
		void ComputeAndFromDecimalString( const std::string& complexString, size_t precision )
		{
			// Parsing Complex String
			bool isNegative = false;
			bool hasFractionalPart = false;
			size_t integerPartEnd = 0;
			size_t fractionalPartStart = 0;

			// Determine if the number is negative
			if (!complexString.empty() && complexString[0] == '-')
			{
				isNegative = true;
				integerPartEnd = 1;
			}

			// Find the end of integer part
			while (integerPartEnd < complexString.size() && std::isdigit(complexString[integerPartEnd]))
			{
				integerPartEnd++;
			}

			// Check if there is a fractional part
			if (integerPartEnd < complexString.size() && complexString[integerPartEnd] == '.')
			{
				hasFractionalPart = true;
				fractionalPartStart = integerPartEnd + 1;
			}

			// Parse integer part
			BigInteger integerPart(complexString.substr(isNegative, integerPartEnd - isNegative));

			// Parse fractional part if exists
			BigInteger fractionalPart = 0;
			BigInteger fractionalMultiplier = 1;
			if (hasFractionalPart)
			{
				for (size_t i = fractionalPartStart; i < complexString.size(); ++i)
				{
					if (!std::isdigit(complexString[i]))
					{
						// Invalid character found in fractional part
						// You may want to handle this error case appropriately
						return;
					}
					fractionalPart = fractionalPart * RADIX + (complexString[i] - '0');
					fractionalMultiplier *= RADIX;
				}
			}

			// Combine integer and fractional parts into the numerator
			numerator = integerPart * fractionalMultiplier + fractionalPart;
			denominator = fractionalMultiplier;

			// Apply the sign
			if (isNegative)
			{
				numerator *= -1;
			}

			// Simplify the fraction to reduce it to its simplest form
			if(SimplifyReduced)
			{
				ReduceSimplify();
			}
		}

		// (A÷B) + (C÷D)
		BigFraction& operator+=( const BigFraction& other )
		{
			if(this->IsNaN() || other.IsNaN())
			{
				this->sign = 1;
				this->numerator = 0;
				this->denominator = 1;
				return *this;
			}

			if (IsZero())
			{
				if (other.IsZero())
				{
					sign = 1;
					numerator = 0;
					denominator = 1;
					return *this;
				}
				else if (other.IsNegative())
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
			else if (other.IsZero())
			{
				return *this;
			}

			// Compute numerator and denominator without changing signs
			BigInteger resultNumerator = sign * numerator * other.denominator + other.sign * other.numerator * denominator;
			BigInteger resultDenominator = denominator * other.denominator;

			// Compute sign based on operands' signs
			sign = (resultNumerator < 0) ^ (resultDenominator < 0) ? -1 : 1;

			// Ensure numerator and denominator are positive
			numerator = resultNumerator.Abs();
			denominator = resultDenominator.Abs();

			// Simplify the fraction to reduce it to its simplest form
			if(SimplifyReduced)
			{
				ReduceSimplify();
			}

			return *this;
		}

		// (A÷B) + (C÷1)
		BigFraction& operator+=( const BigInteger& other )
		{
			if(this->IsNaN())
			{
				this->sign = 1;
				this->numerator = 0;
				this->denominator = 1;
				return *this;
			}

			if (this->IsZero())
			{
				if (other.IsZero())
				{
					// Both are zero, result is zero
					return *this;
				}
				else if (other.IsNegative())
				{
					// Fraction is zero, integer is negative, result is negative fraction
					this->sign = -1;
					this->numerator = other.Abs();
					this->denominator = 1;
					return *this;
				}
				else
				{
					// Fraction is zero, integer is positive, result is positive fraction
					this->sign = 1;
					this->numerator = other;
					this->denominator = 1;
					return *this;
				}
			}
			else if (other == 0)
			{
				// Integer is zero, fraction remains unchanged
				return *this;
			}

			// At this point, we know neither the fraction nor the integer is zero

			// If the fraction and integer have the same sign, we can directly add the integers
			if (this->sign == 1)
			{
				this->numerator += other;
			}
			else
			{
				// If the fraction is negative, we treat it as positive for calculation
				this->numerator -= other;
			}

			if(SimplifyReduced)
			{
				ReduceSimplify();
			}

			return *this;
		}

		// (A÷B) - (C÷D)
		BigFraction& operator-=( const BigFraction& other )
		{
			if(this->IsNaN() || other.IsNaN())
			{
				this->sign = 1;
				this->numerator = 0;
				this->denominator = 1;
				return *this;
			}

			if (this->IsZero())
			{
				if (other.IsZero())
				{
					this->sign = 1;
					this->numerator = 0;
					this->denominator = 1;
					return *this;
				}
				else if (other.IsNegative())
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
			else if (other.IsZero())
			{
				if (this->IsNegative())
				{
					this->sign = -1;
					return *this;
				}
				else
				{
					this->sign = 1;
					return *this;
				}
			}

			// Compute numerator and denominator without changing signs
			BigInteger resultNumerator = sign * numerator * other.denominator - other.sign * other.numerator * denominator;
			BigInteger resultDenominator = denominator * other.denominator;

			// Compute sign based on operands' signs
			sign = (resultNumerator.IsNegative()) ^ (resultDenominator.IsNegative()) ? -1 : 1;

			// Ensure numerator and denominator are positive
			numerator = resultNumerator.Abs();
			denominator = resultDenominator.Abs();

			// Simplify the fraction to reduce it to its simplest form
			if(SimplifyReduced)
			{
				ReduceSimplify();
			}

			return *this;
		}

		// (A÷B) - (C÷1)
		BigFraction& operator-=( const BigInteger& other )
		{
			if (this->IsNaN())
			{
				this->sign = 1; // Assuming 1 represents positive infinity or NaN
				this->numerator = 0;
				this->denominator = 1;
				return *this;
			}

			if (this->IsZero())
			{
				if (other.IsZero())
				{
					// Both fraction and integer are zero, result is zero
					this->sign = 1;
					this->numerator = 0;
					this->denominator = 1;
					return *this;
				}
				else if (other.IsNegative())
				{
					// Fraction is zero, integer is negative, result is a positive integer
					this->sign = 1;
					this->numerator = other.Abs(); // Assuming BigInteger has a unary minus operator
					this->denominator = 1;
					return *this;
				}
				else
				{
					// Fraction is zero, integer is positive, result is a negative integer
					this->sign = -1;
					this->numerator = other;
					this->denominator = 1;
					return *this;
				}
			}
			else if (other.IsZero())
			{
				// Integer is zero, fraction remains unchanged
				return *this;
			}

			// At this point, we know neither the fraction nor the integer is zero

			// If the fraction is positive and the integer is negative, or vice versa, we subtract the absolute value of the integer from the numerator
			if ((this->sign > 0 && other.IsNegative()) || (this->sign < 0 && !other.IsNegative()))
			{
				this->numerator += other.Abs(); // Assuming BigInteger has an Abs method
			}
			else
			{
				// If both are negative or both are positive, we add the integer to the numerator
				this->numerator += other;
			}

			// Simplify the fraction to reduce it to its simplest form
			if(SimplifyReduced)
			{
				ReduceSimplify();
			}

			return *this;
		}

		// (A÷B) × (C÷D)
		BigFraction& operator*=( const BigFraction& other )
		{
			if(this->IsNaN() || other.IsNaN())
			{
				this->sign = 1;
				this->numerator = 0;
				this->denominator = 1;
				return *this;
			}

			if(other.IsZero())
			{
				this->sign = 1;
				this->numerator = 0;
				this->denominator = 1;
				return *this;
			}
			if((other.numerator > 0 &&  other.denominator > 0) && (other.numerator == other.denominator))
			{
				return *this;
			}

			this->numerator *= (other.numerator).Abs();
			this->denominator *= other.denominator.Abs();

			if((this->sign == -1 && other.IsNegative()) || (this->sign == 1 && !other.IsNegative()))
			{
				this->sign = 1;
			}
			else if(this->sign == 1 && other.IsNegative() || (this->sign == -1 && !other.IsNegative()))
			{
				this->sign = -1;
			}

			// Simplify the fraction to reduce it to its simplest form
			if(SimplifyReduced)
			{
				ReduceSimplify();
			}

			return *this;
		}

		// (A÷B) × (C÷1)
		BigFraction& operator*=( const BigInteger& other )
		{
			if(other.IsZero())
			{
				this->sign = 1;
				this->numerator = 0;
				this->denominator = 1;
				return *this;
			}
			if(other == 1)
			{
				return *this;
			}

			this->numerator *= other.Abs();

			if((this->sign == -1 && other.IsNegative()) || (this->sign == 1 && !other.IsNegative()))
			{
				this->sign = 1;
			}
			else if(this->sign == 1 && other.IsNegative() || (this->sign == -1 && !other.IsNegative()))
			{
				this->sign = -1;
			}

			// Simplify the fraction to reduce it to its simplest form
			if(SimplifyReduced)
			{
				ReduceSimplify();
			}

			return *this;
		}

		// (A÷B) ÷ (C÷D)
		BigFraction& operator/=( const BigFraction& other )
		{
			if(this->IsNaN() || other.IsNaN())
			{
				this->sign = 1;
				this->numerator = 0;
				this->denominator = 1;
				return *this;
			}

			// Check for division by zero in the other fraction
			if (other.denominator.IsZero())
			{
				throw std::runtime_error("Division by zero error.");
			}

			// Calculate the new numerator and denominator
			this->numerator = (this->numerator * other.denominator).Abs();
			this->denominator = (this->denominator * other.numerator).Abs();

			// Compute the sign of the result
			this->sign = (this->sign * other.sign) * ((this->sign * other.denominator) > 0 ? 1 : -1);

			// Simplify the fraction to reduce it to its simplest form
			if(SimplifyReduced)
			{
				ReduceSimplify();
			}

			return *this;
		}

		// (A÷B) ÷ (C÷1)
		BigFraction& operator/=( const BigInteger& other )
		{
			if(this->IsNaN())
			{
				this->sign = 1;
				this->numerator = 0;
				this->denominator = 1;
				return *this;
			}

			if (other.IsZero())
			{
				throw std::runtime_error("Division by zero error.");
			}

			// Calculate the new numerator
			this->numerator = (this->numerator).Abs();
			this->denominator *= other.Abs();

			// Compute the sign of the result
			this->sign *= (other > 0) ? 1 : -1;

			// Simplify the fraction to reduce it to its simplest form
			if(SimplifyReduced)
			{
				ReduceSimplify();
			}

			return *this;
		}

		// Compute the square root of the fraction
		BigFraction Sqrt()
		{
			if ( IsNaN() )
			{
				return BigFraction( 0, 1 );
			}
			else if ( IsZero() )
			{
				return BigFraction( 0, 1 );
			}
			else if ( IsNegative() )
			{
				// Square root of a negative number is undefined
				return BigFraction( 0, 1 );
			}
			else
			{
				// Compute square root using the Newton-Raphson method
				BigInteger sqrtNumerator = numerator.Sqrt();
				BigInteger sqrtDenominator = denominator.Sqrt();

				// If the denominator is not a perfect square, take the square root of numerator and denominator separately
				if ( sqrtDenominator * sqrtDenominator != denominator )
				{
					sqrtNumerator *= denominator.Sqrt();
					sqrtDenominator *= denominator.Sqrt();
				}

				// Return the square root of the fraction
				return BigFraction( sqrtNumerator, sqrtDenominator );
			}
		}

		// Compute the reciprocal of the fraction
		BigFraction Reciprocal() const
		{
			if ( IsNaN() )
			{
				return BigFraction( 0, 1 );
			}
			else if ( IsZero() )
			{
				// Reciprocal of zero is undefined
				return BigFraction( 0, 1 );
			}
			else
			{
				return BigFraction( denominator, numerator, sign );
			}
		}

		// Compute the power of the fraction to an integer exponent
		BigFraction Power( const BigInteger& exponent ) const
		{
			if ( exponent.IsZero() )
			{
				return BigFraction( 1, 1 );  // Any number to the power of 0 is 1
			}
			else if ( exponent > 0 )
			{
				// Compute the power as (numerator^exponent) / (denominator^exponent)
				BigInteger numeratorPower = numerator;
				numeratorPower.BigPower( exponent );
				BigInteger denominatorPower = denominator;
				denominatorPower.BigPower( exponent );
				return BigFraction( numeratorPower, denominatorPower );
			}
			else  // exponent < 0
			{
				// Compute the power of the reciprocal and take the reciprocal
				return Reciprocal().Power( exponent.Abs() );
			}
		}

		// Compute the power of the fraction to a fractional exponent
		BigFraction Power( const BigFraction& exponent ) const
		{
			// Check for special cases
			if ( IsNaN() )
			{
				return BigFraction( 0, 1 );
			}
			else if ( IsZero() )
			{
				// Any non-zero number to the power of 0 is 1
				return BigFraction( 1, 1 );
			}
			else if ( exponent.IsZero() )
			{
				// Any number to the power of 0 is 1
				return BigFraction( 1, 1 );
			}
			else if ( numerator.IsZero() )
			{
				// Zero to any non-zero power is zero
				return BigFraction( 0, 1 );
			}
			else if ( !exponent.IsNaN() )
			{
				// If the exponent is an integer, use the integer power function
				return this->Power(exponent.Round());
			}
			else
			{
				//PLEASE FIXME !!! (IS WRONG FORMULA)
				BigInteger numeratorPower = numerator;
				numeratorPower.BigPower( exponent.Round() );
				BigInteger denominatorPower = denominator;
				denominatorPower.BigPower( exponent.Round() );
				return BigFraction( numeratorPower, denominatorPower );
			}
		}

		// Compute the floor of the fraction
		BigInteger Floor() const
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

		// Compute the ceiling of the fraction
		BigInteger Ceil() const
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

		// Compute the rounded value of the fraction
		BigInteger Round() const
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

		static bool IsPerfectPower(const BigInteger& N)
		{
			// Special cases: 0 and 1 are not considered perfect powers
			if (N <= 1)
				return false;

			// Start with b = 2 and iterate until the logarithm of N to the base 2
			for (BigInteger b = 2; b < N.Log2() + 1; ++b)
			{
				// Compute a = N^(1/b)
				BigFraction exponent(1, b, 1);
				BigFraction base(N, 1, 1);
				BigFraction a = base.Power(exponent);

				// Check if a is an integer
				if (!a.IsNaN())
				{
					return true; // N is a perfect power
				}
			}

			return false; // N is not a perfect power
		}

	private:
		//Numerator of a fraction
		BigInteger numerator;
		//Denominator of a fraction
		BigInteger denominator;
		int32_t	sign = 1;  // 1 for positive, -1 for negative

		void ReduceSimplify()
		{
			if(this->numerator.IsZero())
			{
				this->sign = 1;
				this->denominator = 1;
				return;
			}
			if(this->denominator.IsZero())
			{
				std::invalid_argument("Denominator cannot be zero.");
			}
			this->sign = numerator.IsNegative() ? -1 : 1;
			this->numerator = this->numerator.Abs();
			BigInteger gcd = BigInteger::GCD(numerator, denominator);
			this->numerator /= gcd;
			this->denominator /= gcd;
		}
	};
}  // namespace TwilightDream::BigFriction