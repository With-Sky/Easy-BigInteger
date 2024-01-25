//Reference code: https://github.com/NoahBz/Easy-BigInt/blob/main/BigInt.h

#pragma once

#include <iostream>
#include <cmath>
#include <string>
#include <algorithm>
#include <vector>
#include <map>
#include <unordered_map>
#include <thread>
#include <complex>
#include <functional>
#include <type_traits>
#include <cstdint>
#include <intrin.h>
#include <immintrin.h>
#include <random>
#include <bitset>

namespace TwilightDream
{
	struct PrimeNumberTester;
}

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
	inline static const uint16_t EXPONENT_VALUE = 32;
	inline static const uint16_t EXPONENT_FHT_UI16 = 16;
	inline static const uint16_t EXPONENT_VALUE2 = std::numeric_limits<digit_type>::digits / 2;
	inline static const uint16_t EXPONENT = EXPONENT_VALUE < EXPONENT_VALUE2 ? EXPONENT_VALUE : EXPONENT_VALUE2;

	//TODO
	//The limit is now a 64-bit type, try opening up more bits !!!
	using FHT_BITS_TYPE = std::conditional<EXPONENT <= 64, std::uint32_t, std::uint64_t>::type;

	// Set max value of the values[index]
	// 65536, 536870912, 2147483648, 4294967296, 18446744073709551616
	inline static const uint64_t BASE = digit_type(pow(2, EXPONENT));
	
	constexpr int BASE_MULTIPLY_LIMIT = 120; // 0-120 Base
	constexpr int FHT_MULTIPLY_LIMIT = 65536; // 120 - 65536 FHT, 65536 - inf, : Karatsuba

	constexpr int BINARY_SEARCH_DIVISION_LIMIT = 4; //0-4 binary
	constexpr int SHORT_DIVISION_LIMIT = 10;		// 4 - 10 short
	constexpr int DONALD_KNUTH_LONG_DIVISION_LIMIT = 20;// 10 - 20 Knuth, Max of knuth, Min of Newton
	constexpr int MULTI_THREAD_LIMIT = 10000;
	constexpr bool MULTI_THREAD = true;
	constexpr uint64_t STACK_LIMIT = 4096;

	inline std::random_device TrueRandomDevice;
	inline std::default_random_engine RandomEngine( TrueRandomDevice() );
	inline std::uniform_int_distribution<uint64_t> UniformIntDistribution( 0, BASE - 1 );

	namespace HyperIntegerFunctions
	{
		using Float32 = float;
		using Float64 = double;
		using Complex32 = std::complex<Float32>;
		using Complex64 = std::complex<Float64>;

		constexpr Float64 HINT_PI = 3.141592653589793238462643;
		constexpr Float64 HINT_2PI = HINT_PI * 2;
		constexpr size_t  FHT_MAX_LEN = size_t( 1 ) << 18;	// The max length of FHT to avoid the overflow of float64

		template <typename T>
		constexpr T int_floor2( T n )
		{
			constexpr int bits = sizeof( n ) * 8;
			for ( int i = 1; i < bits; i *= 2 )
			{
				n |= ( n >> i );
			}
			return ( n >> 1 ) + 1;
		}

		template <typename T>
		constexpr T int_ceil2( T n )
		{
			constexpr int bits = sizeof( n ) * 8;
			n--;
			for ( int i = 1; i < bits; i *= 2 )
			{
				n |= ( n >> i );
			}
			return n + 1;
		}

		//log2(n) of integer
		template <typename T>
		constexpr int hint_log2( T n )
		{
			constexpr int bits = sizeof( n ) * 8;
			int			  l = -1, r = bits;
			while ( ( l + 1 ) != r )
			{
				int mid = ( l + r ) / 2;
				if ( ( T( 1 ) << mid ) > n )
				{
					r = mid;
				}
				else
				{
					l = mid;
				}
			}
			return l;
		}

		// namespace of fft and fft like algorithm
		namespace Transform
		{
			template <typename T>
			inline void Transform_TwoPoint( T& sum, T& diff )
			{
				T temp0 = sum, temp1 = diff;
				sum = temp0 + temp1;
				diff = temp0 - temp1;
			}
			// Point of unit circle
			template <typename FloatTy>
			inline auto UnitCircleRoot( FloatTy theta )
			{
				return std::polar<FloatTy>( 1.0, theta );
			}

			namespace AlgorithmFHT
			{
				// Look up table of unit roots
				template <typename FloatTy>
				class DynamicTable
				{
				public:
					using Complex = std::complex<FloatTy>;
					using CompVec = std::vector<Complex>;
					
					DynamicTable() {}
					DynamicTable( int log_len, int factor, bool conj = false )
					{
						size_t vec_len = ( 1 << log_len ) / 4;
						table = CompVec( vec_len );
						init( factor, conj );
					}

					void init( int factor, bool conj )
					{
						size_t	len = table.size() * 4;
						FloatTy unity = -HINT_2PI * factor / len;
						if ( conj )
						{
							unity = -unity;
						}
						for ( size_t i = 0; i < table.size(); i++ )
						{
							table[ i ] = UnitCircleRoot( unity * i );
						}
					}

					const Complex& operator[]( size_t n ) const
					{
						return table[ n ];
					}
					Complex& operator[]( size_t n )
					{
						return table[ n ];
					}
					auto get_it( size_t n = 0 )
					{
						return &table[ n ];
					}

				private:
					CompVec table;
				};

				// The Fast Hartley Transform
				// Reference:
				// [1] J¨org Arndt.Matters Computational[M].Heidelberg:Springer Berlin,2011:515-533.
				// https://doi.org/10.1007/978-3-642-14764-7 https://www.jjj.de/fxt/fxtbook.pdf
				// [2] Vlodymyr Myrnyy.A Simple and Efficient FFT Implementation in C++[J/OL].EETimes,2007.
				// https://www.eetimes.com/a-simple-and-efficient-fft-implementation-in-c-part-i/
				// [3] Chinese: RRRR_wys.离散哈特莱变换(DHT)及快速哈特莱变换(FHT)学习.博客园.2018.
				// https://www.cnblogs.com/RRRR-wys/p/10090007.html
				// template class of FHT, using template recursion to speed up code
				template <size_t LEN, typename FloatTy>
				struct FHT
				{
					enum
					{
						fht_len = LEN,
						half_len = LEN / 2,
						quater_len = LEN / 4,
						log_len = hint_log2( fht_len )
					};
					using HalfFHT = FHT<half_len, FloatTy>;
					static DynamicTable<FloatTy> TABLE;
					
					template <typename FloatIt>
					// Decimation in time FHT
					static void DecimationInTime( FloatIt in_out )
					{
						static_assert( std::is_same<typename std::iterator_traits<FloatIt>::value_type, FloatTy>::value, "Must be same as the FHT template float type" );
						HalfFHT::DecimationInTime( in_out );
						HalfFHT::DecimationInTime( in_out + half_len );

						Transform_TwoPoint( in_out[ 0 ], in_out[ half_len ] );
						Transform_TwoPoint( in_out[ quater_len ], in_out[ half_len + quater_len ] );

						auto it0 = in_out + 1, it1 = in_out + half_len - 1;
						auto it2 = in_out + half_len + 1, it3 = in_out + fht_len - 1;
						auto omega_it = TABLE.get_it( 1 );
						for ( ; it0 < it1; ++it0, --it1, ++it2, --it3, omega_it++ )
						{
							auto temp0 = it2[ 0 ], temp1 = it3[ 0 ];
							auto omega = omega_it[ 0 ];
							auto temp2 = temp0 * omega.real() + temp1 * omega.imag();
							auto temp3 = temp0 * omega.imag() - temp1 * omega.real();
							temp0 = it0[ 0 ], temp1 = it1[ 0 ];
							it0[ 0 ] = temp0 + temp2;
							it1[ 0 ] = temp1 + temp3;
							it2[ 0 ] = temp0 - temp2;
							it3[ 0 ] = temp1 - temp3;
						}
					}
					// Decimation in frequency FHT
					template <typename FloatIt>
					static void DecimationInFrequency( FloatIt in_out )
					{
						static_assert( std::is_same<typename std::iterator_traits<FloatIt>::value_type, FloatTy>::value, "Must be same as the FHT template float type" );
						Transform_TwoPoint( in_out[ 0 ], in_out[ half_len ] );
						Transform_TwoPoint( in_out[ quater_len ], in_out[ half_len + quater_len ] );

						auto it0 = in_out + 1, it1 = in_out + half_len - 1;
						auto it2 = in_out + half_len + 1, it3 = in_out + fht_len - 1;
						auto omega_it = TABLE.get_it( 1 );
						for ( ; it0 < it1; ++it0, --it1, ++it2, --it3, omega_it++ )
						{
							auto temp0 = it0[ 0 ], temp1 = it1[ 0 ];
							auto temp2 = it2[ 0 ], temp3 = it3[ 0 ];
							it0[ 0 ] = temp0 + temp2;
							it1[ 0 ] = temp1 + temp3;
							temp0 = temp0 - temp2;
							temp1 = temp1 - temp3;
							auto omega = omega_it[ 0 ];
							it2[ 0 ] = temp0 * omega.real() + temp1 * omega.imag();
							it3[ 0 ] = temp0 * omega.imag() - temp1 * omega.real();
						}

						HalfFHT::DecimationInFrequency( in_out );
						HalfFHT::DecimationInFrequency( in_out + half_len );
					}
				};
				template <size_t LEN, typename FloatTy>
				DynamicTable<FloatTy> FHT<LEN, FloatTy>::TABLE( FHT<LEN, FloatTy>::log_len, 1, true );

				// Template specialization of short FHTs
				template <typename FloatTy>
				struct FHT<0, FloatTy>
				{
					template <typename FloatIt>
					static void DecimationInTime( FloatIt in_out )
					{
					}
					template <typename FloatIt>
					static void DecimationInFrequency( FloatIt in_out )
					{
					}
				};

				template <typename FloatTy>
				struct FHT<1, FloatTy>
				{
					template <typename FloatIt>
					static void DecimationInTime( FloatIt in_out )
					{
					}
					template <typename FloatIt>
					static void DecimationInFrequency( FloatIt in_out )
					{
					}
				};

				template <typename FloatTy>
				struct FHT<2, FloatTy>
				{
					template <typename FloatIt>
					static void DecimationInTime( FloatIt in_out )
					{
						Transform_TwoPoint( in_out[ 0 ], in_out[ 1 ] );
					}
					template <typename FloatIt>
					static void DecimationInFrequency( FloatIt in_out )
					{
						Transform_TwoPoint( in_out[ 0 ], in_out[ 1 ] );
					}
				};

				template <typename FloatTy>
				struct FHT<4, FloatTy>
				{
					template <typename FloatIt>
					static void DecimationInTime( FloatIt in_out )
					{
						auto temp0 = in_out[ 0 ], temp1 = in_out[ 1 ];
						auto temp2 = in_out[ 2 ], temp3 = in_out[ 3 ];
						Transform_TwoPoint( temp0, temp1 );
						Transform_TwoPoint( temp2, temp3 );
						in_out[ 0 ] = temp0 + temp2;
						in_out[ 1 ] = temp1 + temp3;
						in_out[ 2 ] = temp0 - temp2;
						in_out[ 3 ] = temp1 - temp3;
					}
					template <typename FloatIt>
					static void DecimationInFrequency( FloatIt in_out )
					{
						auto temp0 = in_out[ 0 ], temp1 = in_out[ 1 ];
						auto temp2 = in_out[ 2 ], temp3 = in_out[ 3 ];
						Transform_TwoPoint( temp0, temp2 );
						Transform_TwoPoint( temp1, temp3 );
						in_out[ 0 ] = temp0 + temp1;
						in_out[ 1 ] = temp0 - temp1;
						in_out[ 2 ] = temp2 + temp3;
						in_out[ 3 ] = temp2 - temp3;
					}
				};

				template <typename FloatTy>
				struct FHT<8, FloatTy>
				{
					template <typename FloatIt>
					static void DecimationInTime( FloatIt in_out )
					{
						auto temp0 = in_out[ 0 ], temp1 = in_out[ 1 ];
						auto temp2 = in_out[ 2 ], temp3 = in_out[ 3 ];
						auto temp4 = in_out[ 4 ], temp5 = in_out[ 5 ];
						auto temp6 = in_out[ 6 ], temp7 = in_out[ 7 ];
						Transform_TwoPoint( temp0, temp1 );
						Transform_TwoPoint( temp2, temp3 );
						Transform_TwoPoint( temp4, temp5 );
						Transform_TwoPoint( temp6, temp7 );
						Transform_TwoPoint( temp0, temp2 );
						Transform_TwoPoint( temp1, temp3 );
						Transform_TwoPoint( temp4, temp6 );
						Transform_TwoPoint( temp5, temp7 );

						in_out[ 0 ] = temp0 + temp4;
						in_out[ 2 ] = temp2 + temp6;
						in_out[ 4 ] = temp0 - temp4;
						in_out[ 6 ] = temp2 - temp6;
						static constexpr decltype( temp0 ) SQRT_2_2 = 0.70710678118654757;
						temp0 = ( temp5 + temp7 ) * SQRT_2_2;
						temp2 = ( temp5 - temp7 ) * SQRT_2_2;
						in_out[ 1 ] = temp1 + temp0;
						in_out[ 3 ] = temp3 + temp2;
						in_out[ 5 ] = temp1 - temp0;
						in_out[ 7 ] = temp3 - temp2;
					}
					template <typename FloatIt>
					static void DecimationInFrequency( FloatIt in_out )
					{
						auto temp0 = in_out[ 0 ], temp1 = in_out[ 1 ];
						auto temp2 = in_out[ 2 ], temp3 = in_out[ 3 ];
						auto temp4 = in_out[ 4 ], temp5 = in_out[ 5 ];
						auto temp6 = in_out[ 6 ], temp7 = in_out[ 7 ];
						Transform_TwoPoint( temp0, temp4 );
						Transform_TwoPoint( temp1, temp5 );
						Transform_TwoPoint( temp2, temp6 );
						Transform_TwoPoint( temp3, temp7 );
						Transform_TwoPoint( temp0, temp1 );
						Transform_TwoPoint( temp2, temp3 );
						in_out[ 0 ] = temp0 + temp2;
						in_out[ 1 ] = temp1 + temp3;
						in_out[ 2 ] = temp0 - temp2;
						in_out[ 3 ] = temp1 - temp3;
						static constexpr decltype( temp0 ) SQRT_2_2 = 0.70710678118654757;
						temp0 = ( temp5 + temp7 ) * SQRT_2_2;
						temp2 = ( temp5 - temp7 ) * SQRT_2_2;
						Transform_TwoPoint( temp4, temp6 );
						Transform_TwoPoint( temp0, temp2 );
						in_out[ 4 ] = temp4 + temp0;
						in_out[ 5 ] = temp4 - temp0;
						in_out[ 6 ] = temp6 + temp2;
						in_out[ 7 ] = temp6 - temp2;
					}
				};

				// Function to help choose the correct template FHT dit function
				template <size_t LEN = 1>
				inline void fht_dit_template_alt( Float64* input, size_t fht_len )
				{
					if ( fht_len < LEN )
					{
						fht_dit_template_alt<LEN / 2>( input, fht_len );
						return;
					}
					FHT<LEN, Float64>::DecimationInTime( input );
				}
				template <>
				inline void fht_dit_template_alt<0>( Float64* input, size_t fht_len )
				{
				}

				// Function to help choose the correct template FHT dif function
				template <size_t LEN = 1>
				inline void fht_dif_template_alt( Float64* input, size_t fht_len )
				{
					if ( fht_len < LEN )
					{
						fht_dif_template_alt<LEN / 2>( input, fht_len );
						return;
					}
					FHT<LEN, Float64>::DecimationInFrequency( input );
				}
				template <>
				inline void fht_dif_template_alt<0>( Float64* input, size_t fht_len )
				{
				}

				inline auto fht_dit = fht_dit_template_alt<FHT_MAX_LEN>;
				inline auto fht_dif = fht_dif_template_alt<FHT_MAX_LEN>;

				// Use FHT to accelerate convolution of float64 array
				inline void fht_convolution( Float64 fht_ary1[], Float64 fht_ary2[], Float64 out[], size_t fht_len )
				{
					if ( fht_len == 0 )
					{
						return;
					}
					if ( fht_len == 1 )
					{
						out[ 0 ] = fht_ary1[ 0 ] * fht_ary2[ 0 ];
						return;
					}
					
					fht_len = int_floor2( fht_len );
					if ( fht_len > FHT_MAX_LEN )
					{
						std::cout << "FHT len cannot be larger than FHT_MAX_LEN" << std::endl;
						throw("FHT len cannot be larger than FHT_MAX_LEN" );
					}
					fht_dif( fht_ary1, fht_len );
					// When the two float arrays are actually same, execute DIF only once
					if ( fht_ary1 != fht_ary2 )
					{
						fht_dif( fht_ary2, fht_len );
					}
					
					const double inv = 0.5 / fht_len;
					out[ 0 ] = fht_ary1[ 0 ] * fht_ary2[ 0 ] / fht_len;
					out[ 1 ] = fht_ary1[ 1 ] * fht_ary2[ 1 ] / fht_len;
					if ( fht_len == 2 )
					{
						return;
					}
					// Theory of convolution using FHT
					auto temp0 = fht_ary1[ 2 ], temp1 = fht_ary1[ 3 ];
					auto temp2 = fht_ary2[ 2 ], temp3 = fht_ary2[ 3 ];
					Transform_TwoPoint( temp0, temp1 );
					out[ 2 ] = ( temp2 * temp0 + temp3 * temp1 ) * inv;
					out[ 3 ] = ( temp3 * temp0 - temp2 * temp1 ) * inv;
					for ( size_t i = 4; i < fht_len; i *= 2 )
					{
						auto it0 = fht_ary1 + i, it1 = it0 + i - 1;
						auto it2 = fht_ary2 + i, it3 = it2 + i - 1;
						auto it4 = out + i, it5 = it4 + i - 1;
						for ( ; it0 < it1; it0 += 2, it1 -= 2, it2 += 2, it3 -= 2, it4 += 2, it5 -= 2 )
						{
							temp0 = *it0, temp1 = *it1, temp2 = *it2, temp3 = *it3;
							Transform_TwoPoint( temp0, temp1 );
							*it4 = ( temp2 * temp0 + temp3 * temp1 ) * inv;
							*it5 = ( temp3 * temp0 - temp2 * temp1 ) * inv;
							temp0 = *( it1 - 1 ), temp1 = *( it0 + 1 ), temp2 = *( it3 - 1 ), temp3 = *( it2 + 1 );
							Transform_TwoPoint( temp0, temp1 );
							*( it5 - 1 ) = ( temp2 * temp0 + temp3 * temp1 ) * inv;
							*( it4 + 1 ) = ( temp3 * temp0 - temp2 * temp1 ) * inv;
						}
					}
					fht_dit( out, fht_len );
				}
			} // namespace hint_fht
		} // namespace hint_transform

		// FHT multiplication
		template<typename UintTy>
		void FHTMul( UintTy* out, const UintTy* in1, size_t in_len1, const UintTy* in2, size_t in_len2 )
		{
			// Use 16bit binary as an element
			auto   out_16 = reinterpret_cast<uint16_t*>( out );
			auto   in1_16 = reinterpret_cast<const uint16_t*>( in1 );
			auto   in2_16 = reinterpret_cast<const uint16_t*>( in2 );
			size_t in1_len16 = in_len1 * sizeof( UintTy ) / sizeof( uint16_t );
			size_t in2_len16 = in_len2 * sizeof( UintTy ) / sizeof( uint16_t );
			size_t out_len16 = in1_len16 + in2_len16, conv_len = out_len16 - 1, fht_len = int_ceil2( conv_len );

			std::vector<Float64> buffer1( fht_len ), buffer2( fht_len );  // FHT bufffer
			std::copy( in1_16, in1_16 + in1_len16, buffer1.data() );
			std::copy( in2_16, in2_16 + in2_len16, buffer2.data() );

			Transform::AlgorithmFHT::fht_convolution( buffer1.data(), buffer2.data(), buffer1.data(), fht_len );  // FHT convolution

			uint64_t carry = 0;
			for ( size_t i = 0; i < conv_len; i++ )
			{
				carry += uint64_t( buffer1[ i ] + 0.5 );
				out_16[ i ] = carry & UINT16_MAX;
				carry = carry >> 16;
			}
			out_16[ conv_len ] = carry & UINT16_MAX;
		}

		// FHT Square
		template <typename UintTy>
		void FHTSquare( UintTy* out, const UintTy* in, size_t in_len )
		{
			// Use 16bit binary as an element
			auto   out_16 = reinterpret_cast<uint16_t*>( out );
			auto   in_16 = reinterpret_cast<const uint16_t*>( in );
			size_t in_len16 = in_len * sizeof( UintTy ) / sizeof( uint16_t );
			size_t out_len16 = in_len16 * 2, conv_len = out_len16 - 1, fht_len = int_ceil2( conv_len );

			std::vector<Float64> buffer( fht_len );	 // FHT bufffer
			std::copy( in_16, in_16 + in_len16, buffer.data() );

			Transform::AlgorithmFHT::fht_convolution( buffer.data(), buffer.data(), buffer.data(), fht_len );	// 卷积

			uint64_t carry = 0;
			for ( size_t i = 0; i < conv_len; i++ )
			{
				carry += uint64_t( buffer[ i ] + 0.5 );
				out_16[ i ] = carry & UINT16_MAX;
				carry = carry >> 16;
			}
			out_16[ conv_len ] = carry & UINT16_MAX;
		}
	}

	class BigInteger
	{
	private:
		friend struct TwilightDream::PrimeNumberTester;

		void Clean();

		void BitLeftShiftCore( const uint32_t& shift );

		void BitRightShiftCore( const uint32_t& shift );

		BigInteger RightShiftBlock( size_t d ) const;

		BigInteger LeftShiftBlock( size_t d ) const;

		/**
		* Computes the reciprocal of a BigInteger raised to a power n.
		*
		* This function calculates the reciprocal of a BigInteger raised to the power of n, which is
		* equivalent to computing (BASE^n) / (*this), where BASE is the base of the BigInteger
		* representation. The function uses a divide-and-conquer approach to efficiently
		* compute the result, especially for large values of n.
		*
		* @param n The exponent to which the BigInteger should be raised.
		* @return The reciprocal of (BASE^n) / (*this).
		*/
		BigInteger DivisionInvertReciprocal(size_t n) const;

		std::vector<digit_type> values;

	public:
		int8_t sign = 1;

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

		BigInteger( int64_t value );

		explicit BigInteger( const std::string& number_string );

		BigInteger( const std::string& number_string, uint32_t base );

		BigInteger( const BigInteger& num );

		BigInteger( BigInteger&& num ) noexcept;

		BigInteger& operator+=( const BigInteger& other );

		/**
		 * @brief Adds another BigInteger to the current instance.
		 *
		 * This function performs addition of two BigIntegers and updates the current
		 * instance with the result.
		 *
		 * @param other The BigInteger to be added.
		 * @return Reference to the modified current instance after addition.
		 */
		BigInteger& Add( const BigInteger& other );

		BigInteger& operator-=( const BigInteger& other );

		BigInteger& Difference( const BigInteger& other );

		/**
		 * @brief Subtracts another BigInteger from the current instance.
		 *
		 * This function performs subtraction of two BigIntegers and updates the current
		 * instance with the result.
		 *
		 * @param other The BigInteger to be subtracted.
		 * @return Reference to the modified current instance after subtraction.
		 */
		BigInteger& Subtract( const BigInteger& other );

		BigInteger operator++();

		BigInteger operator--();

		BigInteger operator++( int );

		BigInteger operator--( int );

		BigInteger operator-() const;

		BigInteger operator+() const;

		/**
		 * @brief Performs base multiplication of two BigIntegers.
		 *
		 * This function multiplies the current instance by another BigInteger using
		 * base multiplication and updates the current instance with the result.
		 *
		 * @param other The BigInteger to be multiplied.
		 * @return Reference to the modified current instance after multiplication.
		 */
		BigInteger& BaseMultiplication( const BigInteger& other );

		/**
		 * @brief Performs FHT multiplication of two BigIntegers.
		 *
		 * This function multiplies the current instance by another BigInteger using
		 * FHT algorithm and updates the current instance with the result.
		 *
		 * @param other The BigInteger to be multiplied.
		 * @return Reference to the modified current instance after multiplication.
		 */
		BigInteger& FHTMultiplication( const BigInteger& other );

		/**
		 * @brief Squares the current BigInteger by using FHT.
		 *
		 * This function squares the current instance by performing a FHT
		 * algorithm and updates the current instance with the result.
		 *
		 * @return Reference to the modified current instance after squaring.
		 */
		BigInteger& FHTSquare();

		/**
		 * @brief Squares the current BigInteger.
		 *
		 * This function squares the current instance by performing a specialized
		 * multiplication and updates the current instance with the result.
		 *
		 * @return Reference to the modified current instance after squaring.
		 */
		BigInteger& Square();

		BigInteger& MultiplyBase( size_t times );

		void SplitAt( const BigInteger& num, const size_t n, BigInteger& high, BigInteger& low );

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
		BigInteger& KaratsubaMultiplication( const BigInteger& other );

		// result = result.Power(exponent)
		BigInteger& Power( const size_t exponent );

		// result = result.BigPower(exponent)
		BigInteger& BigPower( const BigInteger& exponent );

		BigInteger& Sqrt();

		// a = a^{exponent} mod modulo
		BigInteger& PowerWithModulo(const BigInteger& exponent, const BigInteger& modulo);

		// this = this * other mod modulo
		BigInteger& MultiplyWithModulo(const BigInteger& other, const BigInteger modulo);

		//Use by MontgomeryReduce
		BigInteger& MontDivR( size_t rsize );

		//Use by MontgomeryReduce
		BigInteger& MontModR( size_t rsize );

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
		BigInteger& MontgomeryReduce( size_t rsize, const BigInteger& m, const BigInteger& mprime );

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
		BigInteger& MontgomeryPower( const BigInteger& e, const BigInteger& m, const BigInteger& mprime, const BigInteger r, const size_t rsize );

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
		BigInteger& MontgomeryPower( const BigInteger& e, const BigInteger& m );

		BigInteger& Divide( const BigInteger& other, BigInteger* r = nullptr );

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
		std::pair<BigInteger, BigInteger> BinarySearchDivision(const BigInteger& other);

		BigInteger& ShortDivision( const BigInteger& other, BigInteger* r = nullptr );

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
		BigInteger& DonaldKnuthLongDivision( const BigInteger& other, BigInteger* r = nullptr );

		/**
		* Performs division using Newton's iteration method for BigIntegers.
		*
		* This function implements a division algorithm based on Newton's iteration method, which is an
		* efficient way to calculate the division of two large integers. The method works by
		* approximating the reciprocal of the divisor and iteratively refining the quotient until
		* the remainder is less than the divisor.
		*
		* @param divisor The BigInteger representing the divisor.
		* @return A pair containing the quotient and the remainder after division.
		*/
		std::pair<BigInteger, BigInteger> NewtonIterationDivision(const BigInteger& divisor) const;

		BigInteger& operator*=( const BigInteger& other );

		BigInteger& operator/=( const BigInteger& other );

		BigInteger& operator%=( const BigInteger& other );

		bool IsEven() const;

		bool IsPowerOfTwo() const;

		bool IsZero() const;

		bool IsNegative() const;

		size_t Size() const;

		void CountLeadingZeros(size_t& result) const;

		static uint16_t LeadingZeros( digit_type x );

		size_t BitSize() const;

		bool GetBit( size_t bit_position ) const;

		void SetBit( size_t bit_position );

		void SetBit( bool value, size_t bit_position );

		// Add 'count' number of bits with the specified value to the most significant end of the integer.
		void BigInteger::ExtendedLeadingBits(size_t count, bool value);

		// Remove 'count' number of bits from the most significant end of the integer.
		void BigInteger::SqueezeLeadingBits(size_t count);

		BigInteger& FromInt(uint64_t value);

		int64_t ToInt() const;

		uint64_t ToUnsignedInt() const;

		BigInteger& FromUnsignedInt(uint64_t value);

		BigInteger operator&( const BigInteger& other ) const;

		BigInteger& operator&=( const BigInteger& other );

		BigInteger operator|( const BigInteger& other ) const;

		BigInteger& operator|=( const BigInteger& other );

		BigInteger operator~() const;

		BigInteger operator^( const BigInteger& other ) const;

		BigInteger& operator^=( const BigInteger& other );

		BigInteger& operator<<=( const uint32_t shift );

		BigInteger& operator>>=( const uint32_t shift );

		friend BigInteger operator<<( const BigInteger& lhs, const uint32_t shift )
		{
			BigInteger result( lhs );
			result <<= shift;
			return result;
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
		static BigInteger BitRotateLeft(const BigInteger& bits, uint32_t shift, const uint32_t reference_bit_capacity);

		/**
		 * @brief Performs a right rotation on the binary representation of the BigInteger.
		 *
		 * @param shift The number of positions to rotate the bits to the right.
		 * @param reference_bit_capacity The desired bit capacity for the resulting binary string.
		 * @return copy to the modified BigInteger.
		 */
		static BigInteger BitRotateRight(const BigInteger& bits, uint32_t shift, const uint32_t reference_bit_capacity);

		digit_type& operator[]( const size_t index );

		BigInteger& operator=( const BigInteger& other );

		explicit operator bool()
		{
			return this->IsZero() ? false : true;
		}

		bool operator&&( const BigInteger& other ) const;

		bool operator||( const BigInteger& other ) const;

		bool operator!() const;

		// Absolute value
		BigInteger Abs() const;

		static BigInteger GCD( BigInteger a, BigInteger b );

		static void EGCD( const BigInteger& a, const BigInteger& b, BigInteger* gcd, BigInteger* co1, BigInteger* co2 );

		static BigInteger LCM( const BigInteger& a, const BigInteger& b );

		static BigInteger ModuloInverse( const BigInteger& a, const BigInteger& b );

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
		static BigInteger PollardRho( const BigInteger& n );

		static BigInteger RandomGenerateNBit( size_t n );

		/**
		* Calculates the base-2 logarithm (log2) of a BigInteger using binary search.
		*
		* This function iterates through the bits of the BigInteger from the most significant bit to the least
		* significant bit, checking if the bit is set (1). The position of the first set bit is the log2 of the number.
		*
		* @return The value of log2(n) as an integer, or -1 if n is 0 or 1.
		*/
		BigInteger Log2() const;

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
		static BigInteger Modulo(const BigInteger& number, const BigInteger& modulus);

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
		static BigInteger ModuloArithmetic(ArithmeticMode mode, const BigInteger& a , const BigInteger& b, const BigInteger& modulus);

		std::string ToString() const;

		/**
		* Convert the BigInteger to a binary string representation.
		*
		* @param reference_bit_capacity The desired bit capacity of the resulting binary string.
		* @param have_leading_zeros     Flag indicating whether to include leading zeros in the binary string.
		* @return                       Binary string representation of the BigInteger.
		*/
		std::string ToBinaryString( const uint32_t reference_bit_capacity, bool have_leading_zeros = true ) const;

		/**
		* Convert the BigInteger to a string representation with the specified base.
		*
		* @param base_value The desired base for the string representation.
		* @warning Binary base is unsupported in ToString function.
		* @return String representation of the BigInteger in the specified base.
		*/
		std::string ToString( uint32_t base_value ) const;

		void FromString( const std::string& number_string );

		void FromString( const std::string& number_string, uint32_t base_value );

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

		void Print( bool base10 ) const;

		size_t SipHash(const BigInteger& Integer, std::vector<uint8_t>* keys = nullptr) const;

		friend class HashFunction;

		friend BigInteger operator+( const BigInteger& lhs, const BigInteger& rhs )
		{
			BigInteger result( lhs );
			result += rhs;
			return result;
		}

		friend BigInteger operator-( const BigInteger& lhs, const BigInteger& rhs )
		{
			BigInteger result( lhs );
			result -= rhs;
			return result;
		}

		friend BigInteger operator*( const BigInteger& lhs, const BigInteger& rhs )
		{
			BigInteger result( lhs );
			result *= rhs;
			return result;
		}

		friend BigInteger operator/( const BigInteger& lhs, const BigInteger& rhs )
		{
			BigInteger result( lhs );
			result /= rhs;
			return result;
		}

		friend BigInteger operator%( const BigInteger& lhs, const BigInteger& rhs )
		{
			BigInteger result( lhs );
			result %= rhs;
			return result;
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
	};
	
	class HashFunction
	{
	public:
		// 使用 lambda 表达式来调用自定义的 hash 方法
		size_t operator()(const BigInteger& key) const
		{
			return key.SipHash(key);
		}
	};

}  // namespace TwilightDream::BigInteger
