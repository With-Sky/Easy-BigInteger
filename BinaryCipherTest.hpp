#pragma once

#include "BigInteger.hpp"
#include <array>

#ifndef USE_BIG_INTEGER_ARITHMATIC_VERSION
//#define USE_BIG_INTEGER_ARITHMATIC_VERSION
#endif

struct BinaryCipher
{
	using BigInteger = TwilightDream::BigInteger::BigInteger;

	#ifdef USE_BIG_INTEGER_ARITHMATIC_VERSION
	const BigInteger Number2Power512 = BigInteger(1) << 512;
	#endif	// USE_BIG_INTEGER_ARITHMATIC_VERSION

	const std::bitset<512> Bitset512SizeZero = std::bitset<512>();
	
	//Initial Vector (512 Bit)
	BigInteger Random0, Random1, Random2, Random3, Random4, Random5, Random6, Random7;

	BinaryCipher()
		: Random0( Bitset512SizeZero ), Random1( Bitset512SizeZero ), Random2( Bitset512SizeZero ), Random3( Bitset512SizeZero ), Random4( Bitset512SizeZero ), Random5( Bitset512SizeZero ), Random6( Bitset512SizeZero ), Random7( Bitset512SizeZero )
	{
		
	}

	void InitialWithKey(std::vector<uint8_t> Keys);
	
	BigInteger AdditionBits512( const  BigInteger& a, const  BigInteger& b );

	void KeyExpansion(const BigInteger& BitsKey, BigInteger& BitsExpansionKey);

	void Test()
	{
		std::cout << "################ BinaryCipher Test ################\n" << std::endl;

		const std::bitset<256> Bitset256SizeZero = std::bitset<256>();
		BigInteger BitsKey = BigInteger::BigInteger( Bitset256SizeZero );
		BitsKey.SetBit( 255 );

		BigInteger BitsExpansionKey = BigInteger::BigInteger( Bitset512SizeZero );
		KeyExpansion( BitsKey, BitsExpansionKey );



		BigInteger BitsPlainData = BigInteger::BigInteger( Bitset512SizeZero );

		std::cout << "MatserKey: " << BitsKey.ToBinaryString(512) << std::endl;
		std::cout << "ExpansionKey: " << BitsExpansionKey.ToBinaryString(512) << std::endl;

		std::cout << "Before PlainData: " << BitsPlainData.ToBinaryString(512) << std::endl;

		BigInteger BitsCipherData = BitsPlainData ^ BitsExpansionKey;
		std::cout << "CipherData: " << BitsCipherData.ToBinaryString(512) << std::endl;

		BitsPlainData = BitsCipherData ^ BitsExpansionKey;
		std::cout << "After PlainData: " << BitsPlainData.ToBinaryString(512) << std::endl;

		std::cout << "\n################ BinaryCipher Test End ################" << std::endl;
	}
};

struct BinaryCipherNaive
{
	void KeyExpansion( const std::bitset<256>& BitsKey, std::bitset<512>& BitsExpansionKey );

	//Initial Vector (512 Bit)
	std::bitset<512> RandomBits0, RandomBits1, RandomBits2, RandomBits3, RandomBits4, RandomBits5, RandomBits6, RandomBits7;

	BinaryCipherNaive()
	{
		
	}

	void InitialWithKey(std::vector<uint8_t> Keys);

	void Test()
	{
		std::cout << "################ BinaryCipherNaive Test ################\n" << std::endl;

		std::bitset<256> BitsKey( "1000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000" );
		std::bitset<512> BitsExpansionKey;
		KeyExpansion( BitsKey, BitsExpansionKey );

		std::bitset<512> BitsPlainData;

		std::cout << "MatserKey: " << BitsKey << std::endl;
		std::cout << "ExpansionKey: " << BitsExpansionKey << std::endl;

		std::cout << "Before PlainData: " << BitsPlainData << std::endl;

		std::bitset<512> BitsCipherData = BitsPlainData ^ BitsExpansionKey;
		std::cout << "CipherData: " << BitsCipherData << std::endl;

		BitsPlainData = BitsCipherData ^ BitsExpansionKey;
		std::cout << "After PlainData: " << BitsPlainData << std::endl;

		std::cout << "\n################ BinaryCipherNaive Test End ################" << std::endl;
	}
};