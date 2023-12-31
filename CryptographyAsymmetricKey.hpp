#pragma once
#include "BigInteger.hpp"

/*
	Twilight-Dream
	Cryptography Asymmetric Key Algorithms
*/


namespace TwilightDream::CryptographyAsymmetric
{
	// Optimal asymmetric encryption padding
	struct OAEP
	{

	};

	// Rivest, Adi Shamir, and Leonard Adleman
	class RSA
	{
	private:
		using BigInteger = TwilightDream::BigInteger::BigInteger;
		const BigInteger ONE = BigInteger(1);
		const BigInteger THREE = BigInteger(3);

		struct RSA_AlgorithmNumbers
		{
			using BigInteger = TwilightDream::BigInteger::BigInteger;

			BigInteger EncryptExponent = 0;
			BigInteger DecryptExponent = 0;
			BigInteger AlgorithmModulus = 0;
		};

		BigInteger GeneratePrimeNumber(size_t bit_count, size_t test_count);

	public:

		RSA_AlgorithmNumbers GenerateKeys(size_t bit_count, size_t test_count, bool is_pkcs);

		void Encryption(BigInteger& PlainMessage, const BigInteger& EncryptExponent, const BigInteger& AlgorithmModulus);
		void Decryption(BigInteger& CipherMessage, const BigInteger& DecryptExponent, const BigInteger& AlgorithmModulus);
	};
}
