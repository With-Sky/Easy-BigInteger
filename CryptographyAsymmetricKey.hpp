#pragma once
#include "BigInteger.hpp"
#include <unordered_map>
#include <optional>
#include <future>

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
		const BigInteger ONE = BigInteger( 1 );
		const BigInteger TWO = BigInteger( 2 );
		const BigInteger THREE = BigInteger( 3 );

		struct RSA_AlgorithmNumbers
		{
			using BigInteger = TwilightDream::BigInteger::BigInteger;

			BigInteger EncryptExponent = 0;
			BigInteger DecryptExponent = 0;
			BigInteger AlgorithmModulus = 0;
		};

		size_t bit_count = 0;
		BigInteger MIN = 0;
		BigInteger MAX = 0;
		BigInteger RANGE_INTEGER = 0;

		struct FindPrimeState
		{
			bool IsPrime = false;
			BigInteger Number = BigInteger( 0 );
		};

		void GeneratePrimesInParallelFunction(size_t bit_count, std::vector<FindPrimeState>& primes, size_t primes_index);
		std::optional<FindPrimeState> GeneratePrimesInParallelFunctions(size_t bit_count);
		BigInteger GeneratePrimeNumber( size_t bit_count );

	public:
		/**
		* @brief Generates RSA key pair (public and private keys).
		*
		* @param bit_count The bit length of the RSA modulus.
		* @param is_pkcs A flag indicating whether to use PKCS#1 standard for key generation.
		* @return RSA_AlgorithmNumbers structure containing the generated key pair.
		*/
		RSA_AlgorithmNumbers GenerateKeys( size_t bit_count, bool is_pkcs );

		/**
		* @brief Encrypts a message using RSA.
		*
		* @param PlainMessage The message to be encrypted.
		* @param EncryptExponent The RSA public key exponent for encryption.
		* @param AlgorithmModulus The RSA modulus used for encryption.
		*/
		void Encryption( BigInteger& PlainMessage, const BigInteger& EncryptExponent, const BigInteger& AlgorithmModulus );

		/**
		* @brief Decrypts a message using RSA.
		*
		* @param CipherMessage The encrypted message to be decrypted.
		* @param DecryptExponent The RSA private key exponent for decryption.
		* @param AlgorithmModulus The RSA modulus used for decryption.
		*/
		void Decryption( BigInteger& CipherMessage, const BigInteger& DecryptExponent, const BigInteger& AlgorithmModulus );

		/**
		* @brief Performs self-sanity check for RSA algorithm.
		*
		* Generates multiple key pairs and tests encryption/decryption consistency.
		*
		* @param bit_count The bit length of the RSA modulus.
		* @param test_count The number of iterations for primality testing during key generation.
		* @param rounds The number of iterations for the self-sanity check.
		* @return True if the self-sanity check passes, false otherwise.
		*/
		static bool SelfSanityCheck( size_t bit_count, size_t rounds )
		{
			RSA		   rsa;
			BigInteger OriginalMessage = 0;
			BigInteger EncryptedMessage = 0, DecryptedMessage = 0;
			for ( size_t current_round = 0; current_round < rounds; current_round++ )
			{
				// 生成密钥对
				auto keys = rsa.GenerateKeys( bit_count, true );

				// 加解密测试
				OriginalMessage = BigInteger::RandomGenerateNBit( bit_count / 2 );
				
				EncryptedMessage = OriginalMessage;
				rsa.Encryption( EncryptedMessage, keys.EncryptExponent, keys.AlgorithmModulus );
				
				DecryptedMessage = EncryptedMessage;
				rsa.Decryption( DecryptedMessage, keys.DecryptExponent, keys.AlgorithmModulus );

				// 验证结果
				if ( OriginalMessage != DecryptedMessage )
				{
					std::cerr << "Decryption failed!" << std::endl;
					return false;
				}
			}

			std::cout << "Self sanity check passed." << std::endl;
			return true;
		}
	};
}  // namespace TwilightDream::CryptographyAsymmetric
