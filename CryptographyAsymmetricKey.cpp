#include "CryptographyAsymmetricKey.hpp"

namespace TwilightDream::CryptographyAsymmetric
{
	RSA::BigInteger RSA::GeneratePrimeNumber( size_t bit_count, size_t test_count )
	{
		BigInteger PrimeNumber = 0;
		do
		{
			PrimeNumber = BigInteger::RandomGenerateNBit(bit_count);
		}
		while ( !BigInteger::MillerRabinWithMontgomery(PrimeNumber, test_count) );

		return PrimeNumber;
	}

	RSA::RSA_AlgorithmNumbers RSA::GenerateKeys( size_t bit_count, size_t test_count, bool is_pkcs )
	{
		if(bit_count == 0 || bit_count < 1024)
		{
			throw std::invalid_argument("Invalid bit_count values.");
		}

		if(test_count == 0 || test_count < 10)
		{
			throw std::invalid_argument("Invalid test_count values.");
		}

		// Generate two large prime numbers
		BigInteger PrimeNumberA = GeneratePrimeNumber(bit_count / 2, test_count);
		BigInteger PrimeNumberB = GeneratePrimeNumber(bit_count / 2, test_count);
		
		// Calculate n = p * q
		BigInteger AlgorithmModulus = PrimeNumberA * PrimeNumberB;

		// Calculate totient(n) = phi(n) = (p - 1) * (q - 1)
		BigInteger Totient_PhiFunctionValue = (PrimeNumberA - ONE) * (PrimeNumberB - ONE);

		BigInteger EncryptExponent = 0;
		//Enable security and performance optimization?
		if(is_pkcs)
			EncryptExponent = 65537;
		else
		{
			// Choosing exponents for RSA encryption
			// encryption_exponents = RandomNumberRange(2, 2^{bit\_count} - 1)
			do 
			{
				EncryptExponent = BigInteger::RandomGenerateNBit(bit_count / 2);

				// Ensure that e is odd and greater than 1
				if (EncryptExponent.IsEven() || EncryptExponent == ONE)
				{
					continue;
				}

				// Check if e is not 3
				// encryption_exponents equal 3 is not safe
				if(EncryptExponent == THREE)
				{
					continue;
				}

				// Check if gcd(e, totient(n)) = 1
				if (BigInteger::GCD(EncryptExponent, Totient_PhiFunctionValue) != ONE)
				{
					continue;
				}

				break;
			} while (true);
		}

		// Choosing exponents for RSA decryption
		// Need calculate d such that decryption_exponents = encryption_exponents^{-1} (mod totient(n))
		BigInteger DecryptExponent = BigInteger::ModuloInverse(EncryptExponent, Totient_PhiFunctionValue);

		RSA::RSA_AlgorithmNumbers AlgorithmNumbers = RSA::RSA_AlgorithmNumbers();
		AlgorithmNumbers.EncryptExponent = EncryptExponent;
		AlgorithmNumbers.DecryptExponent = DecryptExponent;
		AlgorithmNumbers.AlgorithmModulus = AlgorithmModulus;

		return AlgorithmNumbers;
	}

	void RSA::Encryption(BigInteger& PlainMessage, const BigInteger& EncryptExponent, const BigInteger& AlgorithmModulus )
	{
		PlainMessage.PowerWithModulo(EncryptExponent, AlgorithmModulus);
	}

	void RSA::Decryption(BigInteger& CipherMessage, const BigInteger& DecryptExponent, const BigInteger& AlgorithmModulus )
	{
		CipherMessage.PowerWithModulo(DecryptExponent, AlgorithmModulus);
	}
}  // namespace TwilightDream::CryptographyAsymmetric
