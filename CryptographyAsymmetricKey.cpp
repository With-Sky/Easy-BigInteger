#include "CryptographyAsymmetricKey.hpp"

namespace TwilightDream::CryptographyAsymmetric
{
	void RSA::GeneratePrimesInParallelFunction(size_t bit_count, std::vector<FindPrimeState>& primes, size_t primes_index)
	{
		BigInteger PrimeNumber = BigInteger::RandomGenerateNBit(bit_count);
		PrimeNumber = MIN + (PrimeNumber % RANGE_INTEGER);
		
		//Prime numbers cannot be found inside an even number except two.
		if(PrimeNumber.IsEven())
		{
			PrimeNumber |= ONE;
		}

		bool IsPrime = BigInteger::IsPrime(PrimeNumber);

		while ( !IsPrime )
		{
			/*
				https://crypto.stackexchange.com/questions/1970/how-are-primes-generated-for-rsa
				Generate a random 512 bit odd number, say p
				Test to see if p is prime; if it is, return p; this is expected to occur after testing about Log(p)/2∼177 candidates
				Otherwise set p=p+2, goto 2
			*/
			if(!IsPrime)
			{
				PrimeNumber += TWO;
			}

			//size_t BitSize = PrimeNumber.BitSize();
			//std::cout << "Generate BigInteger Bit Size Is :" << BitSize << std::endl;
			IsPrime = BigInteger::IsPrime(PrimeNumber);
		}

		primes[primes_index].IsPrime = IsPrime;
		primes[primes_index].Number = PrimeNumber;
	}

	std::optional<RSA::FindPrimeState> RSA::GeneratePrimesInParallelFunctions(size_t bit_count)
	{
		std::vector<std::future<void>> futures;
		const size_t				   max_thread_count = std::thread::hardware_concurrency();
		std::vector<FindPrimeState>	   prime_map( max_thread_count, FindPrimeState() );

		for ( size_t i = 0; i < prime_map.size(); ++i )
		{
			futures.push_back( std::async( std::launch::async, &RSA::GeneratePrimesInParallelFunction, this, bit_count, std::ref( prime_map ), i ) );
		}

		size_t counter = 0; // 初始化为 0，因为你可能需要等待所有线程完成
		while (counter < futures.size())
		{
			for (size_t i = 0; i < futures.size(); ++i)
			{
				auto& future = futures[i];
				auto status = future.wait_for(std::chrono::seconds(1));
				if (status == std::future_status::ready)
				{
					if (i == counter)
					{
						++counter;
					}
				}
				else if (status == std::future_status::timeout)
				{
					std::this_thread::sleep_for(std::chrono::seconds(1));
				}
			}
		}

		for ( const auto& state : prime_map )
		{
			if ( state.IsPrime )
			{
				std::cout << "Generate BigInteger Is Prime: True" << std::endl;
				return state;
			}
			else
			{
				std::cout << "Generate BigInteger Is Prime: False" << std::endl;
			}
		}

		return std::nullopt;  // Return std::nullopt if no prime is found
	}

	RSA::BigInteger RSA::GeneratePrimeNumber( size_t bit_count )
	{
		this->bit_count = bit_count;
		MIN = BigInteger(2).Power(bit_count);
		MAX = BigInteger(2).Power(bit_count + 1) - ONE;
		RANGE_INTEGER = MAX - MIN + ONE;

		if(this->bit_count == 0)
			return BigInteger(0);

		std::optional<RSA::FindPrimeState> state_data = GeneratePrimesInParallelFunctions(bit_count);
		
		if(state_data.has_value())
			return state_data.value().Number;
		else
			return BigInteger(0);
	}

	RSA::RSA_AlgorithmNumbers RSA::GenerateKeys( size_t bit_count, bool is_pkcs )
	{
		/*if(bit_count == 0 || bit_count < 1024)
		{
			throw std::invalid_argument("Invalid bit_count values.");
		}*/

		// Generate two large prime numbers
		BigInteger PrimeNumberA = GeneratePrimeNumber(bit_count / 2);
		BigInteger PrimeNumberB = GeneratePrimeNumber(bit_count / 2);
		
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
		if(PlainMessage >= AlgorithmModulus)
		{
			return;
		}
	
		PlainMessage.PowerWithModulo(EncryptExponent, AlgorithmModulus);
	}

	void RSA::Decryption(BigInteger& CipherMessage, const BigInteger& DecryptExponent, const BigInteger& AlgorithmModulus )
	{
		if(CipherMessage >= AlgorithmModulus)
		{
			return;
		}
	
		CipherMessage.PowerWithModulo(DecryptExponent, AlgorithmModulus);
	}
}  // namespace TwilightDream::CryptographyAsymmetric
