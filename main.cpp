#include "BigInteger.hpp"
#include "BigIntegerTest.hpp"
#include "BigFraction.hpp"
#include "BinaryCipherTest.hpp"
#include "CryptographyAsymmetricKey.hpp"
#include "HardPoly1305.hpp"

#include <iomanip>
#include <cassert>

inline void test_negative_big_fraction()
{
	// Initialize a BigFraction with a negative fraction -3/7
	TwilightDream::BigInteger::BigInteger numerator( 3 );	 // Positive numerator
	TwilightDream::BigInteger::BigInteger denominator( 7 );	 // Positive denominator
	int									  sign = -1;		 // Negative sign for the fraction

	// Create the BigFraction with the negative sign
	TwilightDream::BigFraction::BigFraction negativeFraction( numerator, denominator, sign );

	// Set a very high precision denominator for the final value
	TwilightDream::BigInteger::BigInteger FullPrecision(
		"1000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000", 
		10
	);
	negativeFraction.SetFullPrecision(FullPrecision);

	// Output the fraction to check initialization
	std::cout << "Negative fraction (-3/7): " << negativeFraction << std::endl;

	// Perform operations on the negative fraction
	std::cout << "Negative fraction squared: " << negativeFraction.Power( 2 ) << std::endl;
	std::cout << "Square root of the negative fraction: " << negativeFraction.Sqrt() << std::endl;  // Should handle or raise an error if negative
	std::cout << "Negative fraction multiplied by 2/3: " << ( negativeFraction * TwilightDream::BigFraction::BigFraction( TwilightDream::BigInteger::BigInteger( 2 ), TwilightDream::BigInteger::BigInteger( 3 ), 1 ) ) << std::endl;

	// Check inversion (should handle the sign correctly)
	std::cout << "Reciprocal negative fraction: " << negativeFraction.Reciprocal() << std::endl;
}

auto main( int argument_cout, char* argument_vector[] ) -> int
{
	//TEST PASSED
	/*TwilightDream::BigInteger::BigInteger TestNumber( "123456789", 10 );
	TwilightDream::BigInteger::BigInteger TestNumberPower = TwilightDream::BigInteger::BigInteger( TestNumber ).Power( 14 );
	TwilightDream::BigInteger::ShiftingKthRoot kthRootCalculator( 14 );
	TwilightDream::BigInteger::BigInteger	   TestNumberKthRoot = kthRootCalculator( TestNumberPower );
	if ( TestNumberKthRoot == TestNumber )
	{
		std::cout << "True" << std::endl;
	}
	else
	{
		std::cout << "False" << std::endl;
	}*/

	//TEST PASSED
	/*auto pi = TwilightDream::BigFraction::BigFraction::GenerateSrinivasaRamanujanPI();
	pi.PrecisionMode = TwilightDream::BigFraction::DecimalPrecisionMode::Fixed;
	pi.FixedPrecisionCount = 10000;
	std::cout << pi << std::endl;*/

	//TEST PASSED
	//std::cout << TwilightDream::BigInteger::BigInteger::Factorial(10000).ToString(10) << std::endl;
	
	//test_negative_big_fraction();

	TwilightDream::BigFraction::BigFraction TestBigFraction(12,48);
	TestBigFraction.PrecisionMode = TwilightDream::BigFraction::DecimalPrecisionMode::Fixed;
	TestBigFraction.FixedPrecisionCount = 100;
	std::cout << TestBigFraction.Power(3) << std::endl;
	std::cout << TestBigFraction.Power(TwilightDream::BigFraction::BigFraction(1, 4)) << std::endl;
	std::cout << TestBigFraction.Sqrt() << std::endl;
	std::cout << TestBigFraction.Cbrt() << std::endl;
	std::cout << TestBigFraction.Log() << std::endl;
	std::cout << TestBigFraction.Log10() << std::endl;
	//FIXME
	std::cout << TestBigFraction.NthRoot(4) << std::endl;

	TestBigFraction.PrecisionMode = TwilightDream::BigFraction::DecimalPrecisionMode::Full;
	TwilightDream::BigInteger::BigInteger FullPrecision("1000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000", 10);
	TestBigFraction.SetFullPrecision(FullPrecision);
	std::cout << TestBigFraction.Power(3) << std::endl;
	std::cout << TestBigFraction.Sqrt() << std::endl;
	std::cout << TestBigFraction.Cbrt() << std::endl;
	std::cout << TestBigFraction.Log() << std::endl;
	std::cout << TestBigFraction.Log10() << std::endl;
	//FIXME
	std::cout << TestBigFraction.NthRoot(4) << std::endl;
	
	TwilightDream::BigInteger::Test::test_all();
	TwilightDream::BigInteger::Test::test_custom_data();

	using PrimeNumberTester = TwilightDream::PrimeNumberTester;
	using BigInteger = TwilightDream::BigInteger::BigInteger;

	//BigInteger BigPrime(std::string("16158503035655503650357438344334975980222051334857742016065172713762327569433945446598600705761456731844358980460949009747059779575245460547544076193224141560315438683650498045875098875194826053398028819192033784138396109321309878080919047169238085235290822926018152521443787945770532904303776199561965192760957166694834171210342487393282284747428088017663161029038902829665513096354230157075129296432088558362971801859230928678799175576150822952201848806616643615613562842355410104862578550863465661734839271290328348967522998634176499319107762583194718667771801067716614802322659239302476074096777926805529798115243"), 10);
	//PrimeNumberTester NumberTester;
	//std::cout << std::boolalpha << NumberTester.IsPrime(BigPrime) << std::endl;

	//TEST PASSED
	//TwilightDream::CryptographyAsymmetric::RSA::SelfSanityCheck(4096, 10);

	BinaryCipher BinaryCipherInstance;
	BinaryCipherInstance.Test();
	
	BinaryCipherNaive BinaryCipherNaiveInstance;
	BinaryCipherNaiveInstance.Test();

	//大整数应该不会计算错误。因为已经全部都测试过了。那唯一可能的是我算法写错了，或者还有什么我们的别的不知道的原因
	//Large integers shouldn't be miscalculated. That's because it's all been tested. So the only possibility is that I wrote the algorithm wrong, or something else we don't know.
	//2020-08-25 Updated: Fixed
	test_hard_poly1305();

	return 0;
}
