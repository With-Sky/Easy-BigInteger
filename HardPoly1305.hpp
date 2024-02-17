/*
MIT License

Copyright (c) 2024 Twilight-Dream & With-Sky

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

#if !defined(HARD_POLY1305_WITH_BIGINTEGER)
#define HARD_POLY1305_WITH_BIGINTEGER

#include "BigInteger.hpp"

// 提取子数组函数，类似于 std::span 和 与Python数据切片range
extern std::vector<uint8_t> SubByteArray( const std::vector<uint8_t>& data, ptrdiff_t start, ptrdiff_t end, ptrdiff_t step );

// 将字节数组转换为16进制字符串表示
extern std::string BytesToHexString( const std::vector<uint8_t>& bytes );

// 生成随机字节序列
extern std::vector<uint8_t> generate_random_bytes( size_t size );

class HardPoly1305
{
private:
	using BigSignedInteger = TwilightDream::BigInteger::BigSignedInteger;
	BigSignedInteger p;					  // Poly1305 算法使用的质数
	BigSignedInteger p2;				  // HardPoly1305 算法使用的质数
	BigSignedInteger r = 0;				  // 秘密状态 r
	BigSignedInteger s = 0;				  // 秘密状态 s
	BigSignedInteger clamp_bit_mask = 0;  //比特掩码 - 修剪数值
	BigSignedInteger hash_max_number = BigSignedInteger( 1 ) << 128;

	BigSignedInteger generate_unpredictable_value( const BigSignedInteger& hash_value, const BigSignedInteger& key );


public:
	HardPoly1305() : 
		p( "1361129467683753853853498429727072845819", 10 ),
		p2( "115792089237316195423570985008687907853269984665640564039457584007913129451867", 10 ),
		clamp_bit_mask( "0FFFFFFC0FFFFFFC0FFFFFFC0FFFFFFF", 16 )
	{}

	std::vector<uint8_t> mix_key_and_message( const std::vector<uint8_t>& message, const std::vector<uint8_t>& key );

	std::vector<uint8_t> hard_poly1305_core( const std::vector<uint8_t>& mixed_data, const std::vector<uint8_t>& key );
};

// 测试 HardPoly1305 类
extern void test_hard_poly1305();

#endif	// HARD_POLY1305_WITH_BIGINTEGER