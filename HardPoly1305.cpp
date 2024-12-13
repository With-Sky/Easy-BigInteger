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

#include <iomanip>
#include <cassert>
#include "HardPoly1305.hpp"

std::vector<uint8_t> SubByteArray( const std::vector<uint8_t>& data, ptrdiff_t start, ptrdiff_t end, ptrdiff_t step = 1 )
{
	if ( step == 0 )
	{
		return {};
	}

	ptrdiff_t dataSize = static_cast<ptrdiff_t>( data.size() );

	if ( start < 0 )
	{
		start += dataSize;
	}
	if ( end < 0 )
	{
		end += dataSize;
	}

	start = std::max<ptrdiff_t>( 0, start );
	end = std::min<ptrdiff_t>( end, dataSize );

	if ( start >= end )
	{
		return {};
	}

	std::vector<uint8_t> sub_array;
	if ( step > 0 )
	{
		auto first = data.begin() + start;
		auto last = data.begin() + end;
		std::copy_if( first, last, std::back_inserter( sub_array ), [ n = 0, step ]( const uint8_t& ) mutable { return n++ % step == 0; } );
	}
	else
	{
		auto first = data.rbegin() + ( dataSize - end );
		auto last = data.rbegin() + ( dataSize - start );
		std::copy_if( first, last, std::back_inserter( sub_array ), [ n = 0, step ]( const uint8_t& ) mutable { return n++ % -step == 0; } );
	}

	return sub_array;
}

std::string BytesToHexString( const std::vector<uint8_t>& bytes )
{
	std::ostringstream oss;
	for ( const auto& byte : bytes )
	{
		oss << std::hex << std::setw( 2 ) << std::setfill( '0' ) << static_cast<int>( byte );
	}
	return oss.str();
}

std::vector<uint8_t> generate_random_bytes( size_t size )
{
	std::vector<uint8_t>			bytes( size );
	std::random_device				rd;
	std::mt19937					gen( rd() );
	std::uniform_int_distribution<> dis( 0, 255 );
	for ( size_t i = 0; i < size; ++i )
	{
		bytes[ i ] = static_cast<uint8_t>( dis( gen ) );
	}
	return bytes;
}

HardPoly1305::BigSignedInteger HardPoly1305::generate_unpredictable_value( const BigSignedInteger& hash_value, const BigSignedInteger& key )
{
	//std::cout << "\n--- generate_unpredictable_value ---\n";
	//std::cout << "Initial hash_value: " << hash_value.ToString( 10 ) << '\n';
	//std::cout << "Initial key: " << key.ToString( 10 ) << '\n';

	const BigSignedInteger a = BigSignedInteger( hash_value ) + BigSignedInteger( key );
	const BigSignedInteger b = BigSignedInteger( key ) - BigSignedInteger( hash_value );

	// 复杂公式计算 alpha
	BigSignedInteger alpha = ( ( a * a ) + ( b * b * b ) );
	//std::cout << "Computed alpha: " << alpha.ToString( 10 ) << '\n';

	// 应用 模数 p2
	if ( alpha < 0 )
	{
		//std::cout << "Alpha is negative: " << alpha.ToString( 10 ) << '\n';
		alpha = p2 + alpha; //Adjustment of negative numbers to a positive range [0, p2]
		//std::cout << "Adjusted alpha by adding p2: " << alpha.ToString( 10 ) << '\n';
	}
	if ( alpha >= p2 )
	{
		//std::cout << "Alpha exceeds p2: " << alpha.ToString( 10 ) << '\n';
		alpha = alpha % p2;
		//std::cout << "Alpha reduced modulo p2: " << alpha.ToString( 10 ) << '\n';
	}

	size_t key_bit_length = key.BitLength();
	size_t hash_value_bit_length = hash_value.BitLength();
	//std::cout << "key_bit_length: " << key_bit_length << '\n';
	//std::cout << "hash_value_bit_length: " << hash_value_bit_length << '\n';

	const BigSignedInteger& u0 = a;
	//std::cout << "Computed u0 (hash_value + key): " << u0.ToString( 10 ) << '\n';
	const BigSignedInteger u1 = hash_value - key;
	//std::cout << "Computed u1 (hash_value - key): " << u1.ToString( 10 ) << '\n';
	uint64_t		 left_shift_amount = ( u0 % key_bit_length ).ToUnsignedInt();
	uint64_t		 right_shift_amount = ( u1 % hash_value_bit_length ).ToUnsignedInt();
	BigSignedInteger u2 = key << left_shift_amount;
	//std::cout << "Computed u2 (key << (u0 % key_bit_length)): " << u2.ToString( 10 ) << '\n';
	BigSignedInteger u3 = hash_value >> right_shift_amount;
	//std::cout << "hash_value_bit_length: " << hash_value_bit_length << "\n";
	//std::cout << "right_shift_amount: " << right_shift_amount << "\n";
	//std::cout << "Computed u3 (hash_value >> (u1 % hash_value_bit_length)): " << u3.ToString( 10 ) << '\n';

	BigSignedInteger u4 = hash_value * ( u2 + u3 + alpha );
	//std::cout << "Computed u4 (hash_value * (u2 + u3 + alpha)): " << u4.ToString( 10 ) << '\n';

	u4 = u4 % p;
	//std::cout << "Final u4 modulo p: " << u4.ToString( 10 ) << '\n';
	//std::cout << "--- End of generate_unpredictable_value ---\n";

	// 应用 模数 p 规约大小
	return u4;
}

std::vector<uint8_t> HardPoly1305::mix_key_and_message( const std::vector<uint8_t>& message, const std::vector<uint8_t>& key )
{
	// Mix the message and key
	std::vector<uint8_t> mixed_data( message.size(), 0 );
	size_t				 key_index = 0;

	for ( size_t i = 0; i < message.size(); ++i )
	{
		mixed_data[ i ] = ( message[ i ] + key[ key_index ] ) % 256;
		key_index = ( key_index + 1 ) % key.size();
	}

	if ( mixed_data.size() < 32 )
	{
		//mixed_data = mixed_data concatenation key

		mixed_data.reserve( 32 );
		for ( size_t i = mixed_data.size(); i < key.size(); i++ )
		{
			mixed_data.push_back( key[ i ] );
		}
	}

	return mixed_data;
}

std::vector<uint8_t> HardPoly1305::hard_poly1305_core( const std::vector<uint8_t>& mixed_data, const std::vector<uint8_t>& key )
{
	// 初始化混合消息状态
	BigSignedInteger mixed_number = 0;

	// 获取中间字节部分的密钥
	BigSignedInteger key_number = 0;
	key_number.ImportData( false, SubByteArray( key, 7, 23 ) );
	std::cout << "\n--- hard_poly1305_core ---\n";
	std::cout << "Initial key_number (from key[7:23]): " << key_number.ToString( 10 ) << '\n';
	key_number *= 5;
	std::cout << "After multiplying by 5, key_number: " << key_number.ToString( 10 ) << '\n';

	// 初始化(Aaccumulator) hash_value 为 0
	BigSignedInteger hash_value = 0;

	// 初始化秘密状态 Secret status
	r.ImportData( false, SubByteArray( key, 0, 15 ) );
	s.ImportData( false, SubByteArray( key, 16, 31 ) );
	std::cout << "Initial r: ";
	r.Print( 10 );
	std::cout << "Initial s: ";
	s.Print( 10 );

	// 计算 loop_count
	size_t loop_count = ( mixed_data.size() + 15 ) / 16 + 1;
	//std::cout << "Total loop_count: " << loop_count << '\n';

	BigSignedInteger mixed_data_number = 0;
	mixed_data_number.ImportData( false, mixed_data );
	//std::cout << "Initial mixed_data_number : ";
	mixed_data_number.Print( 10 );
	std::vector<uint8_t> mixed_data_span;

	// HardPoly1305 算法核心计算循环
	for ( size_t i = 1; i <= loop_count - 1; ++i )
	{
		//std::cout << "\n--- Loop " << i << " ---\n";

		// 将 r 进行 clamp
		r &= clamp_bit_mask;
		//std::cout << "Clamped r: " << r.ToString( 10 ) << '\n';

		// 获取 mixed_data 的字节切片
		mixed_data_span = SubByteArray( mixed_data, ( i - 1 ) * 16, i * 16 );

		// 重新分配 bytes 大小并拼接 Byte 0x01
		mixed_data_span.push_back( 0x01 );

		// 打印 mixed_data_span 的16进制表示
		//std::cout << "Mixed data span (hex): " << BytesToHexString( mixed_data_span ) << '\n';

		// 计算 mixed_number
		mixed_number.ImportData( false, mixed_data_span );
		//std::cout << "Computed mixed_number from mixed_data_span: " << mixed_number.ToString( 10 ) << '\n';

		// 更新 hash_value
		hash_value += mixed_number;
		//std::cout << "Updated hash_value after adding mixed_number: " << hash_value.ToString( 10 ) << '\n';

		// 首先 应用 秘密部分 r 进行倍增 然后 应用 模数 p
		hash_value = ( r * hash_value ) % p;
		//std::cout << "Updated hash_value after multiplying by r and mod p: " << hash_value.ToString( 10 ) << '\n';

		// 最后 应用 秘密部分 s
		hash_value += s;
		//std::cout << "Updated hash_value after adding s: " << hash_value.ToString( 10 ) << '\n';

		// 更新秘密状态 r 和 s
		r = generate_unpredictable_value( hash_value, mixed_number );
		s = generate_unpredictable_value( hash_value, key_number );
	}

	// 规约大小 hash = hash (mod 2^128)
	hash_value = hash_value % hash_max_number;
	std::cout << "\nFinal hash_value mod 2^128: " << hash_value.ToString( 10 ) << '\n';
	std::cout << "--- End of hard_poly1305_core ---\n";

	// 返回结果
	std::vector<uint8_t> hash_result;
	bool				 is_negative = false;
	hash_value.ExportData( is_negative, hash_result, 16 );

	std::cout << "Final Tag/Hash Bytes Data:\n";
	for ( const auto& byte : hash_result )
	{
		std::cout << std::hex << std::setw( 2 ) << std::setfill( '0' ) << static_cast<unsigned int>( byte );
	}
	std::cout << '\n';

	return hash_result;
}


void test_hard_poly1305()
{
	using BigInteger = TwilightDream::BigInteger::BigInteger;

	// 创建 HardPoly1305 对象
	HardPoly1305 hard_poly1305;

	// 生成消息和随机密钥
	std::string			 string_message = std::string( "Hello, world!" );
	std::vector<uint8_t> message( string_message.begin(), string_message.end() );
	std::vector<uint8_t> key( 32, 'A' );

	// 计算混合消息数据 将消息与密钥混合
	std::vector<uint8_t> mixed_data = hard_poly1305.mix_key_and_message( message, key );

	auto format_flags = std::cout.flags();

	// 输出密钥、消息和标签
	std::cout << "Message: ";
	for ( const auto& byte : message )
	{
		std::cout << std::hex << std::setw( 2 ) << std::setfill( '0' ) << static_cast<unsigned int>( byte );
	}
	std::cout << '\n';
	std::cout.flags( format_flags );
	std::cout << "message byte length: " << message.size() << std::endl;

	std::cout << "Key: ";
	for ( const auto& byte : key )
	{
		std::cout << std::hex << std::setw( 2 ) << std::setfill( '0' ) << static_cast<unsigned int>( byte );
	}
	std::cout << '\n';
	std::cout.flags( format_flags );
	std::cout << "key byte length: " << key.size() << std::endl;

	std::cout << "MixData: ";
	for ( const auto& byte : mixed_data )
	{
		std::cout << std::hex << std::setw( 2 ) << std::setfill( '0' ) << static_cast<unsigned int>( byte );
	}
	std::cout << '\n';
	std::cout.flags( format_flags );
	std::cout << "mixed_data byte length: " << mixed_data.size() << std::endl;

	// 计算消息的标签
	std::vector<uint8_t> tag = hard_poly1305.hard_poly1305_core( mixed_data, key );

	std::cout.flags( format_flags );
}
