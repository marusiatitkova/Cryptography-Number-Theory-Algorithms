#ifndef LONG_ARITHMETIC_BIGINTEGER_H
#define LONG_ARITHMETIC_BIGINTEGER_H

#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include <sstream>
using namespace std;

namespace LongArithmetic {

    class BigInteger {

    private:
        vector<int> number;
        bool isNegative;
        static const int base = 1000 * 1000 * 1000; //1e9

        void deleteLeadingZeroes();

    public:

        BigInteger();

        BigInteger(string);

        BigInteger(long long);

        BigInteger(const BigInteger &);

        BigInteger &operator = (const BigInteger &);

        BigInteger operator -() const;  // unary minus

        //Add
        BigInteger operator + (const BigInteger &) const;

        BigInteger operator + (long long) const;

        BigInteger &operator += (const BigInteger &);

        BigInteger &operator += (long long);

        //Subtract
        BigInteger operator - (const BigInteger &) const;

        BigInteger operator - (long long) const;

        BigInteger &operator -= (const BigInteger &);

        BigInteger &operator -= (long long);

        //Multiply
        BigInteger operator * (const BigInteger &) const;

        BigInteger operator * (long long) const;

        BigInteger &operator *= (const BigInteger &);

        BigInteger &operator *= (long long);

        //Divide
        BigInteger operator / (const BigInteger &) const;

        BigInteger operator / (long long) const;

        BigInteger &operator /= (const BigInteger &);

        BigInteger &operator /= (long long);

        //Modulo
        BigInteger operator % (const BigInteger &) const;

        BigInteger operator % (long long) const;

        BigInteger &operator %= (const BigInteger &);

        BigInteger &operator %= (long long);

        //Compare
        bool operator < (const BigInteger &) const;

        bool operator > (const BigInteger &) const;

        bool operator <= (const BigInteger &) const;

        bool operator >= (const BigInteger &) const;

        bool operator == (const BigInteger &) const;

        bool operator != (const BigInteger &) const;

        int getSign() const;

        bool is_Negative () const;

        friend BigInteger abs(const BigInteger&);


        //input and output
        friend ostream &operator << (ostream &, const BigInteger &);

        friend istream &operator >> (istream &, BigInteger &);
    };

    //correct math mod
    BigInteger mod(const BigInteger&, const BigInteger&);

    BigInteger pow(const BigInteger &, long long);

    BigInteger pow_mod(const BigInteger &, const BigInteger &, const BigInteger &);

    BigInteger add_mod(const BigInteger &, const BigInteger &, const BigInteger &);

    BigInteger subtract_mod(const BigInteger &, const BigInteger &, const BigInteger &);

    BigInteger multiply_mod(const BigInteger &, const BigInteger &, const BigInteger &);

    BigInteger sqrt(const BigInteger &);

    BigInteger gcd(BigInteger, BigInteger );

    BigInteger gcd_extended_euclid(BigInteger, BigInteger, BigInteger &, BigInteger &);

    BigInteger inverse_mod(BigInteger, BigInteger);

    BigInteger chinese_remainder_theorem(const vector<BigInteger> &, const vector<BigInteger> &);

    const BigInteger zero(0);
}
#endif //LONG_ARITHMETIC_BIGINTEGER_H
