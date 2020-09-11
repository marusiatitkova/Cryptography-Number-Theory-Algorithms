#include "BigInteger.h"

namespace LongArithmetic {

    BigInteger::BigInteger() {
        number.push_back(0);
        isNegative = false;
    }

    BigInteger::BigInteger(long long _number) {
        string s = "";
        if(_number == 0)
            s = "0";
        bool negative = false;
        if (_number < 0) {
            _number = -_number;
            negative = true;
        }

        while (_number > 0)
            s = char(_number % 10 + '0') + s, _number /= 10;

        if (negative) {
            s = '-' + s;
        }
        *this = BigInteger(s);

        //*this = _number;
    }

    BigInteger::BigInteger(const BigInteger &_number) {
        isNegative = _number.isNegative;
        number = _number.number;
    }

    BigInteger::BigInteger(string _number) {
        if (!_number.empty() && _number[0] == '-') {
            isNegative = true;
            _number.erase(_number.begin());
        } else {
            isNegative = false;
        }

        for (int i = (int) _number.length(); i > 0; i -= 9) {
            if (i >= 9)
                number.push_back(stoi(_number.substr(i - 9, 9)));
            else
                number.push_back(stoi(_number.substr(0, i)));
        }

        deleteLeadingZeroes();
    }

    void BigInteger::deleteLeadingZeroes() {
        while (number.size() > 1 && !number.back()) {
            number.pop_back();
        }
    }

    int BigInteger::getSign() const {
        int sign = 1;
        if (isNegative)
            sign = -1;
        return sign;
    }

// Assign value

    BigInteger &BigInteger::operator = (const BigInteger &num) {

        if (this != &num) {
            isNegative = num.isNegative;
            number = num.number;
        }

        return *this;
    }

// Input

    istream &operator >> (istream &in, BigInteger &num) {
        string s;
        in >> s;
        num = BigInteger(s);
        return in;
    }

//Output

    ostream &operator << (ostream &out, const BigInteger &num) {

        if (!num.number.size()) {
            out << 0;
            return out;
        }

        if (num.isNegative) {
            out << '-';
        }

        out << (num.number.empty() ? 0 : num.number.back());

        for (int i = (int) num.number.size() - 2; i >= 0; --i)
            out << setw(9) << setfill('0') << num.number[i];

        return out;
    }

//Unary minus

    BigInteger BigInteger::operator -() const {
        BigInteger ans(*this);
        if (!ans.number.empty() || ans != zero) {
            ans.isNegative ^= true;
        }
        return ans;
    }

// Add

    BigInteger BigInteger::operator + (const BigInteger &b) const {
        BigInteger a(*this);

        if (!a.isNegative && b.isNegative) //+ -
            return a - (-b);
        if (a.isNegative && !b.isNegative) //- +
            return b - (-a);
        if (a.isNegative && b.isNegative) //- -
            return -((-a) + (-b));


        int carry = 0;
        vector<int>::iterator itA = a.number.begin();
        vector<int>::const_iterator itB = b.number.cbegin();
        for (; itA != a.number.end() || itB != b.number.cend() || carry; ++itA) {
            if (itA == a.number.end()) {
                a.number.push_back(0);
                itA = a.number.end() - 1;
            }
            *itA += carry + (itB != b.number.cend() ? *(itB++) : 0);
            carry = *itA >= base;
            if (carry) {
                *itA -= base;
            }
        }

        return a;
    }

    BigInteger BigInteger::operator + (long long num) const {
        return *this + BigInteger(num);
    }

    BigInteger &BigInteger::operator += (const BigInteger &num) {
        return *this = *this + num;
    }

    BigInteger &BigInteger::operator += (long long num) {
        return *this = *this + BigInteger(num);
    }

//Subtract

    BigInteger BigInteger::operator - (const BigInteger &b) const {
        BigInteger a(*this);

        if (!a.isNegative && b.isNegative) //+ -
            return a + (-b);
        if (a.isNegative && !b.isNegative) //- +
            return -((-a) + b);
        if (a.isNegative && b.isNegative) //- -
            return a + (-b);

        if (b > a) {
            return -(b - a);
        }

        int carry = 0;
        int b_size = b.number.size();
        for (int i = 0; i < b_size || carry; ++i) {
            a.number[i] -= carry + (i < b_size ? b.number[i] : 0);
            carry = a.number[i] < 0;
            if (carry)
                a.number[i] += base;
        }

        a.deleteLeadingZeroes();

        return a;
    }

    BigInteger BigInteger::operator - (long long num) const {
        return *this - BigInteger(num);
    }

    BigInteger &BigInteger::operator -= (const BigInteger &num) {
        return *this = *this - num;
    }

    BigInteger &BigInteger::operator -= (long long num) {
        return *this = *this - BigInteger(num);
    }

// Multiply

    BigInteger BigInteger::operator * (const BigInteger &b) const {
        BigInteger res;
        int a_len = number.size();
        int b_len = b.number.size();
        res.isNegative = isNegative ^ b.isNegative;
        res.number.resize(a_len + b_len);

        for (int i = 0; i < a_len; ++i) {
            for (int j = 0, carry = 0; j < b_len || carry; ++j) {
                long long cur = res.number[i + j] + number[i] * 1ll * (j < b_len ? b.number[j] : 0) + carry;
                res.number[i + j] = int(cur % base);
                carry = int(cur / base);
            }
        }

        res.deleteLeadingZeroes();

        return res;
    }

    BigInteger BigInteger::operator * (long long num) const {
        return *this * BigInteger(num);
    }

    BigInteger &BigInteger::operator *= (const BigInteger &num) {
        return *this = *this * num;
    }

    BigInteger &BigInteger::operator *= (long long num) {
        return *this = *this * BigInteger(num);
    }

// Divide

    BigInteger BigInteger::operator / (const BigInteger &num) const {
        BigInteger a(*this);
        a.isNegative = false;
        BigInteger b(num);
        b.isNegative = false;

        BigInteger l(0);
        BigInteger r(a);
        bool sign = isNegative ^ num.isNegative;
        while (r - l > 1) {
            BigInteger mid = (l + r) / 2;
            if (mid * b <= a) {
                l = mid;
            } else {
                r = mid - 1;
            }
        }

        if (r * b <= a) {
            l = r;
        }
        l.isNegative = sign;

        return l;
    }

    BigInteger BigInteger::operator / (long long num) const {
        BigInteger res(*this);

        int carry = 0; // remainder = carry and quotient = number

        for (int i = (int) number.size() - 1; i >= 0; --i) {
            long long cur = res.number[i] + carry * 1ll * base;
            res.number[i] = int(cur / num);
            carry = int(cur % num);
        }

        res.deleteLeadingZeroes();

        return res;
    }

    BigInteger &BigInteger::operator /= (const BigInteger &num) {
        return *this = *this / num;
    }

    BigInteger &BigInteger::operator /= (long long num) {
        return *this = *this / num;
    }

// Modulo

    BigInteger BigInteger::operator % (const BigInteger &num) const {
        return *this - ((*this / num) * num);
    }

    BigInteger BigInteger::operator % (long long num) const {
        BigInteger res(*this);
        int carry = 0; // remainder = carry and quotient = number

        for (int i = (int) number.size() - 1; i >= 0; --i) {
            long long cur = res.number[i] + carry * 1ll * base;
            res.number[i] = int(cur / num);
            carry = int(cur % num);
        }

        res.deleteLeadingZeroes();

        return carry;
    }

    BigInteger &BigInteger::operator %= (const BigInteger &num) {
        return *this = *this % num;
    }

    BigInteger &BigInteger::operator %= (long long num) {
        return *this = *this % num;
    }

// Compare

    bool BigInteger::operator < (const BigInteger &b) const {
        int sign = this -> getSign();
        int b_sign = b.getSign();

        int number_size = number.size();
        int b_size = b.number.size();

        if (sign != b_sign)
            return sign < b_sign;

        if (number_size != b_size)
            return number_size * sign < b_size * b_sign;

        for (int i = number_size - 1; i >= 0; i--) {
            if (number[i] != b.number[i])
                return number[i] * sign < b.number[i] * b_sign;
        }

        return false;
    }

    bool BigInteger::operator > (const BigInteger &num) const {
        return num < *this;
    }

    bool BigInteger::operator <= (const BigInteger &num) const {
        return !(num < *this);
    }

    bool BigInteger::operator >= (const BigInteger &num) const {
        return !(*this < num);
    }

    bool BigInteger::operator == (const BigInteger &num) const {
        return !(*this < num) && !(num < *this);
    }

    bool BigInteger::operator != (const BigInteger &num) const {
        return *this < num || num < *this;
    }

    bool BigInteger::is_Negative() const {
        return isNegative;
    }

    BigInteger abs(const BigInteger & a) {
        BigInteger res(a);
        res.isNegative = false;
        return res;
    }

    BigInteger mod(const BigInteger & a, const BigInteger & modulo) {
        return (modulo + (a % modulo)) % modulo;
    }

    BigInteger pow(const BigInteger &a, long long n) {
        if (n < 0) {
            throw invalid_argument("exponent is negative number");
        }

        if (n == 0) {
            return BigInteger(1);
        }
        if (n == 1) {
            return a;
        }
        if (n % 2 == 1) {
            return pow(a, n - 1) * a;
        }
        BigInteger b = pow(a, n / 2);

        return b * b;
    }

    BigInteger pow_mod(const BigInteger &a, const BigInteger &n, const BigInteger &modulo) {
        if (n == 0) {
            //return mod(BigInteger(1), modulo);
            return BigInteger(1) % modulo;
        }
        if (n == 1) {
            //return mod(a, modulo);
            return a % modulo;
        }
        if (n % 2 == 1) {
            //return mod((pow_mod(a, n - 1, modulo) * a), modulo);
            return (pow_mod(a, n - 1, modulo) * a) % modulo;
        }
        //BigInteger b = mod((pow_mod(a, n / 2, modulo)), modulo);
        BigInteger b = pow_mod(a, n / 2, modulo) % modulo;

        //return mod((b * b), modulo);
        return (b * b) % modulo;
    }

    BigInteger sqrt(const BigInteger &num) {
        if (num.is_Negative()) {
            throw invalid_argument("sqrt of a negative number");
        }

        BigInteger a(num);
        BigInteger b = (a + 1) / 2;
        while (b < a) {
            a = b;
            b = (b + num / b) / 2;
        }
        return a;

    }

    BigInteger gcd(BigInteger a, BigInteger b) {
        a = abs(a);
        b = abs(b);
        return b > zero ? gcd(b, a % b) : a;
    }

    BigInteger add_mod(const BigInteger & a, const BigInteger & b, const BigInteger & modulo) {
        //return mod((mod(a, modulo) + mod(b, modulo)), modulo);
        return ((a % modulo) + (b % modulo)) % modulo;
    }

    BigInteger subtract_mod(const BigInteger & a, const BigInteger & b, const BigInteger & modulo) {
        //return mod((mod(a, modulo) - mod(b, modulo)), modulo);
        return ((a % modulo) - (b % modulo)) % modulo;
    }

    BigInteger multiply_mod(const BigInteger & a, const BigInteger & b, const BigInteger & modulo) {
        //return mod((mod(a, modulo) * mod(b, modulo)), modulo);
        return ((a % modulo) * (b % modulo)) % modulo;
    }

    BigInteger inverse_mod(BigInteger a, BigInteger m) {
        BigInteger x, y;
        BigInteger inv;
        BigInteger g = gcd_extended_euclid(a, m, x, y);
        if (g != BigInteger(1)) {
            //throw invalid_argument("neg doesn't exist");
            return BigInteger(-1);
        } else {
            // m is added to handle negative x
            //inv = mod((mod(x, m) + m), m);
            inv = ((x % m) + m) % m;
            return inv;
        }
    }

    BigInteger gcd_extended_euclid(BigInteger a, BigInteger b, BigInteger &x, BigInteger &y) {
        if(a == zero){
            x = zero;
            y = BigInteger(1);
            return b;
        }
        BigInteger x1, y1;
        BigInteger d = gcd_extended_euclid(b % a, a, x1, y1);
        x = y1 - (b / a) * x1;
        y = x1;
        return d;
    }

    BigInteger chinese_remainder_theorem(const vector<BigInteger> &a, const vector<BigInteger> &r) {
        BigInteger M(1), x(0);
        int len = a.size();
        for (const BigInteger& a_i : a) {
            M *= a_i;
        }

        vector<BigInteger> m(len), rev_m(len);
        for (int i = 0; i < len; ++i) {
            m[i] = M / a[i];
            rev_m[i] = inverse_mod(m[i], a[i]);
            if (rev_m[i] < zero) {
                return BigInteger(-1);
            }
            //x = mod((x + r[i] * m[i] * rev_m[i]), M);
            x = (x + r[i] * m[i] * rev_m[i]) % M;
        }
        return x;
    }


}












