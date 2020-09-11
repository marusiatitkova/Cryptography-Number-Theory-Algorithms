#include <iostream>
#include <map>
#include "Algorithms.h"

bool Algorithms::is_Prime(const BigInteger & n) {
    if (n == 1){
        return false;
    }
    if (n < 100000000) {
        for (long long i = 2; BigInteger(i) * BigInteger(i) <= n; ++i) {
            if (n % i == 0) {
                return false;
            }
        }
        return true;
    }
    return isPrime_MillerRabin(n);
}

BigInteger Algorithms::pollardRho(BigInteger n) {
    BigInteger x = random(1, n - 2);
    BigInteger y = 1; //  F(F(y)) where F(x) = x^2 - 1
    int step = 2; //j (xj)
    BigInteger g = gcd(n, abs(x - y));

    for (int i = 0; g == 1; i++) { // x_i = F^i(x0) y_i = x_2i = F^2i(x0)
        if (i == step ) {
            y = x;
            step = step * 2;
        }
        x = add_mod(x*x, 1, n); // (x * x + 1) % n
        g = gcd(n, abs(x - y));
    }

    return g;
}

BigInteger Algorithms::random(BigInteger minN, BigInteger maxN) {
    return (BigInteger(rand()) % (maxN - minN + 1)) + minN;
}

vector<BigInteger> Algorithms::pollardRhoFactorization(BigInteger n) {
    prime_divisors.clear();

    if (n.is_Negative() || n == zero) {
        throw invalid_argument("number is not positive, is prime, = 0 or = 1");
    }
// if 1 return 1
    factorization(n);

    return prime_divisors;
}

void Algorithms::factorization(BigInteger n) {
    if (is_Prime(n)) {
        prime_divisors.push_back(n);
        return;
    }

    if  (n == 1) {
        return;
    }
    BigInteger divisor = pollardRho(n);

    if (divisor > 1) {
        if (is_Prime(divisor)) {
            prime_divisors.push_back(divisor);
        }
        else
            factorization(divisor);
    }

    BigInteger nextDivisor = n / divisor;
    if (nextDivisor == 1)
        return;
    else {
        factorization(nextDivisor);
    }
}

// a^x = b mod m
// find x
// a^in = b * a^j  i -  giant step; j - baby step
BigInteger Algorithms::discreteLogarithm(BigInteger a, BigInteger b, BigInteger m) {
    BigInteger n = sqrt(m) / 100 + 1;
    BigInteger an = pow_mod(a, n, m);
    BigInteger cur_a = an;
    map<BigInteger, BigInteger> value;
//store values a^ni
    for (int i = 1; BigInteger(i) <= n; i++) {
        value[cur_a] = i;
        cur_a = multiply_mod(cur_a, an, m);
    }
    cur_a = b;
    for (int j = 0; BigInteger(j) < n; j++) {
        if (value.count(cur_a)) {
            BigInteger ans = value[cur_a] * n - j;
            if (ans < m)
                return ans;
        }
        cur_a = multiply_mod(cur_a, a, m);
    }
    return BigInteger(-1);
}

// phi(n) = n * (1 - 1/p1) * ... * (1 - 1/pk) pi - prime numbers
BigInteger Algorithms::eulerFunction(BigInteger n) {
    if (n == 1) {
        return 1;
    }
    auto primes = pollardRhoFactorization(n);
    //delete duplicates
    sort(primes.begin(), primes.end());
    auto last = unique(primes.begin(),primes.end());
    primes.erase(last, primes.end());

    BigInteger res = n;

    for (int i = 0; i < primes.size(); i++) {
        res -= res/primes[i];
    }
    return res;
}

// mu(n) = {−1, 0, 1} depending on the factorization of n into prime factors
// 1 -> 1
// 0 -> has a squared prime factor
// (-1)^k -> n = p1 * .. * pk, where pi - different prime factors
int Algorithms::mobiusFunction(BigInteger n) {
    if (n == BigInteger(1)) {
        return 1;
    }
    auto primes = pollardRhoFactorization(n);
    sort(primes.begin(), primes.end());
    int p = 1;
    for (int i = 1; i <= primes.size(); i++) {
        if(primes[i] == primes[i -1]){
            p++;
        }
    }

    if (p > 1) {
        return 0;
    } else {
        return (int)(pow(-1, primes.size()));
    }
}

// a / p  a - integer, p - prime
// 0 -> a %  p = 0
// 1 -> a is quadratic residue modulo p => x : x^2 = a mod p
// -1 -> else
int Algorithms::legendreSymbol(BigInteger a, BigInteger p) {
    if (a % p == 0) {
        return 0;
    }

    if (pow_mod(a, (p - 1) / 2, p) == 1) {
        return 1;
    }

    return -1;
}

// a / n  = (a / p1) * ... * (a / pk) - legendre symbols
// n = p1 * ... * pk - primes
int Algorithms::jacobiSymbol(BigInteger a, BigInteger n) {
    if (a == BigInteger(1)) {
        return 1;
    }

    auto primes = pollardRhoFactorization(n);
    int symbol = 1;
    for (BigInteger p : primes) {
        symbol *= legendreSymbol(a, p);
        if (symbol == 0) {
            return 0;
        }
    }
    return symbol;
}

bool Algorithms::isPrime_MillerRabin(const BigInteger &n, int k) {
    if (n % 2 == 0) {
        return false;
    }
    BigInteger m = n - 1;
    int t = 0;

    while (m % 2 == 0) {
        m /= 2;
        t++;
    }

    for (int i = 0; i < k; i++) {
        BigInteger a = random(2, n - 2);
        BigInteger u = pow_mod(a, m, n);

        if (u == BigInteger(1) || u == n - 1) {
            continue;
        }

        bool s = false;
        for (int j = 1; j < t; j++) {
            u = multiply_mod(u, u, n);
            if (u == 1) {
                return false;
            }
            if (u == n - 1) {
                s = true;
                break;
            }
        }
        if (!s) {
            return false;
        }
    }
    return true;
}

// x^2 = n mod p
pair<BigInteger, BigInteger> Algorithms::discreteSquareRoot(BigInteger n, BigInteger p) {
    if (legendreSymbol(n, p) == -1) {
        return {BigInteger(-1), BigInteger(-1)};
    }
    BigInteger a;
    for (int i = 0; i < 100; i++) {
        a = random(1, p - 1);
        if (legendreSymbol(a * a - n, p) == -1) {
            break;
        }
    }
    if(legendreSymbol(a * a - n, p) != -1) {
        return {BigInteger(-1), BigInteger(-1)};
    }

    pair<BigInteger, BigInteger> vals{a, 1};
    // w^2 = a * a - n
    // exp = (p + 1)/2
    auto res = powModComplex(vals, (p + 1) / 2, p, a * a - n);

    if (res.second != BigInteger(0)) {
        return {BigInteger(0), BigInteger(0)};
    }

    if (res.first * res.first % p != n) {
        return {BigInteger(0), BigInteger(0)};
    }

    return {res.first, (p - res.first) % p};
}

pair<BigInteger, BigInteger>
Algorithms::multiplyField(pair<BigInteger, BigInteger> a, pair<BigInteger, BigInteger> b, BigInteger m, BigInteger w2) {
    return {(a.first * b.first + a.second * b.second * w2 + m) % m, (a.first * b.second + a.second * b.first + m) % m};
}

pair<BigInteger, BigInteger>
Algorithms::powModComplex(pair<BigInteger, BigInteger> vals, BigInteger n, BigInteger m, BigInteger w2) {
    if (n == 0) {
        return {BigInteger(1) % m, BigInteger(0)};
    }

    if (n == 1) {
        return {vals.first % m, vals.second % m};
    }

    if (n % 2 == 1) {
        auto res = multiplyField(powModComplex(vals, n - 1, m, w2), vals, m, w2);
        return {res.first % m, res.second % m};
    }

    auto ans = powModComplex(vals, n / 2, m, w2);
    ans = multiplyField(ans, ans, m, w2);
    return {ans.first % m, ans.second % m};
}

bool Algorithms::isPrime_SolovayStrassen(const BigInteger &n, int k) {
    if (n < 2)
        return false;
    if (n != 2 && n % 2 == 0)
        return false;

    for (int i = 0; i < k; i++)
    {
        // Generate a random number a
        BigInteger a = random(2, n - 2);
        BigInteger jacobian = (n + jacobiSymbol(a, n)) % n;
        BigInteger mod = pow_mod(a, (n - 1) / 2, n);

        if (jacobian == 0 || mod != jacobian)
            return false;
    }
    return true;
}

pair<ECPoint, BigInteger> Algorithms::generateKeys(EllipticCurve curve) {
    BigInteger n = curve.getOrder();

    BigInteger privateKey = random(1, n -1);
    ECPoint publicKey = curve.multiply(curve.getBasePoint(), privateKey);

    return {publicKey, privateKey};
}

pair<ECPoint, ECPoint> Algorithms::encryption(EllipticCurve curve, BigInteger message) {
    // get public key
    auto keys = generateKeys(curve);
    ECPoint publicKey = keys.first;

    // message to point
    BigInteger M1 = message;
    // y^2 = x^3 + a*x + b
    BigInteger n = (pow(message, 3) + curve.getA() * message + curve.getB()) % curve.getP();
    BigInteger M2 = discreteSquareRoot(n, curve.getOrder()).second;
    ECPoint M = ECPoint(M1, M2);
    ECPoint base = curve.getBasePoint();


    BigInteger k = random(1, curve.getOrder() - 1);

    ECPoint C1 = curve.multiply(base, k);
    ECPoint C2 = curve.add(M, curve.multiply(publicKey, k));

    return {C1, C2};
}

BigInteger Algorithms::decryption(EllipticCurve curve, BigInteger message) {
    // get secret key
    auto  keys = generateKeys(curve);
    BigInteger privateKey = keys.second;

    // get cryptogram
    auto cryptogram = encryption(curve, message);
    ECPoint C1 = cryptogram.first;
    ECPoint C2 = cryptogram.second;

    ECPoint temp = curve.multiply(C1, privateKey);


    temp = curve.neg(temp);
    ECPoint decryptedMessage = curve.add(C2, temp);

    BigInteger decrypt = decryptedMessage.x;
    return decrypt;
}

