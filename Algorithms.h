#ifndef CRYPOGRAPHY_NUMBERTHEORY_ALGORITHMS_ALGORITHMS_H
#define CRYPOGRAPHY_NUMBERTHEORY_ALGORITHMS_ALGORITHMS_H

#include "BigInteger.h"
#include "EllipticCurve.h"
#include <algorithm>

using namespace std;
using namespace LongArithmetic;

class Algorithms {

private:
    vector<BigInteger> prime_divisors;

    BigInteger random(BigInteger, BigInteger);

    void factorization(BigInteger); // recursive function to find all divisors of a number

    bool isPrime_SolovayStrassen(const BigInteger &, int k = 5);

    pair<BigInteger, BigInteger> multiplyField(pair<BigInteger, BigInteger>, pair<BigInteger, BigInteger>, BigInteger, BigInteger);

    pair<BigInteger, BigInteger> powModComplex(pair<BigInteger, BigInteger>, BigInteger, BigInteger, BigInteger);

    pair<ECPoint, BigInteger> generateKeys(EllipticCurve);

public:
    bool is_Prime(const BigInteger &);

    BigInteger pollardRho(BigInteger); // for factorization

    vector<BigInteger> pollardRhoFactorization(BigInteger);

    BigInteger discreteLogarithm(BigInteger, BigInteger, BigInteger);

    BigInteger eulerFunction(BigInteger);

    int mobiusFunction(BigInteger);

    int legendreSymbol(BigInteger, BigInteger);

    int jacobiSymbol(BigInteger, BigInteger);

    bool isPrime_MillerRabin(const BigInteger &, int k = 4);

    pair<BigInteger, BigInteger> discreteSquareRoot(BigInteger, BigInteger);

    pair<ECPoint, ECPoint> encryption(EllipticCurve, BigInteger);

    BigInteger decryption(EllipticCurve, BigInteger);

};


#endif //CRYPOGRAPHY_NUMBERTHEORY_ALGORITHMS_ALGORITHMS_H
