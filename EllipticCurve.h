//
// Created by marus on 12/12/2018.
//

#ifndef CRYPOGRAPHY_NUMBERTHEORY_ALGORITHMS_ELIPTICCURVE_H
#define CRYPOGRAPHY_NUMBERTHEORY_ALGORITHMS_ELIPTICCURVE_H

#include "BigInteger.h"
#include "ECPoint.h"

using namespace std;
using namespace LongArithmetic;

class EllipticCurve {

private:
    // finite field Fp
    BigInteger p;
    // parameters for y^2 = x^3 + a*x + b
    BigInteger a;
    BigInteger b;
    //base point
    ECPoint G = ECPoint(1,1);
    //order of G
    BigInteger n;
    //cofactor
    BigInteger h;
public:

    EllipticCurve() {
        p = BigInteger("6277101735386680763835789423207666416083908700390324961279");
        a = BigInteger(-3);
        b = BigInteger("2455155546008943817740293915197451784769108058161191238065");
        G = ECPoint(BigInteger("60204628237568865675821348058752611191669876636884684818"), BigInteger("174050332293622031404857552280219410364023488927386650641"));
        n = BigInteger("6277101735386680763835789423176059013767194773182842284081");
        h = 1;
    };

    BigInteger getOrder() {
        return n;
    }

    ECPoint getBasePoint() {
        return G;
    }

    BigInteger getA() {
        return a;
    }

    BigInteger getB() {
        return b;
    }

    BigInteger getP() {
        return p;
    }

    ECPoint Double(ECPoint point) {
        ECPoint r = ECPoint(1,1);
        BigInteger L = (BigInteger(3) * point.x * point.x) % p * inverse_mod(BigInteger(2) * point.y, p) % p;
        r.x = (L * L - BigInteger(2) * point.x) % p;
        r.y = (L * (point.x - r.x) - point.y) % p;
        return r;
    }

    ECPoint add(ECPoint point, ECPoint q) {

        if (point.x == q.x && point.y == q.y) return Double(point);

        ECPoint r = ECPoint(1, 1);
        //BigInteger L = (q.y - point.y) / (q.x - point.x);
        BigInteger L = (q.y - point.y) % p * inverse_mod(q.x - point.x, p) % p;
        r.x = (L * L - point.x - q.x) % p;
        r.y = (L * (point.x - r.x) - point.y) % p;
        return r;
    }

    ECPoint multiply(ECPoint point, BigInteger x) {
        ECPoint temp = point;

        while (x != 0) {

            if ((x % 2) != 0)
            {
                if((temp.x == point.x)||(temp.y == point.y))
                    temp = Double(temp);
                else
                    temp = add(temp, point);
                x = x - 1;
            }
            x = x / 2;
            point = Double(point);
        }
        return temp;
    }

    ECPoint neg(ECPoint point) {
        //return (ECPoint){ point.x, -point.y };
        return (ECPoint) {point.x, (-point.y + p) % p};
    }
};

#endif //CRYPOGRAPHY_NUMBERTHEORY_ALGORITHMS_ELIPTICCURVE_H
