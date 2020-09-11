//
// Created by marus on 12/12/2018.
//

#ifndef CRYPOGRAPHY_NUMBERTHEORY_ALGORITHMS_ECPOINT_H
#define CRYPOGRAPHY_NUMBERTHEORY_ALGORITHMS_ECPOINT_H

#include "BigInteger.h"

using namespace std;
using namespace LongArithmetic;

class ECPoint {
public:
    BigInteger x;
    BigInteger y;
    ECPoint(BigInteger _x, BigInteger _y){
        x = _x;
        y = _y;
    };

};


#endif //CRYPOGRAPHY_NUMBERTHEORY_ALGORITHMS_ECPOINT_H
