package main

import (
    "fmt"
    "io"
    "errors"
    "crypto/rsa"
    "crypto/rand"
    "crypto/sha256"
    "math/big"
)

var bigZero = big.NewInt(0)
var bigOne = big.NewInt(1)
var bigTwo = big.NewInt(2)

func GenerateMultiPrimeKey(random io.Reader, nprimes int, bits int) (*rsa.PrivateKey, error) {
    priv := new(rsa.PrivateKey)
    priv.E = 65537
    primes := make([]*big.Int, nprimes)

    if nprimes < 2 {
        return nil, errors.New("crypto/rsa: GenerateMultiPrimeKey: nprimes must be >= 2")
    } else if nprimes == 2 {
        g := new(big.Int)
        p := new(big.Int)
        y := new(big.Int)
        n := new(big.Int).Set(bigOne)
        e := big.NewInt(int64(priv.E))
        priv.D = new(big.Int)
        pminus1 := new(big.Int)
        totient := new(big.Int).Set(bigOne)
        var err error

        for { // pick random prime p of right size, s.t. gcd(e, p-1) = 1
            primes[0], err = rand.Prime(random, bits/2)
            if err != nil {
                return nil, err
            }
            g.GCD(nil, nil, e, pminus1.Sub(primes[0], bigOne))
            if g.Cmp(bigOne) == 0 {
                break
            }
        }
        // fmt.Printf("p: 0x%s\n", primes[0].Text(16))
        primes[1], err = rand.Prime(random, bits/2)
        if p.Mod(primes[1], bigTwo).Cmp(bigZero) == 0 {
            primes[1].Add(primes[1], bigOne)
        }
        // fmt.Printf("q': %v\n", primes[1])
        p.Mul(primes[0], primes[1])
        // fmt.Printf("n': %v\n", p)

        s1 := 4*bits/8
        s2 := 7*bits/8
        for i := 0; i < bits; i++ {
            if s1 <= i && i < s2 {
                n.SetBit(n, i, primes[0].Bit(primes[0].BitLen()/4+i-s1))
            } else {
                n.SetBit(n, i, p.Bit(i))
            }
        }

        primes[1].Div(n, primes[0])
        if p.Mod(primes[1], bigTwo).Cmp(bigZero) == 0 {
            primes[1].Add(primes[1], bigOne)
        }
        // fmt.Printf("q: 0x%s\n", primes[1].Text(16))

        for {
            pminus1.Sub(primes[1], bigOne)
            g.GCD(priv.D, y, e, pminus1)
            //fmt.Printf("g: %v\n", g)
            if g.Cmp(bigOne) != 1 && primes[1].ProbablyPrime(7) {
                break
            }
            y, err = rand.Prime(random, bits/16)
            if err != nil {
                return nil, err
            }
            if p.Mod(y, bigTwo).Cmp(bigZero) != 0 {
                y.Add(y, bigOne)
            }
            primes[1].Xor(primes[1], y)
            n.Mul(primes[0], primes[1])
            //fmt.Printf("q: p Xor %v = %v\n", y, primes[1])
            //fmt.Printf("n: %v\n", n)
        }
        
        pminus1.Sub(primes[0], bigOne)
        totient.Mul(totient, pminus1)
        pminus1.Sub(primes[1], bigOne)
        totient.Mul(totient, pminus1)
        g.GCD(priv.D, y, e, totient)
        if priv.D.Sign() < 0 {
            priv.D.Add(priv.D, totient)
        }
        priv.Primes = primes
        priv.N = n
    } else {
NextSetOfPrimes:
        for {
            todo := bits
            // crypto/rand should set the top two bits in each prime.
            // Thus each prime has the form
            //   p_i = 2^bitlen(p_i) × 0.11... (in base 2).
            // And the product is:
            //   P = 2^todo × α
            // where α is the product of nprimes numbers of the form 0.11...
            //
            // If α < 1/2 (which can happen for nprimes > 2), we need to
            // shift todo to compensate for lost bits: the mean value of 0.11...
            // is 7/8, so todo + shift - nprimes * log2(7/8) ~= bits - 1/2
            // will give good results.
            if nprimes >= 7 {
                todo += (nprimes - 2) / 5
            }
            for i := 0; i < nprimes; i++ {
                var err error
                primes[i], err = rand.Prime(random, todo/(nprimes-i))
                if err != nil {
                    return nil, err
                }
                todo -= primes[i].BitLen()
            }

            // Make sure that primes is pairwise unequal.
            for i, prime := range primes {
                for j := 0; j < i; j++ {
                    if prime.Cmp(primes[j]) == 0 {
                        continue NextSetOfPrimes
                    }
                }
            }

            n := new(big.Int).Set(bigOne)  
            totient := new(big.Int).Set(bigOne)  // φ(n)
            pminus1 := new(big.Int)   
            for _, prime := range primes {
                n.Mul(n, prime)
                pminus1.Sub(prime, bigOne)
                totient.Mul(totient, pminus1)
            }
            if n.BitLen() != bits {
                // This should never happen for nprimes == 2 because
                // crypto/rand should set the top two bits in each prime.
                // For nprimes > 2 we hope it does not happen often.
                continue NextSetOfPrimes
            }

            g := new(big.Int)
            priv.D = new(big.Int)
            y := new(big.Int)
            e := big.NewInt(int64(priv.E))
            g.GCD(priv.D, y, e, totient)

            if g.Cmp(bigOne) == 0 {
                if priv.D.Sign() < 0 {
                    priv.D.Add(priv.D, totient)
                }
                priv.Primes = primes
                priv.N = n

                break
            }
        }
    }

    priv.Precompute()
    return priv, nil
}
func main() {
    bits := 2048
    k, e := GenerateMultiPrimeKey(rand.Reader, 2, bits)
    if e != nil {
        fmt.Println(e)
    }
    p := new(big.Int)
    fmt.Printf("pub N: 0x%s\n", k.N.Text(16))
    fmt.Printf("pub E: %v\n", k.E)
    fmt.Printf("pri D: 0x%s\n", k.D.Text(16))
    secretMessage := []byte("secret")
    label := []byte("orders")
    
    ciphertext, err := rsa.EncryptOAEP(sha256.New(), rand.Reader, &k.PublicKey, secretMessage, label)
    if err != nil {
        fmt.Printf("pub encrypt failed: %v\n", err)
        return
    }
    // fmt.Printf("pub encrypt result: 0x%x\n", string(ciphertext))

    for i:=0; i<bits/8*3; i++ {
        p.SetBit(p, i, k.N.Bit(i+bits/8*4))
    }
    p.Lsh(p, uint(bits/8*1))
    fmt.Printf("key p: 0x%s(%v)\n", k.Primes[0].Text(16), k.Primes[0].ProbablyPrime(15))
    fmt.Printf("key q: 0x%s(%v)\n", k.Primes[1].Text(16), k.Primes[1].ProbablyPrime(15))
    fmt.Printf("recovered p: 0x%s(%v)\n", p.Text(16), p.ProbablyPrime(15))


    plaintext, err := rsa.DecryptOAEP(sha256.New(), rand.Reader, k, ciphertext, label)
    if err != nil {
        fmt.Printf("pri decrypt failed: %v\n", err)
        return
    }

    fmt.Printf("pub encrypt message: %s\n", string(secretMessage))
    fmt.Printf("pri decrypt result: %s\n", string(plaintext))
}



