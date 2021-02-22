# Numerical-Method-for-Comparison-on-Homomorphically-Encrypted-Numbers

Numerical Method for Comparison on Homomorphically Encrypted Numbers [1]

Description:

This paper presents a method that approximately computes the min/max and comparison operations on homomorphically encrypted numbers.

The numbers are encrypted word-wise which the authors show gives much better performance than bit-wise encryption on polynomial operations and comparable performance for comparison operations.

The authors present a way of computing min/max and comparison operations using an iterative approach, the results being a good approximation of the real ones. It must be said that the error is known and can be reduced with a performance compromise.

I have implemented the min/max algorithms and the comparison operation in C++, starting from an existing library (HEAAN) which has been implemented by the authors for a previous paper they wrote, called Homomorphic Encryption for Arithmetic of Approximate Numbers [2].

At the end of the paper, the authors show the results they get by implementing the algorithms using HEAAN[3]. I have done a similar performance evaluation which is presented in the Results excel file.

The source code for the algorithms presented in the paper can be found in HEAAN/src/main.cpp and can be compiled using the makefile in the same folder.


[1] https://eprint.iacr.org/2019/417.pdf
[2] https://eprint.iacr.org/2016/421.pdf
[3] https://github.com/snucrypto/HEAAN 


