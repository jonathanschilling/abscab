# abscab
Accurate Biot-Savart routines with Correct Asymptotic Behaviour

The main article on this software can be found here: [abscab_main.pdf](https://github.com/jonathanschilling/abscab/blob/master/article/abscab_main.pdf)

## Verification Procedure

### Test Points
A set of test points is defined at which the implementations
will be tested against the arbitrary-precision reference will be tested.
These points are computed in Java, since it is a strictly-typed language
and binaries are expected to given platform-independent results.
Nevertheless they are still of finite (64-bit) precision.
The code to generate the test points is [`src/test/java/GenerateTestKnots.java`](https://github.com/jonathanschilling/abscab/blob/master/src/test/java/GenerateTestKnots.java).
