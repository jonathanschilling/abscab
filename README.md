# abscab
Accurate Biot-Savart routines with Correct Asymptotic Behaviour

The main article on this software can be found here: [abscab_main.pdf](https://github.com/jonathanschilling/abscab/blob/master/article/abscab_main.pdf)

## Verification Procedure

### `mpmath` vs. Mathematica

First, the reference implementation using `mpmath` is verified against
an implementation in Mathematica.
For brevity, this is done on a reduced set of test points.

### Test Points
A set of test points is defined at which the implementations
will be tested against the arbitrary-precision reference will be tested.
These points are computed in Java, since it is a strictly-typed language
and binaries are expected to given platform-independent results.
Nevertheless they are still of finite (64-bit) precision.
The code to generate the test points is [`src/test/java/GenerateTestKnots.java`](https://github.com/jonathanschilling/abscab/blob/master/src/test/java/de/labathome/GenerateTestKnots.java).

The test points are saved into text files in [`src/test/resources`](https://github.com/jonathanschilling/abscab/blob/master/src/test/resources):
* [`testPointsRpStraightWireSegment.dat`](https://github.com/jonathanschilling/abscab/blob/master/src/test/resources/testPointsRpStraightWireSegment.dat) contains the value of `r'` at which
  the straight wire segment methods are tested.
* [`testPointsZpStraightWireSegment.dat`](https://github.com/jonathanschilling/abscab/blob/master/src/test/resources/testPointsZpStraightWireSegment.dat) contains the value of `z'` at which
  the straight wire segment methods are tested.
* [`testPointsRpCircularWireLoop.dat`](https://github.com/jonathanschilling/abscab/blob/master/src/test/resources/testPointsRpCircularWireLoop.dat) contains the value of `r'` at which
  the circular wire loop methods are tested.
* [`testPointsZpCircularWireLoop.dat`](https://github.com/jonathanschilling/abscab/blob/master/src/test/resources/testPointsZpCircularWireLoop.dat) contains the value of `z'` at which
  the circular wire loop methods are tested.

Those above files are only provided to have a human-readable equivalent
of the set of test points.
The actual test point data read by the arbitrary-precision software
is in [`testPointsStraightWireSegment.dat`]() for the straight wire segment
and in [`testPointsCircularWireLoop.dat`]() for the circular wire loop.

### Reference Data
The reference data (`A_z` and `B_phi` for a straight wire segment;
`A_phi`, `B_rho` and `B_z` for a circular wire loop) is also available in the folder
[`src/test/resources`](https://github.com/jonathanschilling/abscab/blob/master/src/test/resources).
