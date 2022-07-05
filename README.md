# abscab
Accurate Biot-Savart routines with Correct Asymptotic Behaviour
<img src="abscab_logo.png" alt="ABSCAB logo" width="200" align="right"/>

The main article on this software can be found here: [abscab_main.pdf](article/abscab_main.pdf)

## Implementation

Various implementations are provided in this repository.
Here is an overview:

|   Language   |            main `abscab` file                                  |                unit tests                   | demo code | parallelized |
| ------------ | -------------------------------------------------------------- | ------------------------------------------- | --------- | ------------ |
| Java 8       | [`ABSCAB.java`](src/main/java/de/labathome/abscab/ABSCAB.java) | [`TestABSCAB.java`](src/test/java/de/labathome/abscab/TestABSCAB.java) | [`DemoABSCAB.java`](src/test/java/de/labathome/abscab/DemoABSCAB.java) | :heavy_check_mark: (threads) |
| C99          | [`abscab.h`](src/main/c/abscab.h)                              | [`test_abscab.c`](src/test/c/test_abscab.c) | [`demo_abscab.c`](src/test/c/demo_abscab.c) | :heavy_check_mark: (OpenMP) |
| Fortran 2008 | [`abscab.f08`](src/main/fortran/abscab.f08)                    | [`test_abscab.f08`](src/test/fortran/test_abscab.f08) | [`demo_abscab.f08`](src/test/fortran/demo_abscab.f08) | :heavy_check_mark: (OpenMP) |
| Python 3     | [`abscab.py`](src/main/python/abscab/abscab.py)                | test | demo  | :heavy_multiplication_x: |

## Reference Outputs

<img src="article/img/StraightWireSegment_results.png" alt="A_z and B_phi of Straight Wire Segment: Java vs. reference" height="500"/>

<img src="article/img/CircularWireLoop_A_phi_Java.png" alt="A_phi of Circular Wire Loop: Java vs. reference" width="500"/>
<img src="article/img/CircularWireLoop_B_rho_Java.png" alt="B_rho of Circular Wire Loop: Java vs. reference" width="500"/>
<img src="article/img/CircularWireLoop_B_z_Java.png" alt="B_z of Circular Wire Loop: Java vs. reference" width="500"/>

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
The code to generate the test points is [`src/test/java/GenerateTestKnots.java`](src/test/java/de/labathome/GenerateTestKnots.java).

The test points are saved into text files in [`src/test/resources`](src/test/resources):
* [`testPointsRpStraightWireSegment.dat`](src/test/resources/testPointsRpStraightWireSegment.dat) contains the value of `r'` at which
  the straight wire segment methods are tested.
* [`testPointsZpStraightWireSegment.dat`](src/test/resources/testPointsZpStraightWireSegment.dat) contains the value of `z'` at which
  the straight wire segment methods are tested.
* [`testPointsRpCircularWireLoop.dat`](src/test/resources/testPointsRpCircularWireLoop.dat) contains the value of `r'` at which
  the circular wire loop methods are tested.
* [`testPointsZpCircularWireLoop.dat`](src/test/resources/testPointsZpCircularWireLoop.dat) contains the value of `z'` at which
  the circular wire loop methods are tested.

Those above files are only provided to have a human-readable equivalent
of the set of test points.
The actual test point data read by the arbitrary-precision software
is in [`testPointsStraightWireSegment.dat`](src/test/resources/testPointsStraightWireSegment.dat) for the straight wire segment
and in [`testPointsCircularWireLoop.dat`](src/test/resources/testPointsCircularWireLoop.dat) for the circular wire loop.

### Reference Data
The reference data (`A_z` and `B_phi` for a straight wire segment;
`A_phi`, `B_rho` and `B_z` for a circular wire loop) is also available in the folder
[`src/test/resources`](src/test/resources):
* [`refDataStraightWireSegment.dat`](src/test/resources/refDataStraightWireSegment.dat) was generated by [`computeReferenceStraightWireSegment.py`](src/test/python/computeReferenceStraightWireSegment.py) and contains the quantities `A_z` and `B_phi`
  for a straight wire segment at all test points.
* [`refDataCircularWireLoop.dat`](src/test/resources/refDataCircularWireLoop.dat) was generated by [`computeReferenceCircularWireLoop.py`](src/test/python/computeReferenceCircularWireLoop.py) and contains the quantities `A_phi`, `B_rho` and `B_z`
  for a circular wire loop at all test points.
