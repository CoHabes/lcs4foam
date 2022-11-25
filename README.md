#  lcs4Foam
## In this repository
 - `LCSFunctionObject` :        Implementation of the function object
 - `LCSFunctionObject/libcfd2lcs`: Slightly modified thrid party library of [libcfd2lcs](https://github.com/justin-finn/libcfd2lcs)
 - `LCS_Testcases` : Testcases that show the application of the LCS function object
 - `Allwmake` : Make script for the compilation of the third party library and the function object
 - `Allwclean`: A script to undo the compilation

## Installation
### Prerequisites:

 1. foam-extend 4.1
 2. liblapack
 
### Compilation:
With a working foam-extend 4.1 installation, clone the repository, source your foam environment and build the library:
`./Allwmake`

---
**NOTE**

The default path for the liblapack installation is `/usr/lib/x86_64-linux-gnu/lapack`. If the instllation path differs from this one, the entry for `LAPACK_LIBS` in `LCSFunctionObject/libcfd2lcs/makefiles/Makefile.FOR_OPENFOAM.in` needs to be changed accordingly. 

---
## Using the lcs4foam function object
This function object can be used as every other OpenFOAM function object. For the necessary settings in the controlDict see `LCSFunctionObject/controlDict` or the  tescases in `LCS_Testcases`. 
Depending on the used simulation mesh it might be necessary to use additional meshes for the LCS computations. Therefore, see the testcases
- `LCS_Testcases/LCS_Testcase_cylinder`
- `LCS_Testcases/LCS_Testcase_cylinder_smallLCSMesh`
- `LCS_Testcases/LCS_Testcase_oversetCylinderThreeLevels`