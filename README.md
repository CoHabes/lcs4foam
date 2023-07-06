#  LCS4FOAM

[![DOI](https://zenodo.org/badge/662623039.svg)](https://zenodo.org/badge/latestdoi/662623039)


Developed by the Computational Multiphase Flow research group.

* [Research Group Website](https://www.mathematik.tu-darmstadt.de/cmf/)
* [Report Bug](https://bitbucket.org/C_Habes/lcs4foam/issues?status=new&status=open)
* [Request Feature](https://bitbucket.org/C_Habes/lcs4foam/issues?status=new&status=open)

## About The Project
An OpenFOAM function object to compute finite-time Lyapunov exponent (FTLE) fields during CFD simulation.
This enables the OpenFOAM community to assess the geometry of the material transport in any flow quantitatively on-the-fly using principally any OpenFOAM flow solver.

## In this repository
 - `LCSFunctionObject` :        Implementation of the function object
 - `LCSFunctionObject/libcfd2lcs`: Slightly modified thrid party library of [libcfd2lcs](https://github.com/justin-finn/libcfd2lcs)
 - `LCS_Testcases` : Testcases that show the application of the LCS function object
 - `Allwmake` : Make script for the compilation of the third party library and the function object
 - `Allwclean`: A script to undo the compilation

## Installation
### Prerequisites:

1. Working installation of foam-extend 4.1: For installing foam-extend 4.1 please refer to and follow the [installation instructions](https://openfoamwiki.net/index.php/Installation/Linux/foam-extend-4.1) step-by-step!
    * [Installation/Linux/foam-extend-4.1/Ubuntu](https://openfoamwiki.net/index.php/Installation/Linux/foam-extend-4.1/Ubuntu)
    * [Installation/Linux/foam-extend-4.1/CentOS](https://openfoamwiki.net/index.php/Installation/Linux/foam-extend-4.1/CentOS)
2. liblapack
    * Installable with the package manager of your choice e.g. `sudo apt-get install liblapack-dev`
 
### Compilation:
With a working foam-extend 4.1 installation, clone the repository, source your foam environment and build the library:

```bash
./Allwmake
```

> **_NOTE:_**  The default path for the liblapack installation is `/usr/lib/x86_64-linux-gnu/lapack`. If the instllation path differs from this one, the entry for `LAPACK_LIBS` in `LCSFunctionObject/libcfd2lcs/makefiles/Makefile.FOR_OPENFOAM.in` needs to be changed accordingly.

## Using the lcs4foam function object
This function object can be used as every other OpenFOAM function object. For the necessary settings in the controlDict see `LCSFunctionObject/controlDict` or the  tescases in `LCS_Testcases`.
Depending on the used simulation mesh it might be necessary to use additional meshes for the LCS computations. Therefore, see the testcases

 - `LCS_Testcases/LCS_Testcase_cylinder`
 - `LCS_Testcases/LCS_Testcase_cylinder_smallLCSMesh`
 - `LCS_Testcases/LCS_Testcase_oversetCylinderThreeLevels`

## Contributing

We invite everyone in the FOAM community to collaborate with us and further 
develop LCS4FOAM. For this, start by forking the repository to your own 
account. Then, clone the forked repository to your local machine and create a 
new branch for your changes. Make the necessary modifications, commit your 
changes, and push the branch to your forked repository. Finally, open a pull 
request from your branch to the original repository and provide a clear 
description of your changes. Collaborate with reviewers, address feedback, 
and once approved, your contributions can be merged.

The `master` branch is the corner stone of the development, please branch all of 
your feature/bugFix branches off of it, and "rebase" your branches on it before 
issuing a pull request. When your branch gets merged, it's considered a 
"best-practice" to delete your feature branch and start a fresh one.


## License

Released under the GNU Public License - see code headers for details.


## Contact

[Email @habes](mailto:constantin.habes@tu-darmstadt.de)