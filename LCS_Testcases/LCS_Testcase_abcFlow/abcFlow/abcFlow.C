/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Application
    abcFlow

Description
    Generates an analytical solution for Arnold–Beltrami–Childress flow

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#define ABC_A 0.5  //Foam::sqrt(3.0)
#define ABC_B 0.8  //Foam::sqrt(2.0)
#define ABC_C 0.8  //1.0
#define ABC_D 0.0

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"

#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    dimensionedScalar oneDBm = dimensionedScalar("oneDBm", dimensionSet(0,-1,0,0,0),1.0);
    volScalarField x = U.mesh().C().component(0)*oneDBm;
    volScalarField y = U.mesh().C().component(1)*oneDBm;
    volScalarField z = U.mesh().C().component(2)*oneDBm;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        scalar time = runTime.value();

        dimensionedScalar mDBs = dimensionedScalar("mDBs", dimensionSet(0,1,-1,0,0) ,1.0);

        U = (mDBs*((ABC_A + ABC_D*Foam::sin(time))*sin(z) + ABC_C*cos(y)))  * dimensionedVector(vector(1,0,0)) 
          + (mDBs*(ABC_B*sin(x) + (ABC_A+ABC_D*Foam::sin(time))*cos(z)))    * dimensionedVector(vector(0,1,0)) 
          + (mDBs*(ABC_C*sin(y) + ABC_B*cos(x)))                            * dimensionedVector(vector(0,0,1));


        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                << nl << endl;

    }

    Info<< "end" << endl;

    return 0;
}

// ************************************************************************* //
