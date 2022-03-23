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
    doubleGyre

Description
    Generates an analytical solution for double gyre flow

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#define DG_A 0.1  //Amplitude
#define DG_EPS 0.1
#define DG_OMG 0.62831853071 //Freq = 2*pi/10

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

        scalar aoft = DG_EPS*Foam::sin(DG_OMG*time);
        scalar boft = 1.0 -2.0*DG_EPS*Foam::sin(DG_OMG*time);

        volScalarField fofxt = aoft*pow(x,2.0) + boft*x;
        volScalarField dfdx = 2.0*aoft*x + boft;

        U = (mDBs*-M_PI*DG_A*sin((M_PI*fofxt))*cos((M_PI*y)))    * dimensionedVector(vector(1,0,0)) 
          + (mDBs*M_PI*DG_A*cos((M_PI*fofxt))*sin((M_PI*y))*dfdx)* dimensionedVector(vector(0,1,0)) 
          + (mDBs*0.0)                                           * dimensionedVector(vector(0,0,1));


        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                << nl << endl;

    }

    Info<< "end" << endl;

    return 0;
}

// ************************************************************************* //
