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

\*---------------------------------------------------------------------------*/

#include "LCS.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(LCS, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LCS::LCS
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    name_(name),
    obr_(obr),
    wordData_(dict.lookupOrDefault<word>("wordData", "defaultWord")),
    scalarData_(readScalar(dict.lookup("scalarData"))),
    labelData_(readLabel(dict.lookup("labelData")))
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::LCS::~LCS()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::LCS::start()
{
    // mpi communicator
    MPI_Comm comm = MPI_COMM_WORLD;

    // number of grid points for THIS partition in x=i, y=j and z=k direction
    getNumberOfCellsInDirection();

    // Global offset for these grid points
    int offset[3]= {0,0,0};

    // Allocate space for data
    x_ = (lcsdata_t*)malloc(n_[0]*n_[1]*n_[2]*sizeof(lcsdata_t));
    y_ = (lcsdata_t*)malloc(n_[0]*n_[1]*n_[2]*sizeof(lcsdata_t));
    z_ = (lcsdata_t*)malloc(n_[0]*n_[1]*n_[2]*sizeof(lcsdata_t));
    u_ = (lcsdata_t*)malloc(n_[0]*n_[1]*n_[2]*sizeof(lcsdata_t));
    v_ = (lcsdata_t*)malloc(n_[0]*n_[1]*n_[2]*sizeof(lcsdata_t));
    w_ = (lcsdata_t*)malloc(n_[0]*n_[1]*n_[2]*sizeof(lcsdata_t));
    flag_ = (int*)malloc(n_[0]*n_[1]*n_[2]*sizeof(int));

    // convert OpenFoam mesh to x,y,z arays
    getCellCenterCoords();

    // Initializes the communications and data storage for OpenFoam data input
    cfd2lcs_init_c(comm,n_,offset,x_,y_,z_,flag_);

    // Initialize LCS diagnostics
    initializeLCSDiagnostics();

    // Set CFD2LCS options/parameters
    setLCSoptions();

    return true;
}

void Foam::LCS::read(const dictionary& dict)
{
    dict.readIfPresent("wordData", wordData_);
    dict.lookup("scalarData") >> scalarData_;
    dict.lookup("labelData") >> labelData_;
}


void Foam::LCS::execute()
{
    // Do nothing - only valid on write
}


void Foam::LCS::end()
{
    // Do nothing - only valid on write
}


void Foam::LCS::timeSet()
{
    // Do nothing - only valid on write
}


void Foam::LCS::write()
{
}


void Foam::LCS::getCellCenterCoords()
{
    //TODO
    x_[0]=1.f;
    y_[0]=1.f;
    z_[0]=1.f;
}

void Foam::LCS::getNumberOfCellsInDirection()
{
    //TODO
    n_[0]=1;
    n_[1]=1;
    n_[2]=1;
}

void Foam::LCS::initializeLCSDiagnostics()
{
    //TODO
    lcsdata_t T = 15.0;
    lcsdata_t H = 1.5;
    char labelfwd[LCS_NAMELEN]="fwdFTLE";
    id_fwd = cfd2lcs_diagnostic_init_c(FTLE_FWD,res_,T,H,labelfwd);
}

void Foam::LCS::setLCSoptions()
{
    //TODO
    cfd2lcs_set_option_c("SYNCTIMER",LCS_FALSE);
    cfd2lcs_set_option_c("DEBUG",LCS_FALSE);
    cfd2lcs_set_option_c("WRITE_FLOWMAP",LCS_FALSE);
    cfd2lcs_set_option_c("WRITE_BCFLAG",LCS_FALSE);
    cfd2lcs_set_option_c("INCOMPRESSIBLE",LCS_FALSE);
    cfd2lcs_set_option_c("AUX_GRID",LCS_FALSE);
    cfd2lcs_set_option_c("INTEGRATOR",RK3);
    cfd2lcs_set_option_c("INTERPOLATOR",LINEAR);
    lcsdata_t CFL = 0.9;
    cfd2lcs_set_param_c("CFL", CFL);
}

// ************************************************************************* //
