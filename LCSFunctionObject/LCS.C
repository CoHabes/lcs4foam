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
    obr_(obr)
{
    // Check if the available mesh is an fvMesh, otherwise deactivate
    if (!isA<fvMesh>(obr_))
    {
        active_ = false;
        WarningIn
        (
            "LCS::LCS"
            "("
                "const word&, "
                "const objectRegistry&, "
                "const dictionary&, "
                "const bool"
            ")"
        )   << "No fvMesh available, deactivating " << name_ << nl
            << endl;
    }

    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::LCS::~LCS()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::LCS::start()
{
    if (active_)
    {
        // mpi communicator
        MPI_Comm comm = MPI_COMM_WORLD;

        // number of grid points for THIS partition in x=i, y=j and z=k direction
        // !! for now user has to provide information in dict !!
        // getNumberOfCellsInDirection();

        // Global offset for these grid points
        // parallel computation not supported yet
        int offset[3]= {0,0,0};

        // Allocate space for data used by lcs library
        x_ = static_cast<lcsdata_t*>(malloc(n_[0]*n_[1]*n_[2]*sizeof(lcsdata_t)));
        y_ = static_cast<lcsdata_t*>(malloc(n_[0]*n_[1]*n_[2]*sizeof(lcsdata_t)));
        z_ = static_cast<lcsdata_t*>(malloc(n_[0]*n_[1]*n_[2]*sizeof(lcsdata_t)));
        u_ = static_cast<lcsdata_t*>(malloc(n_[0]*n_[1]*n_[2]*sizeof(lcsdata_t)));
        v_ = static_cast<lcsdata_t*>(malloc(n_[0]*n_[1]*n_[2]*sizeof(lcsdata_t)));
        w_ = static_cast<lcsdata_t*>(malloc(n_[0]*n_[1]*n_[2]*sizeof(lcsdata_t)));
        flag_ = static_cast<int*>(malloc(n_[0]*n_[1]*n_[2]*sizeof(int)));

        // convert OpenFoam mesh cell centers to x,y,z arays 
        // and set boundary type flags for each cell center
        getCellCenterCoords();

        // Initializes the communications and data storage for OpenFoam data input
        cfd2lcs_init_c(comm,&n_[0],offset,x_,y_,z_,flag_);

        // Initialize LCS diagnostics
        initializeLCSDiagnostics();

        // Set CFD2LCS options/parameters
        setLCSoptions();
    }

    return true;
}

void Foam::LCS::read(const dictionary& dict)
{
    if (active_)
    {
        n_ = dict.lookup("n");
        uName_ = dict.lookupOrDefault<word>("velocityName", "U");
    }
}


void Foam::LCS::execute()
{   
    if (active_)
    {
        // Do nothing - only valid on write
    }
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
    const fvMesh& mesh = refCast<const fvMesh>(obr_);

    // Cell centroid coordinates
    const vectorField& centres = mesh.C().internalField();
    // Loop over cell centres
    forAll (centres , celli)
    {
        x_[celli] = mesh.C()[celli].component(0);
        y_[celli] = mesh.C()[celli].component(1);
        z_[celli] = mesh.C()[celli].component(2);
        // set all cells als internal
        flag_[celli] = LCS_INTERNAL;
    }

    // Loop over all boundary patches
    forAll (mesh.boundaryMesh(), patchi)
    {
        // Current poly patch
        const polyPatch& pp = mesh.boundaryMesh()[patchi];
        const word& patchName = pp.name();

        // Velocity field
        const volVectorField& U = mesh.lookupObject<volVectorField>(uName_);
        const fvPatchVectorField& U_p = U.boundaryField()[patchi];


        // determine boundary type
        label boundaryType = LCS_INTERNAL;
        if (pp.type() == "wall")
        {
            if (U_p.type() == "slip")
            {
                boundaryType = LCS_SLIP;
            }else
            {
                boundaryType = LCS_WALL;
            }
        }else if (pp.type() == "empty" || pp.type() == "symmetryPlane" || pp.type() == "wedge")
        {
            boundaryType = LCS_SLIP;
        }else if (pp.type() == "cyclic" || pp.type() == "processor")
        {
            boundaryType = LCS_INTERNAL;
        }
        else
        {
            // TODO: Use U boundary Field to determine boundary field type or use regex for patchnames
            if (patchName == "inlet" || patchName == "Inlet")
            {
                boundaryType = LCS_INFLOW;
            }else if (patchName == "outlet" || patchName == "Outlet")
            {
                boundaryType = LCS_OUTFLOW;
            }
        }
        
        
        // Loop over all faces of boundary patch
        forAll(mesh.boundary()[patchi], facei)
        {
            // Boundary cell index
            const label& bCell = mesh.boundary()[patchi].faceCells()[facei];
            flag_[bCell] = boundaryType;
        }
    }
}

void Foam::LCS::getNumberOfCellsInDirection()
{
    // For now user has to provide information in dict
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
    char option1[] = "SYNCTIMER";
    cfd2lcs_set_option_c(option1,LCS_FALSE);

    char option2[] = "DEBUG";
    cfd2lcs_set_option_c(option2,LCS_FALSE);

    char option3[] = "WRITE_FLOWMAP";
    cfd2lcs_set_option_c(option3,LCS_FALSE);
 
    char option4[] = "WRITE_BCFLAG";
    cfd2lcs_set_option_c(option4,LCS_FALSE);

    char option5[] = "INCOMPRESSIBLE";
    cfd2lcs_set_option_c(option5,LCS_FALSE);

    char option6[] = "AUX_GRID";
    cfd2lcs_set_option_c(option6,LCS_FALSE);

    char option7[] = "INTEGRATOR";
    cfd2lcs_set_option_c(option7,RK3);

    char option8[] = "INTERPOLATOR";
    cfd2lcs_set_option_c(option8,LINEAR);

    char option9[] = "CFL";
    lcsdata_t CFL = 0.9;
    cfd2lcs_set_param_c(option9, CFL);
}

// ************************************************************************* //
