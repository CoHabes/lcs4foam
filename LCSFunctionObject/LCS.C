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
    mesh_(refCast<const fvMesh>(obr)),
    controlDict_(obr.time().controlDict()),
    globalBb_(refCast<const fvMesh>(obr).bounds()),
    localBb_(refCast<const fvMesh>(obr).points(), false)
{
    // Check if the available mesh is an fvMesh, otherwise deactivate
    if (!isA<fvMesh>(mesh_))
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
    Info<< "Start initializing LCS diagnostics" << endl;
    if (active_)
    {
        // Mpi communicator
        MPI_Comm comm = MPI_COMM_WORLD;

        // Compute number of grid points for THIS partition in x=i, y=j and z=k direction
        getNumberOfCellsInDirection();

        // Compute offset of local mesh in respect to the global mesh 
        getOffset();

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
        cfd2lcs_init_c(comm,n_.data(),offset_.data(),x_,y_,z_,flag_);

        // Initialize LCS diagnostics
        initializeLCSDiagnostics();

        // Set CFD2LCS options/parameters
        setLCSoptions();

        // start() has been called
        firstExe_ = false;
    }

    Info<< "Finished initializing LCS diagnostics" << endl;
    return true;
}

void Foam::LCS::read(const dictionary& dict)
{   
    Info<< "Start Reading LCS diagnostic settings" << endl;
    if (active_)
    {
        globalN_ = dict.lookup("n");
        uName_ = dict.lookupOrDefault<word>("velocityName", "U");
        ftleFwd_ = dict.lookupOrDefault<Switch>("ftleFwd", true);
        ftleBkwd_ = dict.lookupOrDefault<Switch>("ftleBkwd", false);

        // check if any of the diagnostics should be executed
        if(!ftleFwd_ & !ftleBkwd_){
        active_ = false;
        WarningIn
        (
            "LCS::read"
            "("
                "const dictionary& dict"
            ")"
        )   << "No FTLE evaluation activated, deactivating " << name_ << nl
            << endl;
        }

        res_ = dict.lookupOrDefault<label>("resolution", 0);
        scalar simulationStartTime = readScalar(controlDict_.lookup("startTime"));
        t0_ = dict.lookupOrDefault<scalar>("lcsStartTime", simulationStartTime);
        scalar simulationEndTime = readScalar(controlDict_.lookup("endTime"));
        scalar lcsEndTime = dict.lookupOrDefault<scalar>("lcsEndTime", simulationEndTime);
        T_ = lcsEndTime - t0_;

        // determine default value for LCS output time interval 
        scalar simulationWriteTimeInterval = 1.0f;
        word writeControl = controlDict_.lookup("writeControl");
        if(writeControl == "runTime" || writeControl == "adjustableRunTime" ){
            simulationWriteTimeInterval = readScalar(controlDict_.lookup("writeInterval"));
        }else if(writeControl == "timeStep"){
            Switch adjustTimeStep = controlDict_.lookupOrDefault<Switch>("adjustTimeStep", false);
            if(adjustTimeStep){
                WarningIn
                (
                    "LCS::read"
                    "("
                        "const dictionary& dict"
                    ")"
                )   << "Adjustable simulation time steps with timeStep writeControl is not support with LCS diagnostics." << nl 
                    << "If no lcsWriteInterval is provided lcsWriteTimeInterval is set to 1s"<< name_ << nl
                    << endl;
            }else{
                scalar simulationWriteStepInterval = readScalar(controlDict_.lookup("writeInterval"));
                scalar simulationDeltaT = readScalar(controlDict_.lookup("deltaT"));
                simulationWriteTimeInterval = simulationDeltaT * simulationWriteStepInterval;
            }
        }else{
            WarningIn
            (
                "LCS::read"
                "("
                    "const dictionary& dict"
                ")"
            )   << "cpuTime and clockTime writeControl is not support with LCS diagnostics." << nl 
                << "If no lcsWriteInterval is provided lcsWriteTimeInterval is set to 1s"<< name_ << nl
                << endl;
        }
        H_ = dict.lookupOrDefault<scalar>("lcsWriteTimeIntervall", simulationWriteTimeInterval);

        // read integrator - Maybe better with selectionTables, but this works for now
        std::map<std::string, int> lcsIntegratorMap { 
            {"euler", EULER}, 
            {"trapezodial", TRAPEZOIDAL}, 
            {"rk2", RK2}, 
            {"rk3", RK3}, 
            {"rk4", RK4} 
        };
        word integrator = dict.lookup("integrator", "rk2");
        auto searchIntegrator = lcsIntegratorMap.find(integrator);
        if (searchIntegrator != lcsIntegratorMap.end()) {
            lcsIntegrator_ = searchIntegrator->second;
        } else {
            WarningIn
            (
                "LCS::read"
                "("
                    "const dictionary& dict"
                ")"
            )   << "LCS integrator " << integrator << "is no valid integration scheme for lcs diagnostic" << nl 
                << "valid types are: euler, trapezodial, rk2, rk3, rk4. " << nl 
                << "Setting LCS integrator to rk2."<< name_ << nl
                << endl;
            lcsIntegrator_ = RK2;
        }


        // read interpolator - Maybe better with selectionTables, but this works for now
        std::map<std::string, int> lcsInterpolatorMap { 
            {"nearestNBR", NEAREST_NBR}, 
            {"linear", LINEAR}, 
            {"quadratic", QUADRATIC}, 
            {"cubic", CUBIC}, 
            {"tse", TSE}, 
            {"tseLimit", TSE_LIMIT}
            };
        word interpolator = dict.lookupOrDefault<word>("interpolator", "linear");
        auto searchInterpolator = lcsInterpolatorMap.find(interpolator);
        if (searchInterpolator != lcsInterpolatorMap.end()) {
            lcsInterpolator_ = searchInterpolator->second;
        } else {
            WarningIn
            (
                "LCS::read"
                "("
                    "const dictionary& dict"
                ")"
            )   << "LCS interpolator " << interpolator << "is no valid interpolation scheme for lcs diagnostic." << nl 
                << "Valid types are: nearestNBR, linear, quadratic, cubic, tse, tseLimit. " << nl 
                << "Setting LCS interpolator to linear."<< name_ << nl
                << endl;
            lcsInterpolator_ = LINEAR;
        }

        // read lcs CFl number
        lcsCFL_ = dict.lookupOrDefault<scalar>("lcsCFL", 0.5f);
        
        // read lcs options
        const dictionary& lcsOptionsDict = dict.subDict("lcsOptions");
        Switch synctimer = lcsOptionsDict.lookupOrDefault<Switch>("synctimer", false);
        optSynctimer_ = (synctimer) ? LCS_TRUE : LCS_FALSE;
        Switch debug = lcsOptionsDict.lookupOrDefault<Switch>("debug", false);
        optDebug_ = (debug) ? LCS_TRUE : LCS_FALSE;
        Switch writeFlowmap = lcsOptionsDict.lookupOrDefault<Switch>("writeFlowmap", false);
        optWriteFlowmap_ = (writeFlowmap) ? LCS_TRUE : LCS_FALSE;
        Switch writeBCFalgs = lcsOptionsDict.lookupOrDefault<Switch>("writeBCFlags", false);
        optWriteBCFalgs_ = (writeBCFalgs) ? LCS_TRUE : LCS_FALSE;
        Switch incompressible = lcsOptionsDict.lookupOrDefault<Switch>("incompressible", false);
        optIncompressible_ = (incompressible) ? LCS_TRUE : LCS_FALSE;
        Switch auxGrid = lcsOptionsDict.lookupOrDefault<Switch>("auxGrid", false);
        optAuxGrid_ = (auxGrid) ? LCS_TRUE : LCS_FALSE;

        Info<< "Finished reading LCS diagnostic settings" << endl;
    }
}

void Foam::LCS::execute()
{   
    if(firstExe_){
        start();
    }

    Info<< "Start executing LCS diagnostics" << endl;
    if (active_)
    {   
        getVelocityField();
        scalar time = mesh_.time().value();
        Info<< "Start updating lcs diagnostics" << endl;
        cfd2lcs_update_c(n_.data(), u_, v_, w_ ,time);
        Info<< "Finished updating lcs diagnostics" << endl;
    }
    Info<< "Finished executing LCS diagnostics" << endl;
}

void Foam::LCS::end()
{
    if (id_fwd_ != -1)
    {
        cfd2lcs_diagnostic_destroy_c(id_fwd_);
    }
    
    if (id_bkwd_ != -1)
    {
        cfd2lcs_diagnostic_destroy_c(id_bkwd_);
    }
}

void Foam::LCS::timeSet()
{
    
}

void Foam::LCS::write()
{
    // 
}

void Foam::LCS::getCellCenterCoords()
{
    Info<< "Start getting cell center coords for LCS diagnostic" << endl;
    // Cell centroid coordinates
    const vectorField& centres = mesh_.C().internalField();
    // Loop over cell centres
    forAll (centres , celli)
    {
        x_[celli] = mesh_.C()[celli].component(0);
        y_[celli] = mesh_.C()[celli].component(1);
        z_[celli] = mesh_.C()[celli].component(2);
        // set all cells als internal
        flag_[celli] = LCS_INTERNAL;
    }

    // Loop over all boundary patches
    forAll (mesh_.boundaryMesh(), patchi)
    {
        // Current poly patch
        const polyPatch& pp = mesh_.boundaryMesh()[patchi];
        const word& patchName = pp.name();

        // Velocity field
        const volVectorField& U = mesh_.lookupObject<volVectorField>(uName_);
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
            boundaryType = LCS_INTERNAL;
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
        forAll(mesh_.boundary()[patchi], facei)
        {
            // Boundary cell index
            const label& bCell = mesh_.boundary()[patchi].faceCells()[facei];
            flag_[bCell] = boundaryType;
        }
    }

    Info<< "Finished getting cell center coords for LCS diagnostic" << endl;
}

void Foam::LCS::getVelocityField()
{
    Info<< "Start getting velocity field for LCS diagnostic" << endl;
    if (mesh_.foundObject<volVectorField>(uName_))
    {
        const volVectorField& u = mesh_.lookupObject<volVectorField>(uName_);
        const vectorField& uIn = u.internalField();

        if (!uIn.empty())
        {
            forAll (uIn, celli)
            {
                u_[celli] = uIn[celli].component(0);
                v_[celli] = uIn[celli].component(1);
                w_[celli] = uIn[celli].component(2);
            }
        }
    }else{
        FatalErrorIn
        (
            "Foam::LCS::getVelocityField()"
        )   << "Velocity field with name " << uName_ << " not found"
            << exit(FatalError);
    }
    Info<< "Finished getting velocity field for LCS diagnostic" << endl;
}

void Foam::LCS::getNumberOfCellsInDirection()
{
    n_.setSize(3);
    // Assuming that LCS mesh has constant cell size along each axis
    n_[0] = round(globalN_[0] * (localBb_.max().component(0) - localBb_.min().component(0)) / (globalBb_.max().component(0) - globalBb_.min().component(0)));
    n_[1] = round(globalN_[1] * (localBb_.max().component(1) - localBb_.min().component(1)) / (globalBb_.max().component(1) - globalBb_.min().component(1)));
    n_[2] = round(globalN_[2] * (localBb_.max().component(2) - localBb_.min().component(2)) / (globalBb_.max().component(2) - globalBb_.min().component(2)));

    Pout<< "local bounding box:" << localBb_ << endl;
    Pout<< "Number of cells in x:" << n_[0] << " y:" << n_[1] << " z:" << n_[2]  << endl;

}

void Foam::LCS::getOffset()
{
    offset_.setSize(3);
    // Assuming that LCS mesh has constant cell size along each axis
    offset_[0] = round(n_[0] * (localBb_.min().component(0) - globalBb_.min().component(0)) / (localBb_.max().component(0) - localBb_.min().component(0)));
    offset_[1] = round(n_[1] * (localBb_.min().component(1) - globalBb_.min().component(1)) / (localBb_.max().component(1) - localBb_.min().component(1)));
    offset_[2] = round(n_[2] * (localBb_.min().component(2) - globalBb_.min().component(2)) / (localBb_.max().component(2) - localBb_.min().component(2)));

    Pout<< "Offset in x:" << offset_[0] << " y:" << offset_[1] << " z:" << offset_[2]  << endl;

}

void Foam::LCS::initializeLCSDiagnostics()
{
    Info<< "Start initializing LCS diagnostic" << endl;
    if(ftleFwd_){
        char labelfwd[LCS_NAMELEN]="fwdFTLE";
        id_fwd_ = cfd2lcs_diagnostic_init_c(FTLE_FWD,res_,T_,H_,labelfwd);
    }

    if(ftleBkwd_){
        char labelbkwd[LCS_NAMELEN]="bkwdFTLE";
        id_bkwd_ = cfd2lcs_diagnostic_init_c(FTLE_BKWD, res_, T_, H_, labelbkwd);
    } 
    Info<< "Finished initializing LCS diagnostic" << endl;
}

void Foam::LCS::setLCSoptions()
{
    Info<< "Start setting LCS diagnostic options" << endl;
    char option1[] = "SYNCTIMER";
    cfd2lcs_set_option_c(option1,optSynctimer_);

    char option2[] = "DEBUG";
    cfd2lcs_set_option_c(option2,optDebug_);

    char option3[] = "WRITE_FLOWMAP";
    cfd2lcs_set_option_c(option3,optWriteFlowmap_);
 
    char option4[] = "WRITE_BCFLAG";
    cfd2lcs_set_option_c(option4,optWriteBCFalgs_);

    char option5[] = "INCOMPRESSIBLE";
    cfd2lcs_set_option_c(option5,optIncompressible_);

    char option6[] = "AUX_GRID";
    cfd2lcs_set_option_c(option6,optAuxGrid_);

    char option7[] = "INTEGRATOR";
    cfd2lcs_set_option_c(option7,lcsIntegrator_);

    char option8[] = "INTERPOLATOR";
    cfd2lcs_set_option_c(option8,lcsInterpolator_);

    char option9[] = "CFL";
    cfd2lcs_set_param_c(option9, lcsCFL_);
    Info<< "Finished setting LCS diagnostic options" << endl;
}

// ************************************************************************* //
