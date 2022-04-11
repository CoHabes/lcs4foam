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

Class
    Foam::LCS

Group

Description
    This function object...

    Example of function object specification:
    \verbatim
    LCS1
    {
        type           LCS;
        functionObjectLibs ("libLCSFunctionObject.so");
        ...
        wordData       someWord;
        scalarData     1.0;
        labelData      1;
    }
    \endverbatim

    \heading Function object usage
    \table
        Property     | Description               | Required | Default value
        type         | type name: LCS | yes      |
        wordData     | some word option...       | no       | defaultWord
        scalarData   | some scalar value...      | yes      |
        labelData    | some label value...       | yes      |
    \endtable

SourceFiles
    LCS.C
    IOLCS.H

\*---------------------------------------------------------------------------*/

#ifndef LCS_H
#define LCS_H

#include "autoPtr.H"
#include "volFields.H"
#include "runTimeSelectionTables.H"
#include "Switch.H"
#include "fvMesh.H"
#include "meshToMesh.H"
#include "fvMeshSubset.H"
#include "cellSet.H"


// include cfd2lcs functionalities
extern "C" {
#include "cfd2lcs_inc_dp.h"
};




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class polyMesh;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                   Class LCS Declaration
\*---------------------------------------------------------------------------*/

class LCS
{
    // Private data

        //- Name of this set of LCS
        word name_;

        //- Refererence to Db/mesh
        const fvMesh& cfdMesh_;

        //- Pointer to const mesh which is used for lcs computations
        autoPtr<const fvMesh> lcsMeshPtr_;

        //- Pointer to subMesh of oversetMesh which is used for lcs computations
        autoPtr<fvMeshSubset> lcsMeshSubset_;

        //- Pointer to mesh to mesh interpolator
        autoPtr<meshToMesh> meshToMeshPtr_;

        //- Pointer to lcs deignated lcs velocity if needed 
        autoPtr<volVectorField> lcsVelFieldPtr_;

        //- Refererence to controlDict
        const dictionary& controlDict_;

        //- global bounding box
        boundBox globalBb_;

        //- local bounding box
        boundBox localBb_;

        //- On/off switch
        bool active_;

        //- Indicator for first execution of function object
        // Workaround because start() somehow does not get called
        bool firstExe_;

        //- Number of grid points in local mesh
        labelList n_;

        //- Local grid offset in respect to global mesh
        labelList offset_;

        //- Grid coordinates and velocities
        lcsdata_t *x_, *y_, *z_, *u_, *v_, *w_;

        //- Boundary condition flag
        label *flag_;

        //- LCS diagnostics ids
        label id_fwd_;
        label id_bkwd_;

        // Read from dictionary
            
            //- If true use specified cellZone of oversetMesh for lcs computation
            Switch isOverset_;

            //- Name of cellZone in oversetMesh which should be used for the lcs computation
            word lcsOversetRegion_;

            //- If true use cfd mesh for lcs computation, if false use extra lcs mesh
            Switch isStaticRectLinear_;

            //- Number of grid points in global mesh
            labelList globalN_;

            //- Name of velocity field, default is "U"
            word uName_;

            //- Resolution for LCS computation
            label res_;

            //- Switches for forward time and backward time diagnostics
            Switch ftleFwd_;
            Switch ftleBkwd_;

            //- LCS start time
            scalar t0_;
            //- LCS diagnostic run time
            scalar T_;

            //- LCS output time interval 
            scalar H_;

            //- LCS options
            struct lcsOpts {
                label synctimer;
                label debug;
                label writeFlowmap;
                label writeBCFalgs;
                label incompressible;
                label auxGrid;
                //- LCS time integration scheme
                label integrator;
                //- LCS interpolation scheme
                label interpolator;
                //- LCS CFL number 
                scalar CFL;

                // Member function to set options
                void setLCSoptions();
            } lcsOpts_;

    // Private Member Functions
        void getBoundBoxes();
        void getCellCenterCoords();
        void getVelocityField();
        void getNumberOfCellsInDirection();
        void getOffset();
        void initializeLCSDiagnostics();
        void createLCSMesh();
        void createLCSSubsetMesh();
        
        //- demand driven mesh to mesh interpolator
        const meshToMesh& meshToMeshInterp();

        //- demand driven volVectorField to designated lcs velocity field
        volVectorField& lcsVelField();


        //- Disallow default bitwise copy construct
        LCS(const LCS&);

        //- Disallow default bitwise assignment
        void operator=(const LCS&);


public:

    //- Runtime type information
    TypeName("LCS");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        LCS
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~LCS();


    // Member Functions

        //- Return name of the LCS
        virtual const word& name() const
        {
            return name_;
        }

        //- start is called at the start of the time-loop
        virtual bool start();

        //- Read the LCS data
        virtual void read(const dictionary&);

        //- Execute
        virtual void execute();

        //- Execute at the final time-loop
        virtual void end();

        //- Called when time was set at the end of the Time::operator++
        virtual void timeSet();

        //- Write the LCS
        virtual void write();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&)
        {}

        //- Update for changes of mesh
        virtual void movePoints(const pointField&)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
