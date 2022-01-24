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

#include "volFields.H"
#include "runTimeSelectionTables.H"
#include "fvMesh.H"

// include cfd2lcs functionalities
#include "cfd2lcs_inc_sp.h"




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

        //- Refererence to Db
        const objectRegistry& obr_;

        //- Number of grid points
        int n_[3];

        //- Resolution for LCS computation
        int res_;

        //- Grid coordinates and velocities
        lcsdata_t *x_, *y_, *z_, *u_, *v_, *w_;

        //- Boundary condition flag
        int *flag_;

        //- LCS diagnostics ids
        int id_fwd = -1;
        int id_bwd = -1;

        // Read from dictionary

            //- word
            word wordData_;

            //- scalar
            scalar scalarData_;

            //- label
            label labelData_;


    // Private Member Functions
        void getCellCenterCoords();
        void getNumberOfCellsInDirection();
        void initializeLCSDiagnostics();
        void setLCSoptions();


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

        //- Execute, currently does nothing
        virtual void execute();

        //- Execute at the final time-loop, currently does nothing
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
