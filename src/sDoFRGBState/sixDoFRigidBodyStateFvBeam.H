/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenFOAM Foundation
    Copyright (C) 2018-2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::functionObjects::sixDoFRigidBodyStateFvBeam

Group
    grpSixDoFRigidBodyFunctionObjects

Description
    Writes the 6-DoF motion state.

    Example of function object specification:
    \verbatim
    sixDoFRigidBodyStateFvBeam
    {
        type           sixDoFRigidBodyStateFvBeam;
        libs           (sixDoFRigidBodyStateFvBeam);
        angleFormat    degrees;
    }
    \endverbatim

Usage
    \table
        Property     | Description                  | Required | Default value
        type         | type name: sixDoFRigidBodyStateFvBeam    | yes |
        angleFormat  | degrees or radians           | no       | radian
    \endtable

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile

SourceFiles
    sixDoFRigidBodyStateFvBeam.C

\*---------------------------------------------------------------------------*/

#ifndef sixDoFRigidBodyStateFvBeam_H
#define sixDoFRigidBodyStateFvBeam_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                    Class sixDoFRigidBodyStateFvBeam Declaration
\*---------------------------------------------------------------------------*/

class sixDoFRigidBodyStateFvBeam
:
    public fvMeshFunctionObject,
    public writeFile
{
    // Private data

        enum class angleTypes
        {
            RADIANS,        //!< Radians
            DEGREES         //!< Degrees
        };

        //- Angle type names
        static const Enum<angleTypes> angleTypeNames_;

        //- Angle format
        angleTypes angleFormat_;


    // Private Member Functions

        //- No copy construct
        sixDoFRigidBodyStateFvBeam(const sixDoFRigidBodyStateFvBeam&) = delete;

        //- No copy assignment
        void operator=(const sixDoFRigidBodyStateFvBeam&) = delete;


protected:

    // Protected Member Functions

        //- Overloaded writeFileHeader from writeFile
        virtual void writeFileHeader(Ostream& os);


public:

    //- Runtime type information
    TypeName("sixDoFRigidBodyStateFvBeam");


    // Constructors

        //- Construct from Time and dictionary
        sixDoFRigidBodyStateFvBeam
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~sixDoFRigidBodyStateFvBeam() = default;


    // Member Functions

        //- Read the sixDoFRigidBodyStateFvBeam data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the sixDoFRigidBodyStateFvBeam
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
