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

\*---------------------------------------------------------------------------*/

#include "sixDoFRigidBodyStateFvBeam.H"
#include "dynamicMotionSolverFvMesh.H"
#include "sixDoFRigidBodyMotionFvBeamSolver.H"
#include "unitConversion.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(sixDoFRigidBodyStateFvBeam, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        sixDoFRigidBodyStateFvBeam,
        dictionary
    );
}
}

const Foam::Enum
<
    Foam::functionObjects::sixDoFRigidBodyStateFvBeam::angleTypes
>
Foam::functionObjects::sixDoFRigidBodyStateFvBeam::angleTypeNames_
({
    { angleTypes::RADIANS, "radians" },
    { angleTypes::DEGREES, "degrees" },
});


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::sixDoFRigidBodyStateFvBeam::writeFileHeader(Ostream& os)
{
    writeHeader(os, "Motion State");
    writeHeaderValue(os, "Angle Units", angleTypeNames_[angleFormat_]);
    writeCommented(os, "Time");

    os  << tab
        << "centreOfRotation" << tab
        << "centreOfMass" << tab
        << "rotation" << tab
        << "velocity" << tab
        << "omega" << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::sixDoFRigidBodyStateFvBeam::sixDoFRigidBodyStateFvBeam
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name, typeName, dict),
    angleFormat_(angleTypes::RADIANS)
{
    read(dict);
    writeFileHeader(file());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::sixDoFRigidBodyStateFvBeam::read(const dictionary& dict)
{
    if (fvMeshFunctionObject::read(dict))
    {
        angleFormat_ =
            angleTypeNames_.getOrDefault
            (
                "angleFormat",
                dict,
                angleTypes::RADIANS
            );

        return true;
    }

    return false;
}


bool Foam::functionObjects::sixDoFRigidBodyStateFvBeam::execute()
{
    return true;
}


bool Foam::functionObjects::sixDoFRigidBodyStateFvBeam::write()
{
    const dynamicMotionSolverFvMesh& mesh =
        refCast<const dynamicMotionSolverFvMesh>(obr_);

    const sixDoFRigidBodyMotionFvBeamSolver& motionSolver_ =
        refCast<const sixDoFRigidBodyMotionFvBeamSolver>(mesh.motion());

    const sixDoFRigidBodyMotionFvBeam& motion = motionSolver_.motion();

    vector rotationAngle
    (
        quaternion(motion.orientation()).eulerAngles(quaternion::XYZ)
    );

    vector angularVelocity(motion.omega());

    switch (angleFormat_)
    {
        case angleTypes::RADIANS:
        {
            // Nothing to do - already in radians
            break;
        }
        case angleTypes::DEGREES:
        {
            rotationAngle.x() = radToDeg(rotationAngle.x());
            rotationAngle.y() = radToDeg(rotationAngle.y());
            rotationAngle.z() = radToDeg(rotationAngle.z());

            angularVelocity.x() = radToDeg(angularVelocity.x());
            angularVelocity.y() = radToDeg(angularVelocity.y());
            angularVelocity.z() = radToDeg(angularVelocity.z());
            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unhandled enumeration " << angleTypeNames_[angleFormat_]
                << abort(FatalError);
        }
    }

    writeCurrentTime(file());
    file()
        << tab
        << motion.centreOfRotation()  << tab
        << motion.centreOfMass()  << tab
        << rotationAngle  << tab
        << motion.v()  << tab
        << angularVelocity << endl;

    return true;
}


// ************************************************************************* //
