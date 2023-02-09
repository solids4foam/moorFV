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

#include "finiteVolumeBeamRestraint.H"
#include "OutputControlDictionary.H" 
#include "addToRunTimeSelectionTable.H"
#include "sixDoFRigidBodyMotionNew.H"
#include "beamModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sixDoFRigidBodyMotionRestraints
{
    defineTypeNameAndDebug(finiteVolumeBeamRestraint, 0);
    addToRunTimeSelectionTable
    (
        sixDoFRigidBodyMotionRestraintNew,
        finiteVolumeBeamRestraint,
        dictionary
    );
};
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionRestraints::finiteVolumeBeamRestraint::finiteVolumeBeamRestraint
(
    const dictionary& sDoFRBMRDict, const Time& time
    //,const IOobject& io   Q.Here?
)
:
    sixDoFRigidBodyMotionRestraintNew(sDoFRBMRDict, time),
    refAttachmentPt_()
{
    read(sDoFRBMRDict);

    Foam::autoPtr<Foam::beamModel> beam =
        Foam::beamModel::New
        (
            const_cast<Time&>(time) , Foam::dynamicFvMesh::defaultRegion
        );

}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionRestraints::finiteVolumeBeamRestraint::~finiteVolumeBeamRestraint()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::sixDoFRigidBodyMotionRestraints::finiteVolumeBeamRestraint::restrain
(
    const sixDoFRigidBodyMotionNew& motion,
    vector& restraintPosition,
    vector& restraintForce,
    vector& restraintMoment
) const
{
    restraintPosition = motion.currentPosition(refAttachmentPt_);

    // volVectorField& DW = const_cast<volVectorField&>
    // (
    //     dict().lookupObject<volVectorField>("DW")
    // );
    //const fvMesh& mesh = DW.mesh();
    // const label patchID =
    //     mesh.boundaryMesh().findPatchID("right"); //manualy
    //     //("name_of_beam_attachment_patch");
    // if (patchID == -1)
    // {
    //     FatalErrorIn("Foam::finiteVolumeBeam::restrainingForce(...)")
    //         << "Attachment patch not found!" << abort(FatalError);
    // }

    if (motion.report())
    {
        Info<< " force " << restraintForce
            << " moment " << restraintMoment
            << endl;
    }
}


bool Foam::sixDoFRigidBodyMotionRestraints::finiteVolumeBeamRestraint::read
(
    const dictionary& sDoFRBMRDict
)
{
    
    sixDoFRigidBodyMotionRestraintNew::read(sDoFRBMRDict);


    sDoFRBMRCoeffs_.lookup("refAttachmentPt") >> refAttachmentPt_;

    return true;
}


void Foam::sixDoFRigidBodyMotionRestraints::finiteVolumeBeamRestraint::write
(
    Ostream& os
) const
{

    os.writeKeyword("refAttachmentPt")
        << refAttachmentPt_ << token::END_STATEMENT << nl;

}

// ************************************************************************* //
