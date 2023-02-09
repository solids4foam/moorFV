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
    const dictionary& sDoFRBMRDict,
    const Time& time
)
:
    sixDoFRigidBodyMotionRestraintNew(sDoFRBMRDict, time),
    beam_(beamModel::New(const_cast<Time&>(time) , "beam")),
    refAttachmentPt_(sDoFRBMRCoeffs_.lookup("refAttachmentPt")),
    patchID_(-1)
{
    // Call base class read function
    read(sDoFRBMRDict);

    // Find attachment patch
    patchID_ =
        beam_->mesh().boundaryMesh().findPatchID
        (
            sDoFRBMRDict.lookup("attachmentPatch")
        );

    if (patchID_ == -1)
    {
        FatalErrorIn("Foam::finiteVolumeBeam::finiteVolumeBeam(...)")
            << "Attachment patch not found!" << abort(FatalError);
    }
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
    // 1. Calculate attachment displacement from the motion object

    restraintPosition = motion.currentPosition(refAttachmentPt_);

    const vector attachmentDisp = vector(0.1, 0.1, 0.1);

    // 2. Set this displacement condition on the attachment patch
    // Note: W is the total displacement

    volVectorField& W = beam_->solutionW();
    W.boundaryField()[patchID_] == attachmentDisp;
    
    // 3. Solve the beam model
    const_cast<finiteVolumeBeamRestraint&>(*this).beam().evolve();

    // 4. Extract the force from the beam attachment
    const surfaceVectorField& Q =
        beam_->mesh().lookupObject<surfaceVectorField>("Q");

    if (Q.boundaryField()[patchID_].size() != 1)
    {
        FatalError
            << "Q.boundaryField()[patchID_].size() != 1: "
            << "this needs to be parallelised" << abort(FatalError);
    }

    const vector attachmentForce = Q.boundaryField()[patchID_][0];

    Info<< "attachment force = " << attachmentForce << endl;
    
    if (motion.report())
    {
        Info<< " force " << restraintForce
            << " moment " << restraintMoment
            << endl;
    }
}


// bool Foam::sixDoFRigidBodyMotionRestraints::finiteVolumeBeamRestraint::read
// (
//     const dictionary& sDoFRBMRDict
// )
// {
    
//     sixDoFRigidBodyMotionRestraintNew::read(sDoFRBMRDict);


//     sDoFRBMRCoeffs_.lookup("refAttachmentPt") >> refAttachmentPt_;

//     return true;
// }


void Foam::sixDoFRigidBodyMotionRestraints::finiteVolumeBeamRestraint::write
(
    Ostream& os
) const
{

    os.writeKeyword("refAttachmentPt")
        << refAttachmentPt_ << token::END_STATEMENT << nl;

}

// ************************************************************************* //
