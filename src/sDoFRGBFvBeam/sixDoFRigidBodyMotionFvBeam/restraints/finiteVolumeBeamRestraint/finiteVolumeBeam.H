/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::sixDoFRigidBodyMotionFvBeamRestraints::finiteVolumeBeam

SourceFiles
    finiteVolumeBeam.C

\*---------------------------------------------------------------------------*/

#ifndef finiteVolumeBeam_H
#define finiteVolumeBeam_H

#include "sixDoFRigidBodyMotionFvBeamRestraint.H"
#include "point.H"
#include "beamModel.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace sixDoFRigidBodyMotionFvBeamRestraints
{

/*---------------------------------------------------------------------------*\
                          Class finiteVolumeBeam Declaration
\*---------------------------------------------------------------------------*/

class finiteVolumeBeam
:
    public sixDoFRigidBodyMotionFvBeamRestraint
{
    // Private data

        //- Pointer to beam model
        mutable autoPtr<Foam::beamModel> beam_;

        //- Restraint state information
        //  This store data to allow consistent restarts
        mutable IOdictionary state_;

        //- Attachment point
        point refAttachmentPt_;

        //- Attachment patch
        word attachmentPatch_;
        word anchorPatch_;

        //- Owner of beam attachment kinematics
        word beamAttachmentKinematics_;

        //- Index of attachment patch
        label patchID_;
        label anchorPatchID_;

        //- Initial W on the patch
        mutable vector initialW_;

        //- Flag to indicate if initial W has been stored
        mutable bool storeInitialW_;

        //- Initial Q ( force ) on the patch
        mutable vector initialQ_;

        //- File for writing attachment force
        mutable autoPtr<OFstream> forceFilePtr_;
        mutable autoPtr<OFstream> anchorForceFilePtr_;
        mutable autoPtr<OFstream> attachmentForceFilePtr_;
        mutable autoPtr<OFstream> displacementFilePtr_;

        //- File for writing the axial (along the line) and shear force at
        //  the anchor and attachment ends
        mutable autoPtr<OFstream> axialForceFilePtr_;

        //- True once beamFoam's rigidBodyEnd has been set up for the
        //  beamFoamCoupled solver
        mutable bool monolithicInitialised_;


    // Private Member Functions

        //- Write the axial and shear force at both ends of the line.
        //  Q on an end face is the whole internal force (axial plus
        //  shear); the axial part is Q along the current line tangent,
        //  positive in tension
        void writeAxialForce(const scalar t) const;


public:

        //- Runtime type information
        TypeName("finiteVolumeBeam");


        // Constructors

        //- Construct from components
        finiteVolumeBeam
        (
            const word& name,
            const dictionary& sDoFRBMRDict,
            const Time& time
        );

        //- Construct and return a clone
        virtual autoPtr<sixDoFRigidBodyMotionFvBeamRestraint> clone() const
        {
            return autoPtr<sixDoFRigidBodyMotionFvBeamRestraint>
            (
                new finiteVolumeBeam(*this)
            );
        }


    //- Destructor
    virtual ~finiteVolumeBeam();


    // Member Functions

        //- Returns non-const access to the beam model
        //const beamModel& beam() const;

        //- Calculate the restraint position, force and moment.
        //  Global reference frame vectors.
        virtual void restrain
        (
            sixDoFRigidBodyMotionFvBeam& motion,
            vector& restraintPosition,
            vector& restraintForce,
            vector& restraintMoment
        ) const;

        //- Monolithic coupling (beamFoamCoupled solver): solve the beam and
        //  the body together in beamFoam for an applied force and moment
        //  (global axes, moment about the centre of mass) and return the
        //  accepted body state. On the first call the body is set up in
        //  beamFoam from the motion; solverDict gives its Newmark and
        //  damping settings. tConstraints and rConstraints project onto the
        //  free translation and rotation directions (global axes)
        RigidBodyEndState solveMonolithic
        (
            const sixDoFRigidBodyMotionFvBeam& motion,
            const vector& force,
            const vector& moment,
            const tensor& tConstraints,
            const tensor& rConstraints,
            const dictionary& solverDict
        ) const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& sDoFRBMRCoeff);

        //- Write
        virtual void write(Ostream&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
