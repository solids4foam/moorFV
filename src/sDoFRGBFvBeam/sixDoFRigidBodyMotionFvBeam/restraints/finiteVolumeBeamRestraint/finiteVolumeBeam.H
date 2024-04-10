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
            const sixDoFRigidBodyMotionFvBeam& motion,
            vector& restraintPosition,
            vector& restraintForce,
            vector& restraintMoment
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
