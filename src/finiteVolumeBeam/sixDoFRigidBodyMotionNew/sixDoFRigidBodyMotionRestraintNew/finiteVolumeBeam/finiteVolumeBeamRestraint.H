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
    Foam::sixDoFRigidBodyMotionRestraints::finiteVolumeBeamRestraint

Description
    sixDoFRigidBodyMotionRestraints model.  Linear spring.

SourceFiles
    finiteVolumeBeamRestraint.C

\*---------------------------------------------------------------------------*/

#ifndef finiteVolumeBeamRestraint_H
#define finiteVolumeBeamRestraint_H

#include "sixDoFRigidBodyMotionRestraintNew.H"
#include "point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace sixDoFRigidBodyMotionRestraints
{

/*---------------------------------------------------------------------------*\
                          Class finiteVolumeBeamRestraint Declaration
\*---------------------------------------------------------------------------*/

class finiteVolumeBeamRestraint
:
    public sixDoFRigidBodyMotionRestraintNew
{
    // Private data

        //- Anchor point, where the spring is attached to an immovable
        //  object
        point anchor_;

        //- Reference point of attachment to the solid body
        point refAttachmentPt_;

        //- Spring stiffness coefficient (N/m)
        scalar stiffness_;

        //- Damping coefficient (Ns/m)
        scalar damping_;

        //- Rest length - length of spring when no forces are applied to it
        scalar restLength_;


public:

    //- Runtime type information
    TypeName("finiteVolumeBeamRestraint");


    // Constructors

        //- Construct from components
        finiteVolumeBeamRestraint
        (
            const dictionary& sDoFRBMRDict, const Time& time
        );

        //- Construct and return a clone
        virtual autoPtr<sixDoFRigidBodyMotionRestraintNew> clone() const
        {
            return autoPtr<sixDoFRigidBodyMotionRestraintNew>
            (
                new finiteVolumeBeamRestraint(*this)
            );
        }


    // Destructor

        virtual ~finiteVolumeBeamRestraint();


    // Member Functions

        //- Calculate the restraint position, force and moment.
        //  Global reference frame vectors.
        virtual void restrain
        (
            const sixDoFRigidBodyMotionNew& motion,
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
