/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::sixDoFRigidBodyMotionRestraints::finiteVolumeBeam

Description
    sixDoFRigidBodyMotionRestraints model.  Finite volume beam solver based on
    Simo–Reissner theory.

    See Bali, S, Tuković, Ž, Cardiff, P, Ivanković, A, Pakrashi, V. A cell-
    centered finite volume formulation of geometrically exact Simo–Reissner
    beams with arbitrary initial curvatures. Int J Numer Methods Eng. 2022;
    123( 17): 3950– 3973. doi:10.1002/nme.6994

SourceFiles
    finiteVolumeBeam.C

Author
    Amirhossein Taran, UCD.
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef finiteVolumeBeam_H
#define finiteVolumeBeam_H

#include "sixDoFRigidBodyMotionRestraint.H"
#include "point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace sixDoFRigidBodyMotionRestraints
{

/*---------------------------------------------------------------------------*\
                          Class finiteVolumeBeam Declaration
\*---------------------------------------------------------------------------*/

class finiteVolumeBeam
:
    public sixDoFRigidBodyMotionRestraint
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
    TypeName("finiteVolumeBeam");


    // Constructors

        //- Construct from components
        finiteVolumeBeam
        (
            const word& name,
            const dictionary& sDoFRBMRDict
        );

        //- Construct and return a clone
        virtual autoPtr<sixDoFRigidBodyMotionRestraint> clone() const
        {
            return autoPtr<sixDoFRigidBodyMotionRestraint>
            (
                new finiteVolumeBeam(*this)
            );
        }


    //- Destructor
    virtual ~finiteVolumeBeam();


    // Member Functions

        //- Calculate the restraint position, force and moment.
        //  Global reference frame vectors.
        virtual void restrain
        (
            const sixDoFRigidBodyMotion& motion,
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
