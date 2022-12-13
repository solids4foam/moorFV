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
    Foam::finiteVolumeBeam

Description
    sixDoFRigidBodyMotionRestraints model.  Finite volume beam solver based on
    Simo–Reissner theory.

    See Bali, S, Tuković, Ž, Cardiff, P, Ivanković, A, Pakrashi, V. A cell-
    centered finite volume formulation of geometrically exact Simo–Reissner
    beams with arbitrary initial curvatures. Int J Numer Methods Eng. 2022;
    123( 17): 3950– 3973. doi:10.1002/nme.6994

Author
    Amirhossein Taran, UCD.
    Philip Cardiff, UCD.

SourceFiles
    finiteVolumeBeam.C

\*---------------------------------------------------------------------------*/

#ifndef finiteVolumeBeam_H
#define finiteVolumeBeam_H

#include "combinedRestraint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class finiteVolumeBeam Declaration
\*---------------------------------------------------------------------------*/

class finiteVolumeBeam
:
    public combinedRestraint
{
    // Private Data

        //- Beam inputs: to be added
        //dimensionedScalar coeff_;


public:

    //- Runtime type information
    TypeName("finiteVolumeBeam");


    // Constructors

        //- Construct from dictionary
        finiteVolumeBeam
        (
            const word& name,
            const dictionary& dict,
            const sixDOFODE& sixDOF
        );

        //- Construct and return a clone
        virtual autoPtr<combinedRestraint> clone() const;


    // Destructor

        virtual ~finiteVolumeBeam();


    // Member Functions

        // Restraint specific functions

            //- Return restraining force (in the global coordinate system)
            virtual vector restrainingForce
            (
                const scalar t,
                const tensor& toRelative,
                const vector& x,
                const vector& u
            ) const;

            //- Return restraining moment (in the local coordinate system)
            virtual vector restrainingMoment
            (
                const scalar t,
                const tensor& toRelative,
	        const vector& omega
            ) const;


    // I-O Functions and Operators

        //- Virtual write function
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
