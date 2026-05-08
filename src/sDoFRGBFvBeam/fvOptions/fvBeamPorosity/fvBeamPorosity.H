/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenCFD Ltd.
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
    Foam::fv::fvBeamPorosity

Group
    grpFvOptionsSources

Description
    Adds an actuator-line momentum source from a beam region to the fluid
    momentum equation.

    The beam region specified by \c beamName must provide the fields
    \c refW, \c W and \c almForce, and the \c fluidCellIDs label list.
    The force is distributed to the fluid cells with a Gaussian kernel of
    width \c epsilon.

Usage
    Example usage:
    \verbatim
        beamPorosity
        {
            type            fvBeamPorosity;
            active          yes;

            fvBeamPorosityCoeffs
            {
                // Optional velocity field name, default is U
                U               U;

                // Optional list of velocity field names
                // UNames       (U);

                // Gaussian kernel width
                epsilon         0.05;

                // Beam region name
                beamName        beam;

                // Fluid density used to scale the momentum source
                rho             1000;
            }
        }
    \endverbatim

Author
    Amirhossein Taran, UCD

SourceFiles
    fvBeamPorosity.C

\*---------------------------------------------------------------------------*/

#ifndef fvBeamPorosity_H
#define fvBeamPorosity_H

#include "fvOption.H"
#include "volFields.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fv
{

/*---------------------------------------------------------------------------*\
                   Class fvBeamPorosity Declaration
\*---------------------------------------------------------------------------*/

class fvBeamPorosity
:
    public option
{
    // Private data

        //- Optional output stream for integrated/source-force diagnostics
        mutable autoPtr<OFstream> forceFilePtr_;


    // Private Member Functions

        //- No copy construct
        fvBeamPorosity(const fvBeamPorosity&) = delete;

        //- No copy assignment
        void operator=(const fvBeamPorosity&) = delete;


protected:

    // Protected data

        //- Local coordinate of the closest point on the nearest beam segment
        //  for each fluid cell. Values are clipped to the range [0, 1].
        scalarField sField_;

        //- Distance from each fluid-cell centre to its nearest beam segment
        scalarField closestBeamCellDist_;

        //- Index of the nearest beam segment associated with each fluid cell
        labelList closestBeamCell_;

        //- Gaussian spreading-kernel width
        scalar eps_;

        //- Name of the beam mesh/object registry used for beam fields
        word beamName_;

        //- Constant fluid density used when assembling the momentum source
        scalar fluidRho_;

        //- Current beam-cell/actuator-point coordinates
        mutable vectorField beamC_;


    // Protected Member Functions

        //- Return the 2-D Gaussian radial kernel value at distance r
        inline scalar eta(const scalar r) const;

        //- Calculate nearest beam segment, segment coordinate and distance
        //  for all fluid cells
        void calculateS();

        //- Spread the beam force to the fluid cells and add it to the
        //  momentum-equation source term
        void applyBeamForce(fvMatrix<vector>& eqn);


public:

    //- Runtime type information
    TypeName("fvBeamPorosity");


    // Constructors

        //- Construct from components
        fvBeamPorosity
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~fvBeamPorosity() = default;


    // Member Functions

        // Add explicit and implicit contributions

            //- Add source contribution to the momentum equation
            virtual void addSup
            (
                fvMatrix<vector>& eqn,
                const label fieldi
            );

            //- Add source contribution to a density-weighted momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const label fieldi
            );


        // IO

            //- Read fvOption coefficients from dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
