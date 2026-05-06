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
    public fv::option
{
    // Private Member Functions
        mutable autoPtr<OFstream> forceFilePtr_;

        //- No copy construct
        fvBeamPorosity(const fvBeamPorosity&) = delete;

        //- No copy assignment
        void operator=(const fvBeamPorosity&) = delete;


protected:

    // Protected data


    // Protected Member Functions
    // Geometric mapping data
    Foam::scalarField sField_;
    Foam::scalarField closestBeamCellDist_;
    Foam::labelList   closestBeamCell_;

    // Gaussian kernel width epsilon
    scalar eps_;

    // Beam mesh name
    word beamName_;

    // Fluid density
    scalar fluidRho_;

    // Beam cell center coordinates
    mutable vectorField beamC_;

    // Protected Member Functions

    // Gaussian kernel
    inline Foam::scalar eta(const scalar r) const;

    // Compute s, r and segment index
    void calculateS();

    // Helper function to apply the force
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

            //- Add implicit contribution to momentum equation
            virtual void addSup
            (
                fvMatrix<vector>& eqn,
                const label fieldi
            );

            //- Add implicit contribution to phase momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const label fieldi
            );


        // IO

            //- Read dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
