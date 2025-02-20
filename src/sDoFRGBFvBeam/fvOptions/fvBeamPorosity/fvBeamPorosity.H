/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 IH-Cantabria
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
    fvOption to subtract the momentum (based on drag) from momentum equation, similar to porosity effects.
    based on mangrove models for porosity.

Usage
    beamPorosity
    {
        type            fvBeamPorosity;
        active          yes;
    }

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

// Forward Declarations
class MangrovesModel;

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

        //- Name of the model
        word modelName_;

        //- Permablity coefficient
        scalar perm_;

        //- Viscosity
        scalar nu_;

        //- Density
        scalar rho_;

        //- penalization factor
        scalar pFactor_;


        //- Return the coefficient for the model
        tmp<volScalarField> coeff(const volVectorField& U, const word& modelName) const;



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
