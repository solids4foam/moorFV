/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "finiteVolumeSixDoFRigidBodyDisplacementPointPatchVectorField.H"
#include "pointPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"
#include "fvMesh.H"
#include "volFields.H"
#include "uniformDimensionedFields.H"
#include "forces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

finiteVolumeSixDoFRigidBodyDisplacementPointPatchVectorField::
finiteVolumeSixDoFRigidBodyDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(p, iF),
    motion_(db().time()),
    initialPoints_(p.localPoints()),
    rhoInf_(1.0),
    rhoName_("rho"),
    lookupGravity_(-1),
    g_(Zero),
    curTimeIndex_(-1)
{}


finiteVolumeSixDoFRigidBodyDisplacementPointPatchVectorField::
finiteVolumeSixDoFRigidBodyDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchField<vector>(p, iF, dict),
    motion_(dict, dict, db().time()),
    rhoInf_(1.0),
    rhoName_(dict.getOrDefault<word>("rho", "rho")),
    lookupGravity_(-1),
    g_(Zero),
    curTimeIndex_(-1)
{
    if (rhoName_ == "rhoInf")
    {
        dict.readEntry("rhoInf", rhoInf_);
    }

    if (dict.readIfPresent("g", g_))
    {
        lookupGravity_ = -2;
    }

    if (!dict.found("value"))
    {
        updateCoeffs();
    }

    if (dict.found("initialPoints"))
    {
        initialPoints_ = vectorField("initialPoints", dict , p.size());
    }
    else
    {
        initialPoints_ = p.localPoints();
    }
}


finiteVolumeSixDoFRigidBodyDisplacementPointPatchVectorField::
finiteVolumeSixDoFRigidBodyDisplacementPointPatchVectorField
(
    const finiteVolumeSixDoFRigidBodyDisplacementPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchField<vector>(ptf, p, iF, mapper),
    motion_(ptf.motion_),
    initialPoints_(ptf.initialPoints_, mapper),
    rhoInf_(ptf.rhoInf_),
    rhoName_(ptf.rhoName_),
    lookupGravity_(ptf.lookupGravity_),
    g_(ptf.g_),
    curTimeIndex_(-1)
{}


finiteVolumeSixDoFRigidBodyDisplacementPointPatchVectorField::
finiteVolumeSixDoFRigidBodyDisplacementPointPatchVectorField
(
    const finiteVolumeSixDoFRigidBodyDisplacementPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(ptf, iF),
    motion_(ptf.motion_),
    initialPoints_(ptf.initialPoints_),
    rhoInf_(ptf.rhoInf_),
    rhoName_(ptf.rhoName_),
    lookupGravity_(ptf.lookupGravity_),
    g_(ptf.g_),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void finiteVolumeSixDoFRigidBodyDisplacementPointPatchVectorField::autoMap
(
    const pointPatchFieldMapper& m
)
{
    fixedValuePointPatchField<vector>::autoMap(m);

    initialPoints_.autoMap(m);
}


void finiteVolumeSixDoFRigidBodyDisplacementPointPatchVectorField::rmap
(
    const pointPatchField<vector>& ptf,
    const labelList& addr
)
{
    const finiteVolumeSixDoFRigidBodyDisplacementPointPatchVectorField& sDoFptf =
        refCast<const finiteVolumeSixDoFRigidBodyDisplacementPointPatchVectorField>(ptf);

    fixedValuePointPatchField<vector>::rmap(sDoFptf, addr);

    initialPoints_.rmap(sDoFptf.initialPoints_, addr);
}


void finiteVolumeSixDoFRigidBodyDisplacementPointPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (lookupGravity_ < 0)
    {
        if (db().time().foundObject<uniformDimensionedVectorField>("g"))
        {
            if (lookupGravity_ == -2)
            {
                FatalErrorInFunction
                    << "Specifying the value of g in this boundary condition "
                    << "when g is available from the database is considered "
                    << "a fatal error to avoid the possibility of inconsistency"
                    << exit(FatalError);
            }
            else
            {
                lookupGravity_ = 1;
            }
        }
        else
        {
            lookupGravity_ = 0;
        }
    }

    const polyMesh& mesh = this->internalField().mesh()();
    const Time& t = mesh.time();
    const pointPatch& ptPatch = this->patch();

    // Store the motion state at the beginning of the time-step
    bool firstIter = false;
    if (curTimeIndex_ != t.timeIndex())
    {
        motion_.newTime();
        curTimeIndex_ = t.timeIndex();
        firstIter = true;
    }

    dictionary forcesDict;

    forcesDict.add("type", functionObjects::forces::typeName);
    forcesDict.add("patches", wordList(1, ptPatch.name()));
    forcesDict.add("rhoInf", rhoInf_);
    forcesDict.add("rho", rhoName_);
    forcesDict.add("CofR", motion_.centreOfRotation());

    functionObjects::forces f("forces", db(), forcesDict);

    f.calcForcesMoments();

    // Get the forces on the patch faces at the current positions

    if (lookupGravity_ == 1)
    {
        uniformDimensionedVectorField g =
            db().time().lookupObject<uniformDimensionedVectorField>("g");

        g_ = g.value();
    }

    // scalar ramp = min(max((t.value() - 5)/10, 0), 1);
    scalar ramp = 1.0;

    motion_.update
    (
        firstIter,
        ramp*(f.forceEff() + motion_.mass()*g_),
        ramp*(f.momentEff() + motion_.mass()*(motion_.momentArm() ^ g_)),
        t.deltaTValue(),
        t.deltaT0Value()
    );

    Field<vector>::operator=
    (
        motion_.transform(initialPoints_) - initialPoints_
    );

    fixedValuePointPatchField<vector>::updateCoeffs();
}


void finiteVolumeSixDoFRigidBodyDisplacementPointPatchVectorField::write(Ostream& os) const
{
    pointPatchField<vector>::write(os);

    os.writeEntry("rho", rhoName_);

    if (rhoName_ == "rhoInf")
    {
        os.writeEntry("rhoInf", rhoInf_);
    }

    if (lookupGravity_ == 0 || lookupGravity_ == -2)
    {
        os.writeEntry("g", g_);
    }

    motion_.write(os);

    initialPoints_.writeEntry("initialPoints", os);

    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePointPatchTypeField
(
    pointPatchVectorField,
    finiteVolumeSixDoFRigidBodyDisplacementPointPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
