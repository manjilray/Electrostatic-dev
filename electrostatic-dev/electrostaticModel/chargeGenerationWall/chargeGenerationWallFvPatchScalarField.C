/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "chargeGenerationWallFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "fvCFD.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "dictionary.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::chargeGenerationWallFvPatchScalarField::
chargeGenerationWallFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),

    mesh_
    (
        this->db()
    ),

    electrostaticProperties_
    (
        mesh_.lookupObject<dictionary>("electrostaticProperties")
    )
{}

Foam::chargeGenerationWallFvPatchScalarField::
chargeGenerationWallFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),

    mesh_
    (
        this->db()
    ),

    electrostaticProperties_
    (
        mesh_.lookupObject<dictionary>("electrostaticProperties")
    ),

    vacuumPermittivity_
    (
        electrostaticProperties_.lookup("vacuumPermittivity")
    ),
    particleRelPermittivity_
    (
        electrostaticProperties_.lookup("particleRelPermittivity")
    ),
    fluidRelPermittivity_
    (
        electrostaticProperties_.lookup("fluidRelPermittivity")
    ),
    saturatedRhoq_
    (
        electrostaticProperties_.lookup("saturatedRhoq")
    ),
    poissonsRatio_
    (
        electrostaticProperties_.lookup("poissonsRatio")
    ),
    youngsModulus_
    (
        electrostaticProperties_.lookup("youngsModulus")
    ),
    chargingEfficiencyWall_
    (
        electrostaticProperties_.lookup("chargingEfficiencyParticle")
    ),
    potentialDifference_
    (
        electrostaticProperties_.lookup("potentialDifference")
    ),
    particleWorkFunction_
    (
        electrostaticProperties_.lookup("particleWorkFunction")
    ),
    particleDiameter_
    (
        electrostaticProperties_.lookup("particleDiameter")
    ),
    particleDensity_
    (
        electrostaticProperties_.lookup("particleDensity")
    ),
    alphaMax_
    (
        electrostaticProperties_.lookup("alphaMax")
    ),
    z0_
    (
        electrostaticProperties_.lookup("z0")
    ),
    ks_
    (
        1.364*(pow(particleDiameter_,2.0))*(pow((particleDensity_*(1.0-poissonsRatio_)/youngsModulus_),2.0/5.0))
    ),
    kqs_
    (
        0.892023*ks_*chargingEfficiencyWall_*vacuumPermittivity_*particleRelPermittivity_*particleDensity_/pow(particleDiameter_,3.0)
    )

{
    this->evaluate();
}

Foam::chargeGenerationWallFvPatchScalarField::
chargeGenerationWallFvPatchScalarField
(
    const chargeGenerationWallFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),

    mesh_
    (
        this->db()
    ),

    electrostaticProperties_
    (
        mesh_.lookupObject<dictionary>("electrostaticProperties")
    )
{}

Foam::chargeGenerationWallFvPatchScalarField::
chargeGenerationWallFvPatchScalarField
(
    const chargeGenerationWallFvPatchScalarField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),

    mesh_
    (
        this->db()
    ),

    electrostaticProperties_
    (
        mesh_.lookupObject<dictionary>("electrostaticProperties")
    )
{}

Foam::chargeGenerationWallFvPatchScalarField::
chargeGenerationWallFvPatchScalarField
(
    const chargeGenerationWallFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),

    mesh_
    (
        this->db()
    ),

    electrostaticProperties_
    (
        mesh_.lookupObject<dictionary>("electrostaticProperties")
    )
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::chargeGenerationWallFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const fvPatch& patch(this->patch());

    // Electric field at the wall boundary
    const fvPatchVectorField& EqsPatchField
    (
        patch.lookupPatchField<volVectorField, vector>("Eq")
    );

    // Effect of electric field on wall charge flux
    scalarField Eqfactor
    (
        (z0_.value()/potentialDifference_.value())*
        (EqsPatchField & patch.Sf())/patch.magSf()
    );

    // Compute particle charge at wall boundary
    scalarField qWall
    (
        saturatedRhoq_.value()*(scalar(1.0)+Eqfactor)
    );

    operator==(qWall);

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::chargeGenerationWallFvPatchScalarField::write(Ostream& os) const
{
    fixedValueFvPatchScalarField::write(os);

    this->writeEntry("value", os);
}


// ************************************************************************* //
namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        chargeGenerationWallFvPatchScalarField
    );
}
// ************************************************************************* //
