/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::hPolymerThermo<EquationOfState>::hPolymerThermo
(
    const EquationOfState& st,
    const scalar cp,
    const scalar hf
)
:
    EquationOfState(st),
    Cp_(cp),
    Hf_(hf)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::hPolymerThermo<EquationOfState>::hPolymerThermo
(
    const word& name,
    const hPolymerThermo& ct
)
:
    EquationOfState(name, ct),
    Cp_(ct.Cp_),
    Hf_(ct.Hf_)
{}


template<class EquationOfState>
inline Foam::autoPtr<Foam::hPolymerThermo<EquationOfState>>
Foam::hPolymerThermo<EquationOfState>::clone() const
{
    return autoPtr<hPolymerThermo<EquationOfState>>
    (
        new hPolymerThermo<EquationOfState>(*this)
    );
}


template<class EquationOfState>
inline Foam::autoPtr<Foam::hPolymerThermo<EquationOfState>>
Foam::hPolymerThermo<EquationOfState>::New(const dictionary& dict)
{
    return autoPtr<hPolymerThermo<EquationOfState>>
    (
        new hPolymerThermo<EquationOfState>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::scalar Foam::hPolymerThermo<EquationOfState>::limit
(
    const scalar T
) const
{
    return T;
}


template<class EquationOfState>
inline Foam::scalar Foam::hPolymerThermo<EquationOfState>::Cp   //- Heat capacity at constant pressure [J/kg/K]
(
    const scalar p,
    const scalar T
) const
{
   // original model
   // return Cp_ + EquationOfState::Cp(p, T);
   
   // new model with cp tabulated:
   #include "cp_Data_Polymer.H"
}


template<class EquationOfState>
inline Foam::scalar Foam::hPolymerThermo<EquationOfState>::Hs //- Sensible enthalpy [J/kg]
(
    const scalar p, const scalar T
) const
{
   // original model
   // return Cp_*(T - Tref_) + Hsref_ + EquationOfState::H(p, T);
   
   // new model:
   
   return Ha(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::hPolymerThermo<EquationOfState>::Hc() const  //- Chemical enthalpy [J/kg]
{
   // original model
   // return Hf_;
   
   return 0.0;
}


template<class EquationOfState>
inline Foam::scalar Foam::hPolymerThermo<EquationOfState>::Ha  //- Absolute Enthalpy [J/kg]
(
    const scalar p, const scalar T
) const
{
   // original model
   // return Hs(p, T) + Hc();
   
   // new model:
   #include "h_Data_Polymer.H"
}


template<class EquationOfState>
inline Foam::scalar Foam::hPolymerThermo<EquationOfState>::S
(
    const scalar p, const scalar T
) const
{
    return Cp_*log(T/Tstd) + EquationOfState::S(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::hPolymerThermo<EquationOfState>::dGdT
(
    const scalar p, const scalar T
) const
{
    return 0;
}


template<class EquationOfState>
inline Foam::scalar Foam::hPolymerThermo<EquationOfState>::dCpdT
(
    const scalar p, const scalar T
) const
{
    return 0;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline void Foam::hPolymerThermo<EquationOfState>::operator+=
(
    const hPolymerThermo<EquationOfState>& ct
)
{
    scalar Y1 = this->Y();

    EquationOfState::operator+=(ct);

    if (mag(this->Y()) > small)
    {
        Y1 /= this->Y();
        scalar Y2 = ct.Y()/this->Y();

        Cp_ = Y1*Cp_ + Y2*ct.Cp_;
        Hf_ = Y1*Hf_ + Y2*ct.Hf_;
    }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::hPolymerThermo<EquationOfState> Foam::operator+
(
    const hPolymerThermo<EquationOfState>& ct1,
    const hPolymerThermo<EquationOfState>& ct2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(ct1)
      + static_cast<const EquationOfState&>(ct2)
    );

    if (mag(eofs.Y()) < small)
    {
        return hPolymerThermo<EquationOfState>
        (
            eofs,
            ct1.Cp_,
            ct1.Hf_
        );
    }
    else
    {
        return hPolymerThermo<EquationOfState>
        (
            eofs,
            ct1.Y()/eofs.Y()*ct1.Cp_
          + ct2.Y()/eofs.Y()*ct2.Cp_,
            ct1.Y()/eofs.Y()*ct1.Hf_
          + ct2.Y()/eofs.Y()*ct2.Hf_
        );
    }
}


template<class EquationOfState>
inline Foam::hPolymerThermo<EquationOfState> Foam::operator*
(
    const scalar s,
    const hPolymerThermo<EquationOfState>& ct
)
{
    return hPolymerThermo<EquationOfState>
    (
        s*static_cast<const EquationOfState&>(ct),
        ct.Cp_,
        ct.Hf_
    );
}


template<class EquationOfState>
inline Foam::hPolymerThermo<EquationOfState> Foam::operator==
(
    const hPolymerThermo<EquationOfState>& ct1,
    const hPolymerThermo<EquationOfState>& ct2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(ct1)
     == static_cast<const EquationOfState&>(ct2)
    );

    return hPolymerThermo<EquationOfState>
    (
        eofs,
        ct2.Y()/eofs.Y()*ct2.Cp_
      - ct1.Y()/eofs.Y()*ct1.Cp_,
        ct2.Y()/eofs.Y()*ct2.Hf_
      - ct1.Y()/eofs.Y()*ct1.Hf_
    );
}


// ************************************************************************* //
