/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::polymerKTransport<Thermo>::polymerKTransport
(
    const Thermo& t,
    const scalar mu,
    const scalar Pr
)
:
    Thermo(t),
    mu_(mu),
    rPr_(1.0/Pr)
{}


template<class Thermo>
inline Foam::polymerKTransport<Thermo>::polymerKTransport
(
    const word& name,
    const polymerKTransport& ct
)
:
    Thermo(name, ct),
    mu_(ct.mu_),
    rPr_(ct.rPr_)
{}


template<class Thermo>
inline Foam::autoPtr<Foam::polymerKTransport<Thermo>>
Foam::polymerKTransport<Thermo>::clone() const
{
    return autoPtr<polymerKTransport<Thermo>>
    (
        new polymerKTransport<Thermo>(*this)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::polymerKTransport<Thermo>>
Foam::polymerKTransport<Thermo>::New
(
    const dictionary& dict
)
{
    return autoPtr<polymerKTransport<Thermo>>
    (
        new polymerKTransport<Thermo>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::polymerKTransport<Thermo>::mu  //- Constant dynamic viscosity [kg/m/s]
(
    const scalar p,
    const scalar T
) const
{
    return mu_;
}


template<class Thermo>
inline Foam::scalar Foam::polymerKTransport<Thermo>::kappa //- Thermal conductivity [W/m/K] --> tabulated
(
    const scalar p,
    const scalar T
) const
{
   // Old version: 
   // return this->Cp(p, T)*mu(p, T)*rPr_;
   
   //Thermal conductivity changed from constant to tabulated data table:
   #include "polymerK_Data_ThermalConduct.H"
}


template<class Thermo>
inline Foam::scalar Foam::polymerKTransport<Thermo>::alphah  //- Thermal diffusivity of enthalpy [kg/m/s]
(
    const scalar p,
    const scalar T
) const
{
// ORIGINAL VERSION Pr constant:
    // return mu(p, T)*rPr_;
    
    // NEW VERSION FOR POLYMER Pr not constant:
    // Pr = mu(p,T)*Cp(p,T)/kappa(p,T)
    // alpha = kappa(p,T)/cp(p,T)*rho(p,T) = kappa(p,T)/Cp(p,T)
    
    return kappa(p,T)/this->Cp(p,T);
    
    // end of alpha version implemented for polymer
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline void Foam::polymerKTransport<Thermo>::operator+=
(
    const polymerKTransport<Thermo>& st
)
{
    scalar Y1 = this->Y();

    Thermo::operator+=(st);

    if (mag(this->Y()) > small)
    {
        Y1 /= this->Y();
        scalar Y2 = st.Y()/this->Y();

        mu_ = Y1*mu_ + Y2*st.mu_;
        rPr_ = 1.0/(Y1/rPr_ + Y2/st.rPr_);
    }
}


template<class Thermo>
inline void Foam::polymerKTransport<Thermo>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::polymerKTransport<Thermo> Foam::operator+
(
    const polymerKTransport<Thermo>& ct1,
    const polymerKTransport<Thermo>& ct2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(ct1) + static_cast<const Thermo&>(ct2)
    );

    if (mag(t.Y()) < small)
    {
        return polymerKTransport<Thermo>
        (
            t,
            0,
            ct1.rPr_
        );
    }
    else
    {
        scalar Y1 = ct1.Y()/t.Y();
        scalar Y2 = ct2.Y()/t.Y();

        return polymerKTransport<Thermo>
        (
            t,
            Y1*ct1.mu_ + Y2*ct2.mu_,
            1.0/(Y1/ct1.rPr_ + Y2/ct2.rPr_)
        );
    }
}


template<class Thermo>
inline Foam::polymerKTransport<Thermo> Foam::operator*
(
    const scalar s,
    const polymerKTransport<Thermo>& ct
)
{
    return polymerKTransport<Thermo>
    (
        s*static_cast<const Thermo&>(ct),
        ct.mu_,
        1.0/ct.rPr_
    );
}


// ************************************************************************* //
