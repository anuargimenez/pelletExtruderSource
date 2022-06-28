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

#include "rhoPolymerEOS.H"
#include "specie.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Specie>
inline Foam::rhoPolymerEOS<Specie>::rhoPolymerEOS
(
    const Specie& sp,
    const scalar b1m,
    const scalar b2m,
    const scalar b3m,
    const scalar b4m,
    const scalar b1s,
    const scalar b2s,
    const scalar b3s,
    const scalar b4s,
    const scalar b5,
    const scalar b6
)
:
    Specie(sp),
    b1m_(b1m),
    b2m_(b2m),
    b3m_(b3m),
    b4m_(b4m),
    b1s_(b1s),
    b2s_(b2s),
    b3s_(b3s),
    b4s_(b4s),
    b5_(b5),
    b6_(b6)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::rhoPolymerEOS<Specie>::rhoPolymerEOS
(
    const word& name,
    const rhoPolymerEOS<Specie>& pg
)
:
    Specie(name, pg),
    b1m_(pg.b1m_),
    b2m_(pg.b2m_),
    b3m_(pg.b3m_),
    b4m_(pg.b4m_),
    b1s_(pg.b1s_),
    b2s_(pg.b2s_),
    b3s_(pg.b3s_),
    b4s_(pg.b4s_),
    b5_(pg.b5_),
    b6_(pg.b6_)
{}


template<class Specie>
inline Foam::autoPtr<Foam::rhoPolymerEOS<Specie>>
Foam::rhoPolymerEOS<Specie>::clone() const
{
    return autoPtr<rhoPolymerEOS<Specie>>(new rhoPolymerEOS<Specie>(*this));
}


template<class Specie>
inline Foam::autoPtr<Foam::rhoPolymerEOS<Specie>> Foam::rhoPolymerEOS<Specie>::New
(
    const dictionary& dict
)
{
    return autoPtr<rhoPolymerEOS<Specie>>(new rhoPolymerEOS<Specie>(dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::scalar Foam::rhoPolymerEOS<Specie>::rho(scalar p, scalar T) const
{
    // OLD MODEL:
    // return p/(this->R()*T);
    
    // TABULATED MODEL:
    // #include "density_Tabulated_Data.H"
    
    // TAIT EOS:
    #include "TaitEOS.H"
}


template<class Specie>
inline Foam::scalar Foam::rhoPolymerEOS<Specie>::H(scalar p, scalar T) const // Enthalpy
{
    // THERMODYNAMICS LAW: H = U + p/rho
    return p/(this->rho(p, T));
}


template<class Specie>
inline Foam::scalar Foam::rhoPolymerEOS<Specie>::Cp(scalar p, scalar T) const
{
    return 0; // Specific Heat Cp tabulated in hPolymer
}


template<class Specie>
inline Foam::scalar Foam::rhoPolymerEOS<Specie>::E(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::rhoPolymerEOS<Specie>::Cv(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::rhoPolymerEOS<Specie>::S(scalar p, scalar T) const
{
    return -this->R()*log(p/Pstd);
}


template<class Specie>
inline Foam::scalar Foam::rhoPolymerEOS<Specie>::psi(scalar p, scalar T) const
{
    // old model, only for ideal gas
    // return 1.0/(this->R()*T);
    
    // new model
    
    return 0.0; 
    
    // this is only true if rho(p,T)=rho(T). However, we won´t use psi models in our model
}


template<class Specie>
inline Foam::scalar Foam::rhoPolymerEOS<Specie>::Z(scalar p, scalar T) const
{
    // old model
    // return 1.0;
    
    // new model
    
    return 0.0; 
    
    // We won´t use it in our model
}


template<class Specie>
inline Foam::scalar Foam::rhoPolymerEOS<Specie>::CpMCv(scalar p, scalar T) const
{
   // old model
   // return this->R();
   
   // new model
   
   return 0.0;

   // cpMcv should not be used with the rhoPolymerEOS model
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::rhoPolymerEOS<Specie>::operator+=(const rhoPolymerEOS<Specie>& pg)
{
    Specie::operator+=(pg);
}


template<class Specie>
inline void Foam::rhoPolymerEOS<Specie>::operator*=(const scalar s)
{
    Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::rhoPolymerEOS<Specie> Foam::operator+
(
    const rhoPolymerEOS<Specie>& pg1,
    const rhoPolymerEOS<Specie>& pg2
)
{
    return rhoPolymerEOS<Specie>
    (
        static_cast<const Specie&>(pg1) + static_cast<const Specie&>(pg2)
    );
}


template<class Specie>
inline Foam::rhoPolymerEOS<Specie> Foam::operator*
(
    const scalar s,
    const rhoPolymerEOS<Specie>& pg
)
{
    return rhoPolymerEOS<Specie>(s*static_cast<const Specie&>(pg));
}


template<class Specie>
inline Foam::rhoPolymerEOS<Specie> Foam::operator==
(
    const rhoPolymerEOS<Specie>& pg1,
    const rhoPolymerEOS<Specie>& pg2
)
{
    return rhoPolymerEOS<Specie>
    (
        static_cast<const Specie&>(pg1) == static_cast<const Specie&>(pg2)
    );
}


// ************************************************************************* //
