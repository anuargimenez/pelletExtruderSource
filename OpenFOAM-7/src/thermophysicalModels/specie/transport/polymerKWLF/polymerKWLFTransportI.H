/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2019 OpenFOAM Foundation
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

#include "specie.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::polymerKWLFTransport<Thermo>::polymerKWLFTransport
(
    const word& name,
    const polymerKWLFTransport& Kwlft
)
:
    Thermo(name, Kwlft),
    mu0_(Kwlft.mu0_),
    Tr_(Kwlft.Tr_),
    Tm_(Kwlft.Tm_),
    TS_(Kwlft.TS_),
    MS_(Kwlft.MS_),
    C1_(Kwlft.C1_),
    C2_(Kwlft.C2_)
{}


template<class Thermo>
inline Foam::autoPtr<Foam::polymerKWLFTransport<Thermo>>
Foam::polymerKWLFTransport<Thermo>::clone() const
{
    return autoPtr<polymerKWLFTransport<Thermo>>
    (
        new polymerKWLFTransport<Thermo>(*this)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::polymerKWLFTransport<Thermo>>
Foam::polymerKWLFTransport<Thermo>::New
(
    const dictionary& dict
)
{
    return autoPtr<polymerKWLFTransport<Thermo>>
    (
        new polymerKWLFTransport<Thermo>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::polymerKWLFTransport<Thermo>::mu
(
    const scalar p,
    const scalar T
) const
{
	const scalar Tg = 373.15; //Tg is an additional parameter for the modeling here, defines where constant zero viscosity is to be calculated. Above Tg zero viscosity follows WLF. Below Tg zero viscosity is constant with same value as Tg in WLF.

	if (((-C1_*(Tm_ - Tr_)/(C2_ + Tm_ - Tr_))) < (0))
	{
		Info << "Warning! Anomaly in the WLF model" << endl;
	}


	if ((T-Tg) < (0))
	{
		
		// if T<Tg=100ºC (transition of MF), constant zero viscosity. Always use Tg=100ºC=373.15 K
		return mu0_*exp((1-this->MF(p, T))*(((-C1_*(Tm_ - Tr_)/(C2_ + Tm_ - Tr_))))+this->MF(p, T)*((-C1_*(Tg - Tr_)/(C2_ + Tg - Tr_))))+10e8*(-this->MF(p, T)+1);	
	}
	else
	{
		// if T>Tg=100ºC (transition of MF), non-constant zero viscosity. Use T.
		return mu0_*exp((1-this->MF(p, T))*(((-C1_*(Tm_ - Tr_)/(C2_ + Tm_ - Tr_))))+this->MF(p, T)*((-C1_*(T - Tr_)/(C2_ + T - Tr_))))+10e8*(-this->MF(p, T)+1);		
	}

}

template<class Thermo>
inline Foam::scalar Foam::polymerKWLFTransport<Thermo>::MF	// Calculation of MF using the Karrenberg approach
(
    const scalar p, const scalar T
) const
{
	return Foam::pow(((Foam::tanh((T-Tm_)*TS_))+1.0)/2.0,MS_);
}


template<class Thermo>
inline Foam::scalar Foam::polymerKWLFTransport<Thermo>::kappa //- Thermal conductivity [W/m/K] --> tabulated
(
    const scalar p, const scalar T
) const
{
   // Old version: 
   // return this->Cp(p, T)*mu(p, T)*rPr_;
   
   //Thermal conductivity changed from constant to tabulated data table:
   #include "polymerKWLF_Data_ThermalConduct.H"
}


template<class Thermo>
inline Foam::scalar Foam::polymerKWLFTransport<Thermo>::alphah //- Thermal diffusivity of enthalpy [kg/m/s]
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
inline void Foam::polymerKWLFTransport<Thermo>::operator+=
(
    const polymerKWLFTransport<Thermo>& Kwlft
)
{
    scalar Y1 = this->Y();

    Thermo::operator+=(Kwlft);

    if (mag(this->Y()) > small)
    {
        Y1 /= this->Y();
        scalar Y2 = Kwlft.Y()/this->Y();

        mu0_ = Y1*mu0_ + Y2*Kwlft.mu0_;
        Tr_ = Y1*Tr_ + Y2*Kwlft.Tr_;
        Tm_ = Y1*Tm_ + Y2*Kwlft.Tm_;
        TS_ = Y1*TS_ + Y2*Kwlft.TS_;
        MS_ = Y1*MS_ + Y2*Kwlft.MS_;
        C1_ = Y1*C1_ + Y2*Kwlft.C1_;
        C2_ = Y1*C2_ + Y2*Kwlft.C2_;
    }
}


template<class Thermo>
inline void Foam::polymerKWLFTransport<Thermo>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::polymerKWLFTransport<Thermo> Foam::operator+
(
    const polymerKWLFTransport<Thermo>& Kwlft1,
    const polymerKWLFTransport<Thermo>& Kwlft2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(Kwlft1) + static_cast<const Thermo&>(Kwlft2)
    );

    if (mag(t.Y()) < small)
    {
        return polymerKWLFTransport<Thermo>
        (
            t,
            0,
            Kwlft1.mu0_,
            Kwlft1.Tr_,
            Kwlft1.C1_,
            Kwlft1.C2_
        );
    }
    else
    {
        scalar Y1 = Kwlft1.Y()/t.Y();
        scalar Y2 = Kwlft2.Y()/t.Y();
        return polymerKWLFTransport<Thermo>
        (
            t,
            Y1*Kwlft1.mu0_ + Y2*Kwlft2.mu0_,
            Y1*Kwlft1.Tr_ + Y2*Kwlft2.Tr_,
            Y1*Kwlft1.C1_ + Y2*Kwlft2.C1_,
            Y1*Kwlft1.C2_ + Y2*Kwlft2.C2_
        );
    }
}


template<class Thermo>
inline Foam::polymerKWLFTransport<Thermo> Foam::operator*
(
    const scalar s,
    const polymerKWLFTransport<Thermo>& Kwlft
)
{
    return polymerKWLFTransport<Thermo>
    (
        s*static_cast<const Thermo&>(Kwlft),
        Kwlft.mu0_,
        Kwlft.Tr_,
        Kwlft.Tm_,
        Kwlft.TS_,
        Kwlft.MS_,
        Kwlft.C1_,
        Kwlft.C2_
    );
}


// ************************************************************************* //
