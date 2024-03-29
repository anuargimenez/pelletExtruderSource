/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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

#include "polymerKWLFTransport.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Thermo>
Foam::scalar Foam::polymerKWLFTransport<Thermo>::readCoeff
(
    const word& coeffName,
    const dictionary& dict
)
{
    return readScalar(dict.subDict("transport").lookup(coeffName));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::polymerKWLFTransport<Thermo>::polymerKWLFTransport(const dictionary& dict)
:
    Thermo(dict),
    mu0_(readCoeff("mu0", dict)),
    Tr_(readCoeff("Tr", dict)),
    Tm_(readCoeff("Tm", dict)),
    TS_(readCoeff("TS", dict)),
    MS_(readCoeff("MS", dict)),
    C1_(readCoeff("C1", dict)),
    C2_(readCoeff("C2", dict))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::polymerKWLFTransport<Thermo>::write(Ostream& os) const
{
    os  << this->specie::name() << endl
        << token::BEGIN_BLOCK  << incrIndent << nl;

    Thermo::write(os);

    dictionary dict("transport");
    dict.add("mu0", mu0_);
    dict.add("Tr", Tr_);
    dict.add("Tm", Tm_);
    dict.add("TS", TS_);
    dict.add("MS", MS_);
    dict.add("C1", C1_);
    dict.add("C2", C2_);

    os  << indent << dict.dictName() << dict
        << decrIndent << token::END_BLOCK << nl;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Thermo>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const polymerKWLFTransport<Thermo>& Kwlft
)
{
    Kwlft.write(os);
    return os;
}


// ************************************************************************* //
