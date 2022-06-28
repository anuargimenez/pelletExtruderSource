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
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
template<class Specie>
Foam::rhoPolymerEOS<Specie>::rhoPolymerEOS(Istream& is)
:
    Specie(is),
    b1m_(readScalar(is)),
    b2m_(readScalar(is)),
    b3m_(readScalar(is)),
    b4m_(readScalar(is)),
    b1s_(readScalar(is)),
    b2s_(readScalar(is)),
    b3s_(readScalar(is)),
    b4s_(readScalar(is)),
    b5_(readScalar(is)),
    b6_(readScalar(is))
{
    is.check("rhoPolymerEOS<Specie>::rhoPolymerEOS(Istream& is)");
}


template<class Specie>
Foam::rhoPolymerEOS<Specie>::rhoPolymerEOS(const dictionary& dict)
:
    Specie(dict),
    b1m_(readScalar(dict.subDict("equationOfState").lookup("b1m"))),
    b2m_(readScalar(dict.subDict("equationOfState").lookup("b2m"))),
    b3m_(readScalar(dict.subDict("equationOfState").lookup("b3m"))),
    b4m_(readScalar(dict.subDict("equationOfState").lookup("b4m"))),
    b1s_(readScalar(dict.subDict("equationOfState").lookup("b1s"))),
    b2s_(readScalar(dict.subDict("equationOfState").lookup("b2s"))),
    b3s_(readScalar(dict.subDict("equationOfState").lookup("b3s"))),
    b4s_(readScalar(dict.subDict("equationOfState").lookup("b4s"))),
    b5_(readScalar(dict.subDict("equationOfState").lookup("b5"))),
    b6_(readScalar(dict.subDict("equationOfState").lookup("b6")))

{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
void Foam::rhoPolymerEOS<Specie>::write(Ostream& os) const
{
    Specie::write(os);

    dictionary dict("equationOfState");
    dict.add("b1m", b1m_);
    dict.add("b2m", b2m_);
    dict.add("b3m", b3m_);
    dict.add("b4m", b4m_);
    dict.add("b1s", b1s_);
    dict.add("b2s", b2s_);
    dict.add("b3s", b3s_);
    dict.add("b4s", b4s_);
    dict.add("b5", b5_);
    dict.add("b6", b6_);



    os  << indent << dict.dictName() << dict;
}

// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Specie>
Foam::Ostream& Foam::operator<<(Ostream& os, const rhoPolymerEOS<Specie>& pg)
{
    os  << static_cast<const Specie&>(pg)
        << token::SPACE << pg.b1m_
        << token::SPACE << pg.b2m_
        << token::SPACE << pg.b3m_
        << token::SPACE << pg.b4m_
        << token::SPACE << pg.b1s_
        << token::SPACE << pg.b2s_
        << token::SPACE << pg.b3s_
        << token::SPACE << pg.b4s_
        << token::SPACE << pg.b5_
        << token::SPACE << pg.b6_;


    os.check("Ostream& operator<<(Ostream&, const rhoPolymerEOS<Specie>&)");
    return os;
}


// ************************************************************************* //
