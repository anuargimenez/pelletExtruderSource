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

Class
    Foam::rhoPolymerEOS

Description
    Perfect gas equation of state.

SourceFiles
    rhoPolymerEOSI.H
    rhoPolymerEOS.C

\*---------------------------------------------------------------------------*/

#ifndef rhoPolymerEOS_H
#define rhoPolymerEOS_H

#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class rhoPolymerEOS;

template<class Specie>
inline rhoPolymerEOS<Specie> operator+
(
    const rhoPolymerEOS<Specie>&,
    const rhoPolymerEOS<Specie>&
);

template<class Specie>
inline rhoPolymerEOS<Specie> operator*
(
    const scalar,
    const rhoPolymerEOS<Specie>&
);

template<class Specie>
inline rhoPolymerEOS<Specie> operator==
(
    const rhoPolymerEOS<Specie>&,
    const rhoPolymerEOS<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const rhoPolymerEOS<Specie>&
);


/*---------------------------------------------------------------------------*\
                           Class rhoPolymerEOS Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class rhoPolymerEOS
:
    public Specie
{
    // Private data

        //- Tait coefficients
        scalar b1m_;
        scalar b2m_;
        scalar b3m_;
        scalar b4m_;
        scalar b1s_;
        scalar b2s_;
        scalar b3s_;
        scalar b4s_;
        scalar b5_;
        scalar b6_;
public:

    // Constructors

	//- Construct from components
        inline rhoPolymerEOS
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

        );
        
        //- Construct from Istream
        rhoPolymerEOS(Istream&);
        
        //- Construct from dictionary
        rhoPolymerEOS(const dictionary& dict);

        //- Construct as named copy
        inline rhoPolymerEOS(const word& name, const rhoPolymerEOS&);

        //- Construct and return a clone
        inline autoPtr<rhoPolymerEOS> clone() const;

        // Selector from dictionary
        inline static autoPtr<rhoPolymerEOS> New(const dictionary& dict);
        
        // Selector from Istream
        inline static autoPtr<rhoPolymerEOS> New(Istream& is);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "rhoPolymerEOS<" + word(Specie::typeName_()) + '>';
        }


        // Fundamental properties

            //- Is the equation of state is incompressible i.e. rho != f(p)
            static const bool incompressible = false;

            //- Is the equation of state is isochoric i.e. rho = const
            static const bool isochoric = false;

            //- Return density [kg/m^3]
            inline scalar rho(scalar p, scalar T) const;

            //- Return enthalpy departure [J/kg]
            inline scalar H(const scalar p, const scalar T) const;

            //- Return Cp departure [J/(kg K]
            inline scalar Cp(scalar p, scalar T) const;

            //- Return internal energy departure [J/kg]
            inline scalar E(const scalar p, const scalar T) const;

            //- Return Cv departure [J/(kg K]
            inline scalar Cv(scalar p, scalar T) const;

            //- Return entropy [J/kg/K]
            inline scalar S(const scalar p, const scalar T) const;

            //- Return compressibility rho/p [s^2/m^2]
            inline scalar psi(scalar p, scalar T) const;

            //- Return compression factor []
            inline scalar Z(scalar p, scalar T) const;

            //- Return (Cp - Cv) [J/(kg K]
            inline scalar CpMCv(scalar p, scalar T) const;


        // IO

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member Operators

        inline void operator+=(const rhoPolymerEOS&);
        inline void operator*=(const scalar);


    // Friend operators

        friend rhoPolymerEOS operator+ <Specie>
        (
            const rhoPolymerEOS&,
            const rhoPolymerEOS&
        );

        friend rhoPolymerEOS operator* <Specie>
        (
            const scalar s,
            const rhoPolymerEOS&
        );

        friend rhoPolymerEOS operator== <Specie>
        (
            const rhoPolymerEOS&,
            const rhoPolymerEOS&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const rhoPolymerEOS&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "rhoPolymerEOSI.H"

#ifdef NoRepository
    #include "rhoPolymerEOS.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
