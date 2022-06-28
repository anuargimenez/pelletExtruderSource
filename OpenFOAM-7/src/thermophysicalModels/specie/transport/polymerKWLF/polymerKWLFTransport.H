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

Class
    Foam::polymerKWLFTransport

Description

    Modifiziertes WLF-Modell basierend auf dem Ansatz nach Gregor Karrenberg mit eigener Modellierung
    Siehe schriftliche Ausarbeitung Masterarbeit Tim Ebinger für die entsprechende Gleichung

    References:
    \verbatim
        Williams, M. L., Landel, R. F., & Ferry, J. D. (1955).
        The temperature dependence of relaxation mechanisms
        in amorphous polymers and other glass-forming liquids.
        Journal of the American Chemical society, 77(14), 3701-3707.
    \endverbatim

SourceFiles
    polymerKWLFTransportI.H
    polymerKWLFTransport.C

\*---------------------------------------------------------------------------*/

#ifndef polymerKWLFTransport_H
#define polymerKWLFTransport_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo> class polymerKWLFTransport;

template<class Thermo>
inline polymerKWLFTransport<Thermo> operator+
(
    const polymerKWLFTransport<Thermo>&,
    const polymerKWLFTransport<Thermo>&
);

template<class Thermo>
inline polymerKWLFTransport<Thermo> operator*
(
    const scalar,
    const polymerKWLFTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const polymerKWLFTransport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                    Class polymerKWLFTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class polymerKWLFTransport
:
    public Thermo
{
    // Private Data

        //- Dynamic viscosity at the reference temperature [Pa.s]
        scalar mu0_;

        //- Reference temperature [T]
        scalar Tr_;

        //-  Melt-Fraction coefficient for Karrenberg MF modeling Tm
        scalar Tm_;

        //-  Melt-Fraction coefficient for Karrenberg MF modeling TS
        scalar TS_;

        //- Melt-Fraction coefficient for Karrenberg MF modeling MS
        scalar MS_;

        //- WLF coefficient 1 []
        scalar C1_;

        //- WLF coefficient 2 [T]
        scalar C2_;



    // Private Member Functions

        //- Read coefficient from dictionary
        scalar readCoeff(const word& coeffName, const dictionary& dict);


public:

    // Constructors

        //- Construct as named copy
        inline polymerKWLFTransport(const word&, const polymerKWLFTransport&);

        //- Construct from dictionary
        polymerKWLFTransport(const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<polymerKWLFTransport> clone() const;

        // Selector from dictionary
        inline static autoPtr<polymerKWLFTransport> New(const dictionary& dict);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "polymerKWLF<" + Thermo::typeName() + '>';
        }

        //- Dynamic viscosity [kg/m/s]
        inline scalar mu(const scalar p, const scalar T) const;

        //- Melt-Fraction
        inline scalar MF(const scalar p, const scalar T) const;

        //- Thermal conductivity [W/m/K]
        inline scalar kappa(const scalar p, const scalar T) const;

        //- Thermal diffusivity of enthalpy [kg/m/s]
        inline scalar alphah(const scalar p, const scalar T) const;

        // Species diffusivity
        // inline scalar D(const scalar p, const scalar T) const;

        //- Write to Ostream
        void write(Ostream& os) const;


    // Member Operators

        inline void operator+=(const polymerKWLFTransport&);
        inline void operator*=(const scalar);


    // Friend operators

        friend polymerKWLFTransport operator+ <Thermo>
        (
            const polymerKWLFTransport&,
            const polymerKWLFTransport&
        );

        friend polymerKWLFTransport operator* <Thermo>
        (
            const scalar,
            const polymerKWLFTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const polymerKWLFTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "polymerKWLFTransportI.H"

#ifdef NoRepository
    #include "polymerKWLFTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
