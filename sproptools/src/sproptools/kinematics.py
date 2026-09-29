#! /usr/bin/env python3

"""
Conversions between different kinematic quantities (momentum, kinetic energy, beta, etc.)
and of differential cosmic-ray fluxes given as a function of one quantity to another.

Example:

    from sproptools.kinematics import Species, KinematicQuantity

    proton = Species.from_name('proton')
    beta = 0.999
    ekn = proton.convert(beta, KinematicQuantity.Beta, KinematicQuantity.KineticEnergyPerNucleon)


"""

from dataclasses import dataclass, astuple
from enum import Enum

import numpy as np

class KinematicQuantity(Enum):
    Momentum = 1
    Rigidity = 2
    Sagitta = 3
    AbsSagitta = 4
    Energy = 5
    KineticEnergy = 6
    EnergyPerNucleon = 7
    KineticEnergyPerNucleon = 8
    Beta = 9
    Gamma = 10
    BetaGamma = 11


def MomentumToEnergy(p, mass):

    return np.sqrt(p * p + mass * mass)


def MomentumToKineticEnergy(p, mass):

    return MomentumToEnergy(p, mass) - mass


def MomentumToRigidity(p, charge):

    if np.allclose(charge, 0.):
        raise ValueError('Charge is zero!')
    return p / np.abs(charge)


def MomentumToAbsSagitta(p, charge):

    R = MomentumToRigidity(p, charge)
    if np.allclose(R, 0.):
        raise ValueError('Rigidity is zero!')

    return 1.0 / R


def MomentumToSagitta(p, charge):

    absS = MomentumToAbsSagitta(p, charge)
    return absS if charge >= 0.0 else -absS


def MomentumToEnergyPerNucleon(p, mass, A):

    if A == 0:
        raise ValueError(f'No atomic mass number given for particle with m={mass:g} GeV in conversion to energy per nucleon.')

    return MomentumToEnergy(p, mass) / A


def MomentumToKineticEnergyPerNucleon(p, mass, A):

    if A == 0:
        raise ValueError(f'No atomic mass number given for particle with m={mass:g} GeV in conversion to Ekn.')

    return MomentumToKineticEnergy(p, mass) / A


def MomentumToBeta(p, mass):

    E = MomentumToEnergy(p, mass)
    if np.allclose(E, 0.):
        raise ValueError('Energy is zero!')

    return p / E


def MomentumToGamma(p, mass):

    if np.allclose(mass, 0.):
        raise ValueError('Mass is zero!')

    return np.sqrt(1.0 + p * p / (mass * mass))


def MomentumToBetaGamma(p, mass):

    if np.allclose(mass, 0.):
        raise ValueError('Mass is zero!')

    return p / mass


def EnergyToMomentum(E, mass):

    #if E < mass:
    #    raise ValueError(f'Invalid energy (lower than the particle mass): E = {E:g}, m = {mass:g}')
    return np.sqrt(E * E - mass * mass)


def KineticEnergyToMomentum(T, mass):

    return EnergyToMomentum(T + mass, mass)


def RigidityToMomentum(R, charge):

    if np.allclose(charge, 0.):
        raise ValueError('Charge is zero!')

    return R * np.abs(charge)


def AbsSagittaToMomentum(S, charge):

    if np.allclose(S, 0.):
        raise ValueError('Sagitta is zero!')

    return RigidityToMomentum(1.0 / S, charge)


def SagittaToMomentum(S, charge):

    return AbsSagittaToMomentum(np.abs(S), charge)


def EnergyPerNucleonToMomentum(EoverA, mass, A):

    if A == 0:
        raise ValueError(f'No atomic mass number given for particle with m={mass:g} GeV in conversion from energy per nucleon.')

    E = EoverA * A
    return EnergyToMomentum(E, mass)


def KineticEnergyPerNucleonToMomentum(ToverA, mass, A):

    if A == 0:
        raise ValueError(f'No atomic mass number given for particle with m={mass:g} GeV in conversion from Ekn.')

    T = ToverA * A
    return KineticEnergyToMomentum(T, mass)


def BetaToMomentum(beta, mass):

    if np.abs(beta) >= 1.0:
        raise ValueError(f'Invalid particle velocity, Beta = {beta:g}')

    return mass * beta / np.sqrt(1.0 - beta * beta)


def GammaToMomentum(gamma, mass):

    if np.abs(gamma) < 1.0:
        raise ValueError(f'Invalid particle Lorentz factor, Gamma = {gamma:g}')

    return mass * np.sqrt(gamma * gamma - 1.0)


def BetaGammaToMomentum(betaGamma, mass):

    return mass * betaGamma




def Convert(value, fromType, toType, mass, charge, atomicMassNumber):

    if fromType == toType:
        return value

    p = ConvertToMomentum(value, fromType, mass, charge, atomicMassNumber)
    return ConvertFromMomentum(p, toType, mass, charge, atomicMassNumber)



def ConvertToMomentum(value, fromType, mass, charge, atomicMassNumber):

    match fromType:
        case KinematicQuantity.Momentum:
            return value
        case KinematicQuantity.Rigidity:
            return RigidityToMomentum(value, charge)
        case KinematicQuantity.AbsSagitta:
            return AbsSagittaToMomentum(value, charge)
        case KinematicQuantity.Sagitta:
            return SagittaToMomentum(value, charge)
        case KinematicQuantity.Energy:
            return EnergyToMomentum(value, mass)
        case KinematicQuantity.KineticEnergy:
            return KineticEnergyToMomentum(value, mass)
        case KinematicQuantity.EnergyPerNucleon:
            return EnergyPerNucleonToMomentum(value, mass, atomicMassNumber)
        case KinematicQuantity.KineticEnergyPerNucleon:
            return KineticEnergyPerNucleonToMomentum(value, mass, atomicMassNumber)
        case KinematicQuantity.Beta:
            return BetaToMomentum(value, mass)
        case KinematicQuantity.Gamma:
            return GammaToMomentum(value, mass)
        case KinematicQuantity.BetaGamma:
            return BetaGammaToMomentum(value, mass)
        case _:
            raise ValueError('No matching kinematic variable type found.')

def ConvertFromMomentum(momentum, toType, mass, charge, atomicMassNumber):

    match toType:
        case KinematicQuantity.Momentum:
            return momentum
        case KinematicQuantity.Rigidity:
            return MomentumToRigidity(momentum, charge)
        case KinematicQuantity.AbsSagitta:
            return MomentumToAbsSagitta(momentum, charge)
        case KinematicQuantity.Sagitta:
            return MomentumToSagitta(momentum, charge)
        case KinematicQuantity.Energy:
            return MomentumToEnergy(momentum, mass)
        case KinematicQuantity.KineticEnergy:
            return MomentumToKineticEnergy(momentum, mass)
        case KinematicQuantity.EnergyPerNucleon:
            return MomentumToEnergyPerNucleon(momentum, mass, atomicMassNumber)
        case KinematicQuantity.KineticEnergyPerNucleon:
            return MomentumToKineticEnergyPerNucleon(momentum, mass, atomicMassNumber)
        case KinematicQuantity.Beta:
            return MomentumToBeta(momentum, mass)
        case KinematicQuantity.Gamma:
            return MomentumToGamma(momentum, mass)
        case KinematicQuantity.BetaGamma:
            return MomentumToBetaGamma(momentum, mass)
        case _:
            raise ValueError('No matching kinematic variable type found.')


def Derivative(x, fromType, toType, mass, charge, atomicMassNumber):

    if fromType == toType:
        return 1.0

    p = ConvertToMomentum(x, fromType, mass, charge, atomicMassNumber)
    b = ConvertFromMomentum(p, KinematicQuantity.Beta, mass, charge, atomicMassNumber)
    E = ConvertFromMomentum(p, KinematicQuantity.Energy, mass, charge, atomicMassNumber)
    c = np.abs(charge)
    m = np.abs(mass)

    A = atomicMassNumber
    if A == 0:
        A = 1

    # compute d(fromType)/dT
    dFromdT = 0.0
    match fromType:
        case KinematicQuantity.Energy | KinematicQuantity.KineticEnergy:
            dFromdT = 1.0
        case KinematicQuantity.Momentum:
            dFromdT = E / p
        case KinematicQuantity.Rigidity:
            dFromdT = 1.0 / c * E / p
        case KinematicQuantity.AbsSagitta:
            dFromdT = c * E / (p * p * p)
        case KinematicQuantity.Sagitta:
            dFromdT = (1.0 if c >= 0.0 else -1.0) * c * E / (p * p * p)
        case KinematicQuantity.EnergyPerNucleon | KinematicQuantity.KineticEnergyPerNucleon:
            dFromdT = 1.0 / A
        case KinematicQuantity.Beta:
            dFromdT = (1.0 - b * b) / (E * b)
        case KinematicQuantity.Gamma:
            dFromdT = 1.0 / m
        case KinematicQuantity.BetaGamma:
            dFromdT = E / (p * m)
        case _:
            raise ValueError('No matching kinematic variable type found.')

    # return d(fromType)/dT * dT/d(toType) = d(fromType)/d(toType)
    match toType:
        case KinematicQuantity.Energy | KinematicQuantity.KineticEnergy:
            return dFromdT
        case KinematicQuantity.Momentum:
            return dFromdT * (p / E)
        case KinematicQuantity.Rigidity:
            return dFromdT * (p / E) * c
        case KinematicQuantity.AbsSagitta:
            return dFromdT * (p * p * p / (E * c))
        case KinematicQuantity.Sagitta:
            return dFromdT * (1.0 if c >= 0.0 else -1.0) * (p * p * p / (E * c))
        case KinematicQuantity.EnergyPerNucleon | KinematicQuantity.KineticEnergyPerNucleon:
            return dFromdT * A
        case KinematicQuantity.Beta:
            return dFromdT * E * b / (1.0 - b * b)
        case KinematicQuantity.Gamma:
            return dFromdT * m
        case KinematicQuantity.BetaGamma:
            return dFromdT * (p * m) / E
        case _:
            raise ValueError('No matching kinematic variable type found.')


def ConvertDifferentialFluxGraph(x, flux, flux_err, fromType, toType, mass, charge, atomicMassNumber):
    """
    Convert graph representing a differential flux, given as a function of one variable to another.

    The proper kinematic derivative will be used in the conversion such that the integral number of events for any given bin is preserved.
    Error bars are scaled such that the relative flux uncertainty is preserved.

    Returns x-values and fluxes if flux_err is None. If flux errors are given, returns x-values, fluxes, and rescaled flux errors.
    """

    x_conv = Convert(x, fromType, toType, mass, charge, atomicMassNumber)
    diff = Derivative(x, fromType, toType, mass, charge, atomicMassNumber)
    flux_conv = diff * flux

    if flux_err is not None:
        flux_err_conv = diff * flux_err
        return x_conv, flux_conv, flux_err_conv

    return x_conv, flux_conv


@dataclass(frozen=True)
class Species:
    mass: float
    charge: int
    atomicMassNumber: int
    name: str

    @classmethod
    def from_name(cls, species):
        match species:
            case 'electron':
                return cls(0.00051099906, -1, 0, species)
            case 'positron':
                return cls(0.00051099906, 1, 0, species)
            case 'proton':
                return cls(0.93827231, 1, 1, species)
            case 'antiproton':
                return cls(0.93827231, -1, 1, species)
            case _:
                raise ValueError(f'Unknown particle species {species}!')
        return cls(0.0, 0.0, 0)

    def convert(self, value, fromType, toType):

        if fromType == toType:
            return value

        p = ConvertToMomentum(value, fromType, *astuple(self)[0:3])
        return ConvertFromMomentum(p, toType, *astuple(self)[0:3])

    def convert_differential_flux_graph(self, x, flux, flux_err, fromType, toType):
        return ConvertDifferentialFluxGraph(x, flux, flux_err, fromType, toType, *astuple(self)[0:3])


if __name__ == '__main__':

    proton = Species.from_name('proton')
    print(proton.convert(0.999, KinematicQuantity.Beta, KinematicQuantity.KineticEnergyPerNucleon))
    print(proton.convert(1.0, KinematicQuantity.Momentum, KinematicQuantity.Rigidity))
    print(Derivative(0.999, KinematicQuantity.Beta, KinematicQuantity.KineticEnergyPerNucleon, *astuple(proton)[0:3]))
