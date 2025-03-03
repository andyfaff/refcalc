#!/usr/bin/python
from __future__ import division
import numpy as np

# import scipy.stats as ss
# import scipy.integrate as spi


def div(d1, d2, L12=2859):
    """
    returns the angular resolutions given a set of collimation conditions
    Parameters:

    d1 - slit 1 opening
    d2 - slit 2 opening
    L12 - distance between slits

    Returns:
        (dtheta, alpha, beta)


    dtheta is the FWHM of the Gaussian approximation to the trapezoidal resolution function
    alpha is the angular divergence of the penumbra
    beta is the angular divergence of the umbra
    When calculating dtheta/theta values for resolution, then dtheta is the value you need to use.
    See equations 11-14 in:
    de Haan, V.-O.; de Blois, J.; van der Ende, P.; Fredrikze, H.; van der Graaf, A.; Schipper, M.; van Well, A. A. & J., v. d. Z. ROG, the neutron reflectometer at IRI Delft Nuclear Instruments and Methods in Physics Research A, 1995, 362, 434-453

    """

    divergence = 0.68 * 0.68 * (d1 * d1 + d2 * d2) / L12 / L12
    alpha = (d1 + d2) / 2.0 / L12
    beta = abs(d1 - d2) / 2.0 / L12
    return np.degrees(np.sqrt(divergence)), np.degrees(alpha), np.degrees(beta)


def qcalc(angle, wavelength):
    """
    calculate Q given angle of incidence and wavelength
    angle - angle of incidence (degrees)
    wavelength -  wavelength of radiation (Angstrom)
    """
    return 4 * np.pi * np.sin(angle * np.pi / 180) / wavelength


def qcalc2(omega, twotheta, phi, wavelength):
    """
    convert angles and wavelength (lambda) to Q vector.
    coordinate system:
    y - along beam direction (in small angle approximation)
    x - transverse to beam direction, in plane of sample
    z - normal to sample plane.
    the xy plane is equivalent to the sample plane.

    angle definitions
    omega - angle of incidence to sample (in xy plane)
    twotheta - angle between xy plane and reflected beam PLUS omega.
    phi - angle between reflected beam and yz plane.
    returns a cartesian vector (Qx, Qy, Qz)
    """

    # convert to radians
    omega = np.radians(omega)
    twotheta = np.radians(twotheta)
    xsi = np.radians(xsi)
    # print "twotheta - omega=", cos(twotheta - omega) * cos(xsi), cos(omega)
    qx = 2 * np.pi / wavelength * np.cos(twotheta - omega) * np.sin(phi)
    qy = (
        2
        * np.pi
        / wavelength
        * (np.cos(twotheta - omega) * np.cos(phi) - np.cos(omega))
    )
    qz = 2 * np.pi / wavelength * (np.sin(twotheta - omega) + np.sin(omega))

    return (qx, qy, qz)


def wavelength(q, angle):
    """
    calculate wavelength given Q vector and angle
    q - wavevector (A^-1)
    angle - angle of incidence (degrees)
    """
    return 4.0 * np.pi * np.sin(angle * np.pi / 180.0) / q


def angle(q, wavelength):
    """
    calculate angle given Q and wavelength
    q - wavevector (A^-1)
    wavelength -  wavelength of radiation (Angstrom)
    """
    return np.asin(q / 4.0 / np.pi * wavelength) * 180 / np.pi


def qcrit(SLD1, SLD2):
    """
    calculate critical Q vector given SLD of super and subphases
    SLD1 - SLD of superphase (10^-6 A^-2)
    SLD2 - SLD of subphase (10^-6 A^-2)
    """
    return np.sqrt(16.0 * np.pi * (SLD2 - SLD1) * 1.0e-6)


def xraylam(energy):
    """
    convert energy (keV) to wavelength (angstrom)
    """
    return 12.398 / energy


def xrayenergy(wavelength):
    """
    convert energy (keV) to wavelength (angstrom)
    """
    return 12.398 / wavelength


# def beamfrac(FWHM, length, angle):
#     '''
#     return the beam fraction intercepted by a sample of length length
#     at sample tilt angle.
#     The beam is assumed to be gaussian, with a FWHM of FWHM.
#     '''
#     height_of_sample = length * np.sin(np.radians(angle))
#     beam_sd = FWHM / 2 / np.sqrt(2 * np.log(2))
#     probability = 2. * (ss.norm.cdf(height_of_sample / 2. / beam_sd) - 0.5)
#     return probability
#
# def beamfrackernel(kernelx, kernely, length, angle):
#     '''
#     return the beam fraction intercepted by a sample of length length
#     at sample tilt angle.
#     The beam has the shape 'kernel', a 2 row array, which gives the PDF for the beam
#     intensity as a function of height. The first row is position, the second row is
#     probability at that position
#     '''
#     height_of_sample = length * np.sin(np.radians(angle))
#     total = spi.simps(kernely, kernelx)
#     lowlimit = np.where(-height_of_sample / 2. >= kernelx)[0][-1]
#     hilimit = np.where(height_of_sample / 2. <= kernelx)[0][0]
#
#     area = spi.simps(kernely[lowlimit: hilimit + 1], kernelx[lowlimit: hilimit + 1])
#     return area / total
