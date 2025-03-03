#!/usr/bin/python
import numpy as np
from . import fminbound
import sys

"""
Script to estimate optimum slit settings.

For a given footprint and angle combination there is a maximum angular resolution that can be achieved whilst keeping the slits equal in size.

Alternatively, for a given angle and a given angular resolution find the optimum slit settings to maximise the footprint/intensity on the sample

See equations 11-14 in:
de Haan, V.-O.; de Blois, J.; van der Ende, P.; Fredrikze, H.; van der Graaf, A.; Schipper, M.; van Well, A. A. & J., v. d. Z. ROG, the neutron reflectometer at IRI Delft Nuclear Instruments and Methods in Physics Research A, 1995, 362, 434-453

Andrew Nelson - 2013

"""


def height_of_beam_after_dx(d1, d2, L12, distance):
    """
    Calculate the widths of beam a given distance away from a collimation slit.

    if distance >= 0, then it's taken to be the distance after d2.
    if distance < 0, then it's taken to be the distance before d1.

    Parameters:
        d1 - opening of first collimation slit
        d2 - opening of second collimation slit
        L12 - distance between first and second collimation slits
        distance - distance from first or last slit to a given position
    Units - equivalent distances (inches, mm, light years)

    Returns:
        (umbra, penumbra)

    """

    alpha = (d1 + d2) / 2.0 / L12
    beta = abs(d1 - d2) / 2.0 / L12
    if distance >= 0:
        return (beta * distance * 2) + d2, (alpha * distance * 2) + d2
    else:
        return (beta * abs(distance) * 2) + d1, (alpha * abs(distance) * 2) + d1


def actual_footprint(d1, d2, L12, L2S, angle):
    """
    Calculate the actual footprint on a sample.
    Parameters:
        d1 - opening of first collimation slit
        d2 - opening of second collimation slit
        L12 - distance between first and second collimation slits
        L2S - distance from second collimation slit to sample

    Returns:
        (umbra_footprint, penumbra_footprint)

    """
    umbra, penumbra = height_of_beam_after_dx(d1, d2, L12, L2S)
    return umbra / np.radians(angle), penumbra / np.radians(angle)


def slitoptimiser(
    footprint,
    resolution,
    angle=1.0,
    L12=2859.5,
    L2S=276,
    LS4=290.5,
    LSD=2500,
    verbose=True,
):
    """
    Optimise slit settings for a given angular resolution, and a given footprint.

    footprint - maximum footprint onto sample (mm)
    resolution - fractional dtheta/theta resolution (FWHM)
    angle      - optional, angle of incidence in degrees

    #slit1-slit2 distance (mm)
    L12 = 2859.5
    #slit2 - sample distance (mm)
    L2S = 276
    #sample - S5 distance (mm)
    LS4 = 290.5
    #sample detector (mm)
    LSD = 2500
    """
    if verbose:
        print("_____________________________________________")
        print("FOOTPRINT calculator - Andrew Nelson 2013")
        print("INPUT")
        print("footprint:", footprint, "mm")
        print("fractional angular resolution (FWHM):", resolution)
        print("theta:", angle, "degrees")

    d1star = lambda d2star: np.sqrt(1 - np.power(d2star, 2))
    L1star = 0.68 * footprint / L12 / resolution

    gseekfun = lambda d2star: np.power(
        (d2star + L2S / L12 * (d2star + d1star(d2star))) - L1star, 2
    )

    res = fminbound.fminbound(gseekfun, 0, 1)

    optimal_d2star = res
    optimal_d1star = d1star(optimal_d2star)
    if optimal_d2star > optimal_d1star:
        # you found a minimum, but this may not be the optimum size of the slits.
        multfactor = 1
        optimal_d2star = 1 / np.sqrt(2)
        optimal_d1star = 1 / np.sqrt(2)
    else:
        multfactor = optimal_d2star / optimal_d1star

    d1 = optimal_d1star * resolution / 0.68 * np.radians(angle) * L12
    d2 = d1 * multfactor

    tmp, height_at_S4 = height_of_beam_after_dx(d1, d2, L12, L2S + LS4)
    tmp, height_at_detector = height_of_beam_after_dx(d1, d2, L12, L2S + LSD)
    tmp, _actual_footprint = actual_footprint(d1, d2, L12, L2S, angle)

    # if verbose:
    # print('\nOUTPUT')
    # if multfactor == 1:
    #     print('Your desired resolution results in a smaller footprint than the sample supports.'
    #     suggested_resolution =  resolution * footprint / _actual_footprint
    #     print 'You can increase flux using a resolution of', suggested_resolution, 'and still keep the same footprint.'
    # print '\nd1', d1, 'mm'
    # print 'd2', d2, 'mm'
    # print '\nfootprint:', _actual_footprint, 'mm'
    # print 'height at S4:', height_at_S4, 'mm'
    # print 'height at detector:', height_at_detector, 'mm'
    # print '\n[d2star', optimal_d2star, ']'
    # print '_____________________________________________'

    return d1, d2


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(
            "Usage:\n\
            python slitoptimiser.py footprint relative_resolution [angle]\n\
            python slitoptimiser.py 50 0.04 [2]"
        )
    else:
        slitoptimiser(*[float(val) for val in sys.argv[1:]])
