import astropy.units.physical
import numpy as np
from astropy.nddata import NDData
import astropy.units as u

class Body:
    mass = 0.0

class SphericalBody(Body):
    radius = 0.0

    def __init__(self, mass, radius):
        self.mass = mass
        self.radius = radius

class Orbit:
    first_body: Body = None
    second_body: Body = None

    def __init__(self, first_body: Body, second_body: Body):
        self.first_body = first_body
        self.second_body = second_body


ANGULAR_MOMENTUM_UNIT: u.Unit = u.kg * u.meter ** 2 / u.second
ANGULAR_MOMENTUM: u.PhysicalType = u.get_physical_type(ANGULAR_MOMENTUM_UNIT)
STANDARD_GRAVITATIONAL_PARAMETER_UNIT: u.Unit = u.meter ** 3 * u.second ** -2
class EllipticOrbit(Orbit):

    total_angular_momentum: NDData = NDData(np.zeros(3) * ANGULAR_MOMENTUM_UNIT)
    total_energy: u.Quantity = u.Quantity(0.0, unit=u.joule)
    semi_major_axis: u.Quantity = u.Quantity(0.0, unit=u.meter)
    eccentricity: u.Quantity = u.Quantity(0.0, unit=u.dimensionless_unscaled)
    semi_lateral_rectum: u.Quantity = u.Quantity(0.0, unit=u.meter)
    gravitational_factor: u.Quantity = u.Quantity(0.0, unit=STANDARD_GRAVITATIONAL_PARAMETER_UNIT)