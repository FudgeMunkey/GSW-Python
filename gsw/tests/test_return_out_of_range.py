from gsw import *
from hypothesis import given, example, settings, strategies as st


# STRATEGIES

# Salinity
global_salinity_st = st.floats(min_value=0, max_value=42, allow_nan=False)
absolute_salinity_st = global_salinity_st  # Paper [0, 42]

# Temperature
global_temperature_st = st.floats(min_value=-2, max_value=40, allow_nan=False)
global_temperature_ice_st = st.floats(min_value=-20, max_value=0, allow_nan=False)
in_situ_temperature_st = global_temperature_st  # Paper? [-2, 40]
conservative_temperature_st = global_temperature_st  # Paper? [-2, 40]

in_situ_temperature_of_ice_st = global_temperature_ice_st  # Guess [-20, -2]

# Pressure
global_pressure_st = st.floats(min_value=0, max_value=11000, allow_nan=False)
sea_pressure_st = global_pressure_st  # Paper [0, 11000]

# Fractions
global_fractions_st = st.floats(min_value=0, max_value=1, allow_nan=False)
saturation_fraction_st = global_fractions_st  # Documentation [0, 1]
mass_fraction_of_ice_st = global_fractions_st  # Documentation [0, 1]

# Others
conductivity_st = st.floats(min_value=3000, max_value=6000, allow_nan=False)  # Paper

# VALID OUTPUT RANGES
practical_salinity_al_range = (0, 500)  # Gaet'ale Pond has highest salinity of 433g/kg
absolute_salinity_final_al_range = practical_salinity_al_range
absolute_salinity_al_range = practical_salinity_al_range
conservative_temperature_al_range = (
    -10,
    60,
)  # Lowest is -2 ish, highest is red sea 54 ish
conservative_temperature_final_al_range = conservative_temperature_al_range
w_ih_final_al_range = None
potential_temperature_ice_al_range = (-30, 10)  # Glaciers can be like -20 ish?


# Utility functions


def check_value_exists(value):
    """Check a single value is not nan/none"""

    if not np.isscalar(value):
        # Should not be checking against list/array/etc
        raise Exception("Value is not scalar")

    if value == None:
        raise Exception("Value is None")
    elif np.isnan(value):
        raise Exception("Value is np.nan")


def check_value_in_range(value, range):
    """Check value is in the specified range"""

    # TODO: Remove this check
    if np.isnan(value):
        # print("Value is nan - skip...")
        return

    if value is None:
        print(f"VALUE AHH: {value}")
        raise Exception("You did not include a value...")

    if not range:
        # TODO: Reraise the error if no range is included
        return
        # raise Exception("You did not include a range...")

    minimum, maximum = range

    if minimum is not None:
        assert minimum <= value

    if maximum:
        assert value <= maximum


@given(
    conductivity_st,
    in_situ_temperature_st,
    sea_pressure_st,
)
def test_SP_from_C(C, t, p):
    practical_salinity_al = SP_from_C(C, t, p)

    check_value_exists(practical_salinity_al)
    check_value_in_range(practical_salinity_al, practical_salinity_al_range)


@given(
    absolute_salinity_st,
    conservative_temperature_st,
    sea_pressure_st,
    mass_fraction_of_ice_st,
    in_situ_temperature_of_ice_st,
)
def test_melting_ice_into_seawater(SA, CT, p, w_Ih, t_Ih):
    (
        absolute_salinity_final_al,
        conservative_temperature_final_al,
        w_ih_final_al,
    ) = melting_ice_into_seawater(SA, CT, p, w_Ih, t_Ih)

    check_value_exists(absolute_salinity_final_al)
    check_value_in_range(absolute_salinity_final_al, absolute_salinity_final_al_range)
    check_value_exists(conservative_temperature_final_al)
    check_value_in_range(
        conservative_temperature_final_al,
        conservative_temperature_final_al_range,
    )
    check_value_exists(w_ih_final_al)
    check_value_in_range(w_ih_final_al, w_ih_final_al_range)


@given(
    in_situ_temperature_st,
    sea_pressure_st,
)
def test_pt0_from_t_ice(t, p):
    potential_temperature_ice_al = pt0_from_t_ice(t, p)

    check_value_exists(potential_temperature_ice_al)
    check_value_in_range(
        potential_temperature_ice_al, potential_temperature_ice_al_range
    )


@given(
    conservative_temperature_st,
    sea_pressure_st,
    saturation_fraction_st,
)
def test_SA_freezing_from_CT(CT, p, saturation_fraction):
    absolute_salinity_al = SA_freezing_from_CT(CT, p, saturation_fraction)

    check_value_exists(absolute_salinity_al)
    check_value_in_range(absolute_salinity_al, absolute_salinity_al_range)


@given(
    conservative_temperature_st,
    sea_pressure_st,
    saturation_fraction_st,
)
def test_SA_freezing_from_CT_poly(CT, p, saturation_fraction):
    absolute_salinity_al = SA_freezing_from_CT_poly(CT, p, saturation_fraction)

    check_value_exists(absolute_salinity_al)
    check_value_in_range(absolute_salinity_al, absolute_salinity_al_range)


@given(
    in_situ_temperature_st,
    sea_pressure_st,
    saturation_fraction_st,
)
def test_SA_freezing_from_t(t, p, saturation_fraction):
    absolute_salinity_al = SA_freezing_from_t(t, p, saturation_fraction)

    check_value_exists(absolute_salinity_al)
    check_value_in_range(absolute_salinity_al, absolute_salinity_al_range)
