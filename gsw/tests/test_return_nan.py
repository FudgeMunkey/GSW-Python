from gsw import *
from hypothesis import given, example, strategies as st

# STRATEGIES

# Salinity
global_salinity_st = st.floats(min_value=0, max_value=42, allow_nan=False)
absolute_salinity_st = global_salinity_st  # Paper [0, 42]
practical_salinity_st = global_salinity_st  # Paper [0, 42]
knudsen_salinity_st = global_salinity_st  # Guess
absolute_salinity_of_sea_ice_st = (
    global_salinity_st  # Guess, probably less than sea water?
)

# Temperature
global_temperature_st = st.floats(min_value=-2, max_value=40, allow_nan=False)
global_temperature_ice_st = st.floats(min_value=-20, max_value=0, allow_nan=False)
in_situ_temperature_st = global_temperature_st  # Paper? [-2, 40]
conservative_temperature_st = global_temperature_st  # Paper? [-2, 40]

in_situ_temperature_of_ice_st = global_temperature_ice_st  # Guess [-20, -2]
in_situ_temperature_of_sea_ice_at_pressure_st = (
    global_temperature_ice_st  # Guess [-20, -2]
)

# Pressure
global_pressure_st = st.floats(min_value=0, max_value=11000, allow_nan=False)
sea_pressure_st = global_pressure_st  # Paper [0, 11000]

# Fractions
global_fractions_st = st.floats(min_value=0, max_value=1, allow_nan=False)
saturation_fraction_st = global_fractions_st  # Documentation [0, 1]
mass_fraction_of_ice_st = global_fractions_st  # Documentation [0, 1]
mass_fraction_of_sea_ice_st = global_fractions_st  # Documentation [0, 1]

# Others
in_situ_seawater_density_rho_st = st.floats(
    min_value=1000, max_value=1060, allow_nan=False
)  # Paper? In situ density?
longitude_st = st.floats(
    min_value=-360, max_value=360, allow_nan=False
)  # Documentation
latitude_st = st.floats(min_value=-90, max_value=90, allow_nan=False)  # Documentation


def check_value_exists(value):
    """Check a single value is not nan/none"""

    if not np.isscalar(value):
        # Should not be checking against list/array/etc
        raise Exception("Value is not scalar")

    if value == None:
        raise Exception("Value is None")
    elif np.isnan(value):
        raise Exception("Value is np.nan")


@given(
    in_situ_seawater_density_rho_st,
    absolute_salinity_st,
    sea_pressure_st,
)
def test_CT_from_rho(rho, SA, p):
    conservative_temperature_al, ct_multiple_al = CT_from_rho(rho, SA, p)

    check_value_exists(conservative_temperature_al)
    check_value_exists(ct_multiple_al)


@given(
    in_situ_seawater_density_rho_st,
    conservative_temperature_st,
    sea_pressure_st,
)
def test_SA_from_rho(rho, CT, p):
    absolute_salinity_al = SA_from_rho(rho, CT, p)

    check_value_exists(absolute_salinity_al)


@given(
    knudsen_salinity_st,
)
def test_SP_from_SK(SK):
    practical_salinity_al = SP_from_SK(SK)

    check_value_exists(practical_salinity_al)


@given(
    absolute_salinity_st,
    conservative_temperature_st,
    sea_pressure_st,
    in_situ_temperature_of_ice_st,
)
def test_ice_fraction_to_freeze_seawater(SA, CT, p, t_Ih):
    sa_freeze_al, ct_freeze_al, w_ih_al = ice_fraction_to_freeze_seawater(
        SA, CT, p, t_Ih
    )

    check_value_exists(sa_freeze_al)
    check_value_exists(ct_freeze_al)
    check_value_exists(w_ih_al)


@given(
    absolute_salinity_st,
    conservative_temperature_st,
    sea_pressure_st,
    in_situ_temperature_of_ice_st,
)
def test_melting_ice_SA_CT_ratio(SA, CT, p, t_Ih):
    melting_ice_sa_ct_ratio_al = melting_ice_SA_CT_ratio(SA, CT, p, t_Ih)

    check_value_exists(melting_ice_sa_ct_ratio_al)


@given(
    absolute_salinity_st,
    conservative_temperature_st,
    sea_pressure_st,
    in_situ_temperature_of_ice_st,
)
def test_melting_ice_SA_CT_ratio_poly(SA, CT, p, t_Ih):
    melting_ice_sa_ct_ratio_al = melting_ice_SA_CT_ratio_poly(SA, CT, p, t_Ih)

    check_value_exists(melting_ice_sa_ct_ratio_al)


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
    check_value_exists(conservative_temperature_final_al)
    check_value_exists(w_ih_final_al)


@given(
    absolute_salinity_st,
    conservative_temperature_st,
    sea_pressure_st,
    absolute_salinity_of_sea_ice_st,
    in_situ_temperature_of_sea_ice_at_pressure_st,
)
def test_melting_seaice_SA_CT_ratio(SA, CT, p, SA_seaice, t_seaice):
    melting_seaice_sa_ct_ratio_al = melting_seaice_SA_CT_ratio(
        SA, CT, p, SA_seaice, t_seaice
    )

    check_value_exists(melting_seaice_sa_ct_ratio_al)


@given(
    absolute_salinity_st,
    conservative_temperature_st,
    sea_pressure_st,
    absolute_salinity_of_sea_ice_st,
    in_situ_temperature_of_sea_ice_at_pressure_st,
)
def test_melting_seaice_SA_CT_ratio_poly(SA, CT, p, SA_seaice, t_seaice):
    melting_seaice_sa_ct_ratio_al = melting_seaice_SA_CT_ratio_poly(
        SA, CT, p, SA_seaice, t_seaice
    )

    check_value_exists(melting_seaice_sa_ct_ratio_al)


@given(
    absolute_salinity_st,
    conservative_temperature_st,
    sea_pressure_st,
    mass_fraction_of_sea_ice_st,
    absolute_salinity_of_sea_ice_st,
    in_situ_temperature_of_sea_ice_at_pressure_st,
)
def test_melting_seaice_into_seawater(SA, CT, p, w_seaice, SA_seaice, t_seaice):
    (
        absolute_salinity_final_al,
        conservative_temperature_final_al,
    ) = melting_seaice_into_seawater(SA, CT, p, w_seaice, SA_seaice, t_seaice)

    check_value_exists(absolute_salinity_final_al)
    check_value_exists(conservative_temperature_final_al)


@given(
    absolute_salinity_st,
    conservative_temperature_st,
    sea_pressure_st,
    absolute_salinity_of_sea_ice_st,
    in_situ_temperature_of_sea_ice_at_pressure_st,
)
def test_seaice_fraction_to_freeze_seawater(SA, CT, p, SA_seaice, t_seaice):
    (
        sa_freeze_al,
        ct_freeze_al,
        w_seaice_al,
    ) = seaice_fraction_to_freeze_seawater(SA, CT, p, SA_seaice, t_seaice)

    check_value_exists(sa_freeze_al)
    check_value_exists(ct_freeze_al)
    check_value_exists(w_seaice_al)


@given(
    conservative_temperature_st,
    sea_pressure_st,
    saturation_fraction_st,
)
def test_SA_freezing_from_CT(CT, p, saturation_fraction):
    absolute_salinity_al = SA_freezing_from_CT(CT, p, saturation_fraction)

    check_value_exists(absolute_salinity_al)


@given(
    conservative_temperature_st,
    sea_pressure_st,
    saturation_fraction_st,
)
def test_SA_freezing_from_CT_poly(CT, p, saturation_fraction):
    absolute_salinity_al = SA_freezing_from_CT_poly(CT, p, saturation_fraction)

    check_value_exists(absolute_salinity_al)


@given(
    in_situ_temperature_st,
    sea_pressure_st,
    saturation_fraction_st,
)
def test_SA_freezing_from_t(t, p, saturation_fraction):
    absolute_salinity_al = SA_freezing_from_t(t, p, saturation_fraction)

    check_value_exists(absolute_salinity_al)


@given(
    in_situ_temperature_st,
    sea_pressure_st,
    saturation_fraction_st,
)
def test_SA_freezing_from_t_poly(t, p, saturation_fraction):
    absolute_salinity_al = SA_freezing_from_t_poly(t, p, saturation_fraction)

    check_value_exists(absolute_salinity_al)


@given(
    absolute_salinity_st,
    conservative_temperature_st,
    saturation_fraction_st,
)
def test_pressure_freezing_CT(SA, CT, saturation_fraction):
    pressure_freezing_al = pressure_freezing_CT(SA, CT, saturation_fraction)

    check_value_exists(pressure_freezing_al)


@given(
    practical_salinity_st,
    longitude_st,
    latitude_st,
)
def test_SA_from_SP_Baltic(SP, lon, lat):
    sa_baltic_al = SA_from_SP_Baltic(SP, lon, lat)

    check_value_exists(sa_baltic_al)


@given(
    absolute_salinity_st,
    longitude_st,
    latitude_st,
)
def test_SP_from_SA_Baltic(SA, lon, lat):
    sp_baltic_al = SP_from_SA_Baltic(SA, lon, lat)

    check_value_exists(sp_baltic_al)