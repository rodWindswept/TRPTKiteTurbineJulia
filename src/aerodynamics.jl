# src/aerodynamics.jl
# Rotor aerodynamic coefficient tables from AeroDyn BEM simulations.
# Source: Rotor_TRTP_Sizing_Iteration2.xlsx — averaged across 4kW, 7kW, 12kW sheets.
# Airfoil: NACA4412, 3 blades, 20° elevation angle, ~2.9° pitch.
#
# cp_at_tsr(λ) — power coefficient Cp as a function of tip speed ratio λ = ω·R/v
# ct_at_tsr(λ) — thrust coefficient CT as a function of tip speed ratio λ
#
# λ=0 anchor: Cp=CT=0 (linear extrapolation from standstill).
# Table range: λ = 0.0 to 8.0 in steps of 0.1 (72 points, averaged from 3 BEM sheets).
# Beyond table: Cp extrapolated (may be negative — freewheeling); CT clamped at CT(8.0).

const BEM_TSR = [
    0.0, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
    2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9,
    3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9,
    4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9,
    5.0, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9,
    6.0, 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9,
    7.0, 7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7, 7.8, 7.9, 8.0
]

# Averaged Cp(λ) from 4kW, 7kW, 12kW BEM sheets.
# Peak Cp ≈ 0.232 at λ ≈ 4.0–4.1. Negative above λ ≈ 6.5 (freewheeling).
const BEM_CP = [
    0.0,       # λ=0.0  (anchor — linear extrapolation to standstill)
    0.014088,  # λ=1.0
    0.017791,  # λ=1.1
    0.022217,  # λ=1.2
    0.027475,  # λ=1.3
    0.033387,  # λ=1.4
    0.039731,  # λ=1.5
    0.046645,  # λ=1.6
    0.054477,  # λ=1.7
    0.063174,  # λ=1.8
    0.072588,  # λ=1.9
    0.082567,  # λ=2.0
    0.092928,  # λ=2.1
    0.103548,  # λ=2.2
    0.114320,  # λ=2.3
    0.125144,  # λ=2.4
    0.135860,  # λ=2.5
    0.146374,  # λ=2.6
    0.156601,  # λ=2.7
    0.166436,  # λ=2.8
    0.175780,  # λ=2.9
    0.184565,  # λ=3.0
    0.192771,  # λ=3.1
    0.200330,  # λ=3.2
    0.207212,  # λ=3.3
    0.213347,  # λ=3.4
    0.218725,  # λ=3.5
    0.223271,  # λ=3.6
    0.226882,  # λ=3.7
    0.229549,  # λ=3.8
    0.231245,  # λ=3.9
    0.231964,  # λ=4.0  ← near-peak
    0.231705,  # λ=4.1  ← peak Cp ≈ 0.232
    0.230432,  # λ=4.2
    0.228213,  # λ=4.3
    0.225001,  # λ=4.4
    0.220875,  # λ=4.5
    0.215827,  # λ=4.6
    0.209899,  # λ=4.7
    0.203098,  # λ=4.8
    0.195505,  # λ=4.9
    0.187135,  # λ=5.0
    0.177985,  # λ=5.1
    0.168125,  # λ=5.2
    0.157610,  # λ=5.3
    0.146471,  # λ=5.4
    0.134739,  # λ=5.5
    0.122447,  # λ=5.6
    0.109604,  # λ=5.7
    0.096227,  # λ=5.8
    0.082335,  # λ=5.9
    0.067953,  # λ=6.0
    0.053087,  # λ=6.1
    0.037777,  # λ=6.2
    0.022030,  # λ=6.3
    0.005833,  # λ=6.4
   -0.010804,  # λ=6.5  ← Cp crosses zero (freewheeling begins)
   -0.027871,  # λ=6.6
   -0.045371,  # λ=6.7
   -0.063319,  # λ=6.8
   -0.081735,  # λ=6.9
   -0.100628,  # λ=7.0
   -0.120022,  # λ=7.1
   -0.139924,  # λ=7.2
   -0.160333,  # λ=7.3
   -0.181271,  # λ=7.4
   -0.202740,  # λ=7.5
   -0.224737,  # λ=7.6
   -0.247279,  # λ=7.7
   -0.270378,  # λ=7.8
   -0.294020,  # λ=7.9
   -0.318214   # λ=8.0
]

# Averaged CT(λ) — increases monotonically with λ (higher speed → more thrust).
# CT(0)=0 (standstill); CT(4.1)≈0.548 at optimal TSR; CT(8)≈0.782.
const BEM_CT = [
    0.0,       # λ=0.0  (anchor)
    0.109112,  # λ=1.0
    0.115347,  # λ=1.1
    0.122075,  # λ=1.2
    0.129817,  # λ=1.3
    0.138704,  # λ=1.4
    0.148395,  # λ=1.5
    0.158711,  # λ=1.6
    0.170230,  # λ=1.7
    0.183324,  # λ=1.8
    0.198242,  # λ=1.9
    0.214285,  # λ=2.0
    0.230933,  # λ=2.1
    0.248135,  # λ=2.2
    0.265742,  # λ=2.3
    0.283610,  # λ=2.4
    0.301605,  # λ=2.5
    0.319587,  # λ=2.6
    0.337549,  # λ=2.7
    0.355340,  # λ=2.8
    0.372851,  # λ=2.9
    0.390012,  # λ=3.0
    0.406851,  # λ=3.1
    0.423314,  # λ=3.2
    0.439494,  # λ=3.3
    0.455257,  # λ=3.4
    0.470488,  # λ=3.5
    0.485291,  # λ=3.6
    0.499204,  # λ=3.7
    0.512387,  # λ=3.8
    0.524933,  # λ=3.9
    0.536767,  # λ=4.0
    0.548029,  # λ=4.1
    0.558593,  # λ=4.2
    0.568706,  # λ=4.3
    0.578145,  # λ=4.4
    0.587116,  # λ=4.5
    0.595494,  # λ=4.6
    0.603501,  # λ=4.7
    0.611133,  # λ=4.8
    0.618378,  # λ=4.9
    0.625253,  # λ=5.0
    0.631727,  # λ=5.1
    0.637943,  # λ=5.2
    0.643866,  # λ=5.3
    0.649590,  # λ=5.4
    0.655083,  # λ=5.5
    0.660429,  # λ=5.6
    0.665593,  # λ=5.7
    0.670665,  # λ=5.8
    0.675650,  # λ=5.9
    0.680583,  # λ=6.0
    0.685427,  # λ=6.1
    0.690242,  # λ=6.2
    0.695251,  # λ=6.3
    0.700027,  # λ=6.4
    0.705058,  # λ=6.5
    0.710135,  # λ=6.6
    0.714991,  # λ=6.7
    0.720098,  # λ=6.8
    0.725312,  # λ=6.9
    0.730379,  # λ=7.0
    0.735431,  # λ=7.1
    0.740637,  # λ=7.2
    0.745730,  # λ=7.3
    0.750807,  # λ=7.4
    0.755968,  # λ=7.5
    0.761131,  # λ=7.6
    0.766170,  # λ=7.7
    0.771345,  # λ=7.8
    0.776569,  # λ=7.9
    0.781613   # λ=8.0
]

"""
    cp_at_tsr(lambda) -> Float64

Return the rotor power coefficient Cp at tip speed ratio `lambda = ω·R/v_hub`.

Interpolated from the averaged AeroDyn BEM table (4kW/7kW/12kW NACA4412 rotors).
- Below λ=0: returns 0.0.
- Within table (0 ≤ λ ≤ 8): linear interpolation.
- Above λ=8: linear extrapolation using the last two table entries.
"""
function cp_at_tsr(lambda::Float64)::Float64
    lambda <= 0.0 && return 0.0
    return _interp_bem(BEM_TSR, BEM_CP, lambda)
end

"""
    ct_at_tsr(lambda) -> Float64

Return the rotor thrust coefficient CT at tip speed ratio `lambda = ω·R/v_hub`.

Interpolated from the averaged AeroDyn BEM table (4kW/7kW/12kW NACA4412 rotors).
- Below λ=0: returns 0.0.
- Within table (0 ≤ λ ≤ 8): linear interpolation.
- Above λ=8: clamped to CT(8) ≈ 0.782 (BEM model less reliable beyond this range).
"""
function ct_at_tsr(lambda::Float64)::Float64
    lambda <= 0.0 && return 0.0
    lambda >= BEM_TSR[end] && return BEM_CT[end]
    return _interp_bem(BEM_TSR, BEM_CT, lambda)
end

# Internal: linear interpolation over a sorted TSR table.
function _interp_bem(tsr_table::Vector{Float64}, coeff_table::Vector{Float64},
                     lambda::Float64)::Float64
    i = searchsortedfirst(tsr_table, lambda)
    i > length(tsr_table) && return coeff_table[end]
    i == 1 && return coeff_table[1]
    t = (lambda - tsr_table[i-1]) / (tsr_table[i] - tsr_table[i-1])
    return coeff_table[i-1] + t * (coeff_table[i] - coeff_table[i-1])
end
