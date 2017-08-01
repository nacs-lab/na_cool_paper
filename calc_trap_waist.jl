#


k_B = 1.38064852e-23
I0 = k_B * 3e-3
m_Na = 23e-3 / 6.02e23
ω2 = 420e3 * 2π
ω3 = 580e3 * 2π
w2 = sqrt(4 * I0 / m_Na / ω2^2)
w3 = sqrt(4 * I0 / m_Na / ω3^2)
@show w2 * 1e6
@show w3 * 1e6
