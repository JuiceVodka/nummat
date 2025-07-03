using dn3

# Run it
mu = 4.0
dt = 1e-5
t_max = 100.0
u = [1., 0.]

T = solve_vanderpol_period(mu, dt, t_max, u)
println("Estimated limit cycle period: ", T)