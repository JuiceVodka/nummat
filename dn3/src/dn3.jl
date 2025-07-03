module dn3

export rk4_step, solve_vanderpol_period

function lc(u, mu)
    y1 = u[1]  # x
    y2 = u[2]  # x'
    return [y2,
            mu * (1 - y1^2) * y2 - y1]
end


"""
u = rk4_step(u, dt, mu, lc)

A function that performs a step of Rk4 method with the given initial state u, step size dt, parameter mu and function fun

"""
function rk4_step(u, dt, mu, fun)
    k1 = fun(u, mu)
    k2 = fun(u .+ 0.5dt .* k1, mu)
    k3 = fun(u .+ 0.5dt .* k2, mu)
    k4 = fun(u .+ dt .* k3, mu)

    return u .+ dt * (k1 .+ 2k2 .+ 2k3 .+ k4) ./ 6
end

"""
    perioda = solve_vanderpol_period(mu, dt, t_max, u)

Function that finds the period of the Van der Pol oscilator with the given mu parameter time step dt, max time t_max and
 initial state u; u = [x, x']
"""
function solve_vanderpol_period(mu::Float64, dt::Float64, t_max::Float64, u, atol=1e-10)
    steps = Int(round(t_max / dt))
    t = 0.0
    prev_y1 = u[1]

    zero_crossings = Float64[]
    prev_period = -1

    for _ in 1:steps
        t += dt
        u = rk4_step(u, dt, mu, lc)
        y1 = u[1]
        y2 = u[2]
        #println(u)
        if prev_y1 < 0 && y1 >= 0 && y2 > 0
            push!(zero_crossings, t)
            if length(zero_crossings) == 2
                if prev_period < 0
                    #println("SAVING FIRST PERIOD")
                    prev_period = zero_crossings[2] - zero_crossings[1]
                    zero_crossings = Float64[]
                    push!(zero_crossings, t)
                elseif abs(prev_period - (zero_crossings[2] - zero_crossings[1])) < atol
                    #println(prev_period, " ", (zero_crossings[2] - zero_crossings[1]))
                    return zero_crossings[2] - zero_crossings[1]
                else
                    #println("NOT YET SAME")
                    prev_period = zero_crossings[2] - zero_crossings[1]
                    zero_crossings = Float64[]
                    push!(zero_crossings, t)
                end
            end
        end

        prev_y1 = y1
    end

    throw("Failed to detect zero crossings to the specified accuracy")
end


end # module dn3


#TODO tests, comments, report