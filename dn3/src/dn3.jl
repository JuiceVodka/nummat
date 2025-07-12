module dn3

export rk4_step, solve_vanderpol_period, lc

"""
    [x', x''] = lc([x, x'])

funkcija ki sestavi sistem diferencialnih enačb prvega reda za zahtevano enačbo limitnega cikla ter začetno stanje [x, x']
"""
function lc(u)
    y1 = u[1]  # x
    y2 = u[2]  # x'
    return [y2,
            4 * (1 - y1^2) * y2 - y1]
end


"""
    u = rk4_step(u, dt, lc)

Funkcija ki naredi en korak RK4 metoda z podanim začetni stanjem u, korakom dt ter funkcijo fun
"""
function rk4_step(u, dt, fun)
    k1 = fun(u)
    k2 = fun(u .+ 0.5dt .* k1)
    k3 = fun(u .+ 0.5dt .* k2)
    k4 = fun(u .+ dt .* k3)

    return u .+ dt * (k1 .+ 2k2 .+ 2k3 .+ k4) ./ 6
end

"""
    perioda = solve_vanderpol_period(dt, t_max, u)

Funkcija ki najde periodo limitnega cikla podane enačbe (Van der Pol oscilator, mu = 4 v našem primeru) z podanim časovnim korakom dt, končnim časom t_max
ter začetnim stanjem u; u = [x, x']
"""
function solve_vanderpol_period(dt::Float64, t_max::Float64, u, atol=1e-10)
    steps = Int(round(t_max / dt))
    t = 0.0
    prev_y1 = u[1]

    zero_crossings = Float64[]
    prev_period = -1

    for _ in 1:steps
        t += dt
        u = rk4_step(u, dt, lc)
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
