using dn3
using Test

@testset "Funkcija lc je implementirana in pravilno sestavi sistem de za podano začetno stanje" begin
    u = [1, 0]
    f = lc(u)
    @test f[1] == 0
    @test f[2] == -1

    u = [1, 1]
    f = lc(u)
    @test f[1] == 1
    @test f[2] == -1

    u = [2, 1]
    f = lc(u)
    @test f[1] == 1
    @test f[2] == -14
end

@testset "Funkcija rk4_step izvede en korak rk4 metode" begin
    fun = lc
    dt = 0.2
    u = [1, 1]
    rez = rk4_step(u, dt, fun)
    @test isapprox(rez, [1.16960830293, 0.6596780627])

    dt = 2
    fun = lc
    u = [1, 1]
    rez = rk4_step(u, dt, fun)
    @test isapprox(rez, [1/3, -1])

    dt = 1
    fun = lc
    u = [1, 1]
    rez = rk4_step(u, dt, fun)
    @test isapprox(rez, [4/3, 7/6])
end


@testset "funkcija solve_vanderpol_period je implementirana in vrne periodo podane funkcije lc" begin
    dt = 0.1
    t_max = 100.
    u = [1, 1]
    period = solve_vanderpol_period(dt, t_max, u)

    @test isapprox(period, 10.2)

    dt = 1e-5
    t_max = 100.
    u = [5, 6]
    period = solve_vanderpol_period(dt, t_max, u)
    aproximation = (3- 2* log1p(1))*4 + (2*pi)/(4^(1/3))
    @test isapprox(period, 10.20352)

    dt = 1e-5
    t_max = 100.
    u = [4, 3]
    period = solve_vanderpol_period(dt, t_max, u)
    aproximation = (3- 2* log1p(1))*4 + (2*pi)/(4^(1/3))
    @test isapprox(period, 10.20352)
end