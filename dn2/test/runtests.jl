using dn2
using Test

@testset "Podatkovni tip Interval je implementiran in deluje" begin
    a = Interval(2, 3)
    @test a.min == 2
    @test a.max == 3
end

@testset "Podatkovni tip Integral je implementiran in deluje" begin
    fun(x) = 2x^2 + 3
    a = Interval(2, 3)
    integral = Integral(fun, a)
    @test integral.int.min == 2
    @test integral.int.max == 3
    @test integral.fun(integral.int.min) == 11
    @test integral.fun(integral.int.max) == 21
end


@testset "Podatkovni tip Gaus_legendre je implementiran in deluje" begin
    g = Gaus_legendre(5)
    @test g.n == 5
end

@testset "funckija integriraj je definirana in deluje za podatkovni tip integral ter Gaus_legendre metodo" begin
    fun(x) = 5
    a = Interval(0, 1)
    metoda = Gaus_legendre(10)
    I = integriraj(Integral(fun, a), metoda)
    @test I == 5

    fun(x) = x
    I = integriraj(Integral(fun, a), metoda)
    @test isapprox(I, 0.5)

    fun(x) = 5
    a = Interval(0, 5)
    metoda = Gaus_legendre(10)
    I = integriraj(Integral(fun, a), metoda)
    @test I == 25

    fun(x) = x
    I = integriraj(Integral(fun, a), metoda)
    @test isapprox(I, 12.5)

    fun(x) = sin(x)
    a = Interval(0, pi)
    metoda = Gaus_legendre(200)
    I = integriraj(Integral(fun, a), metoda)
    @test isapprox(I, 2, atol=1e-10)

    fun(x) = cos(x)
    I = integriraj(Integral(fun, a), metoda)
    @test isapprox(I, 0, atol=1e-10)
end