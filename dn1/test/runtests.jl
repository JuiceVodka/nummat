using Test
using dn1


I = [nothing 4;
    1 3;
    nothing nothing;
    nothing 2]
V = [nothing 12;
    -1 5;
    nothing nothing;
    nothing 3]


@testset "Podatkovni tip RedkaMatrika je definiran" begin
    @test RedkaMatrika(V, I).V == V
    @test RedkaMatrika(V, I).I == I
end

@testset "funkcija size je definirana za redko matriko" begin
    @test size(RedkaMatrika(V, I)) == (4, 4)
end


@testset "funkcija getindex je definirana za Redko matriko in vrne prave vrednosti" begin
    R = RedkaMatrika(V, I)
    @test R[1, 4] == 12
    @test R[3, 3] == 0
    @test R[4, 2] == 3

    @test_throws BoundsError R[0, 3] == 5
    @test_throws BoundsError R[5, 3] == 0
    @test_throws BoundsError R[3, 0] == 5
    @test_throws BoundsError R[3, 5] == 0
end


@testset "Funkcija setindex! je definirana in deluje za RedkoMatriko" begin
    R = RedkaMatrika(V, I)
    R[1, 1] = 30
    @test R[1, 1] == 30 #spremenimo vrednost ki je prej nothing

    R[2, 1] = -15 #spremenimo vrednost ki je prej neka druga vrednost
    @test R[2, 1] == -15

    R[2, 4] = 6 #spremenimo vrednost tako da moramo podaljšati vse vrstice
    @test R[2, 4] == 6
    @test R[2, 3] == 5 #nismo pokvarili prejšnje vrednosti
    print(R.V)

    @test_throws BoundsError R[0, 3] == 5
    @test_throws BoundsError R[5, 3] == 0
    @test_throws BoundsError R[3, 0] == 5
    @test_throws BoundsError R[3, 5] == 0
end


I = [nothing 4;
    1 3;
    nothing nothing;
    nothing 2]
V = [nothing 12;
    -1 5;
    nothing nothing;
    nothing 3]

@testset "Množenje z leve z vektorjem je implementirano in pravilno deluje" begin
    R = RedkaMatrika(V, I)
    x = [1, 1, 1, 1]
    #rezultat : [12, 4, 0, 3]

    @test R*[0, 0, 0, 0] == [0, 0, 0, 0]
    @test R*x == [12, 4, 0, 3]

    x = [1, 2, 3, 4]
    @test R*x == [48, 14, 0, 6]

end

#for SOR to converge matrix has to be diagonaly dominant and SPD
I = [1 4 nothing;
    1 2 3;
    3 nothing nothing;
    2 4 nothing]
V = [15. 12. nothing;
    -1. 10. 5.;
    10. nothing nothing;
    3. 10. nothing]

@testset "Funkcija sor je implementirana in pravilno deluje" begin
    R = RedkaMatrika(V, I)
    println(R * [1, 1, 1, 1])


    b = [27., 14., 10., 13.]
    x0 = [0.1, 1., 0.1, 1.]
    println(R[1, 1], R[2, 2], R[3, 3], R[4, 4])

    x, it = sor(R, b, x0, 0.5)

    @test isapprox(x, [1, 1, 1, 1])

end



#TODO: comments for all the functions, using SOR to solve the problem of graph insertion, report
#za graph insertion porabi  kodo na straneh 76, 77