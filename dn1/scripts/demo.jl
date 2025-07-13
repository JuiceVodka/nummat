using Graphs
using dn1
using Plots

#' ## 18.1.1 SOR iteracija za razpršene matrike

#' Igor Nikolaj Sok, nummat 2025

#' V sklopu te naloge smo morali implementirati podatkovni tip RazpršenaMatrika(V, I), ki hrani razpršeno matriko velikosti nxn s pomočjo matrik V in I velikosti nxm, kjer je m<=n in velja A[i, I[i, j]] = V[i, j]
#' Podatkovni tip je uporaben, ko se srečujemo z matrikami, ki imajo veliko ničelnh elementov. S tem da jih predstavimo na ta način prihranimo veliko prostora v primerjavi s tradicionalno predstavitvijo. V najslabšem možnem primeru (m==n) pa predstavitev z tipom RedkaMatrika zasede dvakrat toliko prostora kot tradicionalna predstavitev. 

#' Drugi del naloge je bila implementacija metode SOR. SOR je iterativna metoda za reševanje linearnega sistema, ki temelji na Gauss-Seidlovi metodi.
#' SOR deluje tako, da v vsakem koraku začetni približek posodobi tako, da sešteje $\omega$ * novi približek ter $(1- \omega)$ * stari približek.
#' Novi približek se po Gaus-Seidlovi metodi izračuna s pomočjo naslednje enačbe:
#' $x_i^{(k)} = \frac{1}{a_{ii}} (b_i - \sum_{j<i} a_{ij} x_j^{(k)} - \sum_{j>i} a_{ij} x_j^{(k-1)})$
#' Metoda SOR se ustavi, ko je $||Ax^{(k)} - b||_{\inf} < \delta$

#' S pomočjo metode SOR ter podatkovnim tipom RedkaMatrika je bilo potrebno rešiti problem vložitev grafa v ravnino ali prostor s fizikalno metodo (opisano v poglavju 6 knjige o numerični matematiki v jeziku julija)
#' Metodo sem preizkusil na grafu krožne lestve. Tri vozlišča sem na začetku fiksiral v ogljišča trikotnika. Rezultat vložitve s SOR metodo je predstavljen na naslednjem grafu.

#helper functions for "vložitev grafa v ravnino ali prostor s fizikalno metodo", adapted from the book
"""
    G = krožna_lestev(n)

Ustvari graf krožna lestev z `2n` točkami.
"""
function krožna_lestev(n)
    G = SimpleGraph(2 * n)
    # prvi cikel
    for i = 1:n-1
        add_edge!(G, i, i + 1)
    end
    add_edge!(G, 1, n)
    # drugi cikel
    for i = n+1:2n-1
        add_edge!(G, i, i + 1)
    end
    add_edge!(G, n + 1, 2n)
    # povezave med obema cikloma
    for i = 1:n
        add_edge!(G, i, i + n)
    end
    return G
end


"""
    A = matrika(G::AbstractGraph, sprem)

Poišči matriko sistema linearnih enačb za vložitev grafa `G` s fizikalno metodo.
Argument `sprem` je vektor vozlišč grafa, ki nimajo določenih koordinat.
Indeksi v matriki `A` ustrezajo vozliščem v istem vrstnem redu,
kot nastopajo v argumentu `sprem`.
"""
function matrika(G::AbstractGraph, sprem)
    # preslikava med vozlišči in indeksi v matriki
    v_to_i = Dict([sprem[i] => i for i in eachindex(sprem)])
    m = length(sprem)
    #A = spzeros(m, m)
    I = Matrix{Union{Int, Nothing}}(undef, m, 1)
    fill!(I, nothing)

    V = Matrix{Union{Float64, Nothing}}(undef, m, 1)
    fill!(V, nothing)
    A = RedkaMatrika(V, I)
    for i = 1:m
        vertex = sprem[i]
        sosedi = neighbors(G, vertex)
        for vertex2 in sosedi
            if haskey(v_to_i, vertex2)
                j = v_to_i[vertex2]
                A[i, j] = 1
            end
        end
        A[i, i] = -length(sosedi)
    end
    return A
end


"""
    b = desne_strani(G::AbstractGraph, sprem, koordinate)

Poišči desne strani sistema linearnih enačb za eno koordinato vložitve grafa `G`
s fizikalno metodo. Argument `sprem` je vektor vozlišč grafa, ki nimajo
določenih koordinat. Argument `koordinate` vsebuje eno koordinato za vsa
vozlišča grafa. Metoda uporabi le koordinato vozlišč, ki so pritrjena.
Indeksi v vektorju `b` ustrezajo vozliščem v istem vrstnem redu,
kot nastopajo v argumentu `sprem`.
"""
function desne_strani(G::AbstractGraph, sprem, koordinate)
    set = Set(sprem)
    m = length(sprem)
    b = zeros(m)
    for i = 1:m
        v = sprem[i]
        for v2 in neighbors(G, v)
            if !(v2 in set) # dodamo le točke, ki so fiksirane
                b[i] -= koordinate[v2]
            end
        end
    end
    return b
end

"""
    vloži!(G::AbstractGraph, fix, točke)

Poišči vložitev grafa `G` v prostor s fizikalno metodo. Argument `fix` vsebuje
vektor vozlišč grafa, ki imajo določene koordinate. Argument `točke` je
začetna vložitev grafa. Koordinate vozlišč, ki niso pritrjena, bodo nadomeščene
z novimi koordinatami.
Metoda ne vrne ničesar, ampak zapiše izračunane koordinate v matriko `točke`.
"""
function vloži!(G::AbstractGraph, fix, točke, omega)
    sprem = setdiff(vertices(G), fix)
    dim, _ = size(točke)
    A = matrika(G, sprem)
    x0 = ones(size(A)[1])
    it = 0
    for k = 1:dim
        b = desne_strani(G, sprem, točke[k, :])
        x = sor(A, b, x0, omega)#x = cg(-A, -b) # matrika A je negativno definitna, mejbi probi z -
        #println(x)
        #println(b)
        #println(length(sprem))
        #println(length(x))
        točke[k, sprem] = x[1]
        x0 = x[1]
        it += x[2]
    end
    return it
end


G = krožna_lestev(10)
tocke = zeros(2, nv(G))
fix = [1, 4, 7]
tocke[:, 1] = [0.0, 0.0]
tocke[:, 4] = [1.0, 0.0]
tocke[:, 7] = [0.5, 1.0]

bestit = 10000
bestomega = 0
its = []
for omega in 0.1:0.1:2
    println(omega)
    global bestit
    global bestomega
    tocke = zeros(2, nv(G))
    fix = [1, 4, 7]
    tocke[:, 1] = [0.0, 0.0]
    tocke[:, 4] = [1.0, 0.0]
    tocke[:, 7] = [0.5, 1.0]
    it = vloži!(G, fix, tocke, omega)
    if it < bestit
        bestit = it
        bestomega = omega
    end
    push!(its, it)
end

it = vloži!(G, fix, tocke, 0.5)
println(it)
println(bestit)
println(bestomega)
#println(tocke)


"""
    narisi_vlozitev(G::AbstractGraph, točke)

funkcija nariše pozicije vozlišč v grafu G z podanimi koordinatami v matriki točke, kjer vsak stolpec predstavlja eno torčko, vrstici pa x in y koordinato
"""
function narisi_vlozitev(G::AbstractGraph, točke)
    x = točke[1, :]
    y = točke[2, :]

    # Draw edges
    for e in edges(G)
        v1 = src(e)
        v2 = dst(e)
        plot!([x[v1], x[v2]], [y[v1], y[v2]], c=:gray, lw=1)
    end

    # Draw vertices
    scatter!(x, y, c=:red, ms=5, label=false)

    # Display the plot
    plot!(aspect_ratio=:equal, legend=false)
end

plot()
narisi_vlozitev(G, tocke)     
#' 

#' Za iskanje optimalnega $\omega$ parametra metode SOR lahko pogledamo, koliko iteracij je potrebno za izračun optimalne vložitve na 10 decimalnih mest natančno. Najmanj iteracij je potrebno pri $\omega = 1.4$ in sicer 98 iteracij
