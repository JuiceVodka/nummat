module dn1

export RedkaMatrika, sor

import Base: getindex, setindex!, size, *, -

"""
    RedkaMatrika(V, I)

podatkovni tip za redke matrike. Polja sta matriki velikosti nxm z neničelnimi elementi, kjer je m<=n ter matrika ki jo opisuje velikosti nxn. 
V mattriki V so zapisane vrednosti neničelnih elementov v posameznih vrsticah. V matriki I so zapisani stolcpi istoležnih vvrednosti v posameznih vrsticah
"""
mutable struct RedkaMatrika
    V
    I
end

"""
    n, n = size(R)

funkcija vrne dimenzije redke matrike R
"""
function size(R::RedkaMatrika)
    return (size(R.V)[1], size(R.V)[1])
end

"""
    negR = -R

funkcija implementira množenje redke matrike z -1
"""
function -(R::RedkaMatrika)
    R.V = -R.V
    return R
end


"""
    a = R[i, i]

pridobivanje vrednosti iz redke matrike R na določenem mestu
"""
function getindex(R::RedkaMatrika, i::Int, j::Int)
    n, _ = size(R)
    if (i > n || j > n || i < 1 || j < 1)
        throw(BoundsError(R, (i, j)))
    end
    

    if (j in R.I[i, :])
        ix = findfirst(==(j), R.I[i, :])
        return R.V[i, ix]
    end

    return 0
end

"""
    R[i, j] = a

funkcija za nastavljanje vrednosti redke matrike R na določenem mestu
"""
function setindex!(R::RedkaMatrika, v, i::Int, j::Int)
    n, m = size(R.I)
    if (i > n || j > n || i < 1 || j < 1)
        throw(BoundsError(R, (i, j)))
    end

    if R[i, j] != 0
        ix = findfirst(==(j), R.I[i, :])
        R.V[i, ix] = v
    else 

        dummy_index = findfirst(==(nothing), R.I[i, :])
        if dummy_index != nothing
            R.I[i, dummy_index] = j
            R.V[i, dummy_index] = v
        else
            R.I = hcat(R.I, fill(nothing, n))
            R.V = hcat(R.V, fill(nothing, n))

            R.I[i, end] = j
            R.V[i, end] = v
        end
    end
end

"""
    x = R*v

funkcija za množenje redke matrike R z vektorjem v z desne
"""
function *(R::RedkaMatrika, v::Vector)
    n, m = size(R.I)
    result = zeros(eltype(v), n)
    for i in 1:n #rows
        for j in 1:m #entries in columns
            col = R.I[i, j]
            if col == nothing
                continue
            end
            val = R.V[i, j]
            result[i] += val * v[col]
            
        end
    end

    return result
end

"""
    x_new = korak_gs(A, x, b)

funkcija izračuna en korak Gaus-Seidlove metode ter vrne nov približek
"""
function korak_gs(A, x, b)
    x_new = copy(x)
    for i in 1:length(x)
        sum_to = 0
        sum_after = 0
        for j in 1:i-1
            sum_to += A[i, j] * x_new[j]
        end
        for j in i+1:length(x)
            sum_after += A[i, j] * x[j]
        end
        #println(1/A[i, i])
        x_new[i] = (1/A[i, i])*(b[i] - sum_to - sum_after)
    end
    return x_new
end

"""
    x, it = sor(A, b, x0, omega)

funkcija izvede sor iteracijo za reševanje sistema Ax = b, z začetnim približkom x0.
 Omega predstavlja parameter, kako močno upoštevamo korak_gs v vsaki iteraciji
"""
function sor(A, b, x0, omega, tol=1e-10)
    x = x0
    it = 0
    #println(A*x - b)
    #println(max(abs(maximum(A*x - b)), abs(minimum(A*x - b))))
    while max(abs(maximum(A*x - b)), abs(minimum(A*x - b))) > tol #infinity norm is max element
        x_new = korak_gs(A, x, b)
        x = omega * x_new + (1-omega) * x
        it += 1
        #println(it, " ", x)
    end

    return x, it
end


end # module dn1
