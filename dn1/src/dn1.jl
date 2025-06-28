module dn1

export RedkaMatrika, sor

import Base: getindex, setindex!, size, *, -


mutable struct RedkaMatrika
    V
    I
end


function size(R::RedkaMatrika)
    return (size(R.V)[1], size(R.V)[1])
end

function -(R::RedkaMatrika)
    R.V = -R.V
    return R
end

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
