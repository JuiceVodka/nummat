module dn2

export Interval, Integral, Gaus_legendre, integriraj


struct Interval
    min
    max
end

struct Integral
    fun
    int::Interval
end

struct Gaus_legendre
    n #st podintervalov
end


function integriraj(integral::Integral, metoda::Gaus_legendre)
    a = integral.int.min
    b = integral.int.max
    h = (b - a)/metoda.n
    #println(a, b, h)

    I = 0

    for i in 0:metoda.n-1
        premik = a + i*h
        x_i = premik + (h/2) * (1- (1/sqrt(3)))
        x_iplusena = premik + (h/2) * (1+ (1/sqrt(3)))
        I += integral.fun(x_i) + integral.fun(x_iplusena)
    end

    I = I * (h/2)
    #println(I)

    return I
end

end # module dn2
