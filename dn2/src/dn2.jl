module dn2

export Interval, Integral, Gaus_legendre, integriraj


"""
    Interval(min, max)

podatkovni tip za interval preko katerega se bo integriralo, podan z spodnjo ter zgornjo mejo
"""
struct Interval
    min
    max
end

"""
    Integral(fun, int::Interval)

podatkovni tip za integral z podano funkcijo, ki jo integriramo ter intervalom na katerem jo integriramo
"""
struct Integral
    fun
    int::Interval
end

"""
    Gaus_legendre(n)

podatkovni tip za gaus-legendre integracijsko pravilo z parametrom n kot številom podintervalov
"""
struct Gaus_legendre
    n #st podintervalov
end

"""
    I = integriraj(integral::Integral, metoda::Gaus_legendre)

funkcija, ki integrira funkcijo integral.fun na intervalu integral.int z gaus_legendrovim integracijskim pravilom z
metoda.n podintervalov
"""
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
