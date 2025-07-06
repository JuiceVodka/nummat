using dn2

f(x) = sin(x) / x

interval = Interval(0, 5)
integral = Integral(f, interval)

atol = 1e-10
max_iter = 100
n = 2
prev_rez = integriraj(integral, Gaus_legendre(2))

#' natančnost na 10 mest dosežemo, ko je razlika med dvema zaporednima približkoma manjša kot 1e-10
for i in 0:max_iter
    global n
    global prev_rez
    n *= 2
    rez = integriraj(integral, Gaus_legendre(n))
    if abs(prev_rez - rez) < atol
        println(rez)
        println("St korakov: ", n)
        break
    end
    println(n, " ", rez)
    println(abs(prev_rez - rez))
    prev_rez = rez
end

#TODO tests, report, comments