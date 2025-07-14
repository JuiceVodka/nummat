using dn2
#' # 18.2.13 Gauss-Legendrove kvadrature
#' Igor Nikolaj Sok, nummat 2025

#' ## Opis problema
#' Potrebno je bilo izpeljati Gaus-Legendrovo integracijsko pravilo na dveh točkah
#' $\int_0^h f(x)\,dx = A f(x_1) + B f(x_2) + R$

#' Za pravilo na dveh točkah potrebujemo eksaktnost za polinom do stopnje 3

#' Kot lokacije vozlišč lahko vzamemo ničle 2. Legendrovega polinoma $P(x) = 0.5 ( 3x^2 - 1)$
#' Ničli sta torej $x1 = \frac{-1}{\sqrt{3}}$ ter $x2 = \frac{1}{\sqrt{3}}$

#' Uteži integriranja dobimo iz eksaktnosti integrala za konstante (na intervalu [-1, 1])
#' $\int_{-1}^1 1dx = 2$
#' $A f(x_1) + B f(x_2) = A + B = 2$
#' Zaradi simetričnosti; $A = B = 1$

#' Pravilo na [-1, 1] je torej enako $\int_{-1}^1 f(x)dx = f(\frac{-1}{\sqrt(3)}) + f(\frac{1}{\sqrt(3)})$

#' Napaka Gaus-Legendrovega pravila na n točkah je:
#' $Rf = \frac{(b-a)^{2n+1} (n!)}{(2n+1)((2n)!)^3}f^{(2n)}(\xi)$

#' V našem primeru je napaka torej enaka $\frac{f^{(4)}(\xi)}{135}$
#' Ko transformiramo na poljuben interval z n podintervali je napaka potem enaka $\frac{n^5 f^{(4)}(\xi)}{4320}$

#' Za transformacijo pravila na interval [0, h] zamenjamo spremenljivko xx v intervalu z $\frac{h}{2}(t+1)$, kjer je t element intervala [-1, 1]
#' dx je torej enako $\frac{h}{2}dt$

#' Integral torej postane: $\int_0^h f(x)dx = \frac{h}{2}(f(\frac{h}{2}(\frac{-1}{\sqrt(3)} + 1)) + f(\frac{h}{2}(\frac{1}{\sqrt(3)} + 1)))$

#' Če ga hočemo iz intervala [0, h] prestaviti na interval [x_i, x_{i+1}] preprosto v funkciji prištejemo $x_i$. Velja da $x_{i+1} - x_i = h$
#' Integral torej postane: 
#' $\int_{x_i}^{x_{i+1}} f(x)dx = \frac{h}{2}(f(x_i + \frac{h}{2}(\frac{-1}{\sqrt(3)} + 1)) + f(x_i + \frac{h}{2}(\frac{1}{\sqrt(3)} + 1)))$

f(x) = sin(x) / x

interval = Interval(0, 5)
integral = Integral(f, interval)

atol = 1e-10
max_iter = 100
n = 2
prev_rez = integriraj(integral, Gaus_legendre(2))

#' V primeru računanja vrednosti integrala $\int_0^5 \frac{sin(x)}{x}dx$ natančnost na 10 mest dosežemo, ko je razlika med dvema zaporednima približkoma pri različni delitni na podintervale manjša kot 1e-10
#' V vsakem koraku število podintervalov podvojimo ter tako dobimo oceno števila potrebnih korakov
#' Ocenjeno število korakov za izračun integrala na 10 mest natančno je 256
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
