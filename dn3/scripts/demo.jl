using dn3
using Plots
#' # 18.3.3 Perioda limitnega cikla
#' Igor Nikolaj Sok, nummat 2025

#' ## Opis problema
#' V kontekstu naloge je bilo potrebno poiskati periodo limitnega cikla za nelinearno diferencialno enačbo drugega reda:
#' $$
#' \ddot{x}(t) - 4(1 - x^2)\dot{x}(t) + x(t) = 0
#' $$

#' Enačba opisuje Van der Polov oscilator s parametrom $\mu = 4$

#' Enačba Van der polovega oscilatorja je sledeča:

#' V kontekstu naloge je bilo potrebno poiskati periodo limitnega cikla za nelinearno diferencialno enačbo drugega reda:
#' $$
#' \frac{d^2x}{dt^2} - \mu (1 - x^2)\frac{dx}{dt} + x = 0
#' $$

#' Limitni cikel predstavlja stabilno periodično rešitev sistema, kamor konvergera trajektorija iz skoraj vseh začetnih pogojev

#' ## Opis Rešitve

#' Enačbo preoblikujemo v sistem prvega reda:

#' $$
#' \begin{cases}
#' \dot{x}_1 = x_2 \\
#' \dot{x}_2 = 4(1 - x_1^2)x_2 - x_1
#' \end{cases}
#' $$

#' Sistem rešimo numerično z uporabo metode Runge-Kutta (RK4) v zaporednih časovnih korakih. Periodo limitnega cikla dobimo tako, da gledamo, kdaj funkcija x(t) prečka ničlo pri čemer je njen odvod pozitiven. Časovna razdalja med dvema zaporednima trenutkoma kjer slednje velja predstavlja periodo limitnega cikla.
#' Ker želimo periodo določiti na 10 decimalnih mest natančno rezultata ne vrnemo že pri prvi dobljeni periodi, temveč nadaljujemo z iskanjem periode. Ko najdemo novo jo primerjamo s prejšnjo. Če se nova periodo od prejšnje razlikuje za manj kot 1e-10 jo lahko vrnemo, saj smo rezultat dobili na 10 decimalnih mest natančno. V nasprotnem primeru shranimo novo periodo ter s postopkom nadaljujemo.


dt = 1e-2
t_max = 100.0
u = [1., 0.]

#' Izračunana perioda limitnega cikla:
T = solve_vanderpol_period(dt, t_max, u)
#'


println("Estimated limit cycle period: ", T)



steps = Int(round(t_max / dt))
traj = zeros(2, steps)
u = [1., 0.]

for i in 1:steps
    global u
    traj[:, i] .= u
    u = rk4_step(u, dt, lc)
end

#' ## Graf odvisnosti x(t) ter x'(t)$

plot(traj[1, :], traj[2, :],
     xlabel="x(t)", ylabel="x'(t)",
     title="Limitni cikel v faznem prostoru",
     legend=false,
     size=(600, 400))
#'


#' ## Graf odvisnosti x(t) ter t  

plot(1:steps, traj[1, :],
     xlabel="t", ylabel="x(t)",
     title="x(t) v odvisnosti od t",
     legend=false,
     size=(600, 400))
#'


#' Na grafih lahko opazimo periodično obnašanje diferencialne enačbe. Na grafu odvisnosti x(t) od t lahko dobro opazimo tudi periodo (razmik med dvema zaporednima oscilacijama)

