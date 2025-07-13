# nummat

## 18.2.13 Gauss-Legendrove kvadrature

Igor Nikolaj Sok

V kontekstu naloge je bilo potrebno izpeljati Gaus-Legendrovo integracijsko pravilo na dveh točkah. vključno z napako Rf.
Integral smo izpeljani tako, da smo ga najprej izpeljali na intervalu [-1, 1] ter ga nato prevedli na integral [0, h] ter kasneje še na interval [$x_i$, $x_{i+1}$]. 
Nato je bilo potrebno še aproksimirati število korakov za izračun integrala $\int_0^5 \frac{sin(x)}{x}dx$ na 10 mest natančno

Da poženemo program, ki število korakov za izračun integrala na 10 mest natančno ter njegove vrednosti v julia REPL poženemo "include("dn3/scripts/dn2.jl")"
Da poženemo teste v pkg načinu v julia REPL poženemo "test dn2"
Za ustvarjanje poročila poženemo "include("dn2/doc/porocilo.jl")"
