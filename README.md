# nummat
Repo for the nummat course at FRI

## 18.1.1 SOR iteracija za razpršene matrike

Igor Nikolaj Sok

V sklopu te naloge je bilo potrebo implementirati podatkovni tip RedkaMatrika, ki dobro deluje za matrike, katerih večin aelementov je nič (bolj prostorsko učinkovite) ter metodo SOR za reševanje linearnih sistemov.
Z implementiranim tipom in metodo je bilo nato potrebno rešiti problem fizikalne vložitve grafov v prostor ali ravnino. Za na problem je bilo potrebno tudi poiskati optimalni $\omega$ parameter metode SOR

Da poženemo program, ki izračuna fizikalno vložitev grafa krožne lestve v ravnino ter izračun optimalnega $\omega$ parametra metode SOR za ta problem v julia REPL poženemo "include("dn1/scripts/dn1.jl")".
Da poženemo teste v pkg načinu v julia REPL poženemo "test dn1".
Za ustvarjanje poročila poženemo "include("dn1/doc/porocilo.jl")"
