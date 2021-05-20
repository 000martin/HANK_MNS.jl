# Replication of McKay, Nakamura, and Steinsson (2016)

This is a replication of ["The Power of Forward Guidance"](https://www.aeaweb.org/articles?id=10.1257/aer.20150063) by Alisdair McKay, Emi Nakamura, and Jón Steinsson (henceforth MNS).

Authors: @mhaense1, @000martin

The project was part of our evaluation for the course [Numerical Methods](https://floswald.github.io/NumericalMethods/) at SciencesPo Paris in Spring 2021. 

The original paper studies the Forward Guidance, an unconventionary monetary tool, which, in standard 
representative agent New Keynesian models, has strong effects on current economic outcomes.

However, as MNS (and we using their model) find, Forward Guidance has considerably less power in models with
heterogenous agents facing unisurable idiosyncratic income risk. 

Conceptually, MNS compute perfect foresight transition paths of their simple Heterogenous Agents New Keynesian (HANK) model following a one-time announcement of an interest change far in the future (20 quarters in the paper).

Our Julia code is able to replicate their main exercises and produce equivalents of Figures 3,4,5 and 6 as well as Table 2 in their paper.
Unfortunalety, it does not (yet) extend to their Zero Lower Bound Analysis featuring time-varying subjective discount factors (Section II.C in the paper).

## How does it work? 

In case you want to simply reproduce these main results, you use our convenience functions, which allow you to obtain the results as simple as follows:

```julia

using HANK_MNS

#Displays our replications of figures 3 and 4
get_figures_3_4()

#Saves a table containing our replication of table 2 
tb2 = get_table_2()

#Displays our replications of figures 5 and 6
#This may take quite a while, since ~20 transition paths must be computed
get_figures_5_6()

```

These convenience functions can also take inputs, allowing to analyze different time horizons or rate changes. See the documentation below for details.

In case one wants to produce a transition path with a custom calibration using the incomplete contracts model, one would proceed as follows:

```julia
using HANK_MNS

#Example: choose a higher aggregate asset level and higher income risk
p = set_parameters(B = 7.8, &sigma; = 0.03^0.5)

#Find the \beta consistent with the target interest rate and asset level in steady state
#and obtain a corresponding steady state structure SS.
#Depending on your calibration, you might have to choose a different betaRange,
#e.g [0.97,0.995] instead of [0.95,0.99]. See function documentation for details. 
b, y, SS = get_steady_state(p,0.9,[0.95,0.99])

#set the correct \beta value
p.&beta; = b

#obtain the transition path for the announcemnt of a 22 quarter ahead interest rate reduction of 60 basis points.
trp = get_transition_full(;TR = 22, RChange = -0.006)

#trp is a structure containing the transition paths for, amongst others, aggregate output, inflation and the wage rate.
#Make a table, plot it,...
```



```@autodocs
Modules = [HANK_MNS]
```


end
