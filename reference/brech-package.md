# brech: Stochastic Modelling of Infectious Disease Outbreaks

The rapid construction and prototyping of stochastic infectious disease
models in small (~100 individuals) to medium (~1000 individuals) size
settings. Includes tools for defining model parameters, rates, and
initial conditions along with model fitting and scenario projections.
This ease of construction does come at the cost of computational
efficiency and so would not be recommended for larger settings or with
more complicated dynamics with multiple rates. Stochastic simulation is
performed using the adaptive tau leaping algorithm Cao (2006)
[doi:10.1063/1.2745299](https://doi.org/10.1063/1.2745299) and model
fitting is implemented using approximate Bayesian computation Kypraios
(2017)
[doi:10.1016/j.mbs.2016.07.001](https://doi.org/10.1016/j.mbs.2016.07.001)
.

## See also

Useful links:

- <https://bccdc-phsa.github.io/brech/>

- <https://github.com/BCCDC-PHSA/brech>

- Report bugs at <https://github.com/BCCDC-PHSA/brech/issues>

## Author

**Maintainer**: Mike Irvine <mike.irvine@bccdc.ca>
([ORCID](https://orcid.org/0000-0003-4785-8998)) \[copyright holder\]
