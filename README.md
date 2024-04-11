esgtoolkit | <a class="github-button" href="https://github.com/Techtonique/esgtoolkit/stargazers" data-color-scheme="no-preference: light; light: light; dark: dark;" data-size="large" aria-label="Star esgtoolkit/esgtoolkit on GitHub">Star</a>
==========

[![Last Commit](https://img.shields.io/github/last-commit/Techtonique/esgtoolkit)](https://github.com/Techtonique/esgtoolkit)
[![HitCount](https://hits.dwyl.com/Techtonique/esgtoolkit.svg?style=flat-square)](http://hits.dwyl.com/Techtonique/esgtoolkit)

**Important:** Starting with 1.0.0, the package is renamed as lower case 'esgtoolkit' (to finally remove all my active packages from CRAN)

A toolkit for Monte Carlo Simulations in Finance, Economics, Insurance, Physics. Multiple simulation models can be created by combining building blocks provided in the package. 

For __more details__, you can read the package  [vignette on 
ResearchGate](https://www.researchgate.net/publication/338549100_esgtoolkit_a_tool_for_stochastic_simulation_v020). Functions' documentation can be found in section 'Reference' of the [website](https://techtonique.github.io/esgtoolkit/). 

# Table of Contents

- [Installation](#Installation)
- [Quickstart](#Quickstart)
- [Contributing](#Contributing)
- [License](#License)


# Installation

- From Github: 

```r
library(devtools)
devtools::install_github("Techtonique/esgtoolkit")
```

- From R universe: 

```r
# Enable universe(s) by techtonique
options(repos = c(
  techtonique = 'https://techtonique.r-universe.dev',
  CRAN = 'https://cloud.r-project.org'))

# Install some packages
install.packages('esgtoolkit')
```

# Quickstart

In addition to the example below, you can read:
- [this blog](https://thierrymoudiki.wordpress.com/)'s archives 
- the functions' examples in section 'Reference' on [the website](https://techtonique.github.io/esgtoolkit/)
contain code + a lot graphs


```r
library(esgtoolkit)

# Geometric Brownian Motion (https://en.wikipedia.org/wiki/Geometric_Brownian_motion)

eps0 <- simshocks(n = 100, horizon = 5, frequency = "quart")
sim.GBM <- simdiff(n = 100, horizon = 5, frequency = "quart",   
               model = "GBM", 
               x0 = 100, theta1 = 0.03, theta2 = 0.1, 
               eps = eps0)
esgplotbands(sim.GBM, xlab = "time", ylab = "values", main = "with esgplotbands")                
matplot(as.vector(time(sim.GBM)), sim.GBM, type = 'l', main = "with matplot")

```


# Contributing

Your contributions are welcome, and valuable. Please, make sure to __read__ the [Code of Conduct](CONTRIBUTING.md) first.

If you're not comfortable with Git/Version Control yet, please use [this form](https://forms.gle/oqvuDU4JQnnmgevx6).

# License

[BSD 3-Clause Clear](https://techtonique.github.io/esgtoolkit/LICENSE-text.html) © Thierry Moudiki, 2014. 


<script async defer src="https://buttons.github.io/buttons.js"></script>
