
<!-- README.md is generated from README.Rmd. Please edit that file -->

# DiscreteGapStatistic

<!-- badges: start -->
<!-- badges: end -->

The goal of DiscreteGapStatistic is to estimate number of clusters from
(multiple choice) categorical response format data given any clustering
algorithm extending the well-known gap statistic using a discrete
distance based approach.

## Installation

You can install the development version of DiscreteGapStatistic from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ecortesgomez/DiscreteGapStatistic")
```

## Example: Math Anxiety Data

Basic usage of `DiscreteGapStatistic` is shown using `likert`’s Math
Anxiety data.

``` r
library(DiscreteGapStatistic)
#> Warning: replacing previous import 'dplyr::recode' by 'likert::recode' when
#> loading 'DiscreteGapStatistic'
```

``` r
library('dplyr')
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
```

First upload libraries and lightly modify the data. Questions are
formatted and the categories are shortened.

``` r
library(likert)
data(mass)

massQ <- colnames(mass)
massQsh <- c('Gender',
             'Q1: Math Interesting', 'Q2: Uptight Math Test', 'Q3: Use Math In Future',
             'Q4: Mind Goes Blank', 'Q5: Math Relates to Life', 'Q6: Ability Solve Math Probls',
             'Q7: Sinking Feeling', 'Q8: Math Challenging', 'Q9: Math Makes Me Nervous',
             'Q10: Take More Math Classes', 'Q11: Math Makes Me Uneasy','Q12: Math Favorite Subj.',
             'Q13: Enjoy Learning Math', 'Q14: Math Makes Feel Confused')

massData <- mass
colnames(massData) <- massQsh
Cats <- setNames(object = c('SD', 'D', 'N', 'A', 'SA'),
                 nm = c('Strongly Disagree', 'Disagree', 
                        'Neutral', 
                        'Agree', 'Strongly Agree') )
massData <- data.frame(Gender = massData$Gender,
                       apply(massData[, -1], 2,
                             function(x) Cats[as.character(x)] %>%
                                factor(levels = Cats)),
                       check.names=FALSE)

rownames(massData) <- c(paste0('F', 1:2), paste0('M', 1),
                        paste0('F', 3:11), paste0('M', 2),
                        paste0('F', 12), paste0('M', 3),
                        paste0('F', 13), paste0('M', 4:5),
                        paste0('F', 14), paste0('M', 6))
```

First we visualize the data using a heatmap similar to the one produced
by `likert::likert.heatmap.plot`.

``` r
likert.heat.plot2(massData[, -1],
                  allLevels = Cats,
                  text.size = 1.5)+
   labs(title = 'Math Anxiety Data Likert Heatmap Summary')+
   theme(axis.text = element_text(size = 5), 
         title = element_text(size = 5))
```

<img src="man/figures/README-LikertHeat-1.png" width="600" style="display: block; margin: auto;" />
