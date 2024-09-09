
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

# Example: Math Anxiety Data

## Summarization and Data Exploration

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

Libraries are uploaded. Questions are lightly formatted and the
categories are shortened.

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

## Likert Heatmap

The dataset is visualized using a heatmap similar to the one produced by
`likert::likert.heatmap.plot`.

``` r
likert.heat.plot2(massData[, -1],
                  allLevels = Cats,
                  text.size = 1.5)+
   labs(title = 'Math Anxiety Data Likert Heatmap Summary')+
   theme(axis.text = element_text(size = 6), 
         title = element_text(size = 6))
```

<img src="man/figures/README-LikertHeat-1.png" width="800" style="display: block; margin: auto;" />

## Distance Matrices

Five categorical distance functions are introduced to quantify
dissimilarities/discrepancies between two categorical vectors (see
function `distancematrix`). The following table describes the available
distances and the names used within the package.

    #> 
    #> Attaching package: 'kableExtra'
    #> The following object is masked from 'package:dplyr':
    #> 
    #>     group_rows

<table class=" lightable-paper" style="font-family: &quot;Arial Narrow&quot;, arial, helvetica, sans-serif; width: auto !important; margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
Distance
</th>
<th style="text-align:left;">
Name
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;border-left:1px solid;">
Hamming
</td>
<td style="text-align:left;font-family: monospace;border-left:1px solid;border-right:1px solid;">
hamming
</td>
</tr>
<tr>
<td style="text-align:left;border-left:1px solid;">
chi-square
</td>
<td style="text-align:left;font-family: monospace;border-left:1px solid;border-right:1px solid;">
chisquare
</td>
</tr>
<tr>
<td style="text-align:left;border-left:1px solid;">
Cramer’s V
</td>
<td style="text-align:left;font-family: monospace;border-left:1px solid;border-right:1px solid;">
cramerV
</td>
</tr>
<tr>
<td style="text-align:left;border-left:1px solid;">
Hellinger
</td>
<td style="text-align:left;font-family: monospace;border-left:1px solid;border-right:1px solid;">
hellinger
</td>
</tr>
<tr>
<td style="text-align:left;border-left:1px solid;">
Bhattacharyya
</td>
<td style="text-align:left;font-family: monospace;border-left:1px solid;border-right:1px solid;">
bhattacharyya
</td>
</tr>
</tbody>
</table>

The resulting distance matrix from a rectangular dataset can easily be
displayed and organized using heatmaps. The following plots visualize
the `massData` using the defined distances with the function
`distanceHeat`. This function only requires the dataset object and the
name of the distance. It can also take advantage of the options and
functionalities available in the `pheatmap` function (see the code
below) from the `pheatmap` R package \[@kolde2019\]. For instance, by
default, the columns are clustered using
`clustering_method = 'complete'` but this parameter can be specified
according to `pheatmap`’s options. Different aesthetic options are
exemplified in the plots below using the introduced distances.

``` r
distanceHeat(x = massData[, -1], distName = 'hamming',
             main = 'Hamming Distance\nMath Anxiety Data', fontsize = 6.5)

distanceHeat(x = massData[, -1], distName = 'cramerV',
             main = "Cramer's V Distance\nMath Anxiety Data", fontsize = 6.5, 
             show_rownames = FALSE, cluster_rows = FALSE, cluster_cols=FALSE)
```

<img src="man/figures/README-DistanceMats1-1.png" width="330" /><img src="man/figures/README-DistanceMats1-2.png" width="330" />

``` r
distanceHeat(x = massData[, -1], distName = 'hellinger',
             main = 'Hellinger Distance\nMath Anxiety Data', fontsize = 6.5, 
             show_rownames = TRUE, border_color = 'black', 
             , cluster_rows = FALSE, cluster_cols=FALSE)

distanceHeat(x = massData[, -1], distName = 'bhattacharyya',
             main = 'Bhattacharyya Distance\nMath Anxiety Data', fontsize = 6.5, 
              show_rownames = TRUE, border_color = 'lightgrey', 
             , cluster_rows = FALSE, cluster_cols=FALSE)
```

<img src="man/figures/README-DistanceMats2-1.png" width="330" /><img src="man/figures/README-DistanceMats2-2.png" width="330" />
