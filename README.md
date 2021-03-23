# dsmart2

## A faster implementation of the dsmart algorithm

Much of the original code can be found at this bitbucket repository: <https://bitbucket.org/brendo1001/dsmart/src/master/>

This has been my attempt at creating a DSMART algorithm that runs faster than the original implementation. I found that this was necessary for large datasets, as the [`raster`](https://github.com/rspatial/raster) package was proving to be insufficient for processing in my own work. There has been a lot of recent work on the [`terra`](https://github.com/rspatial/terra) package, which is the `raster` package replacement and uses C++ computing to significantly increase processing speed. Since the DSMART package used the `raster` package throughout its scripting, it was necessary to go through and rewrite the steps to produce the final outputs using the `terra` package. What was once a process that took a few days with the `raster` implementation of DSMART has turned into a process that takes a few hours with the `terra` implementation.

In order to use a non-default model, all you need to do is define the model in `method.model`. Model specific arguments are then placed in a named list, `args.model`. For example:

    method.model <- "randomForest"
    args.model <- list(mtry = 3, ntree = 100)

These objects are passed into the `disaggregate` function to create the model with your set parameters. A full set of your model specific parameters can be viewed as follows:

    library(mlr3verse)
    lrn <- lrn("classif.randomForest")
    lrn$param_set

    <ParamSet>
                 id    class lower upper nlevels        default value
     1:       ntree ParamInt     1   Inf     Inf            500      
     2:        mtry ParamInt     1   Inf     Inf <NoDefault[3]>      
     3:     replace ParamLgl    NA    NA       2           TRUE      
     4:     classwt ParamUty    NA    NA     Inf                     
     5:      cutoff ParamUty    NA    NA     Inf <NoDefault[3]>      
     6:      strata ParamUty    NA    NA     Inf <NoDefault[3]>      
     7:    sampsize ParamUty    NA    NA     Inf <NoDefault[3]>      
     8:    nodesize ParamInt     1   Inf     Inf              1      
     9:    maxnodes ParamInt     1   Inf     Inf <NoDefault[3]>      
    10:  importance ParamFct    NA    NA       4          FALSE      
    11:    localImp ParamLgl    NA    NA       2          FALSE      
    12:   proximity ParamLgl    NA    NA       2          FALSE      
    13:    oob.prox ParamLgl    NA    NA       2 <NoDefault[3]>      
    14:  norm.votes ParamLgl    NA    NA       2           TRUE      
    15:    do.trace ParamLgl    NA    NA       2          FALSE      
    16: keep.forest ParamLgl    NA    NA       2           TRUE      
    17:  keep.inbag ParamLgl    NA    NA       2          FALSE      
    18: predict.all ParamLgl    NA    NA       2          FALSE      
    19:       nodes ParamLgl    NA    NA       2          FALSE  

In the above example, anything in the `id` column can be placed in the `args.model` list; for example,

    args.model = list(mtry = 3, ntree = 100, importance = "gini")

See `lrn$param_set$levels` for a list of valid entries that may be used, or read the respective model package help page (i.e.: `?randomForest::randomForest`).

Some important notes on the changes made, starting in the `disaggregate()` coding:

1.  `mlr3` coding is used to generate all model types, including the `C5.0` defaults. Model types in this program are only allowed to be classification or probability based. No model resampling is performed here - this is a more recent update. Previously, the `caret` package was used to generate models that were not default, and you could perform model resampling within that coding if you wished. I chose to eliminate resampling here for simplicity, since the idea of using this algorithm is to generate multiple model predictions (realisations) anyways, which then get sent to a summarising function. Additionally, the use of `mlr3` model types allowed for a more simple model prediction script since the folks at `mlr3` generate model predictions all in the same object slots no matter the model type. This worked wonderfully and was easy to provide a simple interface to, much like how the original `rdsmart` package worked well with `caret` models.

2.  The methods for sampling polygons all stayed the same, just using the `terra` package implementation over using `raster` and `sp` packages. While writing this, I found it simpler to roll the `.sampler()` function into the `.getVirtualSamples()` function, as it provided a cleaner comparison to the `.getStrtifiedVirtualSamples()` function. Additionally, rather than using `data.table::rbindlist()` at the end of the sampling, a `lapply()` iteration loop was created where it binds rows to itself through each iteration. This eliminated the use of the `data.table` and `foreach` packages. Messages indicating which polygon is being sampled is given throughout as well, moreso for debugging.

3.  Map predictions are made using the `terra` package. In a previous version of my coding, I implemented a tiling system to load smaller tiles to be predicted using the `stars` package; however, recent updates to the model prediction strategies used in `terra` made my own method obsolete. This also removed the dependency of `stars`. Trained `mlr3` learners (models) are passed to `terra`'s predict function and predictions are carried out across each row. In some model realisations, classes are dropped in order to create models. The missing classes are added back into the model predictions to generate the probabilities of those classes being present (which would be 0). There is an option to set the number of cores for the model prediction, but currently there is a lot of overhead on parallelism in `terra`, which makes it not very efficient to do but can be re-examined in the future ([see this post on GitHub](https://github.com/rspatial/terra/issues/178)).

4.  Since C++ computing is implemented, we no longer have to worry about the number of cores a computer system has. This argument is removed in all contexts.

5.  The `summarise` function uses `sort` and `order` to calculate the most probable classes and determine their respective probabilities. C++ versions of these commands are implemented to significantly reduce processing time. I don't have a lot of experience with C++ so I might not have the most efficient methods coded, but they currently work to mimic the results from base R implementations of the same functions.

6.  Messages pop up during summarizing that indicate what is going on at a given point in time.

Other functions, like the `sort_stack_values()` have not been recreated here, though it wouldn't be difficult to do that either if it is necessary.

I also want to point out that this package has different package dependencies than the original `rdsmart` package. The full list is:

-   `tidyverse`

-   `terra`

-   `mlr3verse`

-   `corrplot`

-   `Rcpp` (only for compiling the `sort_cpp` and `order_cpp` functions)

The `tidyverse` packages that are used throughout include:

-   `dplyr`

-   `tidyr`

-   `magrittr`

-   `tibble`

On the topic of package dependencies, I will also note here that when a model type is chosen, `mlr3` will attempt to load it and set the model arguments (args.model) before generating samples. If the model you are trying to load is not installed, `mlr3` will attempt to install and load it before proceeding. If it fails, it will be on the user to install the modelling package.

If any of the original developers of the `rdsmart()` package get a chance to read/try this, I would greatly appreciate any feedback and criticisms that you may have! I didn't do this to step on any of your toes and I hope this isn't coming off as such. I had been working with the `terra` package recently and thought that this would be a good exercise at trying some new things out.

Finally, I don't know how to write packages properly, so in this set of code here each file is simply `source()`'d in order to carry out the computations. This is fine if it's just me using it, but if this ever gains traction I'll look into packaging this up properly.
