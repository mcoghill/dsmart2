# dsmart2

## A faster implementation of the dsmart algorithm

Much of the original code can be found at this bitbucket repository: <https://bitbucket.org/brendo1001/dsmart/src/master/>

This has been my attempt at creating a DSMART algorithm that runs faster than the original implementation. I found that this was necessary for large datasets, as the [`raster`](https://github.com/rspatial/raster) package was proving to be insufficient for processing in my own work. There has been a lot of recent work on the [`terra`](https://github.com/rspatial/terra) package, which is the `raster` package replacement and uses C++ computing to significantly increase processing speed. Since the DSMART package used the `raster` package throughout its scripting, it was necessary to go through and rewrite the steps to produce the final outputs using the `terra` package. What was once a process that took a few days with the `raster` implementation of DSMART has turned into a process that takes a few hours with the `terra` implementation.

Some important notes on the changes made, starting in the `disaggregate()` coding:

1.  If a repeated cross validation model is used, then the number of model realisations is reduced to 1. This was a liberty that I took upon myself since I figured that in the end it would take a very long time to create multiple repeated cross validation models, and if the models are in an essence "repeated", then I figured that would suffice in place of multiple model realisations. This is likely inaccurate and can be easily removed in lines 205 - 208. This is the only major deviation from the main methodology, all other changes are performed to increase efficiency of processing and to promote the implementation of the `terra` package.

2.  The methods for sampling polygons all stayed the same, just using the `terra` package implementation over using `raster` and `sp` packages. While writing this, I found it simpler to roll the `.sampler()` function into the `.getVirtualSamples()` function, as it provided a cleaner comparison to the `.getStrtifiedVirtualSamples()` function. Additionally, rather than using `data.table::rbindlist()` at the end of the sampling, a `lapply()` iteration loop was created where it binds rows to itself through each iteration. This eliminated the use of the `data.table` and `foreach` packages. Messages indicating which polygon is being sampled is given throughout as well, moreso for debugging.

3.  Map predictions are made using a tiling system (scripts `predict_landscape.R` and `tile_index.R`). In my own objectives, I had over 100 covariate rasters and they were all rather large in size. This proved to be problematic when it came to map predictions, since I could not load all of the covariates into memory to perform the predictions. What happens instead is that the entire study area is broken up into smaller tiles (hard set at 500 x 500 pixels, though this can be adjusted), then model predictions are carried out on the tiles, and finally after the model predictions are finished the tiles are mosaicked back together. This method uses the [`stars`](https://github.com/r-spatial/stars) package to load raster tiles into memory. Depending on if the prediction type was set to "prob" or not, different outputs of the prediction process are curated: "prob" returns a character vector of file paths of the mosaicked rasters, and "raw" will return a mosaicked raster of the raw predictions. Both ways are taken care of when `dsmart()` is called.

4.  Since C++ computing is implemented, we no longer have to worry about the number of cores a computer system has. This argument is removed in all contexts.

5.  An argument in the `predict_landscapte` function is `mask`, which is simply a character vector of the name of the raster covariate layer used to mask the predicted area where necessary. If this argument is left blank (or NULL), it is automatically calculated as the raster covariate layer with the greatest number of data pixels. This argument is followed through in both `disaggregate` and `dsmart` functions, and is automatically calculated in the `disaggregate` function if it is not specified.

6.  General messages about the tiling outputs have been incorporated. Additionally, messages pop up during summarizing that indicate what is going on as well.

Other functions, like the `sort_stack_values()` have not been recreated here, though it wouldn't be difficult to do that either if it is necessary. One thing that I can think of that may be a future issue would be if there were many realisations with large file sizes, then tiling may be necessary to incorporate in order to complete some of the `summarise()` steps, though until it becomes a larger issue I don't really want to incorporate that.

I also want to point out that this package has different package dependencies than the original `rdsmart` package. The full list is:

-   `tidyverse`

-   `terra`

-   `stars`

I also make heavy use of the `sf` package, but that is depended on by the `stars` package when installed and loaded, hence why it is not included in that list. The `tidyverse` packages that are used throughout include:

-   `dplyr`

-   `tidyr`

-   `magrittr`

-   `tibble`

On the topic of package dependencies, I will also note here that either `caret` or `C50` is required for the modelling, though I didn't include it in the above list. It depends on the users input, and if `caret` is used then subsequent modelling packages may be required as well.

If any of the original developers of the `rdsmart()` package get a chance to read/try this, I would greatly appreciate any feedback and criticisms that you may have! I didn't do this to step on any of your toes and I hope this isn't coming off as such. I had been working with the `terra` package recently and thought that this would be a good exercise at trying some new things out.

Finally, I don't know how to write packages properly, so in this set of code here each file is simply `source()`'d in order to carry out the computations. This is fine if it's just me using it, but if this ever gains traction I'll look into packaging this up properly.
