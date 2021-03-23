#' Disaggregating and harmonising soil map units through resampled
#' classification trees
#'
#' This function, together with companion function \code{summarise} implements the
#' DSMART (Disaggregating and harmonising soil map units through resampled
#' classification trees) algorithm as described in Odgers et al. (2014). This is
#' the workhorse function that involves multiple resampling, C5 decision tree
#' model fitting, and subsequent mapping in order to realize potential candidate
#' soil classes within aggregated soil mapping units. There is also added
#' facility to incorporate point observed data into the algorithm too.
#'
#' There is no prescription for the layers that compose \code{covariates} other
#' than that they should represent the \emph{scorpan} factors (McBratney
#' \emph{et al.}, 2003) as effectively as possible. The order of the layers in
#' the SpatRaster is not important. See \code{data(dsT_covariates)} for an
#' example.
#'
#' DSMART assumes, but does not currently check, that \code{covariates},
#' \code{polygons} and any \code{observations} are projected to the same
#' coordinate reference system.
#'
#' DSMART assumes that the soil classes in \code{composition} and those in
#' \code{observations} belong to the same taxonomic level of the same soil
#' classification system. For example if the soil classes in \code{composition}
#' are all soil series, those in \code{observations} should also be series
#' rather than, for example, orders or great groups.
#'
#' \code{dsmart} produces a number of outputs which are saved to subdirectories
#' in the current working directory. The base folder for the outputs is
#' \code{output}. Inside this folder, rasters of the realisations are saved into
#' the \code{realisations} subfolder and the classification models are saved
#' into the \code{models} subfolder.
#'
#' @param covariates A \code{SpatRaster} of \emph{scorpan} environmental
#'   covariates to calibrate the \code{C50} classification trees against. See
#'   \emph{Details} for more information.
#' @param polygons A \code{SpatVector} containing the soil map unit polygons 
#'   that will be disaggregated. The first field of the data frame must be an 
#'   integer that identifies each polygon.
#' @param composition A \code{data.frame} that contains information on the
#'   soil-class composition of each polygon in \code{polygons}. Each row
#'   contains information about one soil class component of one polygon, which
#'   belongs to one soil map unit. First field contains the integer that
#'   identifies the polygon. Second field contains a code that identifies the
#'   soil map unit that the polygon belongs to.
#'
#'   If \code{strata = NULL} (the default), third column contains a code that
#'   identifies the soil class and fourth column contains a number in the range
#'   \code{(0, 100)} that identifies the proportion of the \strong{map unit}
#'   that the soil class corresponds to. See the example data
#'   \code{data(dalrymple_composition)}.
#'
#'   If \code{strata} is a \code{SpatRaster}, third column contains an integer
#'   that identifies the stratum in \code{strata}, fourth column contains a code
#'   that identifies the soil class and fifth column contains a number in the
#'   range \code{(0, 100)} that identifies the proportion of the
#'   \strong{stratum} that the soil class corresponds to.
#' @param rate An integer that identifies the number of virtual samples to draw
#'   from each polygon in each realisation. If \code{method.sample =
#'   "by_polygon"}, the number of samples to draw from each polygon in
#'   \code{polygons}. If \code{method.sample = "by_area"}, the sampling density
#'   in number of samples per square kilometre.
#' @param reals An integer that identifies the number of realisations of the
#'   soil class distribution that DSMART should compute.
#' @param observations \emph{Optional} A \code{data.frame} that contains actual
#'   observations of the soil class at locations across the soil map area. These
#'   data augment the virtual samples and are used in each realisation. Each row
#'   contains information about one soil class observation. First and second
#'   fields contain the \emph{x-} and \emph{y-}components of the observation's
#'   spatial location. Third field is the soil class code. See \emph{Details}.
#' @param method.sample Identifies the sampling method. Valid values are
#'   \code{"by_polygon"} (the default), in which case the same number of samples
#'   are taken from each polygon; or \code{"by_area"}, in which case the number
#'   of samples per polygon depends on the area of the polygon.
#' @param method.allocate Method of allocation of virtual samples to soil
#'   classes. Valid values are \code{"weighted"}, for weighted-random allocation
#'   to a soil class from within the virtual sample's map unit;
#'   \code{"random_mapunit"}, for completely random allocation to a soil class
#'   from within the virtual sample's map unit; and \code{"random_all"}, for
#'   completely random allocation to a soil class from within the entire map
#'   area.
#' @param method.model Method to be used for the classification model. If no
#'   value is passed, a C5.0 decision tree is built. Otherwise, the value must
#'   match a valid 'method' argument in the caret::train function.
#' @param args.model A list of arguments to be passed to the caret::train
#'   function. The list can include a trainControl object, arguments to be
#'   passed directly to the train function and arguments to be passed to the
#'   predictive model.
#' @param strata \emph{optional} An integer-valued \code{SpatRaster} that will
#'   be used to stratify the allocation of virtual samples to soil classes.
#'   Integer values could represent classes of slope position (e.g. crest,
#'   backslope, footslope, etc.) or land use (e.g. cropland, native vegetation,
#'   etc.) or some other variable deemed to be an important discriminator of the
#'   occurrence of soil classes within a map unit.
#' @param outputdir A character string that identifies the location of the main
#'   output directory. The folder \code{output} and its subfolders will be
#'   placed here. Default is the current working directory, \code{getwd()}.
#' @param stub \emph{optional} A character string that identifies a short name
#'   that will be prepended to all output.
#' @param factors A character vector with the names of the covariates that
#'   should be treated as factors.
#' @param type A character vector to specify the type of the predictions. By
#'   default, raw class predictions are used ("response"). If set to "prob", a multiband
#'   SpatRaster with class probabilities will be produced for each realisation.
#'
#' @return A list that contains metadata about the current run of
#'   \code{disaggregate}.
#'
#' @examples
#' # Load datasets
#' data(dalrymple_composition)
#' data(dalrymple_covariates)
#' data(dalrymple_observations)
#' data(dalrymple_polygons)
#'
#' # Run disaggregate without adding observations
#' disaggregate(dalrymple_covariates, dalrymple_polygons, dalrymple_composition,
#'  rate = 15, reals = 10)
#'
#' # Run disaggregate with extra observations
#' disaggregate(dalrymple_covariates, dalrymple_polygons, dalrymple_composition,
#'  observations = dalrymple_observations, rate = 15, reals = 10)
#'
#' @references McBratney, A.B., Mendonca Santos, M. de L., Minasny, B., 2003. On
#'   digital soil mapping. Geoderma 117, 3--52. doi:
#'   \href{http://dx.doi.org/10.1016/S0016-7061(03)00223-4}{10.1016/S0016-7061(03)00223-4}
#'   
#'   Odgers, N.P., McBratney, A.B., Minasny, B., Sun, W., Clifford, D., 2014.
#'   DSMART: An algorithm to spatially disaggregate soil map units, \emph{in:}
#'   Arrouays, D., McKenzie, N.J., Hempel, J.W., Richer de Forges, A.,
#'   McBratney, A.B. (Eds.), GlobalSoilMap: Basis of the Global Spatial Soil
#'   Information System. Taylor & Francis, London, pp. 261--266.
#'
#'   Odgers, N.P., Sun, W., McBratney, A.B., Minasny, B., Clifford, D., 2014.
#'   Disaggregating and harmonising soil map units through resampled
#'   classification trees. Geoderma 214, 91--100. doi:
#'   \href{http://dx.doi.org/10.1016/j.geoderma.2013.09.024}{10.1016/j.geoderma.2013.09.024}
#'
#' @export

disaggregate <- function(
  covariates, polygons, composition, rate = 15, reals = 100, 
  observations = NULL, method.sample = "by_polygon", method.allocate = "weighted", 
  method.model = NULL, args.model = NULL, strata = NULL, outputdir = getwd(), 
  stub = NULL, factors = NULL, type = "response", predict = TRUE) {
  
  # Requires the following packages:
  invisible(sapply(c("tidyverse", "terra", "mlr3verse", "corrplot"), 
                   require, character.only = TRUE))
  
  # Remove terra progress bars
  def_ops <- terra:::spatOptions()$progress
  terraOptions(progress = 0)
  on.exit(terraOptions(progress = def_ops))
  
  # Source external functions used in this script
  source("./R/helpers.R")
  
  # Create list to store output
  output <- base::list()
  
  # Save start time
  output$timing <- base::list(start = Sys.time())
  
  # Check arguments before proceeding
  messages <- c("Attention is required with the following arguments:\n")
  if(!(class(covariates) == "SpatRaster")) {
    messages <- append(messages, "'covariates': Not a valid SpatRaster.\n")
  }
  if(!(class(polygons) == "SpatVector")) {
    messages <- append(messages, "'polygons': Not a valid SpatVector\n")
  }
  if(!(class(composition) == "data.frame")) {
    messages <- append(messages, "'composition': Not a valid data.frame.\n")
  }
  if(rate <= 0) {
    messages <- append(messages, "'n': Value must be greater than 0.\n")
  }
  if(reals <= 0) {
    messages <- append(messages, "'reals': Value be greater than 0.\n")
  }
  if(!(is.null(observations))) {
    if(!(class(observations) == "data.frame")) {
      messages <- append(messages, "'observations': Not a valid data.frame.\n")
    }
  }
  if(!(file.exists(outputdir))) {
    messages <- append(messages, "'outputdir': Output directory does not exist.\n")
  }
  if(!(is.null(strata))) {
    if(!(class(strata) == "SpatRaster")) {
      messages <- append(messages, "'strata': Not a valid SpatRaster\n")
    }
  }
  if(length(messages) > 1) {
    stop(messages)
  }
  
  # If method.model is a character value, enforce proper learner detection
  valid_learners <- mlr3extralearners::list_mlr3learners(filter = list(
    class = "classif", predict_types = type,
    properties = "multiclass"), select = c("name", "id", "required_packages"))
  
  if(is.character(method.model) & length(method.model) == 1) {
    
    filt_learner <- valid_learners[name == method.model]
    if(!nrow(filt_learner)) {
      filt_learner <- valid_learners[id == method.model]
    } else {
      method.model <- paste0("classif.", method.model)
      filt_learner <- valid_learners[id == method.model]
    }
    if(!nrow(filt_learner)) {
      stop("This is not a valid learner for your use. See here for a list of
             \rvalid learner names and ids:
             \rmlr3extralearners::list_mlr3learners(filter = list(
      \r  class = 'classif', predict_types = '", type, "',
      \r  properties = 'multiclass'))")
    }
  } else {
    method.model <- "classif.C50"
  }
  
  # Check that required package is installed and loaded for use
  pkgs <- unlist(valid_learners[id == method.model]$required_packages)
  pkg_require <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
  if(length(pkg_require)) mlr3extralearners::install_learners(method.model)
  invisible(sapply(pkgs, require, character.only = TRUE))
  
  # Check that model arguments will work as expected
  model <- lrn(method.model, predict_type = type)
  model$param_set$values <- args.model
  
  # Set stub to "" if NULL
  if(is.null(stub)) {
    stub <- ""
  }  else if(stub == "") {
    stub <- ""
  }  else if(!(substr(stub, nchar(stub), nchar(stub)) == "_")) {
    stub <- paste0(stub, "_")
  }
  
  # Save function call
  output$call <- base::match.call()
  
  # Save parameters
  output$parameters <- base::list(
    rate = rate, reals = reals, method.sample = method.sample, 
    method.allocate = method.allocate, method.model = method.model, 
    args.model = args.model, stub = stub, factors = factors, type = type)
  
  # Create subdirectories to store results in
  outputdir <- file.path(outputdir)
  ndirs <- length(dir(outputdir, pattern = "output"))
  if(ndirs > 1) {
    last <- max(as.numeric(
      gsub("\\D", "", dir(outputdir, pattern = "output"))), na.rm = TRUE)
    subdir <- paste0("output_", formatC(last, width = 3, format = "d", flag = 0))
    nsubdirs <- length(dir(file.path(outputdir, subdir)))
    if(nsubdirs > 1) 
      subdir <- paste0("output_", formatC(last + 1, width = 3, format = "d", flag = 0))
  } else {
    nsubdirs <- length(dir(file.path(outputdir, "output")))
    subdir <- ifelse(nsubdirs > 1, 
                     paste0("output_", formatC(1, width = 3, format = "d", flag = 0)), 
                     "output")
  }
  dir.create(file.path(outputdir, subdir), showWarnings = FALSE)
  dir.create(file.path(outputdir, subdir, "realisations"), showWarnings = FALSE)
  dir.create(file.path(outputdir, subdir, "models"), showWarnings = FALSE)
  dir.create(file.path(outputdir, subdir, "probabilities"), showWarnings = FALSE)
  
  # Save output locations
  output$locations <- list(
    root = file.path(outputdir, subdir), 
    realisations = c(), 
    models = c())
  
  # Rename composition column names
  if(!(is.null(strata))) {
    names(composition) <- c("poly", "mapunit", "stratum", "soil_class", "proportion")
  } else {
    names(composition) <- c("poly", "mapunit", "soil_class", "proportion")
  }
  
  # Since we only really need the first column of the polygons SpatialPolygonsDataFrame,
  # drop all its other attributes and rename the first column to "poly"
  polygons <- polygons[, 1] %>% stats::setNames("poly")
  
  # Make sure that there are no missing values in map unit composition
  polys_to_remove <- composition %>% stats::complete.cases() %>% 
    magrittr::equals(FALSE) %>% base::which() %>% composition$poly[.] %>% 
    base::unique()
  
  # Remove the polygons that contain missing information from both the polygons
  # SpatVect and the composition data frame
  if(length(polys_to_remove) > 0) {
    polygons <- polygons[!polygons$poly %in% polys_to_remove]
    composition <- composition %>% dplyr::filter(!(poly %in% polys_to_remove))
    
    # Let the user know about what we've done
    warning(base::paste0(
      "The following polygons were removed from further analysis because they have incomplete or undefined map unit compositions: ", 
      base::paste(base::as.character(polys_to_remove), collapse = ", ")))
  }
  
  # Write covariate names to file
  write.table(names(covariates), file.path(outputdir, subdir, 
                                           paste0(stub, "covariate_names.txt")), 
              quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE)
  
  # Get samples for all realisations
  message(Sys.time(), " Generating samples for ", reals, " realisations")
  
  if(is.null(strata)) {
    # Get samples without allocation stratification
    samples <- .getVirtualSamples(
      covariates, polygons, composition, 
      n.realisations = reals, rate = rate, method.sample = method.sample, 
      method.allocate = method.allocate)
  } else {
    samples <- .getStratifiedVirtualSamples(
      covariates, polygons, 
      composition, strata = rast(strata), n.realisations = reals, rate = rate, 
      method.sample = method.sample, method.allocate = method.allocate)
  }
  
  # Make sure that there are no missing values in samples
  samples <- tidyr::drop_na(samples)
  
  # If there are observations, get their covariates
  if(!(is.null(observations))) {
    names(observations) <- c("x", "y", "class")
    observations <- .observations(observations, covariates)
    write.table(observations, file.path(outputdir, subdir, 
                                        paste0(stub, "observations_with_covariates.txt")), 
                sep = ",", quote = FALSE, col.names = TRUE, row.names = FALSE)
  }
  
  # Write samples to text file
  write.table(samples, file.path(outputdir, subdir, paste0(stub, "virtual_samples.txt")), 
              sep = ",", quote = FALSE, col.names = TRUE, row.names = FALSE)
  
  
  # Column number of samples' first covariate
  startcol <- numeric(0)
  if(is.null(strata)) {
    startcol = 8
  } else {
    startcol = 9
  }
  
  # Create correlation matrix - first need variances to see if there are vars
  # with zero variance
  variance <- samples[, startcol:ncol(samples)] %>% 
    dplyr::mutate(across(where(is.factor), as.integer)) %>% 
    apply(2, var, na.rm = TRUE)
  var_valid <- names(variance[variance > 0])
  
  cor_mat <- samples[, startcol:ncol(samples)] %>% 
    dplyr::select(any_of(var_valid)) %>% 
    dplyr::mutate(across(where(is.factor), as.integer)) %>% 
    cor(method = "pearson")
  
  task_cor <- TaskClassif$new(
    id = "correlation", 
    backend = samples[, (startcol - 1):ncol(samples)] %>% 
      dplyr::select(soil_class, all_of(var_valid)) %>% 
      dplyr::mutate(across(where(is.factor) & !all_of("soil_class"), as.integer),
                    soil_class = as.factor(soil_class)),
    target = "soil_class")
  
  cor <- flt("find_correlation", method = "pearson")
  cor$calculate(task_cor)
  cor_mat_out <- as.data.frame(cor_mat) %>% 
    rownames_to_column("feature") %>% 
    merge(as.data.table(cor)) %>% 
    dplyr::mutate(score = 1 - score) %>% 
    dplyr::rename(pearson_score = score)
  
  # Remove correlated variables 
  samples <- samples[, startcol:ncol(samples)] %>% 
    dplyr::select(
      as.data.table(cor) %>% 
        dplyr::filter(score >= (1 - 0.95)) %>% 
        dplyr::pull(feature)) %>% 
    bind_cols(samples[, 1:(startcol - 1)], .)
  
  # Save correlation image and matrix
  cor_out <- function(cor_matrix) {
    jpeg(file.path(outputdir, subdir, paste0(stub, "correlation_matrix.jpg")), 
         width = 10, height = 10, units = "in", pointsize = 12, res = 300)
    corrplot.mixed(cor_mat, lower = "number", upper = "circle", tl.pos = "lt",
                   diag = "u", order = "AOE", lower.col = "black", 
                   title = "Predictor Correlation", cl.cex = 1, number.cex = 0.83,
                   mar = c(0, 0, 2, 0))
    dev.off()
    return(invisible())
  }
  cor_out(cor_mat)
  write.csv(cor_mat_out, file.path(outputdir, subdir, paste0(stub, "correlation_matrix.csv")), 
            row.names = FALSE)
  
  # We submit the target classes to C5.0 as a factor. To do that, we need to 
  # make sure that the factor has the full set of levels, since due to the 
  # randomness of the allocation procedure we can't assume that they are all
  # represented in the samples drawn for a particular realisation. If the levels
  # are not specified correctly and they are not all represented in all sets of
  # samples, it is possible that the integer values that the class predictions
  # are coded to by raster::predict will not refer to the same soil class from
  # realisation to realisation.
  levs <- as.character(unique(composition$soil_class))
  
  if(!(is.null(observations))) {
    # Union map unit levels with observation levels
    levs <- base::union(levs, as.character(unique(observations$soil_class)))
  }
  
  # Sort levels
  levs <- sort(levs)
  
  # Generate lookup table
  lookup <- data.frame(name = levs, code = 1:length(levs), stringsAsFactors = FALSE)
  
  # Write lookup table to file
  write.table(lookup, file.path(outputdir, subdir, paste0(stub, "lookup.txt")), 
              sep = ",", quote = FALSE, col.names = TRUE, row.names = FALSE)
  output$locations$lookup <- file.path(outputdir, subdir, paste0(stub, "lookup.txt"))
  
  # Process realisations
  for (j in 1:reals) {
    message("\n", Sys.time(), " Realisation ", j)
    
    # Extract sample covariates for the current realisation
    s <- samples[which(samples$realisation == j), startcol:ncol(samples)]
    soil_class <- as.character(samples$soil_class[which(samples$realisation == j)])
    
    # Concatenate observations if they exist
    if(!(is.null(observations))) {
      o <- observations %>% dplyr::select(all_of(names(s)))
      s <- rbind(s, o)
      soil_class <- append(soil_class, as.character(observations$soil_class))
    }
    
    # Sort levels and convert soil_class back to factor
    soil_class <- factor(soil_class, levels = levs)
    
    # Convert designated covariates to factors.
    if(!is.null(factors)) {
      s <- s %>% dplyr::mutate(across(factors, as.factor))
    }
    
    # Test if there are levels in soil_class with 0 cases
    # (this is used for a workaround, see below).
    zeroes <- sum(table(soil_class) == 0) > 0
    
    # If there are levels with 0 cases, make a reclassification table for the prediction raster.
    rclt <- cbind(c(1:sum(table(soil_class) != 0)), 
                  c(1:length(table(soil_class)))[table(soil_class) != 0])
    
    # Create mlr3 task from data
    task <- TaskClassif$new(id = "data", backend = cbind(soil_class = soil_class, s), 
                            target = "soil_class")$droplevels()
    
    # Fit model
    model <- model$train(task)
    model_out <- model$model
    out <- utils::capture.output(model_out)
    
    # Save model to text file
    model_dir <- file.path(
      outputdir, subdir, "models", 
      paste0(stub, "model_", formatC(j, width = nchar(reals), format = "d", flag = "0")))
    cat(out, file = paste0(model_dir, ".txt"), sep = "\n", append = FALSE)
    
    # Save model to RDS file
    saveRDS(model_out, file = paste0(model_dir, ".rds"))
    output$locations$models <- c(
      output$locations$models, paste0(model_dir, ".rds")
    )
    
    # Model prediction
    if(predict) {
      # Create predictive function to work with mlr3 learner models for either
      # probability or classification modelling
      pf <- function(model, ...) {
        if(model$predict_type == "prob") {
          p <- model$predict_newdata(...)$data$prob
          if(length(levs) != ncol(p)) {
            missing <- setdiff(levs, colnames(p))
            pm <- matrix(0, ncol = length(missing), nrow = nrow(p), dimnames = list(NULL, missing))
            p <- cbind(p, pm)
            p <- p[, levs]
          }
          p
        } else {
          model$predict_newdata(...)$data$response
        }
      }
      if(type == "prob") {
        r1 <- predict(
          subset(covariates, task$feature_names), model, na.rm = TRUE, fun = pf,
          filename = file.path(
            outputdir, subdir, "realisations",
            paste0(stub, "realisation_", formatC(
              j, width = nchar(reals), format = "d", flag = "0"), ".tif")),
          overwrite = TRUE, wopt = list(names = lookup$name)
        )
      } else {
        # If classification prediction, it will predict for the non-dropped
        # classes. Those classes will be attributed and then an index made
        # afterwards.
        tmp1 <- terra::predict(
          subset(covariates, task$feature_names), model, na.rm = TRUE, fun = pf,
          filename = tempfile(pattern = "spat_", fileext = ".tif"),
          wopt = list(
            datatype = "INT2S",
            names = paste0(stub, "realisation_", formatC(
              j, width = nchar(reals), format = "d", flag = "0"))))
        
        # Convert integer classes to factor levels
        tmp2 <- as.factor(tmp1)
        
        # Get the factor levels and labels from the converted raster
        tmp_labels <- data.frame(levels = levels(tmp2)[[1]]$levels,
                                 code = levels(tmp2)[[1]]$labels)
        
        # Merge the lookup tables so that factor levels show which soil class
        # corresponds to a given level and lookup code
        tmp_lookup <- lookup %>%
          merge(tmp_labels) %>%
          dplyr::select(levels, name, code)
        
        # Write the merged index
        write.csv(tmp_lookup, file.path(
          outputdir, subdir, "realisations",
          paste0(stub, "realisation_", formatC(j, width = nchar(reals),
                                               format = "d", flag = "0"), "_lut.csv")),
          row.names = FALSE)
        
        # Apply the new factor levels to the raster
        levels(tmp2) <- tmp_lookup[, c(1, 2)]
        
        # Write out the completed raster (does so with associated .xml file
        # identifying the class labels)
        r1 <- writeRaster(tmp2, filename = file.path(
          outputdir, subdir, "realisations",
          paste0(stub, "realisation_", formatC(j, width = nchar(reals),
                                               format = "d", flag = "0"), ".tif")),
          overwrite = TRUE, wopt = list(datatype = "INT1U"))
      }
      
      # Save output realisation paths
      output$locations$realisations <- c(
        output$locations$realisations, sources(r1)[, "source"])
      
      # Remove temporary files from system
      suppressWarnings(terra::tmpFiles(old = TRUE, remove = TRUE))
    }
  }
  
  # Save finish time
  output$timing$finish <- Sys.time()
  
  # Return output
  return(output)
}

#END