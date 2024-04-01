# # [SETUP] -----------------------------------------------------------------
# # - Packages ----------------------------------------------------------------
# # CRAN packages
# chr_pkg <- c(
#   'bvls', 'fastglm', 'weights' #Regression models
#   , 'devtools' #GitHub packages (temp)
#   , 'readr' #Read data (temp)
#   , 'dplyr', 'tidyr', 'purrr' #Data wrangling
# )
# 
# # Git packages
# chr_git <- c(
#   'CaoBittencourt' = 'atlas.aeq' #Human capital indispensability coefficient
# )
# 
# # Activate / install CRAN packages
# lapply(
#   chr_pkg
#   , function(pkg){
# 
#     if(!require(pkg, character.only = T)){
# 
#       install.packages(pkg)
# 
#     }
# 
#     require(pkg, character.only = T)
# 
#   }
# )
# 
# # Activate / install Git packages
# Map(
#   function(git, profile){
# 
#     if(!require(git, character.only = T)){
# 
#       install_github(
#         paste0(profile, '/', git)
#         , upgrade = F
#         , force = T
#       )
# 
#     }
# 
#     require(git, character.only = T)
# 
#   }
#   , git = chr_git
#   , profile = names(chr_git)
# )

# [MATCHING FUNCTIONS] -------------------------------------------------------------
# - Regression weights --------------------------------------------
fun_match_weights <- function(
    dbl_var
    , chr_weights = c('unweighted', 'linear', 'quadratic', 'speciality-root', 'attribute-eqvl')
    , dbl_generality = NULL
){
  
  # Argument validation
  stopifnot(
    "'dbl_var' must be a numeric vector or matrix." =
      is.numeric(dbl_var)
  )
  
  stopifnot(
    "'chr_weights' must be either 'unweighted', 'linear', 'quadratic', 'speciality-root', or 'attribute-eqvl'." = 
      any(
        chr_weights == 'unweighted',
        chr_weights == 'linear',
        chr_weights == 'quadratic', 
        chr_weights == 'speciality-root', 
        chr_weights == 'attribute-eqvl'
      )
  )
  
  # Weighting methods
  if(chr_weights == 'linear'){
    # Linear weights
    dbl_var -> dbl_weights
    
  } else if(chr_weights == 'quadratic'){
    # Quadratic weights
    dbl_var ^ 2 -> dbl_weights
    
  } else if(chr_weights == 'speciality-root'){
    # Speciality root weights
    dbl_var ^ (
      1 / (1 - dbl_generality)
    ) -> dbl_weights
    
  } else if(chr_weights == 'attribute-eqvl'){
    # Attribute equivalence weights
    fun_aeq_aequivalence(
      dbl_profile = dbl_var
      , dbl_generality = 
        dbl_generality
      # , dbl_scale_lb = 0
    ) -> dbl_weights
    
  } else {
    # Unweighted
    NULL -> dbl_weights
    
  }
  
  # Output
  return(dbl_weights)
  
}

# - Vectorized regression weights -----------------------------------------
fun_match_vweights <- function(
    df_data_cols
    , chr_weights = NULL
    , dbl_generality = NULL
){
  
  # Arguments validation
  stopifnot(
    "'df_data_cols' must be a data frame." =
      is.data.frame(df_data_cols)
  )
  
  # Data wrangling
  df_data_cols %>% 
    select(where(
      is.numeric
    )) -> df_data_cols
  
  # Map weights function
  if(length(dbl_generality)){
    
    map2(
      .x = df_data_cols
      , .y = dbl_generality
      , .f = fun_match_weights
      , chr_weights = 
        chr_weights
    ) -> list_weights
    
  } else {
    
    map(
      .x = df_data_cols
      , .f = fun_match_weights
      , chr_weights = 
        chr_weights
    ) -> list_weights 
    
  }
  
  if(chr_weights != 'unweighted'){
    
    list_weights %>% 
      as_tibble() ->
      df_weights
    
  } else {
    
    list_weights %>% 
      list_c() -> 
      df_weights
    
  }
  
  rm(list_weights)
  
  # Output
  return(df_weights)
  
}

# - BVLS regression matching ----------------------------------------------
fun_match_bvls <- function(
    df_data_cols
    , dbl_query
    , df_weights = NULL
){
  
  # Arguments validation
  stopifnot(
    "'df_data_cols' must be a data frame." = 
      is.data.frame(df_data_cols)
  )
  
  stopifnot(
    "'dbl_query' must be numeric." =
      all(
        is.numeric(dbl_query)
        , length(dbl_query) ==
          nrow(df_data_cols)
      )
  )
  
  stopifnot(
    "'df_weights' must be either NULL or a numeric data frame." = 
      any(
        all(map_lgl(df_weights, is.numeric))
        , is.null(df_weights)
      )
  )
  
  # BVLS regression
  if(!length(df_weights)){
    
    # Run BVLS regression matching without weights
    map_dbl(
      .x = df_data_cols
      , ~ 
        coef(bvls(
          as.matrix(.x)
          , dbl_query[,]
          , bl = 0
          , bu = 1
        ))
    ) -> dbl_similarity
    
  } else {
    
    # Add weights to regression 
    sqrt(df_weights) ->
      df_weights
    
    df_data_cols *
      df_weights ->
      df_data_cols
    
    # BVLS regression matching with weights
    map2_dbl(
      .x = df_data_cols
      , .y = df_weights
      , ~
        coef(bvls(
          as.matrix(.x)
          , dbl_query[,] * .y
          , bl = 0
          , bu = 1
        ))
    ) -> dbl_similarity
    
  }
  
  # Output
  return(dbl_similarity)
  
}

# - Pearson correlation matching ----------------------------------------------
fun_match_pearson <- function(
    df_data_cols
    , dbl_query
    , df_weights = NULL
){
  
  # Arguments validation
  stopifnot(
    "'df_data_cols' must be a data frame." = 
      is.data.frame(df_data_cols)
  )
  
  stopifnot(
    "'dbl_query' must be numeric." = 
      all(
        is.numeric(dbl_query)
        , length(dbl_query) ==
          nrow(df_data_cols)
      )
  )
  
  stopifnot(
    "'df_weights' must be either NULL or a numeric data frame." = 
      any(
        all(map_lgl(df_weights, is.numeric))
        , is.null(df_weights)
      )
  )
  
  # Pearson correlation
  if(!length(df_weights)){
    
    # Pearson correlation matching without weights
    map_dbl(
      .x = df_data_cols
      , ~ (
        1 +
          wtd.cors(
            dbl_query[,]
            , .x
          )[,]
      ) / 2
    ) -> dbl_similarity
    
  } else {
    
    # Pearson correlation matching with weights
    map2_dbl(
      .x = df_data_cols
      , .y = df_weights
      , ~ (
        1 +
          wtd.cors(
            dbl_query[,]
            , .x
            , weight = .y
          )[,]
      ) / 2
    ) -> dbl_similarity
    
  }
  
  # Output
  return(dbl_similarity)
  
}

# - Spearman correlation matching -------------------------------------------


# - Kendau's tau correlation matching -------------------------------------------


# - Logistic regression matching ------------------------------------------
fun_match_logit <- function(
    df_data_cols
    , dbl_query
    , dbl_scale_ub
    , dbl_scale_lb
    , chr_method = c('logit', 'probit')
    , df_weights = NULL
){
  
  # Arguments validation
  stopifnot(
    "'df_data_cols' must be a data frame." = 
      is.data.frame(df_data_cols)
  )
  
  stopifnot(
    "'dbl_query' must be numeric." = 
      all(
        is.numeric(dbl_query)
        , length(dbl_query) ==
          nrow(df_data_cols)
      )
  )
  
  stopifnot(
    "'dbl_scale_ub' must be numeric." =
      is.numeric(dbl_scale_ub)
  )
  
  stopifnot(
    "'dbl_scale_lb' must be numeric." =
      is.numeric(dbl_scale_lb)
  )
  
  stopifnot(
    "'chr_method' must be either 'logit' or 'probit'." =
      any(
        chr_method == 'logit',
        chr_method == 'probit'
      )
  )
  
  stopifnot(
    "'df_weights' must be either NULL or a numeric data frame." = 
      any(
        all(map_lgl(df_weights, is.numeric))
        , is.null(df_weights)
      )
  )
  
  # Data wrangling
  dbl_scale_ub[[1]] -> dbl_scale_ub
  
  dbl_scale_lb[[1]] -> dbl_scale_lb
  
  as.integer(dbl_scale_ub) -> dbl_scale_ub
  
  as.integer(dbl_scale_lb) -> dbl_scale_lb
  
  chr_method[[1]] -> chr_method
  
  # Convert query to a Bernoulli variable
  list_c(map(
    .x = as.integer(dbl_query[,])
    , ~ rep(
      c(1,0), times = ceiling(c(
        .x, (dbl_scale_ub - dbl_scale_lb) - .x
      ))
    )
  )) -> int_query_bernoulli
  
  rm(dbl_query)
  
  # Convert data to a Bernoulli variable
  map(
    .x = df_data_cols
    , ~
      as.matrix(list_c(map(
        .x = as.integer(.x)
        , ~ rep(
          c(1,0), times = ceiling(c(
            .x, (dbl_scale_ub - dbl_scale_lb) - .x
          )) 
        )
      )))
  ) -> list_data_bernoulli
  
  rm(df_data_cols)
  
  # Logistic regression
  if(!length(df_weights)){
    
    # Run logistic regression matching without weights
    map_dbl(
      .x = list_data_bernoulli
      , ~
        coef(fastglmPure(
          x = .x
          , y = int_query_bernoulli
          , family = binomial(
            link = 'logit'
          )
        ))
    ) -> dbl_similarity
    
    exp(dbl_similarity) /
      (1 + exp(dbl_similarity)) ->
      dbl_similarity
    
  } else {
    
    # Repeat df_weights' rows
    df_weights[rep(
      1:nrow(df_weights)
      , each = 
        dbl_scale_ub - 
        dbl_scale_lb
    ), ] -> df_weights
    
    # Run logistic regression matching with weights
    map2_dbl(
      .x = list_data_bernoulli
      , .y = df_weights
      , ~
        coef(fastglmPure(
          x = .x
          , y = int_query_bernoulli
          , family = binomial(
            link = 'logit'
          ), weights = .y
        ))
    ) -> dbl_similarity
    
    exp(dbl_similarity) /
      (1 + exp(dbl_similarity)) ->
      dbl_similarity
    
  }
  
  # Output
  return(dbl_similarity)
  
}

# # - Euclidean matching ----------------------------------------------------
# fun_match_euclidean <- function(
#     df_data_cols
#     , dbl_query
#     , dbl_scale_ub = 100
#     , dbl_scale_lb = 0
#     , df_weights = NULL
# ){
#   
#   # Arguments validation
#   stopifnot(
#     "'df_data_cols' must be a data frame." = 
#       is.data.frame(df_data_cols)
#   )
#   
#   stopifnot(
#     "'dbl_query' must be numeric." =
#       all(
#         is.numeric(dbl_query)
#         , length(dbl_query) ==
#           nrow(df_data_cols)
#       )
#   )
#   
#   stopifnot(
#     "'df_weights' must be either NULL or a numeric data frame." = 
#       any(
#         all(map_lgl(df_weights, is.numeric))
#         , is.null(df_weights)
#       )
#   )
#   
#   # Data wrangling
#   dbl_scale_ub[[1]] -> dbl_scale_ub
#   
#   dbl_scale_lb[[1]] -> dbl_scale_lb
#   
#   t(df_data_cols) -> df_data_rows
#   
#   if(length(df_weights)){
#     
#     t(df_weights) -> df_weights
#     
#   }
#   
#   as_tibble(t(dbl_query)) -> df_query_rows
#   
#   rm(df_data_cols)
#   rm(dbl_query)
#   
#   # Euclidean distance
#   if(length(df_weights)){
#     
#     # Weighted Euclidean distance
#     df_query_rows[rep(
#       1, nrow(df_data_rows)
#     ), ] -> df_query_rows
#     
#     df_weights * (
#       df_data_rows - 
#         df_query_rows
#     ) ^ 2 -> df_dist
#     
#     rm(df_query_rows)
#     
#     sqrt(rowSums(df_dist)) -> df_dist
#     
#     # Normalize by maximum distance
#     df_dist / sqrt(rowSums(
#       df_weights * 
#         pmax(
#           dbl_scale_ub - df_data_rows,
#           df_data_rows - dbl_scale_lb
#         ) ^ 2
#     )) -> dbl_dist
#     
#   } else {
#     
#     # Unweighted Euclidean distance
#     df_query_rows[rep(
#       1, nrow(df_data_rows)
#     ), ] -> df_query_rows
#     
#     (
#       df_data_rows - 
#         df_query_rows
#     ) ^ 2 -> df_dist
#     
#     rm(df_query_rows)
#     
#     sqrt(rowSums(df_dist)) -> df_dist
#     
#     # Normalize by maximum distance
#     df_dist / sqrt(rowSums(
#       pmax(
#         dbl_scale_ub - df_data_rows,
#         df_data_rows - dbl_scale_lb
#       ) ^ 2
#     )) -> dbl_dist 
#     
#   }
#   
#   rm(df_dist)
#   rm(df_weights)
#   rm(df_data_rows)
#   
#   # Calculate similarity
#   1 - dbl_dist -> dbl_similarity
#   
#   # Output
#   return(dbl_similarity)
#   
# }

# - Euclidean matching ----------------------------------------------------
fun_match_euclidean <- function(
    df_data_cols
    , dbl_query
    , dbl_scale_ub = 100
    , dbl_scale_lb = 0
    , df_weights = NULL
    , lgc_overqualification_sub = T
){
  
  # Arguments validation
  stopifnot(
    "'df_data_cols' must be a data frame." = 
      is.data.frame(df_data_cols)
  )
  
  stopifnot(
    "'dbl_query' must be numeric." =
      all(
        is.numeric(dbl_query)
        , length(dbl_query) ==
          nrow(df_data_cols)
      )
  )
  
  stopifnot(
    "'df_weights' must be either NULL or a numeric data frame." = 
      any(
        all(map_lgl(df_weights, is.numeric))
        , is.null(df_weights)
      )
  )
  
  # Data wrangling
  dbl_scale_ub[[1]] -> dbl_scale_ub
  
  dbl_scale_lb[[1]] -> dbl_scale_lb
  
  t(df_data_cols) -> df_data_rows
  
  if(length(df_weights)){
    
    t(df_weights) -> df_weights
    
  }
  
  as_tibble(t(dbl_query)) -> df_query_rows
  
  rm(df_data_cols)
  rm(dbl_query)
  
  # Euclidean distance
  if(length(df_weights)){
    
    # Weighted Euclidean distance
    df_query_rows[rep(
      1, nrow(df_data_rows)
    ), ] -> df_query_rows
    
    if(lgc_overqualification_sub){
      
      df_weights * (pmax(
        df_query_rows - 
          df_data_rows
        , 0
      )) ^ 2 -> df_dist
      
    } else {
      
      df_weights * (
        df_data_rows - 
          df_query_rows
      ) ^ 2 -> df_dist
      
    }
    
    rm(df_query_rows)
    
    sqrt(rowSums(df_dist)) -> df_dist
    
    # Normalize by maximum distance
    df_dist / sqrt(rowSums(
      df_weights * 
        pmax(
          dbl_scale_ub - df_data_rows,
          df_data_rows - dbl_scale_lb
        ) ^ 2
    )) -> dbl_dist
    
  } else {
    
    # Unweighted Euclidean distance
    df_query_rows[rep(
      1, nrow(df_data_rows)
    ), ] -> df_query_rows
    
    (
      df_data_rows - 
        df_query_rows
    ) ^ 2 -> df_dist
    
    rm(df_query_rows)
    
    sqrt(rowSums(df_dist)) -> df_dist
    
    # Normalize by maximum distance
    df_dist / sqrt(rowSums(
      pmax(
        dbl_scale_ub - df_data_rows,
        df_data_rows - dbl_scale_lb
      ) ^ 2
    )) -> dbl_dist 
    
  }
  
  rm(df_dist)
  rm(df_weights)
  rm(df_data_rows)
  
  # Calculate similarity
  1 - dbl_dist -> dbl_similarity
  
  # Output
  return(dbl_similarity)
  
}

# - Similarity function (col vectors) ---------------------------------------------------
fun_match_similarity_cols <- function(
    df_data_cols
    , dbl_query
    , chr_method = c('bvls', 'logit', 'probit', 'pearson', 'spearman', 'kendau', 'euclidean')
    , df_weights = NULL
    , dbl_scale_ub = 100
    , dbl_scale_lb = 0
    , lgc_sort = F
){
  
  # Argument validation
  stopifnot(
    "'df_data_cols' must be a data frame." =
      is.data.frame(df_data_cols)
  )
  
  stopifnot(
    "'dbl_query' must be numeric." =
      is.numeric(dbl_query)
  )
  
  stopifnot(
    "'chr_method' must be one of the following methods: 'bvls', 'logit', 'probit', 'pearson', 'spearman', 'kendau' or 'euclidean'." =
      any(
        chr_method == 'bvls',
        chr_method == 'logit',
        chr_method == 'probit',
        chr_method == 'pearson',
        chr_method == 'spearman',
        chr_method == 'kendau',
        chr_method == 'euclidean'
      )
  )
  
  stopifnot(
    "'df_weights' must be either NULL or a numeric data frame." = 
      any(
        all(map_lgl(df_weights, is.numeric))
        , is.null(df_weights)
      )
  )
  
  stopifnot(
    "'dbl_scale_ub' must be numeric." =
      is.numeric(dbl_scale_ub)
  )
  
  stopifnot(
    "'dbl_scale_lb' must be numeric." =
      is.numeric(dbl_scale_lb)
  )
  
  # Apply matching method
  if(chr_method == 'bvls'){
    
    # Apply BVLS regression matching
    fun_match_bvls(
      df_data_cols =
        df_data_cols
      , dbl_query = 
        dbl_query
      , df_weights = 
        df_weights
    ) -> dbl_similarity
    
  } else if(chr_method == 'pearson') {
    
    # Apply Pearson correlation matching
    fun_match_pearson(
      df_data_cols = 
        df_data_cols
      , dbl_query = 
        dbl_query
      , df_weights = 
        df_weights
    ) -> dbl_similarity
    
  }  else if(chr_method == 'spearman') {
    
    # Apply Spearman correlation matching
    fun_match_pearson(
      df_data_cols = 
        df_data_cols
      , dbl_query = 
        dbl_query
      , df_weights = 
        df_weights
    ) -> dbl_similarity
    
  }  else if(chr_method == 'kendau') {
    
    # Apply Kendau's tau correlation matching
    fun_match_pearson(
      df_data_cols = 
        df_data_cols
      , dbl_query = 
        dbl_query
      , df_weights = 
        df_weights
    ) -> dbl_similarity
    
  } else if(chr_method == 'euclidean') {
    
    # Apply Euclidean matching
    fun_match_euclidean(
      df_data_cols = 
        df_data_cols
      , dbl_query = 
        dbl_query
      , dbl_scale_ub = 
        dbl_scale_ub
      , dbl_scale_lb = 
        dbl_scale_lb
      , df_weights = 
        df_weights
    ) -> dbl_similarity
    
  } else {
    
    # Apply logistic regression matching
    fun_match_logit(
      df_data_cols = 
        df_data_cols
      , dbl_query = 
        dbl_query
      , dbl_scale_ub = 
        dbl_scale_ub
      , dbl_scale_lb = 
        dbl_scale_lb
      , chr_method = 
        chr_method
      , df_weights = 
        df_weights
    ) -> dbl_similarity
    
  }
  
  # Output
  return(dbl_similarity)
  
}

# - Similarity function (row vectors) ---------------------------------------------------
fun_match_similarity <- function(
    df_data_rows
    , df_query_rows
    , chr_method = c('bvls', 'logit', 'probit', 'pearson', 'spearman', 'kendau', 'euclidean')
    , chr_weights = c('linear', 'quadratic', 'speciality-root', 'attribute-eqvl')
    , dbl_scale_ub = 100
    , dbl_scale_lb = 0
    , chr_id_col = NULL
    , lgc_sort = F
){
  
  # Argument validation
  stopifnot(
    "'df_data_rows' must be a data frame." =
      is.data.frame(df_data_rows)
  )
  
  stopifnot(
    "'df_query_rows' must be a data frame." =
      is.data.frame(df_query_rows)
  )
  
  stopifnot(
    "'chr_method' must be one of the following methods: 'bvls', 'logit', 'probit', 'pearson', 'spearman', 'kendau' or 'euclidean'." =
      any(
        chr_method == 'bvls',
        chr_method == 'logit',
        chr_method == 'probit',
        chr_method == 'pearson',
        chr_method == 'spearman',
        chr_method == 'kendau',
        chr_method == 'euclidean'
      )
  )
  
  stopifnot(
    "'chr_weights' must be either 'unweighted, 'linear', 'quadratic', 'speciality-root', or 'attribute-eqvl'." = 
      any(
        chr_weights == 'unweighted',
        chr_weights == 'linear',
        chr_weights == 'quadratic', 
        chr_weights == 'speciality-root', 
        chr_weights == 'attribute-eqvl'
      )
  )
  
  stopifnot(
    "'dbl_scale_ub' must be numeric." =
      is.numeric(dbl_scale_ub)
  )
  
  stopifnot(
    "'dbl_scale_lb' must be numeric." =
      is.numeric(dbl_scale_lb)
  )
  
  stopifnot(
    "'lgc_sort' must be either TRUE or FALSE." =
      all(
        is.logical(lgc_sort)
        , !is.na(lgc_sort)
      )
  )
  
  stopifnot(
    "'chr_id_col' must be either NULL or a character string." = 
      any(
        is.null(chr_id_col)
        , is.character(chr_id_col)
      )
  )
  
  # Data wrangling
  dbl_scale_ub[[1]] -> dbl_scale_ub
  
  dbl_scale_lb[[1]] -> dbl_scale_lb
  
  chr_method[[1]] -> chr_method
  
  chr_weights[[1]] -> chr_weights
  
  Filter(
    function(x){all(is.numeric(x))}
    , df_query_rows
  ) -> dbl_query
  
  # rm(df_query_rows)
  
  df_data_rows[names(
    dbl_query
  )] -> df_data_cols
  
  # Pivot data
  t(dbl_query) -> 
    dbl_query
  
  as_tibble(t(
    df_data_cols
  )) -> df_data_cols
  
  # ID column
  if(length(chr_id_col)){
    
    chr_id_col[[1]] -> chr_id_col
    
    df_query_rows %>%
      pull(!!sym(
        chr_id_col
      )) -> chr_id_query
    
    df_data_rows %>%
      pull(!!sym(
        chr_id_col
      )) -> chr_id_data
    
    # rm(chr_id_col)
    
  }
  
  # Weights
  if(any(
    chr_weights == 'speciality-root',
    chr_weights == 'attribute-eqvl'
  )){
    # Calculate skill set generality
    # for weighting methods that require generality
    map_dbl(
      df_data_cols
      , fun_gene_generality
      , dbl_scale_lb =
        dbl_scale_lb
      , dbl_scale_ub =
        dbl_scale_ub
    ) -> dbl_generality
    
  } else {
    
    NULL -> dbl_generality
    
  }
  
  # Apply weighting function
  fun_match_vweights(
    df_data_cols
    , chr_weights = 
      chr_weights
    , dbl_generality = 
      dbl_generality
  ) -> df_weights
  
  rm(chr_weights)
  rm(dbl_generality)
  
  # Apply similarity function
  dbl_query %>%
    as_tibble() ->
    df_query_cols
  
  rm(dbl_query)
  
  df_query_cols %>% 
    map(
      ~ fun_match_similarity_cols(
        df_data_cols = df_data_cols
        , dbl_query = as.matrix(.x)
        , chr_method = chr_method
        , df_weights = df_weights
        , dbl_scale_ub = dbl_scale_ub
        , dbl_scale_lb = dbl_scale_lb
      )
    ) -> list_similarity
  
  list_similarity %>% 
    map(
      ~ .x %>%
        set_names(
          chr_id_data
        )
    ) -> list_similarity
  
  rm(chr_method)
  rm(dbl_scale_ub)
  rm(dbl_scale_lb)
  
  # Similarity matrix
  rm(df_query_cols)
  rm(df_data_cols)
  
  list_similarity %>%
    bind_cols() %>%
    as.matrix() ->
    mtx_similarity
  
  if(length(chr_id_col)){
    
    colnames(mtx_similarity) <- chr_id_query
    
    rownames(mtx_similarity) <- chr_id_data
    
    names(list_similarity) <- chr_id_query
    
  }
  
  # Sort
  if(lgc_sort){
    
    list_similarity %>% 
      map(sort, decreasing = T) -> 
      list_similarity
    
  }
  
  # Output
  return(compact(list(
    'list_similarity' = list_similarity
    , 'mtx_similarity' = mtx_similarity
  )))
  
}

# # [TEST] ------------------------------------------------------------------
# # - Data ------------------------------------------------------------------
# df_occupations <- read_csv('/home/Cao/Storage/github/atlas-research/data/occupations/df_occupations_2022.csv')
# 
# df_occupations %>%
#   select(
#     occupation
#     , starts_with('skl_'),
#     , starts_with('abl_'),
#     , starts_with('knw_')
#   ) -> df_occupations
# 
# # - fun_match_similarity -----------------------------------
# fun_match_similarity(
#   df_data_rows =
#     df_occupations %>%
#     filter(
#       occupation %in% c(
#         'Mechanical Engineers',
#         'Physicists',
#         'Credit Analysts',
#         'Dishwashers'
#       )
#     ) %>%
#     slice(
#       2, 3, 1, 4
#     )
#   , df_query_rows =
#     df_occupations %>%
#     filter(
#       occupation %in% c(
#         'Mechanical Engineers',
#         'Physicists',
#         'Credit Analysts',
#         'Dishwashers'
#       )
#     ) %>%
#     slice(
#       2, 3, 1, 4
#     )
#   , chr_method = 'euclidean'
#   # , chr_method = 'bvls'
#   # , chr_method = 'pearson'
#   # , chr_method = 'logit'
#   # , chr_method = 'probit'
# 
#   # , chr_weights = 'unweighted'
#   # , chr_weights = 'linear'
#   # , chr_weights = 'quadratic'
#   # , chr_weights = 'speciality-root'
#   , chr_weights = 'attribute-eqvl'
#   , dbl_scale_ub = 100
#   , dbl_scale_lb = 0
#   , chr_id_col = 'occupation'
#   , lgc_sort = F
# ) -> list_match_euclidean
# 
# list_match_euclidean$
#   mtx_similarity %>%
#   round(2)
# 
