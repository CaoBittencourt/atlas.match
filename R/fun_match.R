# [SETUP] -----------------------------------------------------------------
# - Packages ----------------------------------------------------------------
# CRAN packages
chr_pkg <- c(
  'bvls', 'fastglm', 'weights' #Regression models
  , 'devtools' #GitHub packages (temp)
  , 'readr' #Read data (temp)
  , 'dplyr', 'tidyr', 'purrr' #Data wrangling
)

# Git packages
chr_git <- c(
  'CaoBittencourt' = 'atlas.kind' #Human capital indispensability coefficient
)

# Activate / install CRAN packages
lapply(
  chr_pkg
  , function(pkg){

    if(!require(pkg, character.only = T)){

      install.packages(pkg)

    }

    require(pkg, character.only = T)

  }
)

# Activate / install Git packages
Map(
  function(git, profile){

    if(!require(git, character.only = T)){

      install_github(
        paste0(profile, '/', git)
        , upgrade = F
        , force = T
      )

    }

    require(git, character.only = T)

  }
  , git = chr_git
  , profile = names(chr_git)
)

# [MATCHING FUNCTIONS] -------------------------------------------------------------
# - Regression weights --------------------------------------------
fun_match_weights <- function(dbl_var){
  
  # Argument validation
  stopifnot(
    "'dbl_var' must be a numeric vector or matrix." =
      is.numeric(dbl_var)
  )
  
  # Apply human capital indispensability function
  fun_kind_indispensability(
    dbl_profile = dbl_var,
    dbl_scale_lb = 0
  ) -> dbl_weights
  
  # Apply equivalence function to regression weights
  # fun_eqvl_equivalence(
  #   dbl_var = dbl_var
  #   , dbl_scaling =
  #     dbl_scaling
  # ) -> dbl_weights
  
  # Output
  return(dbl_weights)
  
}

# - Vectorized regression weights -----------------------------------------
fun_match_vweights <- function(df_data_cols){
  
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
  map_df(
    .x = df_data_cols
    , fun_match_weights
  ) -> df_data_cols
  
  # map_df(
  #   .x = df_data_cols
  #   , ~ fun_match_weights(
  #     dbl_var = .x
  #     , dbl_scaling = 
  #       dbl_scaling
  #   )
  # ) -> df_data_cols
  
  # Output
  return(df_data_cols)
  
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

# - KNN/Euclidean matching ----------------------------------------------------
fun_match_knn <- function(
    df_data_cols
    , dbl_query
    , dbl_scale_ub = 100
    , dbl_scale_lb = 0
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
  
  # Data wrangling
  dbl_scale_ub[[1]] -> dbl_scale_ub
  
  dbl_scale_lb[[1]] -> dbl_scale_lb
  
  t(df_data_cols) -> df_data_rows
  
  t(df_weights) -> df_weights
  
  as_tibble(t(dbl_query)) -> df_query_rows
  
  rm(df_data_cols)
  rm(dbl_query)
  
  # Euclidean distance
  df_query_rows[rep(
    1, nrow(df_data_rows)
  ), ] -> df_query_rows
  
  df_data_rows - df_query_rows -> df_dist
  
  rm(df_data_rows)
  rm(df_query_rows)
  
  abs(df_dist) -> df_dist
  
  # Normalize by scale bounds
  df_dist / (
    dbl_scale_ub -
      dbl_scale_lb
  ) -
    dbl_scale_lb / (
      dbl_scale_ub -
        dbl_scale_lb
    ) -> df_dist
  
  # Weigh distances by attribute indispensability
  if(!is.null(df_weights)){
    
    df_dist * df_weights -> df_dist
    
  }
  
  # Calculate similarity
  1 - rowMeans(df_dist) -> 
    dbl_similarity
  
  # Output
  return(dbl_similarity)
  
}

# - Similarity function (col vectors) ---------------------------------------------------
fun_match_similarity_cols <- function(
    df_data_cols
    , dbl_query
    , chr_method = c('bvls', 'logit', 'probit', 'pearson', 'knn')
    , dbl_scale_ub = 100
    , dbl_scale_lb = 0
    , dbl_scaling = 0.25
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
    "'chr_method' must be one of the following methods: 'bvls', 'logit', 'probit', 'pearson', or 'knn'." =
      any(
        chr_method == 'bvls',
        chr_method == 'logit',
        chr_method == 'probit',
        chr_method == 'pearson',
        chr_method == 'knn'
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
  
  # Data wrangling
  dbl_scale_ub[[1]] -> dbl_scale_ub
  
  dbl_scale_lb[[1]] -> dbl_scale_lb
  
  chr_method[[1]] -> chr_method
  
  # Weights
  fun_match_vweights(
    df_data_cols
  ) -> df_weights
  
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
    
  } else if(chr_method == 'knn') {
    
    # Apply KNN/Euclidean matching
    fun_match_knn(
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
    , chr_method = c('bvls', 'logit', 'probit', 'pearson', 'knn')
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
    "'chr_method' must be one of the following methods: 'bvls', 'logit', 'probit', 'pearson', or 'knn'." =
      any(
        chr_method == 'bvls',
        chr_method == 'logit',
        chr_method == 'probit',
        chr_method == 'pearson',
        chr_method == 'knn'
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
  Filter(
    function(x){all(is.numeric(x))}
    , df_query_rows
  ) -> dbl_query
  
  rm(df_query_rows)
  
  df_data_rows[names(
    dbl_query
  )] -> df_data_cols
  
  # Pivot data
  t(dbl_query) -> 
    dbl_query
  
  as_tibble(t(
    df_data_cols
  )) -> df_data_cols
  
  # Apply similarity function
  if(ncol(dbl_query) == 1){
    
    fun_match_similarity_cols(
      df_data_cols = df_data_cols
      , dbl_query = dbl_query
      , chr_method = chr_method
      , dbl_scale_ub = dbl_scale_ub
      , dbl_scale_lb = dbl_scale_lb
    ) -> df_data_rows$similarity
    
    rm(dbl_query)
    rm(chr_method)
    rm(dbl_scale_ub)
    rm(dbl_scale_lb)
    rm(df_data_cols)
    
    # Sort data frame
    if(lgc_sort){
      
      df_data_rows %>%
        arrange(desc(
          similarity
        )) -> df_data_rows
      
    }
    
    list_similarity <- NULL
    mtx_similarity <- NULL
    
  } else {
    
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
          , dbl_scale_ub = dbl_scale_ub
          , dbl_scale_lb = dbl_scale_lb
        )
      ) -> list_similarity
    
    rm(chr_method)
    rm(dbl_scale_ub)
    rm(dbl_scale_lb)
    
    # Similarity matrix
    if(
      all(
        df_query_cols ==
        df_data_cols
      )
    ){
      
      rm(df_query_cols)
      rm(df_data_cols)
      
      list_similarity %>%
        bind_cols() %>%
        as.matrix() ->
        mtx_similarity
      
      if(length(chr_id_col)){
        
        chr_id_col[[1]] -> chr_id_col
        
        df_data_rows %>%
          pull(!!sym(chr_id_col)) ->
          colnames(
            mtx_similarity
          )
        
        colnames(
          mtx_similarity
        ) -> rownames(
          mtx_similarity
        )
        
        colnames(
          mtx_similarity
        ) -> names(
          list_similarity
        )
        
      }
      
    }
    
    df_data_rows <- NULL
    
  }
  
  # Output
  return(compact(list(
    'df_similarity' = df_data_rows
    , 'list_similarity' = list_similarity
    , 'mtx_similarity' = mtx_similarity
  )))
  
}

# # [TEST] ------------------------------------------------------------------
# # - Data ------------------------------------------------------------------
# # Occupations data frame
# df_occupations <- read_csv('C:/Users/Cao/Documents/Github/Atlas-Research-dev/Data/df_occupations_2023_efa.csv')
# df_occupations_old <- read_csv('C:/Users/Cao/Documents/Github/Atlas-Research-dev/Data/df_atlas_complete_equamax_15_factors.csv')
# 
# 
# # My own professional profile
# df_input <- read_csv('https://docs.google.com/spreadsheets/d/e/2PACX-1vSj7u2N59j8MTa7MZqk2Y-VDVWIWEDzAR_0gkb_jB_pBX4sm8yMS1N26ClmY6iWXA/pub?gid=145103706&single=true&output=csv')
# df_input_old <- read_csv('https://docs.google.com/spreadsheets/d/e/2PACX-1vT7Gmo-eVC5C1lksSlefJd8N1lciaAn057hExjeoE5XtbpDvJfz7Bxc_f_hwKSx_A/pub?gid=1515296378&single=true&output=csv')
# 
# # Factor model
# efa_model <- read_rds('C:/Users/Cao/Documents/Github/atlas-research-dev/data/efa/efa_equamax_14factors.rds')
# efa_model_old <- read_rds('C:/Users/Cao/Documents/Github/atlas-research-dev/data/efa/old/2022/efa_model_equamax_15_factors.rds')
# 
# # - Regression weights 1 ----------------------------------------------------
# fun_match_weights(
#   dbl_var = runif(50, 0, 100)
# )
# 
# # - Regression weights 2 --------------------------------------------------
# fun_match_vweights(
#   df_data =
#     df_occupations %>%
#     select(starts_with('item_')) %>%
#     t() %>%
#     as_tibble()
# )
# 
# # - BVLS similarity test ------------------------------------------------------------------
# fun_match_similarity(
#   df_data_rows =
#     df_occupations %>%
#     select(
#       occupation
#       , starts_with('item_')
#     )
#   , chr_method = 'bvls'
#   , df_query_rows =
#     # df_input
#     df_occupations %>%
#     slice_sample(n = 1) %>%
#     select(
#       occupation
#       , starts_with('item_')
#     )
#   , dbl_scale_ub = 100
#   , dbl_scale_lb = 0
#   , lgc_sort = T
# )[[1]] %>%
#   select(
#     occupation,
#     similarity
#   ) %>% 
#   print(
#     n = 100
#   )
# 
# # - Pearson similarity test ------------------------------------------------------------------
# fun_match_similarity(
#   df_data_rows =
#     df_occupations %>%
#     select(
#       occupation
#       , starts_with('item_')
#     )
#   , chr_method = 'pearson'
#   , df_query_rows =
#     # df_input
#     df_occupations %>%
#     slice_sample(n = 1) %>%
#     select(
#       occupation
#       , starts_with('item_')
#     )
#   , dbl_scale_ub = 100
#   , dbl_scale_lb = 0
#   , lgc_sort = T
# )[[1]] %>%
#   select(
#     occupation,
#     similarity
#   ) %>% 
#   print(
#     n = 100
#   )
# 
# # - KNN/Euclidean similarity test ------------------------------------------------------------------
# fun_match_similarity(
#   df_data_rows =
#     df_occupations %>%
#     select(
#       occupation
#       , starts_with('item_')
#     )
#   , chr_method = 'knn'
#   , df_query_rows =
#     df_input
#   # df_occupations %>%
#   # slice_sample(n = 1) %>%
#   # select(
#   #   occupation
#   #   , starts_with('item_')
#   # )
#   , dbl_scale_ub = 100
#   , dbl_scale_lb = 0
#   , lgc_sort = T
# )[[1]] %>%
#   select(
#     occupation,
#     similarity
#   ) %>%
#   print(
#     n = 100
#   )
# 
# # - Logit similarity test ------------------------------------------------------------------
# # New data base
# fun_match_similarity(
#   df_data_rows =
#     df_occupations %>%
#     select(
#       occupation
#       , starts_with('item_')
#     )
#   , chr_method = 'logit'
#   , df_query_rows =
#     df_input
#   # df_occupations %>%
#   # slice_sample(n = 1) %>%
#   # select(
#   #   occupation
#   #   , starts_with('item_')
#   # )
#   , dbl_scale_ub = 100
#   , dbl_scale_lb = 0
#   , lgc_sort = T
# )[[1]] %>%
#   select(
#     occupation,
#     similarity
#   ) %>% 
#   print(
#     n = 100
#   ) -> df_similarity
# 
# # Old data base
# fun_match_similarity(
#   df_data_rows =
#     df_occupations_old %>%
#     select(
#       occupation
#       , ends_with('.l')
#     )
#   , chr_method = 'logit'
#   , df_query_rows =
#     df_input_old
#   , dbl_scale_ub = 100
#   , dbl_scale_lb = 0
#   , lgc_sort = T
# )[[1]] %>%
#   select(
#     occupation,
#     similarity
#   ) %>% 
#   print(
#     n = 100
#   ) -> df_similarity_old
# 
# df_similarity %>% print(n = 20)
# df_similarity_old %>% print(n = 20)
# 
# # - Similarity matrix test ------------------------------------------------
# fun_match_similarity(
#   df_data_rows =
#     df_occupations %>%
#     slice(1:10) %>%
#     select(
#       occupation
#       , starts_with('item_')
#     )
#   , df_query_rows =
#     df_occupations %>%
#     slice_sample(n = 1) %>%
#     select(
#       occupation
#       , starts_with('item_')
#     )
#   , chr_method = 'bvls'
#   , dbl_scale_ub = 100
#   , dbl_scale_lb = 0
#   , chr_id_col =
#     'occupation'
# )
