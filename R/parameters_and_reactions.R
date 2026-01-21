###################################
##                               ##
## Parameter and state creation  ##
##                               ##
###################################

#' load model parameters
#' @param file file path to yaml
#' @export
load_model_params <- function(file=NULL){
  if(is.null(file)){
    file <- system.file("extdata", "params.yaml", package = "brech")
  }
  conditions <- yaml::read_yaml(file)
  params <- conditions$params
  # preserves vector if using age groups
  N <- Reduce('+',conditions$initial_states)

  params$theta <- 1 / params$latent_period
  params$gamma <- 1 / params$infectious_period
  params$case_rate <- 1 / params$case_delay

  # if beta isn't defined then base on R0 and population size
  if(!("beta" %in% names(params))){
    params$beta <- params$r_0 * params$gamma / N
  }
  params$N <- N
  conditions$params <- params
  # init states
  # set initial cases to 0 by convention
  # multiply by another state to retain correct shape
  first_state <- names(conditions$initial_states)[1]
  conditions$initial_states$C <- 0*conditions$initial_states[[first_state]]
  conditions$initial_states$D <- 0*conditions$initial_states[[first_state]]
  return(conditions)
}

#' change vaccination coverage in initial conditions
#' @param sim_scenario output of `load_model_params`
#' @param vaccine_rate numeric proportion vaccinated. Can be a numeric
#'   vector if using age groups. If named then vaccinates according to each
#' age group with format "S_\{age_group\}"
#' @export
vaccinate_initial_conditions <- function(sim_scenario,vaccine_rate){
  age_groups <- names(vaccine_rate)
  initial_states <- sim_scenario$initial_states

  # if not named
  if(is.null(age_groups)){
    S <- (initial_states$S + initial_states$V)
    initial_states$V <- round(S*vaccine_rate)
    initial_states$S <- round(S*(1 - vaccine_rate))
    # add back in discrepancy due to rounding
    initial_states$S <- initial_states$S  + (S - initial_states$V - initial_states$S)

  }else{
    for(age_group in age_groups){
      age_vaccine_rate <- vaccine_rate[[age_group]]
      S <- initial_states[[paste0("S_",age_group)]]
      V <- initial_states[[paste0("V_",age_group)]]
      S <- N <- S + V
      # vaccinate
      V <- round(S*age_vaccine_rate)
      S <- N - V
      # update initial states
      initial_states[[paste0("S_",age_group)]] <- S
      initial_states[[paste0("V_",age_group)]] <- V
    }
  }
  sim_scenario$initial_states <- initial_states
  return(sim_scenario)
}

#' Convert a vector of each state into a labeled state using
#' `age_groups` labels
#' @param sim_scenario output of `load_model_params`
#' @param age_groups vector of age group labels
#' @examples
#' sim_scenario <- list(
#'   "initial_states" = list(
#'     "S" = c(1,2,3),
#'     "I" = c(0,0,1)
#'   )
#' )
#' age_groups <- c("0-5","5-19","19+")
#' create_age_initial_conditions(sim_scenario,age_groups)
#' @export
create_age_initial_conditions <- function(sim_scenario,age_groups){
  initial_states <- sim_scenario$initial_states
  updated_initial_states <- list()
  for(state in names(initial_states)){
    for(i in seq_along(age_groups)){
      age_group <- age_groups[i]
      updated_initial_states[[paste0(state,"_",age_group)]] <- initial_states[[state]][i]
    }
  }
  sim_scenario$initial_states <- updated_initial_states
  return(sim_scenario)
}

#' get reaction list from number of age groups m and
#' contact_matrix
#' @param contact_matrix matrix m x m
#' @export
build_age_reactions <- function(contact_matrix){
  m <- ncol(contact_matrix)
  age_groups <- colnames(contact_matrix)
  if(m != nrow(contact_matrix)){
    stop("contact_matrix needs the same number of rows and columns.")
  }

  reactions <- list()
  for (i in seq_along(age_groups)){
    age_group <- age_groups[i]
    ######################################
    #  detectable infection transitions  #
    ######################################
    transition <- function() {
      trans <- c()
      trans[paste0("S_",age_group)] <- -1
      trans[paste0("E_",age_group)] <- +1
      trans[paste0("D_",age_group)] <- +1
      trans[paste0("INC_",age_group)] <- +1
      trans
    }
    detectable_infection_rate <-
      function(i,age_group){
        force(i)
        force(age_group)
        function(x, p, t) {
          I <- x[paste0("I_", age_groups)]
          S <- x[paste0("S_", age_group)] # just need the ith S component
          lambda <- p$case_ascertainment * p$beta * S * sum(contact_matrix[,i] * (I / p$N))
          lambda  # return i-th force of infection term
        }
      }
    reactions[[paste0("detectable_infection_",age_group)]] <- list(
      "transition"=transition(),
      "rate" = detectable_infection_rate(i,age_group)
    )
    ########################################
    #  undetectable infection transitions  #
    ########################################
    transition <- function() {
      trans <- c()
      trans[paste0("S_",age_group)] <- -1
      trans[paste0("E_",age_group)] <- + 1
      trans[paste0("INC_",age_group)] <- +1
      trans
    }
    undetectable_infection_rate <-
      function(i,age_group){
        force(i)
        force(age_group)
        function(x, p, t) {
          I <- x[paste0("I_", age_groups)]
          S <- x[paste0("S_", age_group)] # just need the ith S component
          lambda <- (1 - p$case_ascertainment) * p$beta * S * sum(contact_matrix[,i] * (I / p$N))
          lambda  # return i-th force of infection term
        }
      }
    reactions[[paste0("undetectable_infection_",age_group)]] <- list(
      "transition"=transition(),
      "rate" = undetectable_infection_rate(i,age_group)
    )
    ############################
    #  incubation transitions  #
    ############################
    transition <- function() {
      trans <- c()
      trans[paste0("E_",age_group)] <- -1
      trans[paste0("I_",age_group)] <- + 1
      trans
    }
    incubation_rate <-
      function(age_group){
        force(age_group)
        function(x, p, t) {
          E <- x[paste0("E_", age_group)] # ith component
          p$theta * E
        }
      }
    reactions[[paste0("incubation_",age_group)]] <- list(
      "transition"=transition(),
      "rate" = incubation_rate(age_group)
    )
    ##############################
    #    recovery transitions    #
    ##############################
    transition <- function() {
      trans <- c()
      trans[paste0("I_",age_group)] <- -1
      trans[paste0("R_",age_group)] <- + 1
      trans
    }
    recovery_rate <-
      function(age_group){
        force(age_group)
        function(x, p, t) {
          I <- x[paste0("I_", age_group)] # ith component
          p$gamma * I
        }
      }
    reactions[[paste0("recovery_",age_group)]] <- list(
      "transition"=transition(),
      "rate" = recovery_rate(age_group)
    )
    ##############################
    #        case detection      #
    ##############################
    transition <- function() {
      trans <- c()
      trans[paste0("C_",age_group)] <- + 1
      trans[paste0("D_",age_group)] <- - 1
      trans
    }
    case_detection_rate <-
      function(age_group){
        force(age_group)
        function(x, p, t) {
          D <- x[paste0("D_", age_group)] # ith component
          p$case_rate * D
        }
      }
    reactions[[paste0("case_detection_",age_group)]] <- list(
      "transition"=transition(),
      "rate" = case_detection_rate(age_group)
    )
    ##############################
    #              PEP           #
    ##############################
    transition <- function() {
      trans <- c()
      trans[paste0("E_",age_group)] <- - 1
      trans[paste0("R_",age_group)] <- + 1
      trans
    }
    pep_administered_rate <-
      function(age_group){
        force(age_group)
        function(x, p, t) {
          E <- x[paste0("E_", age_group)] # ith component
          p$pep_rate * E
        }
      }
    reactions[[paste0("pep_administered_",age_group)]] <- list(
      "transition"=transition(),
      "rate" = pep_administered_rate(age_group)
    )
    ###############################
    #        hospitalization      #
    ###############################
    transition <- function() {
      trans <- c()
      trans[paste0("H_",age_group)] <- + 1
      trans[paste0("I_",age_group)] <- - 1
      trans
    }
    hospitalization_rate <-
      function(i, age_group){
        force(i)
        force(age_group)
        function(x, p, t) {
          I <- x[paste0("I_", age_group)] # ith component
          prop_hosp <- p$hospitalization_proportion[i]
          p$gamma * (prop_hosp/(1-prop_hosp)) * I
        }
      }
    reactions[[paste0("hospitalization_",age_group)]] <- list(
      "transition"=transition(),
      "rate" = hospitalization_rate(i, age_group)
    )

  }

  return(reactions)

}

#' Using the next generation matrix approach, calculates beta based on a fixed
#' R0
#' @param sm Object of class \code{stochastic_model}
#' @param contact_matrix An m x m contact matrix
#' @param population_sizes An m vector of population size
#' @param R0 numeric reproduction number
#' @return numeric
#' @export
calculate_age_structured_beta <- function(sm, contact_matrix, population_sizes, R0){
  # beta calculation based on R0
  recovery_rate <- sm$params$gamma
  adjusted_matrix <- contact_matrix * outer(population_sizes, 1/population_sizes)
  spectral_radius <- max(Mod(eigen(adjusted_matrix)$values))
  beta <- recovery_rate*R0/spectral_radius

  return(beta)
}
