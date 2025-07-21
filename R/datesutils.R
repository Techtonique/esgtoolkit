
make_time_index <- function(start_date, n, frequency) {
  freq_str <- normalize_frequency(frequency)
  
  by <- switch(freq_str,
               "daily" = "day",
               "weekly" = "week",
               "monthly" = "month",
               "quarterly" = "quarter",
               "yearly" = "year",
               stop("Unsupported frequency"))
  
  seq.Date(from = as.Date(start_date), by = by, length.out = n)
}

normalize_frequency <- function(freq) {
  if (is.character(freq)) return(tolower(freq))
  freq_map <- c(`1` = "yearly", `4` = "quarterly", `12` = "monthly",
                `52` = "weekly", `260` = "daily", `365` = "daily")
  freq_map[[as.character(freq)]] %||% "daily"
}

# Converts start_/end_ input to Date or appropriate zoo index type
parse_time_input <- function(x, frequency = NULL) {
  # If x is character -> try as Date
  if (is.character(x)) {
    as.Date(x)
  } else if (is.numeric(x) && length(x) == 2) {
    # x is ts cycle, e.g. c(2025, 4)
    year <- x[1]
    period <- x[2]
    if (is.null(frequency)) stop("frequency must be provided when start_/end_ is a cycle")
    
    freq_str <- normalize_frequency(frequency)
    
    switch(freq_str,
      "yearly" = as.Date(paste0(year, "-01-01")),
      "quarterly" = as.Date(as.yearqtr(year + (period - 1) / 4)),
      "monthly" = as.Date(as.yearmon(year + (period - 1) / 12)),
      "weekly" = as.Date(as.Date(paste0(year, "-01-01")) + (period - 1) * 7),
      "daily" = as.Date(as.Date(paste0(year, "-01-01")) + (period - 1)),
      stop("Unsupported frequency for cycle parsing")
    )
  } else if (inherits(x, "Date")) {
    x
  } else {
    stop("start_ or end_ must be either 'YYYY-MM-DD' string or numeric cycle of length 2")
  }
}


#' Convert Simulation Output to a zoo Time Series
#'
#' Converts the output of `simdiff()` or `simshocks()` (typically a `ts` object or matrix)
#' into a `zoo` time series with proper calendar-based indexing.
#'
#' The function supports flexible time specification using either a `Date` string (e.g. `"2025-01-01"`) 
#' or a cycle notation (e.g. `c(2025, 4)` for 4th quarter of 2025), along with a frequency that can 
#' be numeric (e.g. `12`, `4`, `260`) or string-based (`"monthly"`, `"quarterly"`, etc.).
#'
#' @param sim_data A time series (`ts`), matrix, or vector representing the simulation output.
#' @param start_ Optional. The starting point of the time index, as a `Date` string (`"YYYY-MM-DD"`) 
#'        or a numeric cycle (`c(year, period)`). If `sim_data` is a `ts` object and `start_` is `NULL`, 
#'        the function will attempt to infer the start date automatically.
#' @param frequency Optional. The frequency of the time series. Can be numeric (e.g., `12`, `4`, `260`) 
#'        or character (`"monthly"`, `"quarterly"`, `"daily"`, etc.). Inferred from `sim_data` if missing.
#' @param end_ Optional. A date string or cycle indicating the end of the time series. If provided, it overrides
#'        the number of observations used to construct the time index.
#'
#' @return A `zoo` object with the same data as `sim_data` and a time index based on calendar dates.
#'
#' @examples
#' @importFrom zoo zoo as.yearmon as.yearqtr
#' @export
sim_to_zoo <- function(sim_data, start_ = NULL, frequency = NULL, end_ = NULL) {
  # If sim_data is ts and start_/frequency not given, infer:
  if (inherits(sim_data, "ts")) {
    if (is.null(start_)) start_ <- start(sim_data)
    if (is.null(frequency)) frequency <- frequency(sim_data)
  } else {
    if (is.null(start_) || is.null(frequency)) {
      stop("start_ and frequency must be provided for non-ts data")
    }
  }
  
  # Parse start_ to Date (if cycle, need frequency)
  start_date <- parse_time_input(start_, frequency)
  
  sim_mat <- if (is.ts(sim_data)) {
    as.matrix(sim_data)
  } else if (is.vector(sim_data)) {
    matrix(sim_data, ncol = 1)
  } else {
    sim_data
  }
  
  n <- nrow(sim_mat)
  
  # Optional: parse end_ (if given, override n)
  if (!is.null(end_)) {
    end_date <- parse_time_input(end_, frequency)
    n <- as.integer(as.numeric(difftime(end_date, start_date, units = "days")) / switch(normalize_frequency(frequency),
                                                                                      "daily" = 1,
                                                                                      "weekly" = 7,
                                                                                      "monthly" = 30,
                                                                                      "quarterly" = 91,
                                                                                      "yearly" = 365,
                                                                                      1) ) + 1
    if (n < 1) stop("end_ must be after start_")
  }
  
  dates <- make_time_index(start_date, n, frequency)
  
  zoo::zoo(sim_mat, order.by = dates)
}

