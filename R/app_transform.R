# New variables on the Data page. A transformation is a structured request
# (type, variables, constants), validated against the data and compiled into
# an assignment such as data_1[["log_aadt"]] <- log(data_1[["AADT"]]). User
# text is never parsed as code: variable names become strings, numbers are
# converted with as.numeric(), and category values stay strings.

.app_transform_types <- c(
  "Natural log: log(x)" = "log",
  "Log of x + c: log(x + c)" = "log_c",
  "Square root" = "sqrt",
  "Power: x^p" = "power",
  "Standardize: (x - mean) / sd" = "standardize",
  "Center: x - mean" = "center",
  "Indicator: x compared with a value" = "indicator",
  "Indicator: x equals a category" = "indicator_level",
  "Product: x * y" = "product",
  "Ratio: x / y" = "ratio",
  "Categories: cut x at break points" = "cut",
  "Categorical: factor(x)" = "factor")

.app_operators <- c(">", ">=", "<", "<=", "==", "!=")

# Problems with a transformation request (empty when it can be applied)
check_transform <- function(t, data) {
  p <- character(0)
  add <- function(...) p <<- c(p, paste0(...))
  cols <- names(data)
  if (is.null(t$name) || !nzchar(trimws(t$name))) {
    add("Give the new variable a name.")
  } else if (t$name %in% cols) {
    add("A column named \"", t$name, "\" already exists; choose a new name ",
        "(existing columns are never overwritten, so earlier fits stay ",
        "reproducible).")
  }
  if (!t$type %in% .app_transform_types) {
    add("Unknown transformation.")
    return(p)
  }
  if (is.null(t$x) || !t$x %in% cols) {
    add("Choose the variable to transform.")
    return(p)
  }
  x <- data[[t$x]]
  numeric_types <- c("log", "log_c", "sqrt", "power", "standardize",
                     "center", "indicator", "product", "ratio", "cut")
  if (t$type %in% numeric_types && !is.numeric(x)) {
    add("\"", t$x, "\" is not numeric.")
  }
  xx <- if (is.numeric(x)) x[!is.na(x)] else numeric()
  num <- function(v) suppressWarnings(as.numeric(v))
  if (t$type == "log" && any(xx <= 0)) {
    add("log() needs positive values; \"", t$x, "\" has ", sum(xx <= 0),
        " zero or negative value(s). Use log(x + c).")
  }
  if (t$type == "log_c") {
    c0 <- num(t$constant)
    if (length(c0) != 1 || is.na(c0)) {
      add("Give the constant c as a number.")
    } else if (any(xx + c0 <= 0)) {
      add("x + c must be positive for every row.")
    }
  }
  if (t$type == "sqrt" && any(xx < 0)) {
    add("sqrt() needs nonnegative values.")
  }
  if (t$type == "power" && (length(num(t$constant)) != 1 ||
                            is.na(num(t$constant)))) {
    add("Give the power p as a number.")
  }
  if (t$type == "indicator") {
    if (!t$op %in% .app_operators) add("Unknown comparison.")
    if (length(num(t$value)) != 1 || is.na(num(t$value))) {
      add("Give the comparison value as a number.")
    }
  }
  if (t$type == "indicator_level") {
    if (is.null(t$value) || !nzchar(t$value)) add("Choose a category.")
  }
  if (t$type %in% c("product", "ratio")) {
    if (is.null(t$y) || !t$y %in% cols) {
      add("Choose the second variable.")
    } else if (!is.numeric(data[[t$y]])) {
      add("\"", t$y, "\" is not numeric.")
    } else if (t$type == "ratio" && any(data[[t$y]] == 0, na.rm = TRUE)) {
      add("\"", t$y, "\" has zeros; the ratio would be infinite.")
    }
  }
  if (t$type == "cut") {
    b <- num(strsplit(t$value %||% "", "[,; ]+")[[1]])
    if (!length(b) || anyNA(b)) {
      add("Give the break points as numbers separated by commas.")
    }
  }
  p
}

# The assignment expression of a transformation of data object `obj`
compile_transform <- function(t, obj) {
  d <- .sym(obj)
  col <- function(v) bquote(.(d)[[.(v)]])
  x <- col(t$x)
  num <- function(v) as.numeric(v)
  rhs <- switch(t$type,
    log = bquote(log(.(x))),
    log_c = bquote(log(.(x) + .(num(t$constant)))),
    sqrt = bquote(sqrt(.(x))),
    power = bquote(.(x)^.(num(t$constant))),
    standardize = bquote((.(x) - mean(.(x), na.rm = TRUE)) /
                           stats::sd(.(x), na.rm = TRUE)),
    center = bquote(.(x) - mean(.(x), na.rm = TRUE)),
    indicator = as.call(list(.sym("as.integer"),
                             as.call(list(.sym(t$op), x, num(t$value))))),
    indicator_level = bquote(as.integer(as.character(.(x)) == .(t$value))),
    product = bquote(.(x) * .(col(t$y))),
    ratio = bquote(.(x) / .(col(t$y))),
    cut = {
      b <- sort(unique(num(strsplit(t$value, "[,; ]+")[[1]])))
      bquote(cut(.(x), breaks = .(as.call(c(.sym("c"), -Inf, as.list(b),
                                                 Inf))),
                 right = FALSE))
    },
    factor = bquote(factor(.(x))))
  bquote(.(d)[[.(t$name)]] <- .(rhs))
}

# A short description of a transformation (for the log title)
describe_transform <- function(t) {
  e <- compile_transform(t, "d")[[3]]
  paste0(t$name, " = ", gsub("d\\[\\[\"([^\"]*)\"\\]\\]", "\\1",
                             .deparse(e)))
}

# Apply a transformation to the current data object of the session
app_add_variable <- function(state, t) {
  obj <- state$data_meta$object
  if (is.null(obj)) stop("Load data first.", call. = FALSE)
  env <- state$env
  p <- check_transform(t, get(obj, envir = env))
  if (length(p)) stop(paste(p, collapse = " "), call. = FALSE)
  e <- compile_transform(t, obj)
  run_compiled(list(e), env)
  state$data <- get(obj, envir = env)
  meta <- state$data_meta
  meta$version <- (meta$version %||% 0L) + 1L
  meta$new_variable <- t$name
  state$data_meta <- meta
  app_log(state, paste0("New variable in ", obj, ": ",
                        describe_transform(t)), list(e), "transform",
          t$name)
  invisible(t$name)
}
