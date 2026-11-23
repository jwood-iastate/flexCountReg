# Data page: choose a bundled dataset or upload a CSV file, then inspect it.

mod_data_ui <- function(id) {
  ns <- shiny::NS(id)
  bd <- fcr_bundled_datasets()
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 320,
      shiny::radioButtons(ns("source"), "Data source",
                          c("Example dataset" = "bundled",
                            "Upload a CSV file" = "upload")),
      shiny::conditionalPanel(
        "input.source == 'bundled'", ns = ns,
        shiny::selectInput(ns("bundled"), "Example dataset",
                           stats::setNames(bd$name, bd$title))),
      shiny::conditionalPanel(
        "input.source == 'upload'", ns = ns,
        shiny::fileInput(ns("file"), "CSV file", accept = c(".csv",
                                                             "text/csv")),
        help_text("Only CSV files are read; scripts and saved R objects ",
                  "are not accepted. The file stays on this computer.")),
      shiny::actionButton(ns("use"), "Use dataset",
                          class = "btn-primary"),
      bslib::accordion(
        open = FALSE, class = "mt-3",
        bslib::accordion_panel(
          "Create a variable",
          shiny::selectInput(ns("tf_type"), "Transformation",
                             .app_transform_types),
          shiny::selectInput(ns("tf_x"), "Variable (x)", character()),
          shiny::conditionalPanel(
            "input.tf_type == 'log_c' || input.tf_type == 'power'", ns = ns,
            shiny::textInput(ns("tf_constant"), "Constant (c or p)", "1")),
          shiny::conditionalPanel(
            "input.tf_type == 'indicator'", ns = ns,
            shiny::selectInput(ns("tf_op"), "Comparison", .app_operators),
            shiny::textInput(ns("tf_value"), "Value", "0")),
          shiny::conditionalPanel(
            "input.tf_type == 'indicator_level'", ns = ns,
            shiny::selectInput(ns("tf_level"), "Category", character())),
          shiny::conditionalPanel(
            "input.tf_type == 'product' || input.tf_type == 'ratio'", ns = ns,
            shiny::selectInput(ns("tf_y"), "Second variable (y)",
                               character())),
          shiny::conditionalPanel(
            "input.tf_type == 'cut'", ns = ns,
            shiny::textInput(ns("tf_breaks"), "Break points",
                             placeholder = "e.g. 10, 20, 50")),
          shiny::textInput(ns("tf_name"), "New variable name"),
          help_text("The new column is added to the current dataset and ",
                    "its code to the session script; existing columns are ",
                    "never overwritten."),
          shiny::actionButton(ns("tf_add"), "Add variable",
                              class = "btn-outline-primary")))),
    shiny::uiOutput(ns("about")),
    bslib::card(bslib::card_header("Columns"), fill = FALSE,
                DT::DTOutput(ns("columns"))),
    bslib::card(bslib::card_header("First rows"), fill = FALSE,
                DT::DTOutput(ns("preview"))))
}

# Summary of the columns of a data frame (type, missing values, range)
column_summary <- function(d) {
  data.frame(
    column = names(d),
    type = vapply(d, function(v) class(v)[1], character(1)),
    missing = vapply(d, function(v) sum(is.na(v)), integer(1)),
    distinct = vapply(d, function(v) length(unique(v[!is.na(v)])),
                      integer(1)),
    min = vapply(d, function(v) if (is.numeric(v) && any(!is.na(v)))
      format(signif(min(v, na.rm = TRUE), 4)) else "", character(1)),
    max = vapply(d, function(v) if (is.numeric(v) && any(!is.na(v)))
      format(signif(max(v, na.rm = TRUE), 4)) else "", character(1)),
    count_like = vapply(d, function(v) is.numeric(v) &&
                          all(v[!is.na(v)] >= 0) &&
                          all(abs(v[!is.na(v)] - round(v[!is.na(v)])) <
                                1e-8), logical(1)),
    row.names = NULL, stringsAsFactors = FALSE)
}

mod_data_server <- function(id, state) {
  shiny::moduleServer(id, function(input, output, session) {
    # Variable choices of the "Create a variable" panel
    shiny::observeEvent(state$data_meta, {
      cols <- names(state$data)
      shiny::updateSelectInput(session, "tf_x", choices = cols,
                               selected = intersect(input$tf_x, cols))
      shiny::updateSelectInput(session, "tf_y", choices = cols,
                               selected = intersect(input$tf_y, cols))
    })
    shiny::observe({
      d <- state$data
      x <- input$tf_x
      lev <- if (!is.null(d) && !is.null(x) && x %in% names(d)) {
        v <- d[[x]]
        u <- if (is.factor(v)) levels(v) else sort(unique(as.character(v)))
        utils::head(u, 200)
      } else character()
      shiny::updateSelectInput(session, "tf_level", choices = lev)
    })
    # Suggest a name for the new variable
    shiny::observeEvent(list(input$tf_type, input$tf_x, input$tf_y), {
      x <- input$tf_x %||% ""
      if (!nzchar(x)) return()
      stem <- switch(input$tf_type, log = "log_", log_c = "log_", sqrt = "sqrt_",
                     power = "pow_", standardize = "z_", center = "c_",
                     indicator = "ind_", indicator_level = "is_",
                     product = "prod_", ratio = "ratio_", cut = "cat_",
                     factor = "f_")
      nm <- paste0(stem, make.names(x))
      if (input$tf_type %in% c("product", "ratio") && nzchar(input$tf_y %||% "")) {
        nm <- paste0(nm, "_", make.names(input$tf_y))
      }
      shiny::updateTextInput(session, "tf_name", value = nm)
    })
    shiny::observeEvent(input$tf_add, app_try(what = "Creating the variable", {
      t <- list(type = input$tf_type, x = input$tf_x, y = input$tf_y,
                name = trimws(input$tf_name), op = input$tf_op,
                value = switch(input$tf_type, indicator = input$tf_value,
                               indicator_level = input$tf_level,
                               cut = input$tf_breaks, NULL),
                constant = input$tf_constant)
      app_add_variable(state, t)
      shiny::showNotification(paste0("Added ", t$name, " to ",
                                     state$data_meta$object, "."),
                              type = "message")
    }))

    shiny::observeEvent(input$use, app_try(what = "Loading the data", {
      if (input$source == "bundled") {
        bd <- fcr_bundled_datasets()
        b <- bd[bd$name == input$bundled, ]
        app_load_data(state, "bundled", input$bundled, title = b$title,
                      response = b$response, group = b$group)
      } else {
        f <- input$file
        if (is.null(f)) stop("Choose a CSV file first.", call. = FALSE)
        # Keep a copy of the upload for the session (fileInput's temporary
        # file can be replaced by the next upload)
        keep <- file.path(state$session_dir,
                          paste0("upload-", length(state$uploads) + 1,
                                 ".csv"))
        if (!isTRUE(file.copy(f$datapath, keep, overwrite = TRUE))) {
          stop("Could not keep a copy of the uploaded file.", call. = FALSE)
        }
        app_load_data(state, "upload", f$name, upload_path = keep,
                      title = paste("Uploaded file", f$name))
      }
      m <- state$data_meta
      shiny::showNotification(paste0("Loaded ", m$name, " as ", m$object,
                                     " (", nrow(state$data), " rows)."),
                              type = "message")
    }))

    output$about <- shiny::renderUI({
      m <- state$data_meta
      if (is.null(m)) {
        return(bslib::card(bslib::card_body(
          shiny::p("Choose a dataset and press ", shiny::strong("Use dataset"),
                   ". The abalone example (counts of shell rings) is a ",
                   "good first dataset."))))
      }
      d <- state$data
      bslib::card(bslib::card_header(m$title), bslib::card_body(
        shiny::p(nrow(d), "rows,", ncol(d), "columns,",
                 sum(!stats::complete.cases(d)), "rows with a missing value."),
        if (m$source == "bundled") shiny::p(
          "In the R code this dataset is ", shiny::code(m$object),
          ". Documentation: ", shiny::code(paste0("?flexCountReg::", m$name)),
          ". Suggested response: ", shiny::code(m$response), ".")))
    })
    output$columns <- DT::renderDT({
      shiny::req(state$data)
      DT::datatable(column_summary(state$data), rownames = FALSE,
                    options = list(pageLength = 15, dom = "tip"))
    })
    output$preview <- DT::renderDT({
      shiny::req(state$data)
      DT::datatable(utils::head(state$data, 50), rownames = FALSE,
                    options = list(pageLength = 10, scrollX = TRUE))
    })
  })
}
