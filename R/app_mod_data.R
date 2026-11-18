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
                          class = "btn-primary")),
    shiny::uiOutput(ns("about")),
    bslib::card(bslib::card_header("Columns"),
                DT::DTOutput(ns("columns"))),
    bslib::card(bslib::card_header("First rows"),
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
    shiny::observeEvent(input$use, {
      if (input$source == "bundled") {
        d <- load_bundled_data(input$bundled)
        bd <- fcr_bundled_datasets()
        meta <- list(source = "bundled", name = input$bundled,
                     title = bd$title[bd$name == input$bundled],
                     response = bd$response[bd$name == input$bundled],
                     group = bd$group[bd$name == input$bundled])
      } else {
        f <- input$file
        if (is.null(f)) {
          shiny::showNotification("Choose a CSV file first.", type = "warning")
          return()
        }
        d <- tryCatch(utils::read.csv(f$datapath, check.names = FALSE,
                                      stringsAsFactors = FALSE),
                      error = function(e) {
                        shiny::showNotification(
                          paste("The file could not be read as CSV:",
                                conditionMessage(e)), type = "error")
                        NULL
                      })
        if (is.null(d)) return()
        meta <- list(source = "upload", name = f$name,
                     title = paste("Uploaded file", f$name),
                     response = NULL, group = NA)
      }
      state$data <- d
      state$data_meta <- meta
      shiny::showNotification(paste("Using", meta$name, "-", nrow(d), "rows."),
                              type = "message")
    })

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
          "Documentation: ", shiny::code(paste0("?flexCountReg::", m$name)),
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
