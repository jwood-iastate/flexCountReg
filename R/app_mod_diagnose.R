# Diagnose page: rootogram, randomized quantile residuals and CURE plots of
# the selected fit. Each calculation runs its compiled code (which is added
# to the session script), and the plot is drawn from that same code.

mod_diagnose_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 320,
      shiny::selectInput(ns("fit"), "Model", character()),
      shiny::uiOutput(ns("fit_note")),
      shiny::selectInput(ns("diag"), "Diagnostic",
                         c("Rootogram" = "rootogram",
                           "Quantile residuals" = "quantile",
                           "CURE plot" = "cure")),
      shiny::conditionalPanel(
        "input.diag == 'rootogram'", ns = ns,
        shiny::numericInput(ns("max_count"), "Largest count shown", NA,
                            min = 0, step = 1),
        help_text("Empty: the largest observed count. The last bar is ",
                  "the tail Y > largest count shown.")),
      shiny::conditionalPanel(
        "input.diag == 'cure'", ns = ns,
        shiny::selectInput(ns("cure_var"), "Covariate (empty: fitted values)",
                           character()),
        help_text("The red limits drawn by cureplots are a descriptive ",
                  "reference, not a calibrated band for this model.")),
      shiny::uiOutput(ns("unavailable")),
      shiny::actionButton(ns("go"), "Calculate diagnostic",
                          class = "btn-primary"),
      shiny::downloadButton(ns("download"), "Download plot data"),
      bslib::accordion(
        open = FALSE, class = "mt-3",
        bslib::accordion_panel(
          "Save the plot as JPEG",
          shiny::numericInput(ns("jpg_width"), "Width (inches)", 7,
                              min = 1, max = 30, step = 0.5),
          shiny::numericInput(ns("jpg_height"), "Height (inches)", 5,
                              min = 1, max = 30, step = 0.5),
          shiny::selectInput(ns("jpg_res"), "Resolution (dots per inch)",
                             c("150" = 150, "300 (print)" = 300,
                               "600" = 600, "1200" = 1200),
                             selected = 300),
          help_text("The code that writes the file is added to the session ",
                    "script."),
          shiny::downloadButton(ns("download_jpg"), "Download JPEG")))),
    bslib::card(bslib::card_header(shiny::textOutput(ns("title"))),
                fill = FALSE,
                shiny::plotOutput(ns("plot"), height = "420px"),
                shiny::uiOutput(ns("note"))),
    bslib::card(bslib::card_header("Code that produced it"), fill = FALSE,
                shiny::verbatimTextOutput(ns("code"))),
    bslib::card(bslib::card_header("Plot data"), fill = FALSE,
                DT::DTOutput(ns("table"))))
}

# The data behind a diagnostic object, as a data frame
diagnostic_data <- function(type, value, fit) {
  switch(type,
    rootogram = as.data.frame(unclass(value)[c("count", "observed",
                                               "expected")]),
    quantile = data.frame(fitted = stats::fitted(fit), residual = value),
    cure = as.data.frame(value$data))
}

mod_diagnose_server <- function(id, state) {
  shiny::moduleServer(id, function(input, output, session) {
    current <- shiny::reactiveVal(NULL)
    app_fit_selector(input, session, state, "fit")
    fit_id <- shiny::reactive({
      id <- input$fit
      if (is.null(id) || !nzchar(id)) NULL else id
    })
    output$fit_note <- shiny::renderUI({
      b <- state$fits[[fit_id() %||% ""]]
      help_text(if (is.null(b)) "Estimate a model first (Build page)." else
        paste0("Object ", b$fit_name, "; data ", b$spec$data$object, "."))
    })
    shiny::observe({
      d <- state$data
      num <- names(d)[vapply(d, is.numeric, logical(1))]
      shiny::updateSelectInput(session, "cure_var",
                               choices = c("(fitted values)" = "", num))
    })
    support <- shiny::reactive({
      f <- app_fit_object(state, fit_id() %||% "")
      if (is.null(f)) return(list(ok = FALSE, reason = NULL))
      app_supports(f, input$diag)
    })
    output$unavailable <- shiny::renderUI({
      s <- support()
      if (!isTRUE(s$ok) && !is.null(s$reason)) {
        shiny::div(class = "alert alert-secondary small", s$reason)
      }
    })
    shiny::observeEvent(input$go, app_try(what = "The diagnostic", {
      id <- fit_id()
      if (is.null(id) || is.null(state$fits[[id]])) {
        stop("Estimate a model first.", call. = FALSE)
      }
      s <- support()
      if (!isTRUE(s$ok)) stop(s$reason, call. = FALSE)
      req <- switch(input$diag,
        rootogram = list(type = "rootogram",
                         max_count = if (!is.na(input$max_count))
                           input$max_count),
        quantile = list(type = "quantile"),
        cure = list(type = "cure",
                    var = if (nzchar(input$cure_var %||% "")) input$cure_var))
      res <- app_run_output(state, req, id)
      res$type <- input$diag
      res$data <- diagnostic_data(input$diag, res$value,
                                  app_fit_object(state, id))
      current(res)
    }))
    output$title <- shiny::renderText({
      r <- current()
      if (is.null(r)) return("No diagnostic yet")
      paste0(switch(r$type, rootogram = "Rootogram",
                    quantile = "Randomized quantile residuals",
                    cure = "Cumulative residuals (CURE)"),
             " - ", r$fit$label, " (", r$out, ")")
    })
    output$plot <- shiny::renderPlot({
      r <- current()
      shiny::req(r)
      # Redraw from the logged code, in a child of the session environment
      # (the session's objects are not modified)
      run_compiled(r$exprs, new.env(parent = state$env))
    }, alt = "Diagnostic plot; its data are in the table below.")
    output$code <- shiny::renderText({
      r <- current()
      shiny::req(r)
      paste(vapply(r$exprs, .deparse, character(1)), collapse = "\n")
    })
    output$note <- shiny::renderUI({
      r <- current()
      shiny::req(r)
      if (r$type == "quantile") help_text(
        "Randomized quantile residuals are standard normal when the model ",
        "is correct. For panel models they can be dependent within panels.")
      else if (r$type == "rootogram") help_text(
        "Bars hang from the expected frequencies; bars ending above or below ",
        "zero mark counts the model under- or over-predicts.")
    })
    output$table <- DT::renderDT({
      r <- current()
      shiny::req(r)
      d <- r$data
      num <- vapply(d, is.numeric, logical(1))
      d[num] <- lapply(d[num], signif, 5)
      DT::datatable(d, rownames = FALSE, fillContainer = FALSE,
                    options = list(pageLength = 10))
    })
    output$download <- shiny::downloadHandler(
      filename = function() paste0(current()$out %||% "diagnostic", ".csv"),
      content = function(file) utils::write.csv(current()$data, file,
                                                row.names = FALSE))
    output$download_jpg <- shiny::downloadHandler(
      filename = function() paste0(current()$out %||% "plot", ".jpg"),
      content = function(file) {
        r <- current()
        if (is.null(r)) stop("Calculate a diagnostic first.", call. = FALSE)
        app_save_plot(state, r, file, width = input$jpg_width,
                      height = input$jpg_height,
                      res = as.numeric(input$jpg_res))
      })
    # The button sits in a collapsed panel at first; keep it ready
    shiny::outputOptions(output, "download_jpg", suspendWhenHidden = FALSE)
  })
}
