# Diagnose page: rootogram, randomized quantile residuals and CURE plots of
# the selected fit, each with its data as a table and a download.

mod_diagnose_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 300,
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
      shiny::actionButton(ns("go"), "Calculate diagnostic",
                          class = "btn-primary"),
      shiny::downloadButton(ns("download"), "Download plot data")),
    bslib::card(bslib::card_header(shiny::textOutput(ns("title"))),
                shiny::plotOutput(ns("plot"), height = "420px"),
                shiny::uiOutput(ns("note"))),
    bslib::card(bslib::card_header("Plot data"),
                DT::DTOutput(ns("table"))))
}

mod_diagnose_server <- function(id, state) {
  shiny::moduleServer(id, function(input, output, session) {
    bundle <- shiny::reactive({
      shiny::req(state$active)
      state$fits[[state$active]]
    })
    shiny::observe({
      d <- state$data
      num <- names(d)[vapply(d, is.numeric, logical(1))]
      shiny::updateSelectInput(session, "cure_var",
                               choices = c("(fitted values)" = "", num))
    })
    result <- shiny::eventReactive(input$go, {
      b <- bundle()
      f <- b$fit
      if (input$diag == "rootogram") {
        mc <- if (is.na(input$max_count)) NULL else input$max_count
        rt <- flexCountReg::rootogram(f, max_count = mc, overflow = TRUE,
                                      plot = FALSE)
        list(type = "rootogram", data = as.data.frame(unclass(rt)[1:3]),
             fit = b)
      } else if (input$diag == "quantile") {
        set.seed(1)
        r <- stats::residuals(f, type = "quantile")
        list(type = "quantile", fit = b,
             data = data.frame(fitted = stats::fitted(f), residual = r))
      } else {
        v <- if (nzchar(input$cure_var %||% "")) input$cure_var
        p <- flexCountReg::cureplot(f, indvar = v)
        list(type = "cure", fit = b, plot = p,
             data = as.data.frame(p$data))
      }
    })
    output$title <- shiny::renderText({
      r <- result()
      paste0(switch(r$type, rootogram = "Rootogram",
                    quantile = "Randomized quantile residuals",
                    cure = "Cumulative residuals (CURE)"),
             " - fit: ", r$fit$label)
    })
    output$plot <- shiny::renderPlot({
      r <- result()
      d <- r$data
      if (r$type == "rootogram") {
        se <- sqrt(d$expected)
        so <- sqrt(d$observed)
        graphics::plot(range(d$count) + c(-0.5, 0.5),
                       range(c(0, se, se - so)), type = "n",
                       xlab = "Count (last bar: tail)",
                       ylab = "sqrt(frequency)")
        graphics::rect(d$count - 0.4, se - so, d$count + 0.4, se,
                       col = "grey80", border = "grey30")
        graphics::abline(h = 0, lty = 2)
        graphics::lines(d$count, se, lwd = 2)
        graphics::points(d$count, se, pch = 19)
      } else if (r$type == "quantile") {
        oldpar <- graphics::par(mfrow = c(1, 2))
        on.exit(graphics::par(oldpar))
        graphics::plot(d$fitted, d$residual, xlab = "Fitted mean",
                       ylab = "Quantile residual", pch = 1)
        graphics::abline(h = 0, lty = 2)
        stats::qqnorm(d$residual, main = "Normal Q-Q")
        stats::qqline(d$residual, lty = 2)
      } else {
        print(r$plot)
      }
    }, alt = "Diagnostic plot; its data are in the table below.")
    output$note <- shiny::renderUI({
      r <- result()
      if (r$type == "quantile") shiny::p(class = "small",
        "Randomized quantile residuals are standard normal when the model ",
        "is correct. For panel models they can be dependent within panels.")
      else if (r$type == "rootogram") shiny::p(class = "small",
        "Bars hang from the expected frequencies; bars ending above or below ",
        "zero mark counts the model under- or over-predicts.")
    })
    output$table <- DT::renderDT({
      d <- result()$data
      num <- vapply(d, is.numeric, logical(1))
      d[num] <- lapply(d[num], signif, 5)
      DT::datatable(d, rownames = FALSE, options = list(pageLength = 10))
    })
    output$download <- shiny::downloadHandler(
      filename = function() paste0("diagnostic-", result()$type, ".csv"),
      content = function(file) utils::write.csv(result()$data, file,
                                                row.names = FALSE))
  })
}
