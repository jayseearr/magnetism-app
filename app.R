#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Model functions

# create 2D Ising lattice
ising <- function(N, Nc = NULL) {
    if (!is.null(Nc)) {
        Nr <- N
        Nc <- Nc
    } else {
        Nr <- N
        Nc <- N
    }
    X <- matrix(runif(Nc*Nr) >= 0.5, ncol=Nc, nrow=Nr) * 2.0 - 1
}

# spin for the input Ising lattice
spin <- function(X) {
    sum(X)
}

# energy for only the row/column location
sub_energy <- function(X, x, y, J, B) {
    nr <- nrow(X)
    nc <- ncol(X)
    column <- x
    row <- y
    if (nc > 1) {
        xl <- X[row, ((column - 1) - 1) %% nc + 1]
        xr <- X[row, ((column - 1) + 1) %% nc + 1]
    } else {
        xl <- 0
        xr <- 0
    }
    if (nr > 1) {
        xu <- X[((row - 1) + 1) %% nr + 1, column]
        xd <- X[((row - 1) - 1) %% nr + 1, column]
    } else {
        xu <- 0
        xd <- 0
    }
    J * X[row, column] * (xl + xr + xu + xd) - B * X[row, column]
    
}

# energy for input Ising lattice, with parameters for neighbor interaction (J),
# external field (B)
energy <- function(X, J, B) {
    nr <- nrow(X)
    nc <- ncol(X)
    
    Xu <- rbind(X[nr,], X[1:(nr-1),])
    Xd <- rbind(X[2:nr,], X[1,])
    Xr <- cbind(X[,nc], X[,1:(nc-1)])
    Xl <- cbind(X[,2:nc], X[,1])
    
    # the 0.5 is to avoid double counting nearest neighbors
    0.5 * J * sum(X * (Xr + Xl + Xu + Xd)) - B * spin(X)
    
}

# single metropolis step of the input Ising lattice
# at temperature T, interaction J, and optional previous energy.
metro <- function(X, temperature, J, B, energy = NULL, nsteps = 1, 
                  ntrials = length(X)) {
    if (is.null(energy)) stop("Initial energy is required.")
    if (ntrials < 1) stop("At least 1 trial is required.")
    k <- 1
    p0 <- 0.5  # probability that a dE of 0 will cause a flip (defaults to 1)
    e0 <- energy
    
    # recursive call for multiple steps
    if (nsteps > 1) {
        for (n in 1:nsteps) {
            out <- metro(X, temperature, J, B, energy = energy, nsteps = 1)
            X <- out$X
            energy <- out$energy
        }
        return(out)
    }
    
    nr <- nrow(X)
    nc <- ncol(X)
    row <- round(runif(ntrials, min = 1, max = nrow(X)))
    column <- round(runif(ntrials, min = 1, max = ncol(X)))
    for (n in 1:ntrials) {
        xl <- X[row[n], ((column[n] - 1) - 1) %% nc + 1]
        xr <- X[row[n], ((column[n] - 1) + 1) %% nc + 1]
        xu <- X[((row[n] - 1) + 1) %% nr + 1, column[n]] 
        xd <- X[((row[n] - 1) - 1) %% nr + 1, column[n]]
        dE <- (-2 * J * X[row[n], column[n]] * (xl + xr + xu + xd) 
               - 2 * B * X[row[n], column[n]])
        
        if ( (dE == 0 && p0 >= runif(1)) ||
             (dE < 0) ||
             (temperature > 0 && exp(-dE / k / temperature) >= runif(1))
        ) {
            X[row[n], column[n]] <- -X[row[n], column[n]]
            energy <- energy + dE
        }
    }
    list(X = X, energy = energy, initial.energy = e0)
}

### ----------------------------------------------------------------------------
# UI Function
# Define UI for application that draws a histogram
shinyUI <- fluidPage(
    titlePanel("Ising Lattice Model"),
    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            numericInput("n.spins", "Lattice Size:",
                         min=1, max=1000, value=20),
            actionButton("generate", "Generate New Lattice"),
            
            br(), br(), br(),
            
            radioButtons("interaction", "Spin-Spin Interaction (J)",
                         choiceNames = c("Ferromagnetic", "Antiferromagnetic"),
                         choiceValues = c(-1, 1)),
            br(),
            
            numericInput("temperature", 
                         "Temperature (relative to interaction, J)",
                         min = 0, max = 1000, value = 1, step = 1),
            numericInput("b.field", 
                         "Magnetic field (relative to interaction, J)", 
                         value = 0),
            br(),
            
            actionButton("take.step", "Take a Step"),
            actionButton("take.n.steps", "Take N Steps"),
            actionButton("run", "Run"),
            actionButton("stop", "Stop"),
            numericInput("n.steps", "Number of steps (N)",
                         min=1, max=1000, value=10),
            
            br(),
            sliderInput("time.per.step", "Seconds per step", min = 0.01, max = 1,
                        value = 0.5)
            
            #sliderInput("n.tries.per.step", "Spin flip trials per step",
            #            min = 1, max = 20*20, value = 20*20),
            
        ),
        
        # Plot of the lattice and time series of magnetic moment
        mainPanel(
            plotOutput("lattice.plot",
                       click = "lattice_click",
                       height = 500, width = 500
            ),
            column(4,
                   textOutput("step.info"),
                   br(),
                   htmlOutput("energy.info"),
                   br(),
                   htmlOutput("spin.info"),
                   br()
            ),
            column(8,
                   plotOutput("time.series", height = "250px")
            )
        )
    )
)

### ----------------------------------------------------------------------------
### Server Function
shinyServer <- function(input, output, session) {
    
    state <- reactiveValues(lattice = ising(20, 20),
                            steps.taken = 0,
                            J = -1,
                            B = 0,
                            running = FALSE,
                            dat = NULL,
                            prev.energy = NULL)
    
    metro_step <- function(y, temperature, nsteps, ntrials) {
        X <- y$lattice
        e <- y$prev.energy
        if (is.null(e)) e <- energy(X, J = y$J, B = y$B)
        out <- metro(X = X,
                     temperature = temperature,
                     J = y$J,
                     B = y$B,
                     energy = e,
                     nsteps = nsteps,
                     ntrials = ntrials)
        y$lattice <- out$X
        y$prev.energy <- out$energy
        y$steps.taken <- y$steps.taken + nsteps
        if (is.null(y$dat)) {
            
            y$dat <- data.frame(step = y$steps.taken,
                                moment = sum(y$lattice),
                                energy = out$energy)
        } else {
            y$dat <- rbind(y$dat, c(y$steps.taken,
                                    sum(y$lattice),
                                    out$energy)
                           )
        }
        y
    }
    
    # Change state variables when a UI event occurs
    
    observeEvent(input$interaction,
                 state$J <- as.numeric(input$interaction)
    )
    
    observeEvent(input$generate, {
        state$lattice <- ising(input$n.spins, input$n.spins)
        state$steps.taken <- 0
        # if (input$fix.tries.per.step) {
        #     updateSliderInput(session, inputId = "n.tries.per.step",
        #                       max = input$n.spins^2)
        # } else {
        #     updateSliderInput(session, inputId = "n.tries.per.step",
        #                       max = input$ny.spins * input$nx.spins,
        #                       value = input$ny.spins * input$nx.spins)
        # }
        state$dat <- NULL
    })
    
    tries_per_step <- reactive({
        ifelse(input$fix.tries.per.step, length(state$lattice), input$n.tries.per.step)
    })

    observe({
        state$running
        if (isolate(state$running)) {
            invalidateLater(1000 * input$time.per.step)
            isolate(state <- metro_step(state, input$temperature, 1,
                                        input$n.spins^2)
            )
        }
    })
    
    observeEvent(input$take.step, {
        state <- metro_step(state, input$temperature, 1, input$n.spins^2)
    })
    
    observeEvent(input$take.n.steps, {
        state <- metro_step(state, input$temperature, input$n.steps, 
                            input$n.spins^2)
    })
    
    observeEvent(input$run, {
        state$running <- TRUE
    })
    
    observeEvent(input$stop, {
        state$running <- FALSE
    })
    
    observeEvent(input$lattice_click, {
        col.row <- round(c(input$lattice_click$x, input$lattice_click$y))
        state$lattice[col.row[2], col.row[1]] <- -state$lattice[col.row[2], col.row[1]]
    })
    
    output$lattice.plot <- renderPlot({
        dims <- dim(state$lattice)
        image(1:dims[2], 1:dims[1], t(state$lattice),
              xlab="", ylab="", xaxt="n", yaxt="n", bty="n", asp=1) # asp makes the aspect ratio 1.
    })
    
    output$time.series <- renderPlot({
        if (is.null(state$dat)) {
            #plot(type = "none")
            return(NULL)
        }
        y <- state$dat$moment
        var.color <- "red"
        plot(state$dat$step, y, type="l", col=var.color,
             xlab="step", ylab=input$time.series.var)
    })
    
    output$step.info <- renderText(paste0("Total steps: ", state$steps.taken))
    
    output$energy.info <- renderUI({
        e <- energy(state$lattice, state$J, state$B)
        HTML(paste0("Energy: ", e, "<br>", "  (", e/length(state$lattice), " per spin)"))
    })
    
    output$spin.info <- renderUI({
        mu <- sum(state$lattice)
        HTML(paste0("Moment: ", mu, "<br>", "  (", mu/length(state$lattice), " per spin)"))
    })
    
}

### ----------------------------------------------------------------------------
# Run the application 
shinyApp(ui = shinyUI, server = shinyServer)
