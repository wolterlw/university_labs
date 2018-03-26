library(shiny)
library(ggplot2)
library(pracma)
library(rmutil)
library(deSolve)


ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      body {
        background: #eee;
      }
      h2 {
        text-align: center;
      }
    "))
  ),
  
  tabsetPanel(
    tabPanel(
      "Модель Фергюльста",
      fluidRow(
        
        column(
          width = 12, 
          plotOutput("Fergulst", height = "600px") 
        )
      ),
      fluidRow(
        column(
          width = 4, 
          wellPanel(
            #SD for Fergulst model
            numericInput('r', 'r', 0.1,
                         min = -100, max = 100, step = 0.1),
            numericInput('k', 'k', 100,
                         min = 0, max = 100000),
            numericInput('N0', 'N₀', 10,
                         min = 0, max = 100000)
          )
        )
      )
    ),
    tabPanel(
      "Модель Солоу",
      fluidRow(
        column(
          width = 12, 
          plotOutput("Solow", height = "600px")
        )
      ),
      fluidRow(
        column(
          width = 4, 
          wellPanel(
            numericInput('s', 's', 0.2, min = -100, max = 100, step = 0.1),
            numericInput('k0', 'k₀', 0.8, min = 0, max = 100000),
            numericInput('a', 'a', 2.5, min = 0, max = 100000),
            numericInput('m', 'm', 0.1, min = 0, max = 100000),
            numericInput('alpha', 'α', 0.3, min = 0, max = 100000),
            numericInput('q', 'q', 0.1, min = 0, max = 100000)
          )
        )
      )
    ),
    tabPanel(
      "Рівняння вимушених коливань",
      fluidRow(
        column(12, plotOutput("Fluct", height = "600px"))
      ),
      fluidRow(
        column(12, plotOutput("Fluct2", height = "600px"))
      ),
      fluidRow(
        column(
          width = 4, 
          wellPanel(
            numericInput('delta', 'δ', 0.1, min = 0, max = 1, step = 0.1),
            numericInput('w0', 'w₀', 2, min = 0, max = 100000),
            numericInput('f0', 'f₀', 0, min = 0, max = 100000),
            numericInput('w', 'w', 1, min = 0, max = 100000),
            numericInput('x0_der', 'x₀ похідна', 2, min = 0, max = 100000),
            numericInput('x0_fluct', 'x₀', 3, min = 0, max = 100000)
          )
        )
      )
    ),
    tabPanel(
      "Система «хижак-жертва»",
      fluidRow(
        column(12,  plotOutput("predator_prey", height = "600px"))
      ),
        fluidRow(
        column(12, plotOutput("predator_prey2", height = "600px"))
      ),
      fluidRow(
        column(
          width = 4, 
          wellPanel(
            numericInput('alpha_x', 'αₓ', 0.05, min = 0, max = 1, step = 0.1),
            numericInput('beta_x', 'βₓ', 2, min = 0, max = 100000),
            numericInput('alpha_y', 'αy', 5, min = 0, max = 100000),
            numericInput('beta_y', 'βy', 1, min = 0, max = 100000),
            numericInput('y0', 'y₀', 5, min = 0, max = 100000),
            numericInput('x0_prey', 'x₀', 40, min = 0, max = 100000)
          )
        )
      )
    )
  )
);
 
server <- function(input, output) {
    
    output$Fergulst <- renderPlot({
      Fergulst <- function(N, t) input$r*N*(input$k - N) 
      x_vec <- seq(0, 10, 0.01) 
      solvFirst <- runge.kutta(Fergulst, input$N0 ,x_vec) 
      solvFirst[solvFirst<0] <- 0
      
      ggplot(data.frame(time = x_vec, population = solvFirst), aes(x = time, y = population)) + 
        geom_line(colour = 'red')
      
    })
    
    output$Solow <- renderPlot({
      Slow <- function(k, t) {
        input$s*input$a*k^input$alpha - (input$m + input$q)*k
      }
      
      x_vec2 <- seq(0, 100, 1) 
      solvFirst <- runge.kutta(Slow, input$k0 ,x_vec2) 
      solvFirst[solvFirst<0] <- 0
      
      ggplot(data.frame(time = x_vec2, capitallabor = solvFirst), aes(x = time, y = capitallabor)) + 
        geom_line(colour = 'red') 
      
    })
    
    output$Fluct <- renderPlot({
      f0 <- input$f0
      w <- input$w
      delta <- input$delta
      
      Fluct <- function(t, x)
        as.matrix(c(input$f0*cos(input$w*t) - 2*input$delta*x[1] - input$w0^2*x[2], x[1]))
      t0 <- 0; tf <- 20
      x0 <- as.matrix(c(input$x0_der, input$x0_fluct))
      sol <- ode23(Fluct, t0, tf, x0)
      
      plot(c(0, 20), c(-3, 3), type = "n",
           xlab = "час", ylab = "", main = "")
      lines(sol$t, sol$y[, 2], col = "red")
      grid()
      
    })
    output$Fluct2 <- renderPlot({
      # Fluct model (fazovoe)
      f0 <- input$f0
      w <- input$w
      delta <- input$delta
      
      Fluct <- function(t, x)
        as.matrix(c(input$f0*cos(input$w*t) - 2*input$delta*x[1] - input$w0^2*x[2], x[1]))
      t0 <- 0; tf <- 20
      x0 <- as.matrix(c(input$x0_der, input$x0_fluct))
      sol <- ode23(Fluct, t0, tf, x0)
      
      plot(c(-5, 5), c(-5, 5), type = "n",
           xlab = "час", ylab = "", main = "")
      lines(sol$y[, 1], sol$y[, 2], col = "red")
      grid()
      
    })
    
    output$predator_prey <- renderPlot({
      LotVmod <- function (Time, State, Pars) {
        with(as.list(c(State, Pars)), {
          dx = x*(alpha*y - beta)
          dy = y*(gamma - delta*x)
          return(list(c(dx, dy)))
        })
      }
      
      Pars <- c(alpha = input$alpha_x, beta = input$beta_x, gamma = input$alpha_y, delta = input$beta_y)
      State <- c(x = input$x0_prey, y = input$y0)
      Time <- seq(0, 100, by = 0.01)
      
      out <- as.data.frame(ode(func = LotVmod, y = State, parms = Pars, times = Time))
      
      matplot(out[,-1], type = "l", xlab = "час", ylab = "популяція")
      
      legend("topright", c("Хижаки", "Жертви"), lty = c(1,2), col = c(1,2), box.lwd = 0)
    })
    
    output$predator_prey2 <- renderPlot({
      LotVmod <- function (Time, State, Pars) {
        with(as.list(c(State, Pars)), {
          dx = x*(alpha*y - beta)
          dy = y*(gamma - delta*x)
          return(list(c(dx, dy)))
        })
      }
      
      Pars <- c(alpha = input$alpha_x, beta = input$beta_x, gamma = input$alpha_y, delta = input$beta_y)
      State <- c(x = input$x0_prey, y = input$y0)
      Time <- seq(0, 100, by = 0.01)
      
      out <- as.data.frame(ode(func = LotVmod, y = State, parms = Pars, times = Time))
      
      plot(out[,2], out[,3], type = "p", xlab = "хижаки", ylab = "жертви")
    })
  }
  
shinyApp(ui = ui, server = server)