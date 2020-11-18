ui <- fluidPage(
  
  # App title ----
  titlePanel("PRYNT Application"),
  mainPanel(hr(),"en 2 phrases",br(),br(),br() ),
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    # sidebarPanel(
      verticalLayout(
        wellPanel(
          fluidRow(
            column(2),
            column(6,
              textAreaInput("caption", "List of deregulated proteins (one of each line)", "", width = "600px",height = "200px"),
              fluidRow(
                column(3,"Not found : "),
                column(6,
                verbatimTextOutput("value"),tags$head(tags$style("#value{color: red;font-size: 20px;font-style: bold;}")),
                )
              ),hr(),
              fluidRow(
                column(3),
                column(3,
                  actionButton("confirmdatabutton",h4("Run PRYNT"),align="center", width ='300px')
                )
              )
            )
          )
        )
    
    
),
    
    # Main panel for displaying outputs ----
    mainPanel(
      conditionalPanel(condition ="input.confirmdatabutton!=0" ,
                 withSpinner(ui_element =dataTableOutput('table') ,color="#0dc5c1",type = 8,)

                        )
      # dataTableOutput("table")
    )
  )
)