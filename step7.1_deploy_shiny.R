library(rsconnect)

rsconnect::setAccountInfo(name='ln-shiny-cancer',
                          token='9B847A81B8B91B4B628E99EB776281A0',
                          secret='r5q2uyjMHc6pfLCk1AhRiwZ6xOE4s960ldC+Rzin')

rsconnect::deployApp('./liling_bucm.phewas_MR/')
