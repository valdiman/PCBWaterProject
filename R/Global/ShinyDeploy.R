
install.packages('rsconnect')

rsconnect::setAccountInfo(name='valdiman', token='0F7C710D1AF5503E411EF82C60E7F695',
                          secret='eecFcwuY6FM97Qvet85zDZJtjF9pth/JsTk6OhxQ')

library(rsconnect)

rsconnect::deployApp()
