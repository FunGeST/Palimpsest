#' "%!in%" 
#'
#' The opposite of "%in%"
#' @keywords Signatures
#' @export
#' @examples
#' days <- c("monday","tuesday","wednesday","thursday","friday","saturday","sunday")
#' weekdays <- days[days %!in% c("saturday", "sunday")]

'%!in%' <- function(x,y)!('%in%'(x,y))
