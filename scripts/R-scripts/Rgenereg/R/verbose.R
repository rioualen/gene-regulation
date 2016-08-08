#' @title Display messages at a given verbosity level
#'
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#'
#' @description Display messages depending on user-defined verbosity level
#'
#' @details
#' First version: 2015-03.
#' Last modification: 2015-03.
#'
#' @param verbosity   Level of verbosity above which the message should be printed.
#' @param print.date=TRUE   Print date and time
#'
#' @examples
#'
#' verbosity <- 1 ## Define level of verbosity
#'
#' ## This message will be printed because the level is <= verbosity
#' verbose("This is printed", 1)
#'
#' ## This message will not be printed because the verbosity is inferior to the specified level
#' verbose("This is not printed", 2)
#'
#' @export
verbose <- function(message.content,
                    level=1,
                    print.date=TRUE) {
  if (!exists("verbosity")) {
    verbosity <- 1
  }
  if (verbosity >= level) {
    if (print.date) {
      message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\t", message.content)
    } else {
      message(message.content)
    }
  }
}
