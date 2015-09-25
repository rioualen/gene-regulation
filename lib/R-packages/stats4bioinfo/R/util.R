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

# #' @title Shade the area below a function, between two limits.
# #'
# #' @author Amalhelu
# #' 
# #' @source http://stackoverflow.com/questions/6786982/shaded-area-under-two-curves-using-r
# #'
# #' @description Sades the area under a curve, by plotting a filled polygon.
# #'
# #' @details
# #' Imported from an answer by Amalhelu on stackoverflow.
# #'
# #' @param fun   function
# #' @param from  min X value for the shaded area
# #' @param to  max X value for the shaded area
# #' @param length=100 number of points for the polygon
# #' @param col="gray" Shading color, passed to polygon() 
# #' @param border  Border color. By default the border is not drawed (NA), to let the initial curve unchanged.
# #' @param ... Additional parameters are passed to polygon
# #'
# #' @examples
# #'
# #' ## Shade the tails of a Normal distribution
# #' y <- function(x)sapply(x, function(xt)dnorm(xt,mean=2,sd=2))
# #' curve(y,from=-6,to=10, n=1000, col="darkblue")
# #' shadeUnderCurve(y,from=6, to=10, length=100, col="#DDBBFF", border="red")
# #' 
# #' ## A more complex example from Amalhelu
# #' ## Define two gaussian functions with different means and standard deviations
# #' y1 <- function(x)sapply(x, function(xt)dnorm(xt,mean=0,sd=1))
# #' y2 <- function(x)sapply(x, function(xt)dnorm(xt,mean=3,sd=2))
# #' 
# #' ## Define a function below both curves
# #' my.fun <- function(x){sapply(x, function(xt)min(y1(xt), y2(xt)))}
# #' 
# #' ## Plot the curves + area under the minimum.
# #' ## Note that we first shade the area, and then plot the curve to fix some overlap due to bitmap resolution
# #' plot(y1, -10, 10, col="darkred", type="n")
# #' shadeUnderCurve(my.fun, -10, 10, length=1000, col="#DDBBFF")
# #' curve(y2, add=TRUE, col="darkblue")
# #' curve(y1, add=TRUE, col="darkred")
# #'
# #' @export
# shadeUnderCurve <- function(fun, from, to, length=100, col="gray", border=NA, ...){
#   xvals <- seq(from, to, length=length)
#   dvals <- match.fun(fun)(xvals)
#   polygon(c(xvals,rev(xvals)),c(rep(0,length),rev(dvals)),col=col, border=border,...)
# }
# 


#' @title Plot a curve with a shaded area.
#'
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#'
#' @description Plot the curve for a given function, and shade the area 
#' below it between two given limits.
#'
#' @details
#' First version: 2015-03
#' Last modification: 2015-03
#'
#' @param fun   function
#' @param xlim  min and max X value to plot the curve
#' @param shade.from=xlim[1]  min X value for the shaded area
#' @param shade.to=xlim[2]    max X value for the shaded area
#' @param length=1000 number of points for the polygon
#' @param shade.col="gray" Shading color, passed to polygon() 
#' @param shade.border  Border color. By default the border is not drawed (NA), to let the initial curve unchanged.
#' @param draw.curve = TRUE If FALSE, only plot the shaded area
#' @param add=FALSE    Add the shaded area to existing plot
#' @param ... Additional parameters are passed to fun via lapply
#'
#' @examples
#'
#' ## Shade the interval (-1,1) below Student density function
#' z <- 3.5
#' shadeArea(dt, xlim=c(-5,+5), shade.from=-1, shade.to=1, df=4, curve.col="darkred", shade.col="#DDBBEE", shade.border="red")
#' 
#' @export
shadeArea <- function(fun, 
                      xlim, 
                      shade.from=xlim[1],
                      shade.to=xlim[2],
                      length = 1000,
                      curve.col="blue",
                      shade.col="gray", 
                      shade.border=NA, 
                      draw.curve=TRUE,
                      add=FALSE,
                      ...) {
  
  x.values <- seq(xlim[1], xlim[2], length=length)
  y.values <- lapply(x.values, fun, ...)

  if (draw.curve) {  
    if (add) {
      lines(x.values,y.values, type='l', col=curve.col)
    } else {
      plot(x.values,y.values, type='l', col=curve.col)
    }
  } else if (!add) {
    plot(x.values,y.values, type='n')
  }
  
  x.to.shade <- x.values [x.values <= shade.to & x.values >= shade.from]
  y.to.shade <- y.values [x.values <= shade.to & x.values >= shade.from]
  polygon(c(x.to.shade,rev(x.to.shade)),c(rep(0,length(x.to.shade)),rev(y.to.shade)),col=shade.col, border=shade.border)
}

#' @title Plot a distribution and shade one or two tails.
#'
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#'
#' @description Plot a distribution and shade one or two tails.
#'
#' @details
#' First version: 2015-03
#' Last modification: 2015-03
#'
#' @param fun   function
#' @param q     quantile above (below) which the curve should be shaded on the upper (lower) tail.
#' @param tails="both" Tails under which the curve should be shaded. Accepted values: "both", "lower", "upper".
#' @param ... Additional parameters are passed to shadeArea()
#'
#' @examples
#'
#' ## Shade the interval (-1,1) below Student density function
#' shadeTails(dt, xlim=c(-5,5), q=2,df=3, shade.col="pink", tails="both",
#'            curve.col="darkblue",shade.border="darkblue")
#' 
#' @export
shadeTails <- function(fun,
                       xlim,
                       q,
                       tails="both",
                       add=FALSE,
                       ...) {

  
  if (tails%in% c("upper", "both")) {
    verbose("Shading upper tail")
    shadeArea(fun, xlim, shade.from=q, add=add, ...)
  }
  
  if (tails%in% c("lower", "both")) {
    verbose("Shading lower tail")
    shadeArea(fun, xlim, shade.to=-q, add=TRUE, draw.curve=FALSE,...)
  }
  
  
}
