#------------------------------------
# New S3 methods defined in `cpca`
#------------------------------------

#' @rdname comprcompClass
#' @export varplot
varplot <- function(object, ...) UseMethod("varplot")

#' @rdname comprcompClass
#' @export compvar
compvar <- function(object, ...) UseMethod("compvar")


#' @rdname comprcompClass
#' @export compscore
compscore <- function(object, ...) UseMethod("compscore")

#------------------------------------
# S3 methods exported from `pls`
#------------------------------------

# @rdname comprcompClass
# @export scoreplot
#scores <- function(object, ...) UseMethod("scores")

# @rdname comprcompClass
# @export scoreplot
#scoreplot <- function(object, ...) UseMethod("scoreplot")

