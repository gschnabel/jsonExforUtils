
#' Parse Reaction Expression
#'
#' @param reacStr reactiong string as given in EXFOR, i.e. of form (TARGET(PROJ,PROCESS)RESIDUAL,...)
#' @param pos position in \code{reacStr} from where the parser should start reading.
#'
#' @return Returns a tree realized as recursive list.
#'         A leaf node is of the structure as explained in \link{parseReacStr}.
#'         A non-leaf node represents an arithmetic expression and contains the following elements:
#'         \describe{
#'           \item{start}{start position of expression in \code{reacStr}}
#'           \item{end}{end position of expression in \code{reacStr}}
#'           \item{reac}{the matched part of the reaction string}
#'           \item{type}{contains the string "expr"}
#'           \item{ops}{arithmetic operators separating the \code{nodes}}
#'           \item{nodes}{sub-expressions being part of this expression (list of lists)}
#'         }
#' @export
#'
parseReacExpr <- function(reacStr, pos=1, quiet=TRUE) {

  if (length(reacStr)!=1)
    stop("reaction string must be vector of length 1")
  if (!is.character(reacStr))
    stop("reaction string must be of type character")

  startPos <- pos
  # parse atomic reaction string
  tmpStr <- substring(reacStr, pos)
  if (grepl("^\\( *[^(]", tmpStr))
    return(parseReacStr(reacStr, pos, quiet))

  pos <- pos + 1
  tmpStr <- substring(reacStr, pos)
  if (!grepl("^ *\\(", tmpStr))
    return(debugReacInfo("expected bracket",
                         pos, reacStr, quiet))

  # parse recursive structure
  opVec <- character(0)
  tokenInfoList <- list()
  enteredLoop <- FALSE
  while (grepl("^ *[-+*/=]? *\\(", tmpStr)) {
    bracketPos <- as.numeric(regexpr("\\(", tmpStr)) + pos - 1
    if (bracketPos > pos) {
      curOp <- substring(reacStr, pos, bracketPos - 1)
      curOp <- gsub(" +", "", curOp)
    } else {
      curOp <- ""
    }
    curTokenInfo <- parseReacExpr(reacStr, bracketPos, quiet)
    if (is.null(curTokenInfo)) return(NULL)
    pos <- curTokenInfo$end + 1
    tmpStr <- substring(reacStr, pos)
    opVec <- c(opVec, curOp)
    tokenInfoList <- c(tokenInfoList, list(curTokenInfo))

    if (!enteredLoop && !isTRUE(curOp==""))
      return(debugReacInfo("no arithmetic operator is permitted here", pos, reacStr, quiet))
    enteredLoop <- TRUE
  }
  if (!enteredLoop)
    return(debugReacInfo("expected opening bracket", pos, reacStr, quiet))

  stopifnot(all(curOp[-1]!=""))
  tmpStr <- substring(reacStr,pos)
  if (!grepl("^ *\\)", tmpStr))
    return(debugReacInfo("closing bracket missing", pos, reacStr, quiet))
  opVec <- opVec[-1]
  stopifnot(all(opVec %in% c("+","-")) ||
              all(opVec %in% c("*")) ||
              all(opVec %in% c("/")) ||
              opVec=="=")

  matchedReacStr <- substring(reacStr, startPos, pos)

  tokenInfo <- list(
    start = startPos,
    end = pos,
    reac =  matchedReacStr,
    type = "expr",
    ops = opVec,
    nodes = tokenInfoList
  )
  return(tokenInfo)
}




#' Parse Reaction String
#'
#' @param reacStr the reaction string
#' @param pos position in \code{reacStr} where parsing should start
#'
#' @return A list with the following components:
#'         \describe{
#'           \item{start}{start position of reaction string}
#'           \item{end}{end position of reaction string}
#'           \item{reac}{the matched reaction string}
#'           \item{type}{contains the string "reac"}
#'           \item{projectile}{projectile}
#'           \item{process}{usually gives ejectile information}
#'           \item{target}{list describing target with elements \code{A, Z, sym}}
#'           \item{residual}{list describing residual nucleus with elements \code{A, Z, sym, meta}}
#'           \item{quantspec}{string with additional flags indicating among other things observable type, e.g. DA or SIG}
#'         }
#' @export
#'
parseReacStr <- function(reacStr, pos=1, quiet=TRUE) {

  if (!grepl("^ *\\(", reacStr))
    return(debugReacInfo("reaction strings have to be enclosed in brackets",
                         pos, reacStr, quiet))

  tmpStr <- substring(reacStr, pos)
  startIdx <- as.numeric(regexpr("\\(", tmpStr)) + pos - 1
  pos <- startIdx + 1

  ############################
  # extract target info
  ############################
  tmpStr <- substring(reacStr, pos)
  pat <- "^([0-9]+)-([A-Z]+)-([0-9]+)-?(M[1-9]?|G)?\\("
  regRes <- regexec(pat, tmpStr, perl=TRUE)
  if (isTRUE(regRes[[1]] == -1))
    return(debugReacInfo("target isotope specification wrong",
                         pos, reacStr, quiet))
  matchPos <- as.integer(regRes[[1]][1])
  matchLen <- attr(regRes[[1]], "match.length")[1]
  pos <- pos + matchPos + matchLen - 2
  regVals <- regmatches(tmpStr, regRes)[[1]]
  if (length(regVals) != 5)
    stop("stupid programmer! Why did you think length(regVals)==5")

  sysTarget <- list(
    Z = as.integer(regVals[2]),
    A = as.integer(regVals[4]),
    sym = regVals[3],
    meta = if (regVals[5]=="") NULL else regVals[5]
  )
  stopifnot(!is.na(sysTarget$A) && !is.na(sysTarget$Z))

  ############################
  # extract scattering process
  ############################
  tmpStr <- substring(reacStr, pos)
  pat <- "\\(([^,()]+),([^,()]+)\\)"
  regRes <- regexec(pat, tmpStr)
  if (isTRUE(regRes[[1]] == -1))
    return(debugReacInfo("incident particle/emitted particle specification wrong",
                       pos, reacStr, quiet))
  matchPos <- as.integer(regRes[[1]][1])
  matchLen <- attr(regRes[[1]], "match.length")[1]

  regVals <- regmatches(tmpStr, regRes)[[1]]
  sysProjectile <- regVals[2]
  sysProcess <- regVals[3]
  pos <- pos + matchPos + matchLen - 1

  ############################
  # extract residual nucleus
  ############################
  tmpStr <- substring(reacStr, pos)
  if (charAt(1,tmpStr)==",")
  {
    sysResidual <- NULL
    pos <- pos + 1
  }
  else if (grepl("^MASS,", tmpStr))
  {
    sysResidual <- "MASS"
    pos <- pos + 5
  }
  else
  {
    pat <- "^(?<Z>[0-9]+)-(?<sym>[A-Z]+)-(?<A>[0-9]+)-?(?<meta>(?:(?:[ML][1-9]?|G)[+/]?)+)?,"
    regRes <- regexpr2(pat, tmpStr)
    if (isTRUE(regRes == -1))
      return(debugReacInfo("residual isotope specification wrong",
                           pos, reacStr, quiet))
    matchPos <- as.integer(regRes)
    matchLen <- attr(regRes, "match.length")
    pos <- pos + matchPos + matchLen - 1  # also skip the subsequent comma
    sysResidual <- as.list(attr(regRes, "capture.strings"))
    names(sysResidual) <- attr(regRes, "capture.names")
    sysResidual$A <- as.integer(sysResidual$A)
    sysResidual$Z <- as.integer(sysResidual$Z)
    if (sysResidual$meta == "") sysResidual$meta <- NULL
  }

  ############################
  # extract additional reaction flags
  ############################
  tmpStr <- substring(reacStr, pos)
  pat <- "^([,A-Z/0-9]+)\\)"
  regRes <- regexec(pat, tmpStr)
  if (isTRUE(regRes[[1]] == -1))
    return(debugReacInfo("problem with additional reaction flags",
                         pos, reacStr, quiet))

  matchPos <- as.integer(regRes[[1]][1])
  matchLen <- attr(regRes[[1]], "match.length")[1]

  regVals <- regmatches(tmpStr, regRes)[[1]]
  reacModifiers <- regVals[2]
  pos <- pos + matchPos + matchLen - 1
  endPos <- pos - 1

  tokenInfo <- list(
    start = startIdx,
    end = endPos,
    reac = substring(reacStr, startIdx, endPos),
    type = "reac",

    target = sysTarget,
    residual = sysResidual,
    projectile = sysProjectile,
    process = sysProcess,
    quantspec = reacModifiers
  )
  return(tokenInfo)
}



#' Convert Reaction Structure to Reaction String
#'
#' @param reacStruc A nested list as returned by \link{parseReacExpr}
#'
#' @return Returns a reaction string according to the EXFOR format
#' @export
#'
reacStrucToStr <- function(reacStruc) {

  if (reacStruc$type == "reac") {

    targetStr <- with(reacStruc, paste0(target$Z, "-", target$sym, "-", target$A))
    if (is.null(reacStruc$residual)) residualStr <- ""
    else residualStr <- with(reacStruc, paste0(residual$Z, "-", residual$sym, "-", residual$A))
    if (!is.null(reacStruc$residual$meta))
      residualStr <- paste0(residualStr, "-", reacStruc$residual$meta)

    paste0("(", targetStr, "(", reacStruc$projectile, ",", reacStruc$process, ")",
           residualStr, ",", reacStruc$quantspec, ")")
  }
  else
  {
    stopifnot(reacStruc$type == "expr")
    curStr <- paste0("(", reacStrucToStr(reacStruc$nodes[[1]]))
    for (i in seq_along(reacStruc$ops)) {
      curStr <- paste0(curStr, " ", reacStruc$ops[i], " ",
                       reacStrucToStr(reacStruc$nodes[[i+1]]))
    }
    curStr <- paste0(curStr,")")
    curStr
  }
}


regexpr2 <- function(pat, strs) {

  rexRes <- regexpr(pat, strs, perl=TRUE)
  captureNames <- attr(rexRes, "capture.names")
  captureStart <- attr(rexRes, "capture.start")
  captureLength <- attr(rexRes, "capture.length")
  captureEnd <- captureStart + captureLength - 1
  namedMatches <- matrix(NA_character_,
                         ncol = length(captureNames),
                         nrow = length(strs))
  for (i in seq_along(strs))
    if (rexRes[i] != -1)
      namedMatches[i,] <- substring(strs[i], captureStart[i,], captureEnd[i,])
  colnames(namedMatches) <- captureNames
  attr(rexRes, "capture.strings") <- namedMatches
  rexRes
}



debugReacInfo <- function(errmsg, pos, reacStr, quiet) {

  if (!quiet)
    cat(paste0(errmsg, "\n",
               "reaction string: ", reacStr, "\n",
               "problem near position ", pos,
               " --- ", substring(reacStr, pos, pos+15), "\n"))
  NULL
}



