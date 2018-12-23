## All functions of the correction of SG_Seq

convertToSGFeatures2 <- function (x, coerce = FALSE, merge = FALSE) 
{
  if (!is(x, "TxFeatures")) {
    stop("x must be a TxFeatures object")
  }
  if (length(x) == 0) {
    return(SGFeatures())
  }
  if (coerce) {
    features <- granges(x)
    mcols(features)$type <- as.character(type(x))
    splice5p <- mcols(features)$type %in% c("I", "L")
    splice3p <- mcols(features)$type %in% c("I", "F")
    splice5p[mcols(features)$type == "J"] <- NA
    splice3p[mcols(features)$type == "J"] <- NA
    mcols(features)$type[mcols(features)$type != "J"] <- "E"
    mcols(features)$splice5p <- splice5p
    mcols(features)$splice3p <- splice3p
    mcols(features)$txName <- txName(x)
    mcols(features)$geneName <- geneName(x)
  }
  else {
    features <- processFeatures2(x, merge = FALSE)
  }
  features <- SGSeq:::addFeatureID(features)
  features <- SGSeq:::addGeneID(features)
  features <- SGSeq:::SGFeatures(features)
  if (!coerce) {
    features <- annotate2(features, x)
  }
  return(features)
}

processFeatures2 <- function (features, coerce =F, merge = F) 
{
  junctions <- granges(features)[type(features) == "J"]
  junctions_D <- flank(junctions, -1, TRUE)
  junctions_A <- flank(junctions, -1, FALSE)
  mcols(junctions)$type <- rep("J", length(junctions))
  if (is(features, "TxFeatures")) {
    exons <- features[type(features) %in% c("I", "F", "L", 
                                            "U")]
    exons_D <- flank(features[type(features) %in% c("I", 
                                                    "F")], -1, FALSE)
    exons_A <- flank(features[type(features) %in% c("I", 
                                                    "L")], -1, TRUE)
  }
  else if (is(features, "SGFeatures")) {
    exons <- features[type(features) == "E"]
    exons_D <- flank(features[splice3p(features)], -1, FALSE)
    exons_A <- flank(features[splice5p(features)], -1, TRUE)
  }
  exons <- granges(exons)
  exons_D <- granges(exons_D)
  exons_A <- granges(exons_A)
  D <- unique(c(junctions_D, exons_D))
  mcols(D)$type <- rep("D", length(D))
  A <- unique(c(junctions_A, exons_A))
  mcols(A)$type <- rep("A", length(A))
  splicesites <- c(D, A)
  other <- c(junctions, splicesites)
  exons <- disjoin(exons)
  exons_start <- flank(exons, -1, TRUE)
  exons_end <- flank(exons, -1, FALSE)
  i_q <- which(!exons_end %over% splicesites)
  i_s <- which(!exons_start %over% splicesites)
  ol <- findOverlaps(suppressWarnings(flank(exons[i_q], 1, 
                                            FALSE)), exons_start[i_s])
  if ((length(ol) > 0)) {
    qH <- i_q[queryHits(ol)]
    sH <- i_s[subjectHits(ol)]
    i_to_be_merged <- union(qH, sH)
    d <- data.frame(from = qH, to = sH)
    v <- data.frame(name = i_to_be_merged)
    g <- graph.data.frame(d = d, directed = TRUE, vertices = v)
    k <- clusters(g)$membership
    exons_to_be_merged <- split(exons[i_to_be_merged], k)
    exons_merged <- unlist(reduce(exons_to_be_merged))
    if (length(exons_to_be_merged) != length(exons_merged)) {
      stop("cannot merge non-adjacent exons")
    }
    if (merge) exons <- c(exons[-i_to_be_merged], exons_merged)
  }
  exons_start <- flank(exons, -1, TRUE)
  exons_end <- flank(exons, -1, FALSE)
  splice5p <- rep(FALSE, length(exons))
  i_spliced <- unique(queryHits(findOverlaps(exons_start, 
                                             A)))
  i_adjacent <- unique(queryHits(findOverlaps(suppressWarnings(flank(exons, 
                                                                     1, TRUE)), exons)))
  splice5p[setdiff(i_spliced, i_adjacent)] <- TRUE
  splice3p <- rep(FALSE, length(exons))
  i_spliced <- unique(queryHits(findOverlaps(exons_end, D)))
  i_adjacent <- unique(queryHits(findOverlaps(suppressWarnings(flank(exons, 
                                                                     1, FALSE)), exons)))
  splice3p[setdiff(i_spliced, i_adjacent)] <- TRUE
  mcols(exons)$type <- rep("E", length(exons))
  mcols(exons)$splice5p <- splice5p
  mcols(exons)$splice3p <- splice3p
  mcols(other)$splice5p <- rep(NA, length(other))
  mcols(other)$splice3p <- rep(NA, length(other))
  features <- setNames(c(exons, other), NULL)
  features <- sort(features)
  return(features)
}

annotate2 <- function (query, subject) 
{
  #
  #query <- features
  #subject <- a
  #
  if (!is(subject, "TxFeatures")) {
    stop("subject must be a TxFeatures object")
  }
  if (is(query, "SGFeatures")) {
    query <- annotateFeatures2(query, subject)
  }
  else if (is(query, "SGVariants")) {
    query <- updateObject(query, verbose = TRUE)
    query_class <- class(query)
    query_mcols <- mcols(query)
    query_unlisted <- unlist(query, use.names = FALSE)
    extended <- addDummySpliceSites(query_unlisted)
    extended <- annotate(extended, subject)
    i <- match(featureID(query_unlisted), featureID(extended))
    query_unlisted <- extended[i]
    query <- relist(query_unlisted, query)
    mcols(query) <- query_mcols
    query <- new(query_class, query)
    query <- annotatePaths(query)
  }
  else if (is(query, "Counts")) {
    rd <- rowRanges(query)
    rd <- annotate(rd, subject)
    rowRanges(query) <- rd
  }
  return(query)
}

annotateFeatures2 <- function (query, subject) 
{
  #
  #query <- features
  #subject <- a
  #
  i <- which(type(subject) %in% c("F", "L"))
  if (length(i) > 0) {
    subject <- c(subject[-i], mergeExonsTerminal2(subject[i], 
                                                  1))
  }
  if (is(query, "TxFeatures")) {
    hits <- matchTxFeatures(query, subject)
  }
  else if (is(query, "SGFeatures")) {
    hits <- SGSeq:::matchSGFeatures(query, subject)
  }
  qH <- queryHits(hits)
  sH <- subjectHits(hits)
  for (option in c("tx", "gene")) {
    q_id <- factor(slot(query, "featureID"))
    s_ann <- slot(subject, paste0(option, "Name"))
    id_ann <- SGSeq:::splitCharacterList(s_ann[sH], q_id[qH])
    q_ann <- setNames(id_ann[match(q_id, names(id_ann))], 
                      NULL)
    slot(query, paste0(option, "Name")) <- q_ann
  }
  if (is(query, "SGFeatures")) {
    query2 <- SGSeq:::propagateAnnotation(query)
  }
  return(query)
}

mergeExonsTerminal2 <- function (features, min_n_sample = 1) 
{
  #
  #features <- subject[i]
  #min_n_sample <- 1
  #
  index <- which(type(features) %in% c("F", "L"))
  if (length(index) > 0) {
    features <- features[index]
    splicesite <- SGSeq:::feature2name(features, collapse_terminal = F) #here is where is merged the starts and ends.
    # collapse_terminal = TRUE y ahora es FALSE
    splicesite_n <- table(splicesite)
    i <- which(splicesite %in% names(which(splicesite_n >= 
                                             min_n_sample)))
    features <- features[i]
    splicesite <- splicesite[i]
    splicesite <- factor(splicesite)
    splicesite_i <- split(seq_along(features), splicesite)
    splicesite_w <- split(width(features), splicesite)
    splicesite_i <- mapply(function(i, w) {
      i[which.max(w)]
    }, i = splicesite_i, w = splicesite_w, SIMPLIFY = TRUE)
    exons <- features[splicesite_i]
    for (ann in c("txName", "geneName")) {
      exons_ann <- SGSeq:::splitCharacterList(slot(features, ann), 
                                              splicesite)
      slot(exons, ann) <- setNames(exons_ann, NULL)
    }
  }
  else {
    si <- seqinfo(features)
    exons <- TxFeatures()
    seqinfo(exons) <- si
  }
  return(exons)
}