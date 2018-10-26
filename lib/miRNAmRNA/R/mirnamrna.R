##' gtTable extracts all the targets and their p-values from a GT object.
##'
##' Details follow.
##' 
##' @title Extract gtTable
##' @param object gt-object
##' @return gt-table
##' @author Maarten van Iterson, Sander Bervoets
##' @export
gtTable <- function(object)
{
  # Estrae per ogni covariate la correlazione con la response
  test <- function(set) object@functions$test(set, calculateP=TRUE)

  # size(object):  number of covariates
  leaves <- t(sapply(1:size(object), function(i) test(i)))
            
  ##calculate zscores  
  zsc <- (leaves[,"S"]  - leaves[,"ES"]) / leaves[,"sdS"]
            
  ##association
  positive <- object@functions$positive()
	
  tbl <- data.frame(Pvalue=leaves[,"p"], Association=positive, Weights = weights(object), zscores=zsc)
  
  rownames(tbl) <- object@functions$cov.names(1:size(object))	
  tbl
}

##' toTable generates table from gtrun result list
##'
##' Details follow
##' @title toTable
##' @param results gtrun results 
##' @param method see p.adjust methods
##' @param alpha singificance level
##' @param level local or global
##' @param org either Hs or Mm are implemented yet
##' @return data.frame
##' @author Maarten van Iterson
##' @export
toTable <- function(results, method="BH", alpha=0.05, level=c("global", "local"), org=c("Hs", "Mm", "None"))
  {
    if(org == "Hs")
      require(org.Hs.eg.db)
    else  if(org == "Mm")
      require(org.Mm.eg.db)
    
    mirs <- results$mirs
    targets <- results$targets
    table <- data.frame()
    cnames <- c("mRNA", "miRNA", "Pvalue")
    level <- match.arg(level)
    
    if(length(method)==2)
      {
        method.global <- method[1]
        method.local <- method[2]
      }
    else
      method.global <- method.local <- method

    if(length(alpha) == 2)
      {
        alpha.global <- alpha[1]
        alpha.local <- alpha[2]
      }
    else
      alpha.global <- alpha.local <- alpha
    
    if(!is.data.frame(mirs))
      {
        cnames <- c("miRNA", "mRNA", "Pvalue")
        tmp <- mirs
        mirs <- targets
        targets <- tmp
      }
    mirs <- mirs[p.adjust(mirs$Pvalue, method=method.global) < alpha.global,]
    for(row in rownames(mirs))
      {
        trgts <- targets[[row]]
        if(level=="local")
          {
            if(sum(p.adjust(trgts$Pvalue, method=method.local) < alpha.local) > 0)
              trgts <- trgts[p.adjust(trgts$Pvalue, method=method.local) < alpha.local, ]              
          }
        if(nrow(trgts) == 0)
          next
        
        rows <- cbind(rownames(trgts),
                      row,
                      mirs[rownames(mirs) == row, "Pvalue"],
                      trgts
                      )
        table <- rbind(table, rows)           
      }
    
    colnames(table)[1:3] <- cnames
    rownames(table) <- 1:nrow(table)
    if(org == "Hs")
      table$Symbol <- unlist(mget(as.character(table$mRNA), org.Hs.egSYMBOL))
    else if( org == "Mm")
      table$Symbol <- unlist(mget(as.character(table$mRNA), org.Mm.egSYMBOL))
    table
  }


toTable.TBLAB <- function(results, method="BH", alpha=0.05, level=c("global", "local"))
{
    mirs <- results$mirs
    targets <- results$targets
    table <- data.frame()
    cnames <- c("mRNA", "miRNA", "Pvalue")
    level <- match.arg(level)
    
    if(length(method)==2)
    {
        method.global <- method[1]
        method.local <- method[2]
    }
    else
        method.global <- method.local <- method
    
    if(length(alpha) == 2)
    {
        alpha.global <- alpha[1]
        alpha.local <- alpha[2]
    }
    else
        alpha.global <- alpha.local <- alpha
    
    if(!is.data.frame(mirs))
    {
        cnames <- c("miRNA", "mRNA", "Pvalue")
        tmp <- mirs
        mirs <- targets
        targets <- tmp
    }
    mirs <- mirs[p.adjust(mirs$Pvalue, method=method.global) < alpha.global,]
    for(row in rownames(mirs))
    {
        trgts <- targets[[row]]
        if(level=="local")
        {
            if(sum(p.adjust(trgts$Pvalue, method=method.local) < alpha.local) > 0)
                trgts <- trgts[p.adjust(trgts$Pvalue, method=method.local) < alpha.local, ]              
        }
        if(nrow(trgts) == 0)
            next
        
        rows <- cbind(rownames(trgts),
                      row,
                      mirs[rownames(mirs) == row, "Pvalue"],
                      trgts
        )
        table <- rbind(table, rows)           
    }
    
    colnames(table)[1:3] <- cnames
    rownames(table) <- 1:nrow(table)
    # if(org == "Hs")
    #     table$Symbol <- unlist(mget(as.character(table$mRNA), org.Hs.egSYMBOL))
    # else if( org == "Mm")
    #     table$Symbol <- unlist(mget(as.character(table$mRNA), org.Mm.egSYMBOL))
    table = table[, c(2, 1, 3:7)]
    colnames(table)[c(3,4)] = c("Pvalue_global", "Pvalue_pair")
    table
}


#   for the reversed model
toTabler.TBLAB <- function(results, method="BH", alpha=0.05, level=c("global", "local"))
{
    mrnas <- results$mrnas
    regulators <- results$regulators
    table <- data.frame()
    cnames <- c("miRNA", "mRNA", "Pvalue")
    level <- match.arg(level)
    
    if(length(method)==2)
    {
        method.global <- method[1]
        method.local <- method[2]
    }
    else
        method.global <- method.local <- method
    
    if(length(alpha) == 2)
    {
        alpha.global <- alpha[1]
        alpha.local <- alpha[2]
    }
    else
        alpha.global <- alpha.local <- alpha
    
    if(!is.data.frame(mrnas))
    {
        cnames <- c("mRNA", "miRNA", "Pvalue")
        tmp <- mrnas
        mrnas <- regulators
        regulators <- tmp
    }
    mrnas <- mrnas[p.adjust(mrnas$Pvalue, method=method.global) < alpha.global,]
    for(row in rownames(mrnas))
    {
        regs <- regulators[[row]]
        if(level=="local")
        {
            if(sum(p.adjust(regs$Pvalue, method=method.local) < alpha.local) > 0)
                regs <- regs[p.adjust(regs$Pvalue, method=method.local) < alpha.local, ]              
        }
        if(nrow(regs) == 0)
            next
        
        rows <- cbind(rownames(regs),
                      row,
                      mrnas[rownames(mrnas) == row, "Pvalue"],
                      regs
        )
        table <- rbind(table, rows)           
    }
    
    colnames(table)[1:3] <- cnames
    rownames(table) <- 1:nrow(table)
    colnames(table)[c(3,4)] = c("Pvalue_global", "Pvalue_pair")
    
    # if(org == "Hs")
    #     table$Symbol <- unlist(mget(as.character(table$mRNA), org.Hs.egSYMBOL))
    # else if( org == "Mm")
    #     table$Symbol <- unlist(mget(as.character(table$mRNA), org.Mm.egSYMBOL))
    table
}



##' rungt is a wrapper for the gt() function
##'
##' Details follow.
##' 
##' @title Run the global test on the list of mirs
##' @param mirs list of mirs
##' @param X mRNA expression
##' @param Y miRNA expression
##' @param path path to database
##' @param dbName database name
##' @param tables prediction databases
##' @param numOverlapping number of at least overlapping targets between databases
##' @param top number of significant targets returned -1 is all 
##' @return list of microRNAs and targets
##' @author Maarten van Iterson, Sander Bervoets
##' @export
rungt <- function (mirs, X, Y, path, dbName, tables, numOverlapping, top = -1) 
{
    targets <- list()
    micrornas <- matrix(NA, nrow = length(mirs), ncol = 2)

    #   Riempie 
    #   1.  targets     (lista) con elementi etichettati con mirs[i] che sono le LISTE di mRNA regolati
    #   2.  micrornas   (array) con righe 
    for (i in 1:length(mirs)) {
        
        #   Esistono predizioni per il miRNA mirs[i] ?
        ovl <- overlap(path, dbName, tables = tables, mirs[i], numOverlapping, full = FALSE)
        if (is.null(ovl)) 
            next # No, salta questo miRNA: rimane NA nella riga micrornas[i]
        
        #   Fissato il mirna mirs[i] e ottenute le predizioni ovl
        #   sX  prendiamo soltanto i valori di espressione per questi ultimi.
        #   sy  prendiamo i valori di espressione per IL miRNA mirs[i].
        sX <- X[ which(rownames(X) %in% rownames(ovl)), , drop = FALSE]
        sy <- Y[ which(rownames(Y) %in% mirs[i]), ]
        
        #   sy: response
        #   sX: miRNA targets x samples; t(sX): mRNAs are the covariates
        obj <- gt(sy, t(sX), directional = 1e-06)        

        tbl <- gtTable(obj)
        tbl <- tbl[order(tbl$Association, tbl$Pvalue), ]
        if (top == -1) 
            targets[[mirs[i]]] <- tbl
        else targets[[mirs[i]]] <- tbl[1:max(c(top, nrow(tbl))),]
        
        micrornas[i, 1] <- p.value(obj)     # pvalue
        micrornas[i, 2] <- size(obj)        # numero di target mRNAs regolati
    }
    
    #   E alla fine crea un unico oggetto LISTA con $mirs, $targets
    colnames(micrornas) <- c("Pvalue", "targets")
    rownames(micrornas) <- mirs

    micrornas <- data.frame(na.omit(micrornas))
    micrornas <- micrornas[order(micrornas$Pvalue), ]
    list(mirs = micrornas, targets = targets)

    }


#
#   Modello ORIGINALE:  models miRNA expression as a function of predicted target mRNAs
#
rungt.TBLAB <- function (mirs, X, Y, path, dbName, credentials=NULL, putatives=NULL, debug=NULL) 
{
    targets <- list()
    micrornas <- matrix(NA, nrow = length(mirs), ncol = 2)
    skipped = c();
    
    #   Riempie 
    #   1.  targets     (lista) con elementi etichettati con mirs[i] che sono le LISTE di mRNA regolati
    #   2.  micrornas   (array) con righe 
    
    # con = openDBConnection(credentials, dbName)       # only MySQL
    
    for (i in 1:length(mirs)) {
        
        #   Esistono predizioni per il miRNA mirs[i] ?
        #ovl <- getPutativeTargets.TBLAB.mysql(con, mirs[i])                # uses MySQL
        #ovl <- getPutativeTargets.TBLAB.sqlite3(path, dbName, mirs[i])     # uses sqlite3

        ovl <- getPutativeTargets.TBLAB.PI(PI, mirs[i])                    # uses putatives boolean interaction matrix
        if (is.null(ovl) || length(ovl) == 0) {
          if (!is.null(debug)) {
            cat(sprintf("ovl is NULL!  i=%d, mirna=%s\n", i, mirs[i]))
          }
          next # No, salta questo miRNA: rimane NA nella riga micrornas[i]
        }
        
        #   Fissato il mirna mirs[i] e ottenute le predizioni ovl
        #   sX  prendiamo soltanto i valori di espressione per questi ultimi.
        #   sy  prendiamo i valori di espressione per IL miRNA mirs[i].
        sX <- as.matrix(X[ which(rownames(X) %in% ovl), , drop = FALSE])
        sy <- as.matrix(Y[ which(rownames(Y) %in% mirs[i]), ])

        if (all(sy == 0)) {
            cat(sprintf("-- SKIPPED: %s, all the expression values in the sample were 0.\n", mirs[i]))
            skipped = c(skipped, mirs[i])
            #next
        }
        
        #   Affinché gt funzioni, la riga sy **DEVE** essere un tipo vector "numeric" "doublE"
        sy.names = colnames(sy)     
        sy = as.vector(sy)
        names(sy) = sy.names
        
        #   sy: response
        #   sX: miRNA targets x samples; t(sX): mRNAs are the covariates
        obj <- gt(sy, t(sX), directional = 1e-06)        
        
        tbl <- gtTable(obj)
        tbl <- tbl[order(tbl$Association, tbl$Pvalue), ]

        # if (top == -1) 
        #     targets[[mirs[i]]] <- tbl
        # else targets[[mirs[i]]] <- tbl[1:max(c(top, nrow(tbl))),]
        targets[[mirs[i]]] <- tbl
                
        micrornas[i, 1] <- p.value(obj)     # pvalue
        micrornas[i, 2] <- size(obj)        # numero di target mRNAs regolati
    }
    
    # closeDBConnection(con)    only MySQL
    
    #   E alla fine crea un unico oggetto LISTA con $mirs, $targets
    colnames(micrornas) <- c("Pvalue", "targets")
    rownames(micrornas) <- mirs
    
    micrornas <- data.frame(na.omit(micrornas))
    micrornas <- micrornas[order(micrornas$Pvalue), ]
    list(mirs = micrornas, targets = targets, skipped=skipped)
    
}


#
#   Modello REVERSED:   models mRNA expression as a function of miRNA (predicted regulators) expression
#                       E' tutto SPECULARE rispetto al modello diretto.
#
runrgt.TBLAB <- function (trgs, X, Y, path, dbName, credentials=NULL, putatives=NULL, debug=NULL) 
{
    regulators <- list()
    mrnas <- matrix(NA, nrow = length(trgs), ncol = 2)
    skipped = c();
    
    #   Riempie 
    #   1.  targets     (lista) con elementi etichettati con mirs[i] che sono le LISTE di mRNA regolati
    #   2.  micrornas   (array) con righe 
    
    # con = openDBConnection(credentials, dbName)       # only MySQL
    
    ovl = NULL
    sX = NULL
    sy = NULL
    
    for (i in 1:length(trgs)) {
        
        #   Esistono predizioni per il miRNA mirs[i] ?
        #ovl <- getPutativeTargets.TBLAB.mysql(con, mirs[i])                # uses MySQL
        #ovl <- getPutativeTargets.TBLAB.sqlite3(path, dbName, mirs[i])     # uses sqlite3
        
        ovl <- getPutativeRegulators.TBLAB.PI(PI, trgs[i])                  # uses putatives boolean interaction matrix
        if (is.null(ovl) || length(ovl) == 0) {
            if (!is.null(debug)) {
              cat(sprintf("- mRNA: %s\n", trgs[i]))
              cat(sprintf("ovl is NULL!  i=%d, trgs=%s\n", i, trgs[i]))
            }
            next # No, salta questo miRNA: rimane NA nella riga micrornas[i]
        }

        #   Fissato il mirna mirs[i] e ottenute le predizioni ovl:
        #   sX  prendiamo soltanto i valori di espressione per questi ultimi.
        #   sy  prendiamo i valori di espressione per IL miRNA mirs[i].
        sX <- as.matrix(Y[ which(rownames(Y) %in% ovl), , drop = FALSE])
        sy <- as.matrix(X[ which(rownames(X) %in% trgs[i]), ])
        
        if (all(sy == 0)) {
            cat(sprintf("-- SKIPPED: %s, all the expression values in the sample were 0.\n", trgs[i]))
            skipped = c(skipped, trgs[i])
            next
        }   
        
        #   Affinché gt funzioni, la riga sy **DEVE** essere un tipo vector "numeric" "doublE"
        sy.names = colnames(sy)     
        sy = as.vector(sy)
        names(sy) = sy.names
        
        #   sy: response
        #   sX: miRNA targets x samples; t(sX): mRNAs are the covariates
        obj <- gt(sy, t(sX), directional = 1e-06)        
        
        tbl <- gtTable(obj)
        tbl <- tbl[order(tbl$Association, tbl$Pvalue), ]
        
        # if (top == -1) 
        #     targets[[mirs[i]]] <- tbl
        # else targets[[mirs[i]]] <- tbl[1:max(c(top, nrow(tbl))),]
        regulators[[trgs[i]]] <- tbl
        
        mrnas[i, 1] <- p.value(obj)     # pvalue
        mrnas[i, 2] <- size(obj)        # numero di regulators (miRNA) che regolano questo mrnas
    }
    
    # closeDBConnection(con)    only MySQL
    
    #   E alla fine crea un unico oggetto LISTA con $mirs, $targets
    colnames(mrnas) <- c("Pvalue", "targets")
    rownames(mrnas) <- trgs
    
    mrnas <- data.frame(na.omit(mrnas))
    mrnas <- mrnas[order(mrnas$Pvalue), ]
    list(mrnas = mrnas, regulators = regulators, skipped=skipped)
    
}


##' reversed version rungt is a wrapper for the gt() function
##'
##' Details follow.
##' 
##' @title Run the global test on the list of mirs
##' @param targets list of targets
##' @param Y miRNA expression
##' @param X mRNA expression
##' @param A data.frame containing predicted microRNA target pairs
##' @return list of microRNAs and targets
##' @author Maarten van Iterson
##' @export
runrgt <- function (targets, Y, X, A) 
{
    mirs <- list()
    Targets <- matrix(NA, nrow = length(targets), ncol = 2)
    
    for (i in 1:length(targets)) {

        miRNAs <- subset(A, mRNA == targets[i])$miRNA
        
        sx <- X[which(rownames(X) == targets[i]),]
        sY <- Y[which(rownames(Y) %in% miRNAs), , drop = FALSE]

        if(nrow(sY) == 0)
          next
               
        obj <- gt(sx, t(sY), directional = 1e-06)
        
        tbl <- gtTable(obj)
        tbl <- tbl[order(tbl$Association, tbl$Pvalue), ]
        mirs[[targets[i]]] <- tbl
        Targets[i, 1] <- p.value(obj)
        Targets[i, 2] <- size(obj)
    }
    colnames(Targets) <- c("Pvalue", "mirs")
    rownames(Targets) <- targets
    Targets <- data.frame(na.omit(Targets))
    Targets <- Targets[order(Targets$Pvalue), ]
    list(targets = Targets, mirs = mirs)
}

