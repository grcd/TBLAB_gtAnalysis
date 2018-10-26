
#  ----------------------------------------------------------------------
#
#  CODICE MODIFICATO PER TBLAB 
#
#  ----------------------------------------------------------------------

#  ----------------------------------------------------------------------
#  FUCTION: runrgt.TBLAB
#
#  runrgt.TBLAB è un wrapper per la funzione gt() del package globaltest.
#  E' la versione 'reversed' (r) di rungt.TBLAB e calcola quello che 
#  nell'articolo di Van Iteerson è chiamato "reverse model".
#
#  Nel reverse model, per ogni target (mRNA) si prendono i suoi regolatori 
#  putativi (miRNA) e si fitta il modello.  
#
#  Si ottengono due informazioni :
#
#    1. "Qual è il target più fortemente regolato dai suoi regolatori putativi?"
#       I target vengono ordinati per p-value globale che indica quanto è più
#       significativa la regolazione dall'INSIEME dei regolatori putativi.
#
#    2. "Preso un target, chi sono i regolatori principali ?"
#       Per ogni target abbiamo la lista di miRNA putativi, ordinati per p-value
#       locale che indica quanto è più significativa la regolazione da parte del
#       singolo miRNA.
#
#  NOTE SULL'IMPLEMENTAZIONE
#    Qui manipoliamo il modello.
#    Al livello del modello, si ha sempre:
#
#    Y = response (output dependent variable)
#    X = covaratiates (input indipendent variables)
#
#   Dunque, nel modello reversed, si ha che 
#   
#   Y[i, :]:                 espressione i-esimo mRNA target
#   X[c(x1, x2, ..., x_k)]:  espressione regolatori miRNA putativi per il target
#

runrgt.TBLAB <- function (trgs, X, Y, path, dbName, 
                          credentials=NULL, 
                          putatives=NULL) {
  # Models target mRNA expression as a function of putative miRNA expression
  #
  # Args:
  #   trgs:       mRNA target list 
  #   X:          miRNA expression matrix (rows=miRNAs, cols=samples)
  #   Y:          mRNA expression matrix  (rows=mRNAS, cols=samples)
  #   putatives:  putative interaction (PI) boolean incidence matrix
  #               (rows=miRNA, cols=mRNA)
  #   path:       not used 
  #   dbName:     not used
  #   credentials:not used
  #
  # Returns:
  #   A list object o containing:
  #     o$mrnas:   mRNA targets ordered by global p-values
  #     o$regulators: list obj where each element corresponds to a specific target.
  #                   For example:  o$regulators$`ENSG00000110713`  contains the list
  #                                 of putative regolators (miRNAs) for ENSG00000110713
  #                                 ordered by:
  #                                   1. Association: 0=negative, 1=positive
  #                                   2. p-value
  #     o$skipped: skipped targets (null expression across all samples) 

  regulators <- list()
  mrnas <- matrix(NA, nrow = length(trgs), ncol = 2)
  skipped = c();

  sy = NULL  # espressione target
  ovl = NULL # lista di regolatori putativi per il dato target
  sX = NULL  # sottomatrice espressione miRNA regolatori putativi per il target
  
  for (i in 1:length(trgs)) {
    
    cat(sprintf("runrgt.TBLAB [DEBUG]: mRNA: %s\n", trgs[i]))
    
    #  Esistono predizioni per il miRNA mirs[i] ?
    ovl <- getPutativeRegulators.TBLAB.PI(PI, trgs[i])
    if (is.null(ovl) || length(ovl) == 0) {
      cat(sprintf("runrgt.TBLAB [DEBUG]: ovl is NULL, i=%d, trgs=%s\n", i, trgs[i]))
      next   # No, salta questo mRNA: rimane NA nella riga mrnas[i]
    }
    
    # Fissato il target (mrnas[i]) e ottenuti i regolatori putativi (ovl):
    #   sX  prendiamo soltanto i valori di espressione per i regolatori ovl
    #   sy  prendiamo i valori di espressione per IL target mrnas[i].
    sX <- as.matrix(Y[ which(rownames(Y) %in% ovl), , drop = FALSE])
    sy <- as.matrix(X[ which(rownames(X) %in% trgs[i]), ])
    
    # Se il target ha espressione nulla su tutti i campioni, lo saltiamo lo
    # aggiungiamo alla lisa degli skipped.
    if (all(sy == 0)) {
      cat(sprintf("runrgt.TBLAB [DEBUG]: SKIPPED: %s, expression values across all the samples were 0.\n", trgs[i]))
      skipped = c(skipped, trgs[i])
      next
    }   
    
    # Affinché gt funzioni, la riga sy **DEVE** essere un tipo vector "numeric" "doublE"
    sy.names = colnames(sy)     
    sy = as.vector(sy)
    names(sy) = sy.names
    
    # Calcolo gt
    obj <- gt(sy, t(sX), directional = 1e-06)        
    
    # I risultati di gt per il target mrnas[i] vengono formattati dalla funzione
    # gtTable ed inseriti in tbl. 
    # tbl contiene le associazioni mRNA-miRNA con i p-value locali.
    tbl <- gtTable(obj)
    tbl <- tbl[order(tbl$Association, tbl$Pvalue), ]
    regulators[[trgs[i]]] <- tbl

    # Infine il p-value globale
    mrnas[i, 1] <- p.value(obj)     # pvalue globale 
    mrnas[i, 2] <- size(obj)        # numero di regulatori (miRNA) che regolano questo mrnas

  }
  
  # Alla fine, creo un unico list object con $mrnas, $regulators e $skipped
  colnames(mrnas) <- c("Pvalue", "targets")
  rownames(mrnas) <- trgs
  
  mrnas <- data.frame(na.omit(mrnas))
  mrnas <- mrnas[order(mrnas$Pvalue), ]
  list(mrnas = mrnas, regulators = regulators, skipped=skipped)
  
}


#  ----------------------------------------------------------------------
#  FUCTION: rungt.TBLAB
#
#  rungt.TBLAB è un wrapper per la funzione gt() del package globaltest.
#  E' la versione 'direct' di rungt.TBLAB e calcola quello che 
#  nell'articolo di Van Iteerson è chiamato "direct model".
#
#  E' esattamente la versione speculare di runrgt.TBLAB, fornita per completezza.

rungt.TBLAB <- function (mirs, X, Y, path, dbName, 
                         credentials=NULL, 
                         putatives=NULL) {
  # Models miRNA expression as a function of putatives mRNA expression 
  #
  # Args:
  #   mirs:       list of miRNAs 
  #   X:          mRNA expression matrix (rows=mRNAs, cols=samples)
  #   Y:          miRNA expression matrix  (rows=miRNAs, cols=samples)
  #   putatives:  putative interaction (PI) boolean incidence matrix
  #               (rows=miRNA, cols=mRNA)
  #   path:       not used 
  #   dbName:     not used
  #   credentials:not used
  #
  # Returns:
  #   A list object o containing:
  #     o$mirs:    miRNAs ordered by global p-values
  #     o$targets: list obj where each element corresponds to a specific mRNAs.
  #                For example:  o$targets$`ENSG00000110713`  contains the list
  #                              of putative regolators (miRNAs) ordered by 
  #                                1. Association: 0=negative, 1=positive
  #                                2. p-value
  #     o$skipped: skipped targets (null expression across all samples) 
  
  targets <- list()
  micrornas <- matrix(NA, nrow = length(mirs), ncol = 2)
  skipped = c();
  
  for (i in 1:length(mirs)) {
    
    #   Esistono predizioni per il miRNA mirs[i] ?
    ovl <- getPutativeTargets.TBLAB.PI(PI, mirs[i])                     # uses putatives boolean interaction matrix
    if (is.null(ovl) || length(ovl) == 0) 
      next # No, salta questo miRNA: rimane NA nella riga micrornas[i]
    
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
  
  #   E alla fine crea un unico oggetto LISTA con $mirs, $targets
  colnames(micrornas) <- c("Pvalue", "targets")
  rownames(micrornas) <- mirs
  
  micrornas <- data.frame(na.omit(micrornas))
  micrornas <- micrornas[order(micrornas$Pvalue), ]
  list(mirs = micrornas, targets = targets, skipped=skipped)
  
}


#  ----------------------------------------------------------------------
#  FUCTION: toTabler.TBLAB
#
#  Prende in input la tabella di risultati raw prodotti da runrgt.TABLE e
#  aggiusta i p-value usando FDR, secondo la soglia fornita.
#
#  Restituisce una tabella in cui ogni riga è una associazione miRNA-mRNA
#  con il corrispondente p-value globale e locale.
#

toTabler.TBLAB <- function(results, method="BH", alpha=0.05, level=c("global", "local")) {
  # Adjusts runrgt.TBLAB p-values for a given FDR threshold.
  #
  # Args:
  #   results:  output from runrgt.TBLAB
  #   method:   FDR adjusting method
  #   alpha:    FDR threshold
  #   level:    at which level(s) apply the FDR correction?
  #
  # Returns:
  #   A table containing miRNA-mRNA pairs with corresponding global and local
  #   adjusted p-values.
  
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


#  ----------------------------------------------------------------------
#  FUCTION: toTable.TBLAB
#
#  versione speculare a toTabler.TBLAB, funziona con risultati dati da
#  rungt.TBLAB (direct model).

toTable.TBLAB <- function(results, method="BH", alpha=0.05, level=c("global", "local"))
{
  # Adjusts rungt.TBLAB p-values for a given FDR threshold.
  #
  # Args:
  #   results:  output from rungt.TBLAB
  #   method:   FDR adjusting method
  #   alpha:    FDR threshold
  #   level:    at which level(s) apply the FDR correction?
  #
  # Returns:
  #   A table containing miRNA-mRNA pairs with corresponding global and local
  #   adjusted p-values.
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


#  ----------------------------------------------------------------------
#  CODICE ORIGINALE (Van Iterson)
#  ----------------------------------------------------------------------

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
