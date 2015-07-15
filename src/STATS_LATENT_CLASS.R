#Licensed Materials - Property of IBM
#IBM SPSS Products: Statistics General
#(c) Copyright IBM Corp. 2011, 2014
#US Government Users Restricted Rights - Use, duplication or disclosure 
#restricted by GSA ADP Schedule Contract with IBM Corp.

# Author: JKP, IBM SPSS
# history
# 10-May-2011 - original version
# 28-Jun-2012 - added gettext calls
# 11-Feb-2013 - add support for scoring the input cases
# 01-Nov-2013 - adjustments for plot changes

helptext="The STATS LATENT CLASS command requires the R Integration Plug-in
and the R poLCA package.

STATS LATENT CLASS  MANIFESTVARS=variable list 
[COVARIATES=variable list] CLASSES=n 
[/OPTIONS [REPETITIONS=m]  [MAXITER=value] [TOL=value]
[/OUTPUT [GRAPHS={YES*|NO}]
[/SAVE PREDCELLDS=datasetname] [SCOREDS=datasetname IDVAR=varname]
[MISSING={INCLUDE*|LISTWISE}]].
[/HELP]

Example:
STATS LATENT CLASS MANIFESTVARS = a b c CLASSES=3
/OPTIONS REPETITIONS = 5.

MANIFESTVARS and COVARIATES list the categorical and covariate variables.
The categorical variables must have integer values starting at 1.
If covariates are listed, latent class regression models are estimated.  MANIFESTVARS
is required, but COVARIATES is optional.

CLASSES is the number of latent classes

In order to minimize the chance of finding only a local maximum of the likelihood, the
procedure can reestimate at different random starting values.  REPETITIONS, default 10,
specifies how many sets to try
MAXITER, default 1000, specifies the maximum number of iterations
TOL specifies the convergence criterion.

PREDCELLS can specify a new dataset to be created holding the classification probabilities -
observed and expected

SCOREDS can specify a new dataset to be created holding the scores for the model.
If SCORES is specified, IDVAR is required and will be included in the dataset in order
to facilitate merging.
The scores dataset has variables, the ID, MaxProb, and MaxClass followed by
probability variables for each class.

MISSING=INCLUDE uses all available nonmissing data.  LISTWISE does listwise deletion of
cases.  Missing covariate values always cause the case to be deleted.

STATS LATENT CLASS /HELP prints this information and does nothing else.

"

latent <- function(manifestvars, covariates = NULL, classes = 2, repetitions = 10, predcellds = NULL,
missingopt = "include", graphs = TRUE, maxiter = 1000, tol =1e-8, scoreds=NULL, idvar=NULL, ignore=NULL) {

    domain<-"STATS_LATENT_CLASS"
    setuplocalization(domain)
    
    gtxt <- function(...) {
        return(gettext(...,domain=domain))
    }

    gtxtf <- function(...) {
        return(gettextf(...,domain=domain))
    }

    tryCatch(library(poLCA), error=function(e) {
        stop(gtxt("The R poLCA package is required but could not be loaded."),call. = FALSE)
    })

    verbose = FALSE  #suppress print output
        # null items get ignored in the input list - covariates and idvar could be NULL
    nummanifest = length(manifestvars)
    numcov = length(covariates)
    dta<-spssdata.GetDataFromSPSS(c(manifestvars, covariates, idvar), missingValueToNA = TRUE)
    numcases = nrow(dta)   # before pruning missing
    if (missingopt != "include") {
        dta = dta[complete.cases(dta), ]
    } else if (numcov > 0) {
        dta = dta[apply(!is.na(dta[(nummanifest+1) : (nummanifest+numcov)]), 1, all),]
    }
  
    if (nrow(dta) == 0)
        stop(gtxt("There are no cases to analyze"), call.=FALSE)
    if (nummanifest <2)
        stop(gtxt("Too few manifest variables specified"), call. = FALSE)
    if (!is.null(scoreds) && is.null(idvar))
        stop(gtxt("An ID variable is required when saving scores"), .call = FALSE)
    StartProcedure(gtxt("Latent Class Analysis"), "STATS LATENT CLASS") 

    if (missingopt == "include") {
        missingopt = FALSE
    } else {
        missingopt = TRUE
    }
    # build formula for model
    pp=paste(manifestvars, collapse = ",")
    if (is.null(covariates)) 
    cov <- "1"
    else
    cov <- paste(covariates, collapse = "+")
    f=paste(paste("cbind(", pp, sep = ""), ")~", sep = "")
    ff = paste(f, cov, sep="")
    res = poLCA(as.formula(ff), dta, nclass=classes, maxiter=maxiter, graphs=FALSE, tol=tol, na.rm=missingopt,
                nrep=repetitions, verbose=verbose)
    if (!is.null(res)) {
        # fit statistics
        lbls = c(gtxt("Number of Cases"), gtxt("Number of Complete Cases"), gtxt("Number of Parameters Estimated"),
            gtxt("Residual D.F"), gtxt("Maximum Log-Likelihood"), 
            gtxtf("AIC(%s)", classes), 
            gtxtf("BIC(%s)", classes), 
            gtxtf("LR/Deviance(%s)", classes), 
            gtxtf("Chi-squared(%s)", classes),
            gtxt("Number of repetitions"))

        vals = c(res$N, res$Nobs, res$npar,
            res$resid.df, res$llik,
            res$aic, res$bic, res$Gsq, res$Chisq, repetitions)
        fits=data.frame(vals, row.names = lbls)
        spsspivottable.Display(fits, title = gtxtf("Fit Statistics for %s Latent Classes", classes),
            collabels = c(gtxt("Statistic")),
            templateName = "LCAFITSTATISTICS",
            outline = gtxt("LCA Fit Statistics"))

        # estimated class conditional probabilities
        maxval = max(dta[,1:nummanifest], na.rm = TRUE)
        probtbl = spss.BasePivotTable(gtxt("Estimated Class Conditional Probabilities"),
            "LCAESTPROBS",)
        coldim = BasePivotTable.Append(probtbl, Dimension.Place.column, gtxt("Probabilities"))
        outerrowdim = BasePivotTable.Append(probtbl, Dimension.Place.row, gtxt("Variable"))
        innerrowdim = BasePivotTable.Append(probtbl, Dimension.Place.row, gtxt("Class"))
        collist = lapply(as.character(1:maxval), spss.CellText.String)
        BasePivotTable.SetCategories(probtbl, coldim, collist)

        varlist = lapply(manifestvars, spss.CellText.String)
        BasePivotTable.SetCategories(probtbl, outerrowdim, varlist)

        classlist = lapply(as.character(1:classes), spss.CellText.String)
        BasePivotTable.SetCategories(probtbl, innerrowdim, classlist)

        for (v in 1:nummanifest) {
            vmax = max(dta[,v], na.rm = TRUE)
            for (c in 1:classes) {
                thisrow = lapply(res$probs[[v]][c,], spss.CellText.Number)
                rowlen = length(thisrow)
                if (rowlen < maxval) {
                    thisrow = c(thisrow, lapply(rep("NA", maxval - rowlen), spss.CellText.String))
                }
                BasePivotTable.SetCellsByRow(probtbl, list(varlist[[v]], classlist[[c]]), thisrow)
            }
        }


        if (graphs) {
          tryCatch(plot(res), 
            error = function(e) {
              tryCatch((exists("plot.poLCA") && poLCA:plot.poLCA(res)) || poLCA:plot.poLCA(res),
                error=function(e) {}, warning = function(w) {})
            }
          )
            # The plot function is not exported prior to R 2.12
            #tryCatch((exists("plot.poLCA") && poLCA:plot.poLCA(res)) || poLCA:plot.poLCA(res),
            #    error=function(e) {}, warning = function(w) {})
        }

        # regressions
        if (!is.null(covariates)) {
            for (c in 1:(classes - 1)) {
                tstat = res$coeff[,c] / res$coeff.se[,c]
                sig = 2*(1 - pt(abs(tstat), res$resid.df))
                df = data.frame(res$coeff[,c], res$coeff.se[,c], tstat, sig)
                spsspivottable.Display(df, title=gtxtf("Estimates for %s Latent Classes: %s / 1", classes, c+1),
                    rowdim=gtxt("Variables"), collabels=c(gtxt("Coefficient"), gtxt("Std. Error"), gtxt("T Statistic"), gtxt("Sig.")),
                    templateName="LCAREGRESSION",
                    outline = gtxtf("Estimates for Latent Classes - %s", c+1))
            }
            # plot each covariate effect

            # extract covariates - need data.frame wrapper in case just one column
            covs = data.frame(dta[,c((nummanifest + 1):(nummanifest + numcov))])
            means = colMeans(covs, na.rm = TRUE)
            for (cv in 1:numcov) {
                vmax = max(covs[,cv], na.rm = TRUE)
                vmin = min(covs[,cv], na.rm = TRUE)
                m = matrix(nrow = vmax - vmin + 1, ncol = numcov + 1)
                m[,1] = 1
                for (i in 1:numcov) {
                    if (i == cv)
                    m[,i + 1] = c(vmin:vmax)
                    else
                    m[,i + 1] = means[[i]]
                }

                expb = exp(m %*% res$coeff)
                matplot(c(vmin:vmax), (cbind(1, expb) / (1+rowSums(expb))), type = "l", lwd = 3,
                    xlab = covariates[[cv]],
                    ylim = c(0:1),  ylab = gtxt("Probability of Latent Class Membership"))
            }
        }

        # latent class sizes
        spsspivottable.Display(data.frame(res$P), title = gtxt("Latent Class Proportions"),
            rowdim=gtxt("Class"), collabels = c(gtxt("Proportion")),
            templateName = "LCACLASSSIZES",
            outline = gtxt("Latent Class Proportions"))
    }
    spsspkg.EndProcedure()
        #Save predcells as a dataset
    if (!is.null(predcellds) && !is.null(res)) { 
        vnames = list()
        numv = length(res$predcell)
        i = 0
        for (v in names(res$predcell)) {
            i = i+1
            if (i < numv) {
                vnames[[i]] = c(v, "", 0, "F8.0", "nominal")
            }
            else {
                vnames[[i]] = c(v, "", 0, "F8.4", "scale")
            }
        }
        vdict <-spssdictionary.CreateSPSSDictionary(vnames)
        tryCatch({
        spssdictionary.SetDictionaryToSPSS(predcellds, vdict)
        spssdata.SetDataToSPSS(predcellds, res$predcell)
        }, error=function(e) {print(sprintf(gtxt("Could not create dataset: %s.  The dataset name must not already be in use"),
            predcellds))})
        if (is.null(scoreds)) 
        spssdictionary.EndDataStep()
    }

    if (!is.null(res) && !is.null(scoreds)) {
        f <- function(v) {
            # compute max, and index of max
            return(c(max(v), which.max(v)))
        }
        #                                                    ID             MaxProb, MaxClass           all scores
        scoredf = data.frame(dta[, length(dta)], t(apply(res$posterior, 1, f)), res$posterior, stringsAsFactors=FALSE)
        classnames = paste("Class", 1:classes,sep="_")  
        names(scoredf) = c("ID", "MaxProb", "MaxClass", classnames)

      # create SPSS dataset
      vardict = spssdictionary.GetDictionaryFromSPSS()
      idvartype = vardict[3, match(idvar, vardict[1,])]
      if (idvartype == 0) {
          idvarformat = "F10.0"
          } else {
              idvarformat = paste("A", idvartype, sep="")
          }
      vnames = list(c("ID", gtxt("Case Id"), idvartype, idvarformat, "nominal"), c("MaxProb", gtxt("Maximum Probability"), 0, "F8.5", "scale"),
                    c("MaxClass", gtxt("Highest Probability Class"), 0, "F6.0", "nominal"))
      i = 0
      for (v in classnames) {  # class names are just numbers
        i = i+1
        vnames[[3 + i]] = c(v, gtxt("Class Probability"), 0, "F8.4", "scale")
      }
      vnames[[1]][1] = makeidname(idvar, vnames)  # ensure ID name is unique
      vdict <-spssdictionary.CreateSPSSDictionary(vnames)
      tryCatch({
      spssdictionary.SetDictionaryToSPSS(scoreds, vdict)
      spssdata.SetDataToSPSS(scoreds, scoredf)},
      error=function(e) {print(sprintf(gtxt("Could not create dataset: %s.  The dataset name must not already be in use"),
        scoreds))})
      spssdictionary.EndDataStep()
  }

    # clean up workspace
    res <- tryCatch(rm(list=ls()), warning = function(e) {return(NULL)})
}

# override for api to account for extra parameter in V19 and beyond
StartProcedure <- function(procname, omsid) {
    if (substr(spsspkg.GetSPSSVersion(),1, 2) >= 19) {
        spsspkg.StartProcedure(procname, omsid)
    }
    else {
        spsspkg.StartProcedure(omsid)
    }
}

makeidname = function(idvar, vnames) {
    # Make a unique name for id variable
    
    # extract the variable name from list of proto-dictionary elements.  ID standin assumed to be first
    inuse = tolower(sapply(vnames, "[[", 1) [-1])  # extract the variable names from list of proto-dictionary elements
    idvarlc = tolower(idvar)
    if (!(idvarlc %in% inuse)) {
        return(idvar)
    } else {
        return(paste(idvar,"", sep="_"))
    }
}

    
# a functionally similar function is exported in R2.12 but not in 2.10 :-(
poLCAplot <- function(x) {
    # plot LCA results
    colkts <- sapply(x$probs, ncol)
    classes <- length(x$P)
    if (max(colkts) == 2) {
        poLCA.makeplot.dich(x$probs, x$P, x$y, NULL)
        } else {
            layout(matrix(seq(1, (classes + 1)), classes + 1,1),heights = c(rep(5, classes),1))
            for (cl in 1:classes) {
                sharelabel = paste(gtxt("Class "), cl, gtxt(": population share = "), round(x$P[cl], 3), sep = "")
                poLCA.makeplot.poly(x$probs, cl, x$y, colkts, sharelabel)
            }
        }
    par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)
}
setuplocalization = function(domain) {
    # find and bind translation file names
    # domain is the root name of the extension command .R file, e.g., "SPSSINC_BREUSCH_PAGAN"
    # This would be bound to root location/SPSSINC_BREUSCH_PAGAN/lang

    fpath = Find(file.exists, file.path(.libPaths(), paste(domain, ".R", sep="")))
    bindtextdomain(domain, file.path(dirname(fpath), domain, "lang"))
} 

Run<-function(args){
    
    cmdname = args[[1]]
    args <- args[[2]]
    oobj<-spsspkg.Syntax(templ=list(
        spsspkg.Template("MANIFESTVARS", subc="",  ktype="existingvarlist", var="manifestvars", islist=TRUE),
        spsspkg.Template("COVARIATES", subc="",  ktype="existingvarlist", var="covariates", islist=TRUE),
        spsspkg.Template("CLASSES", subc="",  ktype="int", var="classes", vallist=list(1)),
        spsspkg.Template("REPETITIONS", subc="OPTIONS",  ktype="int", var="repetitions", vallist=list(1), islist=FALSE),
        spsspkg.Template("MISSING", subc="OPTIONS", ktype="str", var="missingopt", 
                         vallist=list("include", "listwise")),
        spsspkg.Template("MAXITER", subc="OPTIONS", ktype="int", var="maxiter", vallist=list(1)),
        spsspkg.Template("TOL", subc="OPTIONS", ktype="float", var="tol", vallist=list(0.)),
        spsspkg.Template("GRAPHS", subc="OUTPUT", ktype="bool", var="graphs", islist=FALSE),
        spsspkg.Template("PREDCELLDS", subc="SAVE", ktype="literal", var="predcellds", islist=FALSE),
        spsspkg.Template("SCOREDS", subc="SAVE", ktype="literal", var="scoreds", islist=FALSE),
        spsspkg.Template("IDVAR", subc="SAVE", ktype="existingvarlist", var="idvar", islist=FALSE),
        spsspkg.Template("FILLER", subc="SAVE", ktype="str", var="ignore", islist=FALSE)
    ))        
    if ("HELP" %in% attr(args,"names")) {
        #writeLines(helptext)
        helper(cmdname)
    } else {
        res <- spsspkg.processcmd(oobj,args,"latent")
    }
}

helper = function(cmdname) {
    # find the html help file and display in the default browser
    # cmdname may have blanks that need to be converted to _ to match the file
    
    fn = gsub(" ", "_", cmdname, fixed=TRUE)
    thefile = Find(file.exists, file.path(.libPaths(), fn, "markdown.html"))
    if (is.null(thefile)) {
        print("Help file not found")
    } else {
        browseURL(paste("file://", thefile, sep=""))
    }
}
if (exists("spsspkg.helper")) {
assign("helper", spsspkg.helper)
}