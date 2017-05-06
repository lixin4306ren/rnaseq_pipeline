ea.browse <- function(res, nr.show=-1, graph.view=NULL, html.only=FALSE)
{
    method <- ifelse( is(res$method, "character"), res$method, NA)
    eset <- res$eset
    alpha <- res$alpha
    gs <- res$gs

    # create out dir
    out.dir <- config.ebrowser("OUTDIR.DEFAULT")
    if(!file.exists(out.dir)) dir.create(out.dir)
    rep.dir <- file.path(out.dir, "reports")
    
    # how many gene sets to show in the output?
    if(nr.show < 1) nr.show <- res$nr.sigs
    if(nr.show > nrow(res$res.tbl)) nr.show <- nrow(res$res.tbl)

    # add description & nr.genes per gene set
    res <- res$res.tbl[seq_len(nr.show),]

    # expecting transition from gene set lists to collections in the near future
    # gsc <- res$gsc
    gsc <- gs.list.2.gs.coll(gs[res[,1]])
    res[,1] <- sapply(res[,1], function(s) unlist(strsplit(s, "_"))[1])
    
    is.kegg <- is(collectionType(gsc[[1]]), "KEGGCollection")
    is.go <- is(collectionType(gsc[[1]]), "GOCollection")

    gs.title <- sapply(gsc, description)
    nr.genes <- sapply(gsc, function(g) length(geneIds(g)))

    cnames <- c(colnames(res)[1], "TITLE")
    resn <- DataFrame(res[,1], gs.title)
    if(!("NR.GENES" %in% colnames(res)))
    {
        cnames <- c(cnames, "NR.GENES")
        resn <- DataFrame(resn, nr.genes)
    }
    cnames <- c(cnames, colnames(res)[2:ncol(res)])
    resn <- DataFrame(resn, res[,2:ncol(res)]) 
    colnames(resn) <- cnames
    res <- resn

    # make gene pages
    # TODO: ensure in sbea and nbea that we are running only
    # on intersecting genes between gs and eset
    im <- incidence(gsc)
    org <- organism(gsc[[1]])
    if(org == "") org <- annotation(eset)
    if(!length(org)) stop("Organism annotation not found!\n", 
        "Organism under study must be annotated via annotation(eset)")

    message("Creating gene report ...")
    eset <- eset[colnames(im),]
    fDat <- fData(eset)[,sapply(c("FC.COL", "ADJP.COL"), config.ebrowser)]
    gt <- suppressMessages(gene.table(im, org, fcs=fDat))
    gn.cols <- sapply(c("SYM.COL", "GN.COL"), config.ebrowser)
    fData(eset)[,gn.cols] <- gt[,gn.cols] 
    gt.reps <- sapply(gsc, function(s) gene.report(s, gt, out.dir)) 
    
    # make gene set page
    # (1) link gene report
    link <- paste0(out.dir,"/reports/",names(gsc), ".html")
    res[,"NR.GENES"] <- hwrite(res[,"NR.GENES"], link=link, table = FALSE)
   
    # set view: volcano, heatmap 
    message("Creating set view ...")
    out.prefix <- file.path(rep.dir, names(gsc))
    names(out.prefix) <- names(gsc)
    vcol <- sapply(gsc, function(s) 
        view.set(eset[geneIds(s),], out.prefix[setName(s)]))
    vcol<-paste0(out.dir,"/reports/",vcol)
    vcol <- hwriteImage(sub("sview.html", "volc.png", vcol),
        link=vcol, table = FALSE, height=50, width=50, target="_blank")
    res <- DataFrame(res, vcol)
    colnames(res)[ncol(res)] <- "SET.VIEW" 

    # path view: kegg maps
    if(is.kegg)
    { 
        message("Creating kegg view ...")
        vcol <- sapply(gsc, function(s) 
            view.path(setName(s), eset[geneIds(s),], out.prefix[setName(s)]))
        vcol<-paste0(out.dir,"/reports/",vcol)
        vcol <- hwriteImage(sub("kview.html", "kpath.png", vcol),
            link=vcol, table = FALSE, height=50, width=50, target="_blank")
        res <- DataFrame(res, vcol)
        colnames(res)[ncol(res)] <- "PATH.VIEW"
    }

    # graph view: ggea graphs
    if(!is.null(graph.view)) 
    {
        message("Creating graph view ...")
        vcol <- sapply(gsc, function(s) view.graph(eset[geneIds(s),], query.grn(
            geneIds(s), graph.view, index=FALSE), alpha, out.prefix[setName(s)]))
        vcol <- hwriteImage(sub("html$", "png", vcol),
            link=vcol, table = FALSE, height=50, width=50, target="_blank")
        res <- DataFrame(res, vcol)
        colnames(res)[ncol(res)] <- "GRAPH.VIEW"
    }
 
    # (2) link KEGG / GO 
    link <- NULL
    GS.COL <- config.ebrowser("GS.COL")
    if(is.kegg) link <- sapply(gsc, function(s) 
        get.html.of.marked.pathway(setName(s), geneIds(s)[fDat[geneIds(s),2] < alpha]))
    else if(is.go) link <- paste0(config.ebrowser("GO.SHOW.URL"), res[,GS.COL])
    if(!is.null(link)) res[,GS.COL] <- 
        hwrite(res[,GS.COL], link=link, table=FALSE)
    if(is.kegg) tmp_post<-"kegg"
    else if(is.go) tmp_post<-"go"
    htmlRep <- HTMLReport(shortName=paste(tmp_post,method,sep="_"),
        title=paste(toupper(method), config.ebrowser("RESULT.TITLE"), sep=" - "),
        basePath="./", reportDirectory="")
    res <- as.data.frame(res)
    publish(res, htmlRep) 
    rep <- finish(htmlRep)
    if(!html.only) if (interactive()) browseURL(rep)
    return(res)
}

gs.list.2.gs.coll <- function(gs.list)
{
    gs.type <- auto.detect.gs.type(names(gs.list)[1])
    ctype <- paste0(gs.type, "Collection")
    gs.list.new <- sapply(names(gs.list), 
        function(s){
                spl <- unlist(strsplit(s,"_")) 
                args <- list()
                if(gs.type != "Computed") args <- list(spl[1])
                descr <- ""
                if(length(spl) > 1) 
                    descr <- paste(spl[2:length(spl)],collapse=" ")
                org <- ifelse(gs.type == "KEGG", substring(spl[1], 1, 3), "")
                
                sname <- gsub("[<>:\\?\\|\"\\*\\/]", "", spl[1]) 

                gset <- GeneSet(setName=sname, 
                        geneIds=gs.list[[s]],
                        collectionType=do.call(ctype, args),
                        shortDescription=descr,
                        organism=org)
                return(gset)
        })
    gsc <- GeneSetCollection(gs.list.new) 
    return(gsc)
}

auto.detect.gs.type <- function(gs.id)
{
    if(substring(gs.id, 1, 3) == "GO:") return("GO")
    else if(grepl("^[a-z]{3}[0-9]{5}", gs.id)) return("KEGG")
    else return("Computed") 
}

gene.table <- function(im, org, fcs=NULL, grn=NULL)#, context="")
{
    # load org pkg
    org.pkg <- .org2pkg(org)
    .isAvailable(org.pkg)
    org.pkg <- get(org.pkg)

    # (1) gene identification 
    EZ.COL <- config.ebrowser("EZ.COL")
    gn.cols <- sapply(c("SYM.COL", "GN.COL"), config.ebrowser)
    gt <- select(org.pkg, keys=colnames(im), columns=gn.cols, keytype=EZ.COL)

    # (2) fold change
    if(!is.null(fcs))
    {
        fcs[,1] <- round(fcs[,1], digits=2)
        fcs[,2] <- signif(fcs[,2], digits=2)
        gt <- cbind(gt, fcs) 
    }

    # (3) interactions
    if(!is.null(grn))
    {
        ias.per.gene <- sapply(colnames(im), 
            function(gene) grn[grn[,1] == gene | grn[,2] == gene,,drop=FALSE])
        nr.ias.per.gene <- sapply(ias.per.gene, nrow)
        gt <- cbind(gt, nr.ias.per.gene)
        colnames(gt)[ncol(gt)] <- "INTACTS"
    }

    # (4) nr.sets
#    gene.occ.freq <- colSums(im)
#    gt <- cbind(gt, gene.occ.freq)
#    colnames(gt)[ncol(gt)] <- "SETS"
#
#    # (5) pubmed
#    pmids <- mapIds(org.pkg, keys=colnames(im), 
#        column="PMID", keytype=EZ.COL, multiVals="list")
#
#    # context search ?
#    # if(context != "") pmids <- search.abstracts(pmids, context=context)
#  
#    nr.articles <- sapply(pmids, length)
#    pubmed.urls <- sapply(pmids, function(p) paste0(
#        config.ebrowser("PUBMED.URL"), paste(p, collapse=",")), USE.NAMES=FALSE)
#    articles <- paste(nr.articles, pubmed.urls)
#
#    gt <- cbind(gt, articles)
#    colnames(gt)[ncol(gt)] <- config.ebrowser("PMID.COL")
    return(gt)
}

.org2pkg <- function(org, type=c("OrgDb", "TxDb", "BSgenome"))
{
    type <- match.arg(type)

    SPECIES <- rbind(
        c("anopheles", "Anopheles gambiae", "Ag", "aga", "anoGam", "7165"),
        c("arabidopsis", "Arabidopsis thaliana", "At", "ath", NA, "3702"),
        c("bovine", "Bos taurus", "Bt", "bta", "bosTau", "9913"),
        c("canine", "Canis familiaris", "Cf", "cfa", "canFam", "9615"),
        c("chicken", "Gallus gallus", "Gg", "gga", "galGal", "9031"), 
        c("chimp", "Pan troglodytes", "Pt", "ptr", "PanTro", "9598"),
        c("ecoliK12", "Escherichia coli K12", "EcK12", "eco", NA, "562,83333,511145"), 
        c("ecoliSakai", "Escherichia coli Sakai", "EcSakai", "ecs", NA, "83334"),
        c("fly", "Drosophila melanogaster", "Dm", "dme", "dm", "7227"),
        c("human", "Homo sapiens", "Hs", "hsa", "hg", "9606"),
        c("malaria", "Plasmodium falciparum", "Pf", "pfa", NA, "5833"),
        c("mouse", "Mus musculus", "Mm", "mmu", "mm", "10090"),
        c("pig", "Sus scrofa", "Ss", "ssc", "susScr", "9823"),
        c("rat", "Rattus norvegicus", "Rn", "rno", "rn", "10116"), 
        c("rhesus", "Macaca mulatta", "Mmu", "mcc", "rheMac", "9544"),  
        c("worm", "Caenorhabditis elegans", "Ce", "cel", "ce", "6239"),
        c("xenopus", "Xenopus laevis", "Xl", "xla", "NA", "8355"),
        c("yeast", "Saccharomyces cerevisiae", "Sc", "sce", "sacCer", "4932,559292"),
        c("zebrafish", "Danio rerio", "Dr", "dre", "danRer", "7955")
    )
    colnames(SPECIES) <- c("common", "tax", "bioc", "kegg", "ucsc", "ncbi")
    

    # org specification via 
    # (a) 3-letter code, e.g. 'hsa' 
    # (b) genome assembly, e.g. 'hg38'
    is.genome <- sub("[0-9]+$", "", org) %in% SPECIES[,"ucsc"]
    if(is.genome)
    {
        ucsc.id <- org
        i <- grep(sub("[0-9]+$", "", org), SPECIES[,"ucsc"]) 
        bioc.id <- SPECIES[i, "bioc"]
    }
    else
    {
        ind <- apply(SPECIES, 1, function(r) org %in% r)
        if(any(ind)) i <- which(ind)[1]
        else stop(paste0("unrecognized organism ID \'", org, "\'"))
        bioc.id <- SPECIES[i, "bioc"]
        ucsc.id <- SPECIES[i, "ucsc"]
    }

    # TxDB, BSgenome, or OrgDB package?
    if(type %in% c("TxDb", "BSgenome"))
    {
        pkg.string <- paste0("^", type, ".", bioc.id, "[a-z]+.UCSC.", ucsc.id)
        pkg <- grep(pkg.string, .availableOrgPkgs(type), value=TRUE)
        if(length(pkg) == 0)
            pkg <- grep(pkg.string, .availableOrgPkgs(type, local=FALSE), value=TRUE)
        if(length(pkg) == 0)
            stop(paste("No corresponding", type, "package for", org))
        else if(length(pkg) > 1)
        {
            message("Found several genome assemblies")
            sapply(pkg, function(p) 
                message(paste(which(pkg==p), p, sep=": ")))
            n <- readline(paste0("Choose assembly (1-", length(pkg),") : "))
            pkg <- pkg[as.integer(n)]

            #message("Found several genome assemblies")
            #message(paste("Using latest:", pkg))
            #ver <- sapply(pkg, 
            #    function(p)
            #    {
            #        spl <- unlist(strsplit(p, "\\."))
            #        ind <- length(spl)
            #        if(type == "TxDb") ind <- ind - 1
            #        ass <- spl[ind]
            #        ver <- sub("^[a-zA-Z]+", "", ass)
            #        return(as.integer(ver))
            #    })
            #pkg <- pkg[which.max(ver)]
        }
    }
    else
    {
        id.type <- .getOrgIdType(bioc.id)
        pkg <- paste("org", bioc.id, id.type, "db", sep=".")
    }
    return(pkg)
}

.getOrgIdType <- function(org)
{
    it <- "eg"
    if(org == "At") it <- "tair"
    else if(org == "Pf") it <- "plasmo"
    else if(org == "Sc") it <- "sgd"
    return(it)
}   

.availableOrgPkgs <- function(type=c("OrgDb", "TxDb", "BSgenome"), local=TRUE)
{
    if(local) pkgs <- .packages(all.available=TRUE)
    else pkgs <- available.packages(paste0("http://bioconductor.org/",
        "packages/release/data/annotation/src/contrib"))[, "Package"]
    
    type <- match.arg(type)
    org.string <- "^org.[A-z][a-z]+.[a-z]+.db$"
    if(type == "TxDb") 
        org.string <- "^TxDb.[A-Z][a-z]+.UCSC.[a-z]{2}[A-Za-z]*[0-9]{1,3}.[a-z]{3,5}Gene$"
    else if(type == "BSgenome") 
        org.string <- "^BSgenome.[A-Z][a-z]+.UCSC.[a-z]{2}[A-Za-z]*[0-9]{1,3}$"
    org.pkgs <- grep(org.string, pkgs, value=TRUE)
    names(org.pkgs) <- NULL 
    return(org.pkgs)
}

.isAvailable <- function(pkg, type="annotation")
{
    if(!(pkg %in% .packages(all.available=TRUE)))
    {   
        message(paste0("Corresponding ", type,  " package not found: ", 
            pkg, "\nMake sure that you have it installed."))
        choice <- readline("Install it now? (y/n): ")
        if(choice == "y")
        {   
            biocLite <- NULL
            source("http://bioconductor.org/biocLite.R")
            biocLite(pkg, suppressUpdates=TRUE, suppressAutoUpdate=TRUE)
        }   
        else stop(paste("Package", pkg, "is not available"))
    }   
    require(pkg, character.only = TRUE)
}

gene.report <- function(s, gt, out.dir)
{
    htmlRep <- HTMLReport(basePath=out.dir, reportDirectory="reports",
        shortName=setName(s), title=paste(setName(s), "Gene Report", sep=": "))
    publish(gt[geneIds(s),], htmlRep, .modifyDF=list(ncbi.gene.link))#, pubmed.link),
        #colClasses = c(rep("sort-string-robust", 3), rep("sort-num-robust", ncol(gt)-3 )))
    rep <- finish(htmlRep)
    return(rep)
}

ncbi.gene.link <- function(object, ...)
{
    EZ.COL <- config.ebrowser("EZ.COL")
    col <- as.character(object[,EZ.COL])
    link <- paste0(config.ebrowser("GENE.URL"), col)
    object[,EZ.COL] <- hwrite(col, link=link, table=FALSE)
    return(object)
}

view.set <- function(eset, out.prefix)
{
    out.files <- paste(out.prefix, 
        c("sview.html", "volc.png", "hmap.png"), sep="_" )
    
    if(!file.exists(out.files[1]))
    {
        ##
        # 1 PLOT: heatmap & volcano
        ##
        volc.html <- make.volc.html(eset, out.files[2])
        hmap.html <- make.hmap.html(eset, out.files[3])
        cont <- make.view(volc.html, hmap.html) 
        cat(cont, file=out.files[1])
    }
    views <- basename(out.files[1])
    return(views)

    # flat file set report
    # Do we need this anymore?
#    report.file <- sub("pdf$", "txt", plot.file)
#    if(!file.exists(report.file))
#    {
#        fDat <- signif(fDat[order(fDat[,2]),], digits=2)
#        fDat <- cbind(rownames(fDat), fDat)
#        colnames(fDat)[1] <- "GENE"
#        write.table(fDat, file=report.file, row.names=FALSE, quote=FALSE, sep="\t")
#    }
#    views <- c(views, basename(report.file))
#    return(rev(views))
}

view.path <- function(s, eset, out.prefix)
{
    org <- substring(s, 1, 3)
    pwy.id <- sub("^[a-z]{3}", "", s)
    fc <- fData(eset)[,config.ebrowser("FC.COL")]
    gn.cols <- sapply(c("SYM.COL", "GN.COL"), config.ebrowser)
    gnam <- apply(fData(eset)[,gn.cols], 1, paste, collapse=": ")
    names(fc) <- names(gnam) <- featureNames(eset)

    out.files <- paste(out.prefix, 
        c("kview.html", "kpath.png", "kgraph.png"), sep="_")
    
    if(!file.exists(out.files[1]) && pwy.id != "01100")
    {
        ##
        # 1 PLOT: kegg.native & kegg.graph
        ##
        kpath.html <- make.kpath.html(fc, pwy.id, org, out.files[2])
        kgraph.html <- make.kgraph.html(fc, gnam, pwy.id, org, out.files[3])
        cont <- make.view(kgraph.html, kpath.html, gene.html.pos="topright") 
        cat(cont, file=out.files[1])
    }
    
    views <- basename(out.files[1])
    return(views)
} 

view.graph <- function(eset, sgrn, alpha, out.prefix)
{
    out.files <- paste0(out.prefix, "_gview.", c("html", "png", "txt"))
    
    if(!file.exists(out.files[1]))
    {
        ##
        # 1 PLOT: ggea.graph
        ##
        ggraph.html <- make.ggraph.html(eset, sgrn, alpha, out.files[2])
        void.html <- "void.html"
        cat(hmakeTag('html'), file=file.path(dirname(out.prefix),void.html))
        cont <- make.view(ggraph.html, void.html)   
        cat(cont, file=out.files[1])
    }
    
    views <- basename(out.files[1])
    return(views)
}

make.hmap.html <- function(eset, img.file)
{
    width <- config.ebrowser("PLOT.WIDTH") 
    height <- config.ebrowser("PLOT.HEIGHT")

    # 1: make the plot
    # (a) heatmap 1: all features vs all samples
    expr <- exprs(eset)
    rownames(expr) <- fData(eset)[,config.ebrowser("SYM.COL")]
    png(img.file, width=width, height=height)
    exprs.heatmap(expr=expr, grp=pData(eset)[,config.ebrowser("GRP.COL")])
    dev.off()
    img.tag <- hwriteImage(basename(img.file))

    # (b) heatmap 2: most signif features
    max.row <- 40
    fc <- abs(fData(eset)[,config.ebrowser("FC.COL")])
    p <- -log(fData(eset)[,config.ebrowser("ADJP.COL")], base=10)
    ind <- (fc >= 1) | (p >= 1)
    eset <- eset[ind,]
    if(nrow(eset) > 1)
    {
        if(nrow(eset) > max.row)
        {
            fc <- fc[ind]
            p <- p[ind]
            score <- sqrt(fc^2 + p^2)
            eset <- eset[order(score, decreasing=TRUE)[seq_len(max.row)],]
            #    # select most variable samples
            #if(ncol(eset) > max.col)
            #{
            #    svar <- apply(exprs(eset), 2, var)
            #    eset <- eset[,order(svar, decreasing=TRUE)[seq_len(max.col)]]
            #}
        }
        expr <- exprs(eset)
        rownames(expr) <- fData(eset)[,config.ebrowser("SYM.COL")]
        img2 <- sub(".png$", "2.png", img.file)
        png(img2, width=width, height=height)
        exprs.heatmap(expr=expr, grp=pData(eset)[,config.ebrowser("GRP.COL")])
        dev.off()
        img.tag <- paste0(img.tag, hwriteImage(basename(img2)))
    }
    
    # 2: make the html
    hmap.html <- sub("png$", "html", img.file)
    cont <- hmakeTag('html', hmakeTag('body', img.tag))
    cat(cont, file=hmap.html)
    return(basename(hmap.html))
}

make.volc.html <- function(eset, img.file)
{    
    width <- config.ebrowser("PLOT.WIDTH") 
    height <- config.ebrowser("PLOT.HEIGHT")

    fc <- fData(eset)[,config.ebrowser("FC.COL")]
    p <- fData(eset)[,config.ebrowser("ADJP.COL")]
    # 1: make the plot
    png(img.file, width=width, height=height)
        volcano(fc, p)
        mai.px <- par('mai') * 75
        usr.grd <- par('usr')
    dev.off()
    p <- -log(p, base=10)
    # 2: make the html

    # INFER LINK COORDS
    #
    # mai: a numerical vector of the form ‘c(bottom, left, top, right)
    #       which gives the margin size specified in inches
    #
    # default plot: (75px == 1 inch)
    #                       btm     left    top     right
    #       par('mai') :    1.02    0.82    0.82    0.42
    #           px:         76.5    61.5    61.5    31.5
    
    # actual plot region:
    #   - origin(x0,y0) on the top left is at (61.5, 61.5)
    #   - width: 500 - 61.5 - 31.5 = 407
    #   - height: 500 - 76.5 - 61.5 = 362
    
    plot.width <- width - sum(mai.px[c(2,4)])
    plot.height <- height - sum(mai.px[c(1,3)])
    
    # re-scale user grid
    min.x <- min(0, usr.grd[1])
    fc <- fc + abs(min.x)
    max.x <- usr.grd[2] + abs(min.x)
    min.y <- min(0, usr.grd[3])
    p <- p + abs(min.y)
    max.y <- usr.grd[4] + abs(min.y)
    
    y.scalef <- plot.height / max.y
    x.scalef <- plot.width / max.x
    
    cx <- fc * x.scalef + mai.px[2]
    cy <- p * y.scalef + mai.px[1]
    
    # turn y-coords around, (0,0) is topleft and not btmleft
    cy <- height - cy
    coord <- cbind(cx, cy, cx+10, cy-10)
        
    # volcano html
    volc.html <- sub("png$", "html", img.file) 
    con <- file(volc.html, open="w")
    
    gn.cols <- sapply(c("SYM.COL", "GN.COL"), config.ebrowser)
    titles <- apply(fData(eset)[,gn.cols], 1, paste, collapse=": ") 
    refs <- paste0(config.ebrowser("GENE.URL"), featureNames(eset))
    geneplotter::imageMap(coord, con, 
        list(HREF=refs, TITLE=titles, TARGET=rep("gene", nrow(coord))), 
        basename(img.file))    
    close(con)
    return(basename(volc.html))
}

make.kpath.html <- function(fc, pwy.id, org, img.file)
{
    # 1: make the plot
    # run pathview for getting native overplotted image
    out.dir <- dirname(img.file)
    #suppressWarnings(suppressMessages(
    pathview(gene.data=fc, 
        pathway.id=pwy.id, species=org, kegg.dir=out.dir, out.suffix="kpath")
    # ))
    pv.out <- file.path(getwd(), paste0(org, pwy.id, ".kpath.png"))
    file.rename(from=pv.out, to=img.file)
    
    # 2: make the html
    kpath.html <- sub("png$", "html", img.file)
    cont <- hmakeTag('html', hmakeTag('body', hwriteImage(basename(img.file))))
    cat(cont, file=kpath.html)
    return(basename(kpath.html))
}

make.kgraph.html <- function(fc, gname, pwy.id, org, img.file)
{
    width <- config.ebrowser("PLOT.WIDTH") 
    height <- config.ebrowser("PLOT.HEIGHT")

    # 1: make the plot
    # run pathview2 for getting the kegg graph
    out.dir <- dirname(img.file)
    png(img.file, width=width, height=height)
    par(mai=rep(0,4))
    gr <- 
        suppressWarnings(suppressMessages(
        pathview2(gene.data=fc, 
        pathway.id=pwy.id, species=org, kegg.dir=out.dir)
         ))
    dev.off()   
 
    # 2: make the html
    kgraph.html <- sub("png$", "html", img.file)
    if(is(gr, "graph"))
    {
        nd <- nodeRenderInfo(gr)$kegg.ids
        nam <- sapply(names(nd), function(n) 
            ifelse(nd[[n]][1] %in% names(gname), gname[nd[[n]][1]], nodeRenderInfo(gr)$label[[n]]))
        names(nam) <- names(nd)
        kstr <- sapply(nd, function(n) 
            paste(paste(org, n, sep=":"), collapse="+"), USE.NAMES=FALSE)
        con <- file(kgraph.html, open="w")
        refs <- paste0(config.ebrowser("KEGG.GENE.URL"), kstr)
        biocGraph::imageMap(gr, con=con,
            tags=list(HREF=refs, TITLE = nam, TARGET = rep("gene", length(nd))),
            imgname=basename(img.file), width=width, height=height)    
        close(con)
    }
    else cat(hmakeTag('html'), file=kgraph.html)
    return(basename(kgraph.html))
}

make.ggraph.html <- function(eset, sgrn, alpha, img.file)
{
    width <- config.ebrowser("PLOT.WIDTH") 
    height <- config.ebrowser("PLOT.HEIGHT")

    sggea.graph <- NULL
    if(nrow(sgrn) > 0)
        sggea.graph <- construct.ggea.graph(grn=sgrn, eset=eset, alpha=alpha)
    
    # ggea graph png
    png(img.file, width=width, height=height)
    par(mai=rep(0,4))
    if(!is.null(sggea.graph))
        sggea.graph <- plot.ggea.graph(sggea.graph,
            show.scores=(numEdges(sggea.graph) < config.ebrowser("NR.SHOW")))
    else plot(NA, axes=FALSE, xlim=c(0,1), ylim=c(0,1),
        ylab="", xlab="", main="No edges in network for this set!")
    dev.off()

    # txt report
    report.file <- sub("png$", "txt", img.file)
    if(!is.null(sggea.graph))
    {
            consistency <- sggea.graph@renderInfo@edges$label
            cons.tbl <- cbind(names(consistency), consistency)
            cons.tbl <- cons.tbl[order(as.numeric(consistency), decreasing=TRUE),]
            colnames(cons.tbl) <- c("EDGE", "CONSISTENCY")
            write.table(cons.tbl,
                file=report.file, row.names=FALSE, quote=FALSE, sep="\t")
    }
    else cat("No edges in network for this set!", file=report.file)

    ggraph.html <- sub("view.png$", "graph.html", img.file)
    if(!is.null(sggea.graph))
    {
        gn.cols <- sapply(c("SYM.COL", "GN.COL"), config.ebrowser)
        gnam <- apply(fData(eset)[,gn.cols], 1, paste, collapse=": ")
        names(gnam) <- featureNames(eset)

        # image map
        nd <- nodes(sggea.graph)
        nam <- gnam[nd]
        con <- file(ggraph.html, open="w")
        refs <- paste0(config.ebrowser("GENE.URL"),  nd)
        biocGraph::imageMap(sggea.graph, con=con,
            tags=list(HREF=refs, TITLE = nam, TARGET = rep("gene", length(nd))),
            imgname=basename(img.file), width=width, height=height)    
        close(con)
    }
    else cat(hmakeTag('html'), file=ggraph.html)
    return(basename(ggraph.html))
}

make.view <- function(html1, html2, gene.html.pos=c("bottom", "topright"))
{
    gene.html.pos <- match.arg(gene.html.pos)
    s <- unlist(strsplit(html1, "_"))[1]
    
    head <- hmakeTag("head", hmakeTag("title", s))
    html1.tag <- hmakeTag("frame", name="volc", src=html1)
    html2.tag <- hmakeTag("frame", name="hmap", src=html2)
    html3.tag <- hmakeTag("frame", name="gene", scrolling="auto")
    
    if(gene.html.pos == "topright")
    {
        bkp <- html3.tag
        html3.tag <- html2.tag
        html2.tag <- bkp
    }

    f1.tag <- hmakeTag("frameset", 
        paste0(sub("</frame>", "", c(html1.tag, html2.tag)), collapse=""), 
        cols=paste0(config.ebrowser("PLOT.WIDTH") + 30, ",*"), border=0)
    html3.tag <- sub("</frame>", "", html3.tag)
    f2.tag <- hmakeTag("frameset", paste(f1.tag, html3.tag),
        rows=paste0(config.ebrowser("PLOT.HEIGHT") + 30, ",*"), border=0)
    cont <- hmakeTag("html", paste0(head,f2.tag))
    return(cont)
}

get.html.of.marked.pathway <- function(pwy, oids)
{
    pwy <- sub("^path:", "", pwy)
    oids <- gsub("[a-z]{3}:", "", oids)
    coids <- paste(oids, collapse="+")
    request <- pwy
    if(nchar(coids)) request <- paste(request, coids, sep="+")
    return(paste0(config.ebrowser("KEGG.SHOW.URL"), request))
}

pathview2 <- function(
    gene.data, pathway.id, species = "hsa", kegg.dir=".",
    gene.idtype="entrez", gene.annotpkg=NULL, min.nnodes=3, 
    map.null=TRUE, map.symbol=TRUE, node.sum="sum", limit=1, bins=10,
    both.dirs=TRUE, low="green", mid="gray", high="red",
    na.col="transparent", afactor=1, text.width = 15
)
{
    gd.names=names(gene.data)
    ng=length(gene.data)
    nsamp.g=1
    gene.idtype=toupper(gene.idtype)

    bods <- get(data("bods", package="pathview"))
    gene.idtype.list <- get(data("gene.idtype.list", package="pathview"))
    species.data=kegg.species.code(species, na.rm=T, code.only=FALSE)
  
    species=species.data["kegg.code"]
    entrez.gnodes=species.data["entrez.gnodes"]==1
    if(is.na(species.data["ncbi.geneid"])){
        if(!is.na(species.data["kegg.geneid"])){
            msg.fmt="Only native KEGG gene ID is supported for this species,\nmake sure it looks like \"%s\"!"
            msg=sprintf(msg.fmt, species.data["kegg.geneid"])
            message("Note: ", msg)
        } else{
        stop("This species is not annotated in KEGG!")
        }
    }
    if(is.null(gene.annotpkg)) gene.annotpkg=bods[match(species, bods[,3]),1]
    if(length(grep("ENTREZ|KEGG", gene.idtype))<1){
        if(is.na(gene.annotpkg)) stop("No proper gene annotation package available!")
        if(!gene.idtype %in% gene.idtype.list) stop("Wrong input gene ID type!")
        gene.idmap=id2eg(gd.names, category=gene.idtype, pkg.name=gene.annotpkg)
        gene.data=mol.sum(gene.data, gene.idmap)
        gene.idtype="ENTREZ"
    }
    if(gene.idtype=="ENTREZ" & !entrez.gnodes){
        message("Info: Getting gene ID data from KEGG...")
        gene.idmap=keggConv("ncbi-geneid", species)
        message("Info: Done with data retrieval!")
        kegg.ids=gsub(paste(species, ":", sep=""), "", names(gene.idmap))
        ncbi.ids=gsub("ncbi-geneid:", "", gene.idmap)
        gene.idmap=cbind(ncbi.ids, kegg.ids)
        gene.data=mol.sum(gene.data, gene.idmap)
        gene.idtype="KEGG"
    }

    #parse  
    warn.fmt="Parsing %s file failed, please check the file!"
    pathway.name = paste(species, pathway.id, sep = "")
    kfiles=list.files(path=kegg.dir, pattern="[.]xml|[.]png")
    tfiles=paste(pathway.name, c("xml","png"), sep=".")
    if(!all(tfiles %in% kfiles)){
        dstatus=download.kegg(pathway.id = pathway.id, species = species, kegg.dir=kegg.dir)
        if(dstatus=="failed") {
        warn.fmt="Failed to download KEGG xml/png files, %s skipped!"
        warn.msg=sprintf(warn.fmt, pathway.name)
        message("Warning: ", warn.msg)
        return(invisible(0))
        }
    }
    
    xml.file <- file.path(kegg.dir, paste(pathway.name, "xml", sep="."))
    gR1=try(pathview:::parseKGML2Graph2(xml.file, genes=F, expand=FALSE, split.group=FALSE), silent=TRUE)
    node.data=try(node.info(gR1), silent=T)
    if(class(node.data)=="try-error"){
      warn.msg=sprintf(warn.fmt, xml.file)
      message("Warning: ", warn.msg)
      return(invisible(0))
    }

    gene.node.type="gene"
    plot.data.gene=node.map(gene.data, node.data, 
        node.types=gene.node.type, node.sum=node.sum, entrez.gnodes=entrez.gnodes)
    kng=plot.data.gene$kegg.names
    kng.char=gsub("[0-9]", "", unlist(kng))
    if(any(kng.char>"")) entrez.gnodes=FALSE
    if(map.symbol & entrez.gnodes) {
              if(is.na(gene.annotpkg)) {
                warn.fmt="No annotation package for the species %s, gene symbols not mapped!"
                warn.msg=sprintf(warn.fmt, species)
                message("Warning: ", warn.msg)
              } else {
                plot.data.gene$labels <- eg2id(as.character(
                    plot.data.gene$kegg.names), category="SYMBOL", pkg.name=gene.annotpkg)[,2]
                mapped.gnodes <- rownames(plot.data.gene)
                node.data$labels[mapped.gnodes] <- plot.data.gene$labels
              }
    }
    cols.ts.gene <- node.color(plot.data.gene, limit, bins, both.dirs=both.dirs,
        discrete=FALSE, low=low, mid=mid, high=high, na.col=na.col)
           
    #group nodes mapping and merge
    grp.idx=node.data$size>1
    if(sum(grp.idx)>0){
        sub2grp=cbind(unlist(node.data$component[grp.idx], use.names=F), rep(names(grp.idx)[grp.idx], node.data$size[grp.idx]))
        du.idx=duplicated(sub2grp[,1])
        if(sum(du.idx)>0){
            du.rn=sub2grp[,1] %in% sub2grp[du.idx,1]
            message("Warning: reconcile groups sharing member nodes!")
            print(sub2grp[du.rn,])
            du.grps=sub2grp[du.idx,]
            rn=which(du.idx)
            for(r in rn){
                comps=node.data$component[[sub2grp[r,2]]]
                comps=comps[comps!=sub2grp[r,1]]
                node.data$component[[sub2grp[r,2]]]=comps
                node.data$size[sub2grp[r,2]]=node.data$size[sub2grp[r,2]]-1
            }
            sub2grp=sub2grp[!du.idx,]
        }
        rownames(sub2grp)=sub2grp[,1]
        for(gn in names(grp.idx)[grp.idx])
            gR1=combineKEGGnodes(node.data$component[[gn]], gR1, gn)
    } else sub2grp=NULL
    nNames=nodes(gR1)
    nSizes=node.data$size[nNames]

    # GENES ONLY!
    nNames <- nNames[node.data$type[nNames] == "gene"]
    gR1 <- subKEGGgraph(nNames, gR1)

    #unconnected nodes processing
    deg=degree(gR1)
    deg=deg$inDegree+deg$outDegree
    if(sum(deg<1)>0){
        gR2=subKEGGgraph(nNames[deg>0], gR1)
        nNames=nNames[deg>0]
        nSizes=nSizes[deg>0]
        if(!is.null(sub2grp))
            sub.idx=sub2grp[,1] %in% nNames |sub2grp[,2] %in% nNames
        else sub.idx=0
    } else {
        gR2=gR1
        if(!is.null(sub2grp)) sub.idx=rep(T, nrow(sub2grp))
        else sub.idx=0
    }

    if(length(nNames)<2){
        msg=sprintf("%s not rendered, 0 or 1 connected nodes!\nTry \"kegg.native=T\" instead!", pathway.name)
        message("Note: ", msg)
        return(list())
    }
 
     
    #give up the KEGG positions, use graphviz layout
    #general attributes
    attrs=list()
    attrs$graph$rankdir="LR"
    attrs$node <- list(fixedsize=FALSE)

    #node attributes
    ntype=node.data$type[nNames]
    cpd.idx=ntype=="compound"
    map.idx=ntype=="map"
    rect.idx=!(cpd.idx|map.idx)
    nAttr=list()
    nAttr$label=rep('', length(nNames))
    shapes=node.data$shape[nNames]
    if(any(cpd.idx)) shapes[cpd.idx]="ellipse"
    if(any(map.idx)) shapes[map.idx]="plaintext"
    nAttr$shape=shapes
    nAttr$height=.75*17/46*nSizes*afactor
    nAttr$width=rep(.75, length(nNames))*afactor
    
    if(any(cpd.idx)){
        nAttr$height[cpd.idx]=nAttr$height[cpd.idx]*1.5
        nAttr$width[cpd.idx]=nAttr$width[cpd.idx]*1.5
    }
    if(any(map.idx)){
        nAttr$height[map.idx]=nAttr$height[map.idx]*1.5
        nAttr$width[map.idx]=nAttr$width[map.idx]*2
    }
    nAttr<- lapply(nAttr, function(x){ names(x) <- nNames; return(x) })

    na.col=rgb(t(col2rgb(na.col)), maxColorValue=255)
    fillcol=rep(na.col, length(nNames))
    names(fillcol)=nNames

    #edge attributes
    subdisplay <- subtypeDisplay(gR2)
    if(length(subdisplay)<1) eAttrs=list() else{
        KEGGEdgeSubtype <- get(data("KEGGEdgeSubtype", package="pathview"))
        na.rn=apply(subdisplay, 2, function(x) sum(is.na(x))==7)
        if(sum(na.rn)>0) subdisplay[,na.rn]=KEGGEdgeSubtype[KEGGEdgeSubtype[,1]=="others",rownames(subdisplay)]
        eLabel <- subdisplay["label", ]
        eCol <- subdisplay["color", ]
        eTextCol <- subdisplay["fontcolor", ]
        eLty <- subdisplay["style", ]
        eArrowhead <- subdisplay["arrowhead", ]
        if (ncol(subdisplay) == 1) {
            tmp <- colnames(subdisplay)[1]
            names(eLabel) <- names(eCol) <- names(eTextCol) <- tmp
            names(eLty) <- names(eArrowhead) <- tmp
        }
        eAttrs <- list(lty = eLty, col = eCol, textCol = eTextCol, 
                 label = eLabel, arrowhead = eArrowhead)
    }
  
    gR2.layout=gR2
    edgeRenderInfo(gR2.layout)=eAttrs
    layoutType= "dot"
    gR2.layout <- Rgraphviz::layoutGraph(gR2.layout, attrs = attrs, nodeAttrs=nAttr, layoutType=layoutType)
    edgeRenderInfo(gR2.layout)=eAttrs
    nri=nodeRenderInfo(gR2.layout)
    loc=list(x=nri$nodeX, y=nri$nodeY)
    if(sum(rect.idx)>0){
        w.unit=min(nri$lWidth[rect.idx])
        h.unit=min(nri$height[rect.idx])
    }
    cni=nSizes>1
    if(sum(cni)>0){
        xloc=rep(loc[[1]][cni], nSizes[cni])
        sn.y=unlist(sapply(nSizes[cni], function(x) seq(-(x-1)/2, (x-1)/2,1)),use.names =F)
        yloc=rep(loc[[2]][cni], nSizes[cni])+h.unit*sn.y
    } else xloc=yloc=NULL
    xloc.nd=c(xloc,loc[[1]][nSizes==1 & rect.idx])
    yloc.nd=c(yloc,loc[[2]][nSizes==1 & rect.idx])
    labs=node.data$labels
    labs[nNames[map.idx]]=sapply(labs[nNames[map.idx]],wordwrap,width=text.width, break.word=F)
    labs[nNames[cpd.idx]]=sapply(labs[nNames[cpd.idx]],wordwrap,width=text.width, break.word=T)

    cols.ts.gene=cbind(cols.ts.gene)
    nc.gene=max(ncol(cols.ts.gene),0)
    pn.suffix=colnames(cols.ts.gene)

    nidx.gene=which(nNames %in% rownames(cols.ts.gene))
    cidx.gene=match(nNames[nidx.gene], rownames(cols.ts.gene))
    sci.gene=match(sub2grp[sub.idx,1], rownames(cols.ts.gene))
    sci.node=match(sub2grp[sub.idx,1], nNames)

    #initialize node colors, independent of user data
    if(sum(rect.idx)>0){
        cn.col=rep(NA, sum(sub.idx))
        cn.col=fillcol[sci.node]
        names(cn.col)=sub2grp[sub.idx,1]
        rect.col=c(cn.col,fillcol[nSizes==1 & rect.idx])
        rect.col[rect.col==na.col]=NA
        rect.col=matrix(rect.col, ncol=1)
    }
    if(sum(cpd.idx)>0){
        ell.col=fillcol[cpd.idx]
        ell.col[ell.col==na.col]=NA
        ell.col=matrix(ell.col, ncol=1)
        w.e=min(nri$lWidth[cpd.idx])
        h.e=min(nri$height[cpd.idx])
        xloc.e=loc[[1]][cpd.idx]
        yloc.e=loc[[2]][cpd.idx]
    }

    fillcol=rep(na.col, length(nNames))
    names(fillcol)=nNames

    if(!is.null(cols.ts.gene) & sum(rect.idx)>0){
        fillcol[nidx.gene]=cols.ts.gene[cidx.gene,]
        cn.col=matrix(NA, nrow=sum(sub.idx), ncol=nc.gene)
        cn.col[]=cols.ts.gene[sci.gene,]
        rownames(cn.col)=sub2grp[sub.idx,1]
        if(nc.gene>1) rect.col=rbind(cn.col,fillcol[nSizes==1 & rect.idx])
        else rect.col=c(cn.col,fillcol[nSizes==1 & rect.idx])
        rect.col[rect.col==na.col]=NA
    }
  
    lab.names <- labs[nNames[!cpd.idx]]
    names(lab.names) <- nodes(gR2.layout)
    nodeRenderInfo(gR2.layout)$label <- lab.names
    nodeRenderInfo(gR2.layout)$fill <- fillcol
    
    # make edge colors
    eri <- edgeRenderInfo(gR2.layout)
    etype <- ifelse(eri$arrowhead == "normal", 
        1, ifelse(eri$arrowhead == "tee", -1, 0))
    ind <- etype %in% c(1,-1)
    nd <- unique(plot.data.gene[,c("labels","mol.data")])
    mol.data <- as.vector(nd[,"mol.data"])
    mol.data <- sapply(mol.data, function(i) ifelse(i > 1, 1, ifelse(i < -1, -1, i)))
    names(mol.data) <- as.vector(nd[,"labels"])
    nmol.data <- sapply(nNames, function(n) mol.data[lab.names[[n]]], USE.NAMES=FALSE)
    names(nmol.data) <- nNames

    from <- nmol.data[eri$enamesFrom[ind]]
    to <- nmol.data[eri$enamesTo[ind]]
    grn <- cbind(from, to, etype[ind])
    edge.cons <- apply(grn, 1, 
        function(x) if(any(is.na(x))) return(0) else return(is.consistent(x)))
    ecol <- determine.edge.color(edge.cons)
    edgeRenderInfo(gR2.layout)$col[ind] <- ecol 

    rg <- Rgraphviz::renderGraph(gR2.layout)

    nodeRenderInfo(rg)$kegg.ids <- 
        node.data$kegg.names[nodes(gR2.layout)]

    return(invisible(rg))
}


is.consistent <- function(grn.rel)
{
    act.cons <- mean(abs(grn.rel[1:2]))
    if(length(grn.rel) == 2) return(act.cons)   
    if(sum(sign(grn.rel[1:2])) == 0) act.cons <- -act.cons
    return( ifelse(grn.rel[3] == 1, act.cons, -act.cons) ) 
}

determine.edge.color <- function(edge.cons) 
    ifelse(edge.cons < 0, rgb(0,0,abs(edge.cons)), rgb(abs(edge.cons),0,0))

##
## determine.node.color
##
## function returns a color depending on fc sign and significance
##
determine.node.color <- function(node.info)
{
    color <- rgb(0,0,0)
    if(any(is.na(node.info))) return(color)

    # intensity of color is 1 - P,
    # i.e. the more significant, the more intense
    intens <- 1 - node.info[2]

    # fold change positive, i.e. upregulated -> green
    if(sign(node.info[1]) == 1) color <- rgb(0, intens, 0)
    # fold change negative, i.e. downregulated -> red
    else if(sign(node.info[1]) == -1) color <- rgb(intens, 0 , 0)

    return(color)
}
