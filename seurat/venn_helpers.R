# Venn diagram helpers (two and three sets) based on the VennDiagram package.
library(gridExtra)
library(grid)

plotVenn2 <- function(list.items, categories, title="",
                      fill.colors=c("darkgreen",'deepskyblue'),
                      cat.colors=c("darkgreen",'deepskyblue3')){

  library(VennDiagram)

  a1 <- unique(list.items[[1]])
  a2 <- unique(list.items[[2]])

  x <- list('area1'=length( a1 ),
            'area2'=length( a2 ),
            'cross.area'=length( intersect( a1, a2 ))
  )

  if(missing(categories)){
    ns <- paste0("\n(n=",lengths(list(a1,a2)),")")
    categories <- paste0(names(list.items), ns)
  }

  cat("Length area1 (",categories[1],"):", x$area1,"\n")
  cat("Length area2 (",categories[2],"):", x$area2)

  #grid.newpage()
  vd <- draw.pairwise.venn(area1 = x$area1,
                           area2 = x$area2,
                           cross.area = x$cross.area,
                           category = categories,
                           fill=fill.colors,
                           euler.d = T, scaled = F,
                           print.mode = c('raw','percent'),
                           alpha=c(.6,.6),
                           label.col='white',
                           fontfamily = 'Helvetica',
                           cat.fontfamily = 'Helvetica',
                           cex = rep(1.8, 3),
                           cat.col = cat.colors,
                           cat.cex = rep(1.2, 2), # rep(0.95, 2)
                           margin = 0.09,
                           cat.pos = c(-20,20),
                           #cat.pos = c(-40,40),
                           #cat.pos = c(-60,60),
                           cat.dist = rep(0.05, 2),
                           ind = F) #cat.dist = 0.1

  grid.arrange(gTree(children=vd), top=title)

}


calculate.triple.venn.overlap <- function(a, b, c){

  in.123 <- sort(intersect(a, intersect(b,c)))
  in.12 <- sort(intersect(a, b))
  in.13 <- sort(intersect(a, c))
  in.23 <- sort(intersect(b, c))

  list("area1"=length(a),
       "area2"=length(b),
       "area3"=length(c),
       "n12"=length(in.12),
       "n23"=length(in.23),
       "n13"=length(in.13),
       "n123"=length(in.123),
       "in.123"=in.123,
       "in.12"=in.12,
       "in.13"=in.13,
       "in.23"=in.23)

}

plotVenn3 <- function(list.items, categories, title="", label.color='white',
                      fill.colors=c("darkgreen",'deepskyblue1','darkred'),
                      cat.colors=c("darkgreen",'deepskyblue3','darkred')){

  library(VennDiagram)

  x <- calculate.triple.venn.overlap(a = list.items[[1]], b = list.items[[2]], c = list.items[[3]])

  if(missing(categories)){
    ns <- paste0("\n(n=",lengths(list.items),")")
    categories <- paste0(names(list.items), ns)
  }

  vd <- draw.triple.venn(area1 = x$area1, area2 = x$area2, area3 = x$area3, n12 = x$n12, n23 = x$n23, n13 = x$n13, n123 = x$n123,
                         category = categories,
                         fill=fill.colors,
                         euler.d = T,
                         scaled = T,
                         print.mode = c('raw','percent'),
                         sigdigs = 1,
                         alpha=c(.6,.6,.6),
                         label.col=label.color,
                         fontfamily = 'Helvetica',
                         cat.fontfamily = 'Helvetica',
                         cex = rep(1.2, 7),
                         cat.col = cat.colors,
                         cat.cex = rep(0.95, 3), # rep(1.05, 3)
                         margin = 0.09, # 0.051,
                         cat.dist = rep(0.1, 3), ind = F) #cat.dist = 0.1

  grid.arrange(gTree(children=vd), top=title)

}
