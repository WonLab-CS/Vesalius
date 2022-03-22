################################################################################
###############################   Vesalius      ################################
################################################################################

#----------------------------/Vesalius Objects/--------------------------------#



setClassUnion("mat",c("matrix","dgCMatrix"))
setClassUnion("ter",c("data.frame","NULL"))


### Need to clean these up
### Increase the robustness of this
setClass("vesaliusObject",
         slots=list(tiles="data.frame",
                    embeddings = "matrix",
                    territories = "ter",
                    counts  = "mat"),
         prototype=c(tiles = data.frame(),
                     embeddings = list()),
        validity= function(object){
        if(class(object@tiles) != "data.frame"){
            stop("Unsupported Tile Format")
        }
        if(class(object@embeddings)[1]!= "matrix"){
            stop("Unsupported Embeddings Format")
        }
        if(class(object@territories)!= "data.frame" & class(object@territories)!= "NULL"){
            stop("Unsupported territories Format")
        }
        if(class(object@counts)!= "matrix" & class(object@counts)!= "dgCMatrix"){
            stop("Unsupported Count Format")
        }
        return(TRUE)
    }

)

vesaliusObject <- function(tiles=NULL,
                           embeddings=NULL,
                           territories=NULL,
                           counts=NULL){

    ves <- new("vesaliusObject",
               tiles = tiles,
               embeddings = embeddings,
               territories = territories,
               counts = counts)
    return(ves)
}






setMethod("show",
    signature = "vesaliusObject",
    definition = function(object){
      .simpleBar(TRUE)
      ### To be modified later
      cat("Vesalius Object containing:\n")
      ntiles <- length(unique(object@tiles$barcodes))
      cat(ntiles,"tiles \n")

      ## This will change if we add multiple embedding type
      embeds <- ncol(object@embeddings)
      cat(embeds, "embeddings\n")

      ##If there are territories
      if(!is.null(object@territories)){
          ter <-length(unique(object@territories$territory))
          cat(ter, "territories\n")
      }
      .simpleBar(TRUE)
    }

)
