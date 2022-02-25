library(tidyverse)
library(rjson)
getwd()

masterTable <- read_tsv("../master_table_from_henry_20210603.txt")


formatFunction <- function(recodeJsonString, originalCodingJsonString) {
  if (!is.na(originalCodingJsonString) 
      && !is.null(originalCodingJsonString) 
      && originalCodingJsonString != "") {
    
    originalCoding <- unlist(fromJSON(originalCodingJsonString))
    
    recoding <- as.numeric(names(originalCoding))
    names(recoding) <- names(originalCoding)
    
    if (!is.na(recodeJsonString)
        && !is.null(recodeJsonString)
        && recodeJsonString != "") {

      recodingAsStrings <- unlist(fromJSON(recodeJsonString))

      recoding[names(recodingAsStrings)] <- as.numeric(recodingAsStrings)

    }
    recodingIndices <- recoding[names(originalCoding)]
    names(recodingIndices) <- originalCoding
    
    recodingIndices <- recodingIndices[order(recodingIndices)]
    
    recodingIndices[names(recodingIndices)] <- as.character(recodingIndices)
    
    recodingIndices[recodingIndices == "NaN"] <- "Excluded"

    return(paste(names(recodingIndices), recodingIndices, sep = " = ", collapse = ", "))
    
  }
}

masterTable$newValueCodingColumn <- NA_character_

for (rowname in rownames(masterTable)) {
  masterTable[rowname, "newValueCodingColumn"] <- formatFunction(masterTable[rowname, "total_recode_json", drop = T], masterTable[rowname, "answer_options_en_json", drop = T])
}


outputTable <- masterTable[c("from_googledrive_label_en", "Question en", "newValueCodingColumn", "model_type", "Significant in heatmap")]

write.table(outputTable, "../output_table_20210603.txt", sep = "\t", col.names = T, row.names = F, quote = F)
