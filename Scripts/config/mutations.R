ToMAF=function (ann, Center = NULL, refBuild = "hg19",
          table = "refGene", ens2hugo = TRUE, basename = NULL, 
          MAFobj = FALSE, sampleAnno = NULL) 
{
  require(data.table)
  start_time = proc.time()

  ann = as.data.table(ann)
  if (!"Tumor_Sample_Barcode" %in% colnames(ann)) {
    colnames(ann)[which(colnames(ann) == "sample")] = "Tumor_Sample_Barcode"
    cat("--Tumor_Sample_Barcode column not found. Creating sample IDs from filenames\n")
  }

  essential.col = c("Chr", "Start", "End", "Ref", "Alt", "Func.refGene", 
                    "Gene.refGene", "GeneDetail.refGene", "ExonicFunc.refGene", 
                    "AAChange.refGene")
  
  for (i in 1:length(essential.col)) {
    colId = suppressWarnings(grep(pattern = paste0("^", 
                                                   essential.col[i], "$"), x = colnames(ann), ignore.case = TRUE))
    if (length(colId) == 1) {
      colnames(ann)[colId] = essential.col[i]
    }
  }
  
  if (length(essential.col[!essential.col %in% colnames(ann)]) > 
      0) {
    message("Available fields:")
    print(colnames(ann))
    message(paste0("Missing required field in input file: "))
    print(essential.col[!essential.col %in% colnames(ann)])
    stop()
  }
  if (is.null(Center)) {
    Center = NA
  }
  ann[, `:=`(Hugo_Symbol, unlist(data.table::tstrsplit(Gene.refGene, 
                                                       split = ";", keep = 1)))]
  
  annovar_values = c(exonic = "RNA", splicing = "Splice_Site", 
                     ncRNA = "RNA", UTR5 = "5'UTR", UTR3 = "3'UTR", intronic = "Intron", 
                     upstream = "5'Flank", downstream = "3'Flank", intergenic = "IGR", 
                     `frameshift insertion` = "Frame_Shift_Ins", `frameshift deletion` = "Frame_Shift_Del", 
                     `frameshift block substitution` = "Frameshift_INDEL", 
                     `frameshift substitution` = "Frameshift_INDEL", stopgain = "Nonsense_Mutation", 
                     stoploss = "Nonstop_Mutation", `nonframeshift insertion` = "In_Frame_Ins", 
                     `nonframeshift deletion` = "In_Frame_Del", `nonframeshift block substitution` = "Inframe_INDEL", 
                     `nonframeshift substitution` = "Inframe_INDEL", `nonsynonymous SNV` = "Missense_Mutation", 
                     `synonymous SNV` = "Silent", unknown = "Unknown", ncRNA_exonic = "RNA", 
                     ncRNA_intronic = "RNA", ncRNA_UTR3 = "RNA", ncRNA_UTR5 = "RNA", 
                     ncRNA = "RNA", ncRNA_splicing = "RNA")
  ann_exonic = ann[Func.refGene %in% "exonic"]
  ann_res = ann[!Func.refGene %in% "exonic"]
  
  if (nrow(ann_exonic) == 0 & nrow(ann_res) == 0) {
    stop("No suitable exonic or intronic variants found!")
  }
  if (nrow(ann_exonic) > 0) {
    cat("-Processing Exonic variants\n")
    ann_exonic[, `:=`(Func.refGene, data.table::tstrsplit(x = as.character(ann_exonic$Func.refGene), 
                                                          split = ";", keep = 1))]
    cat("--Adding Variant_Classification\n")
    ann_exonic[, `:=`(Variant_Classification, annovar_values[ExonicFunc.refGene])]
    cat("--Parsing aa-change\n")
    aa_change = unlist(data.table::tstrsplit(x = as.character(ann_exonic$AAChange.refGene), 
                                             split = ",", fixed = TRUE, keep = 1))
    aa_tbl = lapply(aa_change, function(x) {
      x = unlist(strsplit(x = x, split = ":", fixed = TRUE))
      if (length(x) == 5) {
        tx = x[2]
        exon = x[3]
        txChange = x[4]
        aaChange = x[5]
      }
      else {
        tx = NA
        exon = NA
        txChange = NA
        aaChange = NA
      }
      data.table::data.table(tx, exon, txChange, aaChange)
    })
    aa_tbl = data.table::rbindlist(l = aa_tbl)
    if (length(aa_change) != nrow(ann_exonic)) {
      stop("Something went wrong parsing aa-change")
    }
    ann_exonic = cbind(ann_exonic, aa_tbl)
  }
  if (nrow(ann_res) > 0) {
    cat("-Processing Non-exonic variants\n")
    ann_res[, `:=`(Func.refGene, data.table::tstrsplit(x = as.character(ann_res$Func.refGene), 
                                                       split = ";", keep = 1))]
    cat("--Adding Variant_Classification\n")
    ann_res[, `:=`(Variant_Classification, annovar_values[Func.refGene])]
    ann = data.table::rbindlist(l = list(ann_exonic, ann_res), 
                                use.names = TRUE, fill = TRUE)
  }
  else {
    ann = ann_exonic
  }
  cat("-Adding Variant_Type\n")
  ann[, `:=`(ref_alt_diff, nchar(Ref) - nchar(Alt))]
  ann$Variant_Type = apply(ann[, .(Ref, Alt)], 1, function(x) {
    # xx = which(x == "-")
    # if (length(xx) == 0) {
    #   if (any(nchar(x) > 1)) {
    #     return("MNP")
    #   }
    #   else {
    #     return("SNP")
    #   }
    # }
    # else if (names(xx) == "Ref") {
    #   return("INS")
    # }
    # else if (names(xx) == "Alt") {
    #   return("DEL")
    # }
    # else {
    #   return(NA)
    # }
    
    if(nchar(x[1])>nchar(x[2])) {
      return("DEL")
    } else if (nchar(x[1])<nchar(x[2])) {
      return("INS")
    } else {
      return("SNP")
    }
    
  })
  ann_mnps = ann[Variant_Type %in% "MNP"]
  if (nrow(ann_mnps) > 0) {
    ann = ann[!Variant_Type %in% "MNP"]
    ann_mnps[, `:=`(Variant_Classification, "Unknown")]
    ann = rbind(ann, ann_mnps)
    rm(ann_mnps)
  }
  ann_indel = ann[Variant_Classification %in% c("Frameshift_INDEL", 
                                                "Inframe_INDEL")]
  if (nrow(ann_indel) > 0) {
    cat("-Fixing ambiguous INDEL annotations\n")
    ann = ann[!Variant_Classification %in% c("Frameshift_INDEL", 
                                             "Inframe_INDEL")]
    vc_fixed = lapply(1:nrow(ann_indel), function(i) {
      x = ann_indel[i, Variant_Classification]
      if (x == "Frameshift_INDEL") {
        if (ann_indel[i, ref_alt_diff] > 0) {
          return("Frame_Shift_Del")
        }
        else {
          return("Frame_Shift_Ins")
        }
      }
      else if (x == "Inframe_INDEL") {
        if (ann_indel[i, ref_alt_diff] > 0) {
          return("In_Frame_Del")
        }
        else {
          return("In_Frame_Ins")
        }
      }
      else {
        return(x)
      }
    })
    ann_indel[, `:=`(Variant_Classification, unlist(vc_fixed))]
    ann = rbind(ann, ann_indel)
  }
  if (table == "ensGene") {
    if (ens2hugo) {
      ens = system.file("extdata", "ensGenes.txt.gz", 
                        package = "maftools")
      cat("-Converting Ensemble Gene IDs into HGNC gene symbols\n")
      ens = data.table::fread(file = ens, sep = "\t", 
                              stringsAsFactors = FALSE)
      ann = merge(ann, ens, by.x = "Hugo_Symbol", by.y = "ens_id", 
                  all.x = TRUE)
      ann[, `:=`(ens_id, Hugo_Symbol)]
      ann[, `:=`(Hugo_Symbol, hgnc_symbol)]
      ann[, `:=`(Entrez_Gene_Id, Entrez)]
      cat("--Done. Original ensemble gene IDs are preserved under field name ens_id\n")
    }
  }
  ann[, `:=`(ref_alt_diff, NULL)]
  colnames(ann)[which(colnames(ann) %in% c("Chr", "Start", 
                                           "End", "Ref", "Alt"))] = c("Chromosome", "Start_Position", 
                                                                      "End_Position", "Reference_Allele", "Tumor_Seq_Allele2")
  ord1 = colnames(x = ann)[colnames(x = ann) %in% c("Hugo_Symbol", 
                                                    "Chromosome", "Start_Position", "End_Position", "Reference_Allele", 
                                                    "Tumor_Seq_Allele2", "Variant_Classification", "Variant_Type", 
                                                    "Tumor_Sample_Barcode", "tx", "exon", "txChange", "aaChange")]
  ord2 = colnames(x = ann)[!colnames(x = ann) %in% c("Hugo_Symbol", 
                                                     "Chromosome", "Start_Position", "End_Position", "Reference_Allele", 
                                                     "Tumor_Seq_Allele2", "Variant_Classification", "Variant_Type", 
                                                     "Tumor_Sample_Barcode", "tx", "exon", "txChange", "aaChange")]
  ann = ann[, c(ord1, ord2), with = FALSE]
  if (!is.null(basename)) {
    data.table::fwrite(x = ann, file = paste(basename, "maf", 
                                             sep = "."), sep = "\t")
  }
  cat("Finished in", data.table::timetaken(start_time), "\n")
  if (MAFobj) {
    m = read.maf(maf = ann, clinicalData = sampleAnno, verbose = FALSE)
    return(m)
  }
  else {
    return(ann)
  }
}


