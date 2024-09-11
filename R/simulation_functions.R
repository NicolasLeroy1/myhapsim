#' Haplotype Generation for Genitors
#'
#' This function generates haplotypes for a specified genitor from a VCF window.
#'
#' @param vcf_window A matrix representing the VCF data for a given window.
#' @param genitor_id The ID of the genitor for which haplotypes are generated.
#' 
#' @return A vector of haplotypes for the specified genitor.
haplotyping_genitor = function(vcf_window, genitor_id) {
  haps = sapply(vcf_window[genitor_id,],function(cell)strsplit(cell,"|")[[1]][c(1,3)])
  apply(haps,1,function(cell)Reduce(paste0,cell))
}

#' Convert VCF Window to Haplotype Sequences
#'
#' This function converts a VCF window into an array of haplotypic sequences, two for each individual.
#'
#' @param vcf_window A matrix representing the VCF data for a given window.
#' 
#' @return A matrix of haplotypic sequences for each individual in the VCF window.
haplotyping = function(vcf_window) {
  haps = apply(vcf_window,c(1,2),function(cell)strsplit(cell,"|")[[1]][c(1,3)])
  apply(haps,c(1,2),function(cell)Reduce(paste0,cell))
}

#' Convert Haplotype Names to Character Vectors
#'
#' This function converts haplotype names into vectors of characters.
#'
#' @param haps_names A vector of haplotype names.
#' 
#' @return A list of character vectors representing the haplotype sequences.
as_vector_list = function(haps_names) {
  lapply(haps_names, function(hap) strsplit(hap, "")[[1]])
}

#' Generate Haplotype Table from VCF Window
#'
#' This function generates a table of haplotypic sequences and their occurrences from a VCF window.
#'
#' @param vcf_window A matrix representing the VCF data for a given window.
#' 
#' @return A sorted table of haplotypic sequences and their frequencies.
get_haps_table = function(vcf_window) {
  haps = haplotyping(vcf_window)
  sort(table(haps),decreasing=TRUE)
}

#' Compute Distance Matrix Between Haplotypes
#'
#' This function computes a distance matrix between haplotypes based on the table of haplotypic sequences.
#'
#' @param haps_table A table of haplotypic sequences.
#' 
#' @return A matrix representing the pairwise distances between haplotypes.
get_distance_matrix_haps = function(haps_table) {
  vector_list = as_vector_list(names(haps_table))
  n = length(vector_list)
  mat = matrix(NA, n, n)
  for(i in 1:n) {
    for(j in 1:n) {
      mat[i,j] = sum(vector_list[[i]] != vector_list[[j]])
    }
  }
  mat
}

#' Compute Distance Matrix Between Parent and Child Haplotypes
#'
#' This function computes a distance matrix between parent and child haplotypes.
#'
#' @param genitor_haps_vector_list A list of character vectors representing the haplotypes of the parents.
#' @param children_hap_vector_list A list of character vectors representing the haplotypes of the children.
#' 
#' @return A matrix representing the distances between parent and child haplotypes.
parents_children_distance_matrix = function(genitor_haps_vector_list, children_hap_vector_list) {
  n = length(genitor_haps_vector_list)
  p = length(children_hap_vector_list)
  mat = matrix(NA, n, p)
  for(i in 1:n) {
    for(j in 1:p) {
      mat[i,j] = sum(genitor_haps_vector_list[[i]] != children_hap_vector_list[[j]])
    }
  }
  mat
}

#' Classify Haplotypes Based on Crossbreeding Table
#'
#' This function classifies haplotypes of children based on a crossbreeding table and VCF window.
#'
#' @param vcf_window A matrix representing the VCF data for a given window.
#' @param table_croisement A table representing the crossbreeding relationships.
#' @param phenot A dataframe containing phenotype data, including parent information.
#' 
#' @return A list containing classified haplotypes for children, and haplotypes for Y, L, and D genitors.
classify_haps <- function(vcf_window, table_croisement,phenot) {
  Y_mother_id = colnames(table_croisement)[4]
  L_mother_id = colnames(table_croisement)[1:3]
  D_father_id = rownames(table_croisement)
  classified_children_haps = haplotyping(vcf_window[phenot$MATRICULE,])
  for(father_id in rownames(table_croisement)) {
    for(mother_id in colnames(table_croisement)) {
      if(table_croisement[father_id,mother_id] > 0) {
        children_id = phenot$MATRICULE[(phenot$PERE==father_id) & (phenot$MERE==mother_id)]
        vcf_children = vcf_window[children_id,]

        haps_father = haplotyping_genitor(vcf_window,father_id)
        haps_mother = haplotyping_genitor(vcf_window,mother_id)
        haps_children = haplotyping(vcf_children)

        genitor_hap_table = sort(c(table(haps_father),table(haps_mother)),decreasing=TRUE)
        children_hap_table = sort(table(haps_children),decreasing=TRUE)
        genitor_haps_vector_list = as_vector_list(names(genitor_hap_table))
        children_haps_vector_list = as_vector_list(names(children_hap_table))
        distance_matrix = parents_children_distance_matrix(genitor_haps_vector_list, children_haps_vector_list)
        p = length(children_hap_table)
        for(j in (p + 1 - (1:p))) {
          i = which.min(distance_matrix[, j])
          new_hap = names(genitor_hap_table)[i]
          haps_children[ haps_children == names(children_hap_table)[j] ] = new_hap
        }
        classified_children_haps[,children_id] = haps_children
      }
    }
  }
  Y_haps = haplotyping_genitor(vcf_window,Y_mother_id)
  L_haps = haplotyping(vcf_window[L_mother_id,])
  D_haps = haplotyping(vcf_window[D_father_id,])
  list(classified_children_haps=classified_children_haps,Y_haps=Y_haps,L_haps=L_haps,D_haps=D_haps)
}
         
#' Simulate QTL Effects on Haplotypes
#'
#' This function simulates the effects of quantitative trait loci (QTL) on haplotypes.
#'
#' @param n_qtl Number of QTLs to simulate.
#' @param phenot A dataframe containing phenotype data.
#' @param chr_vcf_list A list of VCF matrices for different chromosomes.
#' @param total_effect_size The total effect size to distribute among the QTLs.
#' @param radii A list of radii to use when simulating QTL effects.
#' 
#' @return A list containing simulated QTL effects and haplotype classifications.
haplo_qtl_generation = function(n_qtl,phenot,chr_vcf_list,total_effect_size=1,radii=list(15)) {
  table_croisement = table(phenot$PERE,phenot$MERE)
  effect_size = rep(1,n_qtl)
  effect_size = effect_size*total_effect_size/sum(effect_size)
  lapply(1:n_qtl,function(qtl) {
    # Generate random QTL position and radius
    chromosome = sample(1:length(chr_vcf_list),1)
    radius = sample(radii,1)[[1]]
    position = sample((radius+1):(ncol(chr_vcf_list[[chromosome]])-radius),1)
    window = ((position-radius):(position + radius))

    # Classify haplotypes at this position and radius
    vcf_window = chr_vcf_list[[chromosome]][,window]
    haps_classification = classify_haps(vcf_window,table_croisement,phenot)
    haps = haps_classification$classified_children_haps
    Y_haps = haps_classification$Y_haps
    L_haps = haps_classification$L_haps
    D_haps = haps_classification$D_haps
    hap_origin = names(table(c(Y_haps,L_haps,D_haps)))
    hap_origin = hap_origin[hap_origin %in% haps]

    # Generate additive effect linked to haplotype classification
    y_effect = rep(0,ncol(haps))
    hap_number = sample(length(hap_origin)%/%2,1)
    effect_size_hap = rep(effect_size[qtl]/hap_number,hap_number)
    qtl_haps = sample(hap_origin,hap_number)
    for(i in 1:hap_number) {
      qtl_hap = qtl_haps[i]
      pop1 = (haps[1,] == qtl_hap)
      pop2 = (haps[2,] == qtl_hap)
      if(sum(pop1)==0){hap_effect = sqrt(effect_size_hap[i]) * pop2/sqrt(sum(pop2))}
      else if(sum(pop2)==0){hap_effect = sqrt(effect_size_hap[i]) * pop1/sqrt(sum(pop1))}
      else {hap_effect = sqrt(effect_size_hap[i]) * (pop1/sqrt(sum(pop1)) + pop2/sqrt(sum(pop2)))}
      y_effect = y_effect + hap_effect
    }
    y_effect = y_effect - mean(y_effect)
    y_effect = y_effect * sqrt(effect_size[qtl])/sqrt(mean(y_effect^2))
    # Return list of results
    list(chromosome = chromosome, position = position,radius = radius,haps = haps,Y_haps=Y_haps,D_haps=D_haps,L_haps=L_haps,y_effect=y_effect)
  })
}

#' Second QTL Generation Function
#'
#' This function generates quantitative trait loci (QTL) effects for haplotypes, considering multiple genitors and SNP effects.
#'
#' @param phenot A dataframe containing phenotype data.
#' @param chr_vcf_list A list of VCF matrices for different chromosomes.
#' @param total_effect_size The total effect size to distribute among the QTLs. Defaults to 1.
#' @param radii A list of radii to use when simulating QTL effects. Defaults to list(15).
#'
#' @return A list of results for each QTL, including chromosome, position, radius, haplotypes, and effect sizes.
second_haplo_qtl_generation = function(phenot,chr_vcf_list,total_effect_size=1,radii=list(15)) {
  table_croisement = table(phenot$PERE,phenot$MERE)
  effect_size = lapply(c("Y","L","D","SNP"),function(genitor){
    list(sizes = (1:3)*total_effect_size/(24),genitor=genitor)
    })
  lapply(effect_size,function(effect) {
    if(effect$genitor=="SNP"){
      lapply(effect$sizes,function(size){
        snp = rep(0,3)
        while(var(snp)==0){
          chromosome = sample(1:length(chr_vcf_list),1)
          position = sample(1:ncol(chr_vcf_list[[chromosome]]),1)
          snp = chr_vcf_list[[chromosome]][,position][phenot$MATRICULE]
          snp = unlist(lapply(snp,function(cell)sum(as.numeric(strsplit(cell,"|")[[1]][c(1,3)]))))
        }
        y_effect = snp - mean(snp)
        y_effect = y_effect * sqrt(size)/sqrt(mean(y_effect^2))
        return(list(chromosome = chromosome, position = position,radius = 0,haps = NULL,Y_haps=NULL,D_haps=NULL,L_haps=NULL,y_effect=y_effect,size=size,genitor=effect$genitor))
      })
    }
    else{
      lapply(effect$sizes,function(size){
        # Generate random QTL position and radius
        chromosome = sample(1:length(chr_vcf_list),1)
        radius = sample(radii,1)[[1]]
        position = sample((radius+1):(ncol(chr_vcf_list[[chromosome]])-radius),1)
        window = ((position-radius):(position + radius))

        # Classify haplotypes at this position and radius
        vcf_window = chr_vcf_list[[chromosome]][,window]
        haps_classification = classify_haps(vcf_window,table_croisement,phenot)
        haps = haps_classification$classified_children_haps
        Y_haps = names(table(haps_classification$Y_haps))
        Y_haps = Y_haps[Y_haps %in% haps]
        L_haps = names(table(haps_classification$L_haps))
        L_haps = L_haps[L_haps %in% haps]
        D_haps = names(table(haps_classification$D_haps))
        D_haps = D_haps[D_haps %in% haps]

        # Generate additive effect linked to haplotype classification

        qtl_hap = sample(get(paste0(effect$genitor,"_haps")),1)
        pop1 = (haps[1,] == qtl_hap)
        pop2 = (haps[2,] == qtl_hap)
        hap_effect = 1*pop1 + 1*pop2
        while(var(hap_effect)==0) {
          warning("QTL haplotype is not present in the population")
          qtl_hap = sample(get(paste0(effect$genitor,"_haps")),1)
          pop1 = (haps[1,] == qtl_hap)
          pop2 = (haps[2,] == qtl_hap)
          hap_effect = 1*pop1 + 1*pop2
        }
        y_effect = hap_effect
        y_effect = y_effect - mean(y_effect)
        y_effect = y_effect * sqrt(size)/sqrt(mean(y_effect^2))
        # Return list of results
        list(chromosome = chromosome, position = position,radius = radius,haps = haps,Y_haps=Y_haps,D_haps=D_haps,L_haps=L_haps,y_effect=y_effect,size=size,genitor=effect$genitor)
    })
    }
  })
}

#' Third QTL Generation Function with Chromosome Selection
#'
#' This function generates QTL effects for haplotypes, allowing for the selection of a subset of chromosomes.
#'
#' @param n_qtl Number of QTLs to simulate.
#' @param phenot A dataframe containing phenotype data.
#' @param chr_vcf_list A list of VCF matrices for different chromosomes.
#' @param n_chr Number of chromosomes to sample for QTL generation.
#' @param total_effect_size The total effect size to distribute among the QTLs. Defaults to 1.
#' @param radii A list of radii to use when simulating QTL effects. Defaults to list(15).
#'
#' @return A list of results for each QTL, including chromosome, position, radius, haplotypes, and effect sizes.
third_haplo_qtl_generation = function(n_qtl,phenot,chr_vcf_list,n_chr,total_effect_size=1,radii=list(15)) {
  table_croisement = table(phenot$PERE,phenot$MERE)
  effect_size = rep(1,n_qtl)
  effect_size = effect_size*total_effect_size/sum(effect_size)
  possible_chr = sample(1:length(chr_vcf_list),n_chr)
  lapply(1:n_qtl,function(qtl) {
    # Generate random QTL position and radius
    chromosome = sample(possible_chr,1)
    radius = sample(radii,1)[[1]]
    position = sample((radius+1):(ncol(chr_vcf_list[[chromosome]])-radius),1)
    window = ((position-radius):(position + radius))

    # Classify haplotypes at this position and radius
    vcf_window = chr_vcf_list[[chromosome]][,window]
    haps_classification = classify_haps(vcf_window,table_croisement,phenot)
    haps = haps_classification$classified_children_haps
    Y_haps = haps_classification$Y_haps
    L_haps = haps_classification$L_haps
    D_haps = haps_classification$D_haps
    hap_origin = names(table(c(Y_haps,L_haps,D_haps)))
    hap_origin = hap_origin[hap_origin %in% haps]

    # Generate additive effect linked to haplotype classification
    y_effect = rep(0,ncol(haps))
    hap_number = sample(length(hap_origin)%/%2,1)
    effect_size_hap = rep(effect_size[qtl]/hap_number,hap_number)
    qtl_haps = sample(hap_origin,hap_number)
    for(i in 1:hap_number) {
      qtl_hap = qtl_haps[i]
      pop1 = (haps[1,] == qtl_hap)
      pop2 = (haps[2,] == qtl_hap)
      if(sum(pop1)==0){hap_effect = sqrt(effect_size_hap[i]) * pop2/sqrt(sum(pop2))}
      else if(sum(pop2)==0){hap_effect = sqrt(effect_size_hap[i]) * pop1/sqrt(sum(pop1))}
      else {hap_effect = sqrt(effect_size_hap[i]) * (pop1/sqrt(sum(pop1)) + pop2/sqrt(sum(pop2)))}
      y_effect = y_effect + hap_effect
    }
    y_effect = y_effect - mean(y_effect)
    y_effect = y_effect * sqrt(effect_size[qtl])/sqrt(mean(y_effect^2))
    # Return list of results
    list(chromosome = chromosome, position = position,radius = radius,haps = haps,Y_haps=Y_haps,D_haps=D_haps,L_haps=L_haps,y_effect=y_effect)
  })
}

#' Simulation Function for Predicting QTL Effects
#'
#' This function simulates multiple scenarios of QTL effects, applying specified methods and generating results for different signal-to-noise ratios.
#'
#' @param methods A vector of methods to apply for prediction.
#' @param sim_method The simulation method to use for generating similarity matrices.
#' @param matrix_interval The interval size for cutting SNP matrices.
#' @param signal_sizes A vector of signal sizes to use for simulations.
#' @param n_qtl Number of QTLs to simulate.
#' @param n_sample Number of samples to generate.
#' @param svd_inertias A vector of singular value decomposition (SVD) inertias to use for generating similarity matrices.
#' @param chr_list A list of chromosome data.
#' @param chr_vcf_list A list of VCF matrices for different chromosomes.
#' @param phenot A dataframe containing phenotype data.
#' @param position_list A list of positions to use for cutting SNP matrices.
#' @param radii A list of radii to use when simulating QTL effects. Defaults to list(15).
#'
#' @return A list containing the simulation results, including predictions for different methods and signal sizes, and the corresponding QTL effects.
simulation_function = function(methods,sim_method,matrix_interval,signal_sizes,n_qtl,n_sample,svd_inertias,chr_list,chr_vcf_list,phenot,position_list,radii=list(15)) {
  X_list = lmrbay::cut_snp_matrices(chr_list,position_list,interval=matrix_interval,keep_position=TRUE)
  matrix_position_list = X_list$position_list
  X_list = X_list$X_list
  svd_list = lapply(svd_inertias,function(svd_inertia){
    temp = lmrbay::get_snp_similarity_matrices(X_list,svd_inertia,method=sim_method)
    V_list = lapply(temp,function(chr){
      lapply(chr,function(svd)svd$V)
    })
    D_list = lapply(temp,function(chr){
      lapply(chr,function(svd)svd$D)
    })
    list(V_list = V_list,D_list = D_list,svd_inertia=svd_inertia)
  })
  samples = lapply(1:n_sample,function(i){
    print(paste0("sample",i))
    qtl_list = haplo_qtl_generation(n_qtl,phenot,chr_vcf_list,radii=radii)
    total_effect = Reduce("+",lapply(qtl_list,function(qtl)qtl$y_effect))
    epsilon = rnorm(length(total_effect))
    epsilon = epsilon - mean(epsilon)
    epsilon = epsilon/sqrt(mean(epsilon^2))
    result = list()
    for(signal_size in signal_sizes) {
      noise_size = (1-signal_size)/signal_size
      y = total_effect+epsilon*sqrt(noise_size)*sqrt(mean(total_effect^2))
      for(svd in svd_list) {
        for(method in methods) {
          pred_name = paste(method,"signal",signal_size,"inertia",svd$svd_inertia,sep="_")
          result[[pred_name]] = get(paste0("rhm_",method))(y,svd$V_list,svd$D_list)
        }
      }
    }
    list(result=result,qtl_list=qtl_list,epsilon=epsilon,effect=total_effect)
  })
  chr_start = cumsum(c(0,sapply(position_list,function(chr)max(chr))))
  for(chr in 1:length(position_list)) {
    position_list[[chr]] = position_list[[chr]] + chr_start[chr]
    matrix_position_list[[chr]] = matrix_position_list[[chr]] + chr_start[chr]
  }
  return(list(samples=samples,position_list=position_list,matrix_position_list=matrix_position_list,chr_start=chr_start))
}

#' Second Simulation Function with Enhanced QTL Generation
#'
#' This function simulates multiple scenarios of QTL effects, using an enhanced QTL generation method and applying specified prediction methods.
#'
#' @param methods A vector of methods to apply for prediction.
#' @param sim_method The simulation method to use for generating similarity matrices.
#' @param matrix_interval The interval size for cutting SNP matrices.
#' @param signal_sizes A vector of signal sizes to use for simulations.
#' @param n_qtl Number of QTLs to simulate.
#' @param n_sample Number of samples to generate.
#' @param svd_inertias A vector of singular value decomposition (SVD) inertias to use for generating similarity matrices.
#' @param chr_list A list of chromosome data.
#' @param chr_vcf_list A list of VCF matrices for different chromosomes.
#' @param phenot A dataframe containing phenotype data.
#' @param position_list A list of positions to use for cutting SNP matrices.
#'
#' @return A list containing the simulation results, including predictions for different methods and signal sizes, and the corresponding QTL effects.
second_simulation_function = function(methods,sim_method,matrix_interval,signal_sizes,n_qtl,n_sample,svd_inertias,chr_list,chr_vcf_list,phenot,position_list) {
  X_list = lmrbay::cut_snp_matrices(chr_list,position_list,interval=matrix_interval,keep_position=TRUE)
  matrix_position_list = X_list$position_list
  X_list = X_list$X_list
  svd_list = lapply(svd_inertias,function(svd_inertia){
    temp = lmrbay::get_snp_similarity_matrices(X_list,svd_inertia,method=sim_method)
    V_list = lapply(temp,function(chr){
      lapply(chr,function(svd)svd$V)
    })
    D_list = lapply(temp,function(chr){
      lapply(chr,function(svd)svd$D)
    })
    list(V_list = V_list,D_list = D_list,svd_inertia=svd_inertia)
  })
  samples = lapply(1:n_sample,function(i){
    print(paste0("sample",i))
    qtl_list = second_haplo_qtl_generation(phenot,chr_vcf_list)
    total_effect = Reduce("+",lapply(qtl_list,function(origin){
      Reduce("+",lapply(origin,function(qtl)qtl$y_effect))
    }))
    epsilon = rnorm(length(total_effect))
    epsilon = epsilon - mean(epsilon)
    epsilon = epsilon/sqrt(mean(epsilon^2))
    result = list()
    for(signal_size in signal_sizes) {
      noise_size = (1-signal_size)/(signal_size)
      y = total_effect+epsilon*sqrt(noise_size)*sqrt(mean(total_effect^2))
      for(svd in svd_list) {
        for(method in methods) {
          pred_name = paste(method,"signal",signal_size,"inertia",svd$svd_inertia,sep="_")
          result[[pred_name]] = get(paste0("rhm_",method))(y,svd$V_list,svd$D_list)
        }
      }
    }
    list(result=result,qtl_list=qtl_list,epsilon=epsilon,effect=total_effect)
  })
  chr_start = cumsum(c(0,sapply(position_list,function(chr)max(chr))))
  for(chr in 1:length(position_list)) {
    position_list[[chr]] = position_list[[chr]] + chr_start[chr]
    matrix_position_list[[chr]] = matrix_position_list[[chr]] + chr_start[chr]
  }
  return(list(samples=samples,position_list=position_list,matrix_position_list=matrix_position_list,chr_start=chr_start))
}

#' Third Simulation Function with Selected Chromosomes
#'
#' This function simulates multiple scenarios of QTL effects, allowing for the selection of specific chromosomes and applying specified prediction methods.
#'
#' @param methods A vector of methods to apply for prediction.
#' @param sim_method The simulation method to use for generating similarity matrices.
#' @param matrix_interval The interval size for cutting SNP matrices.
#' @param signal_sizes A vector of signal sizes to use for simulations.
#' @param n_qtl Number of QTLs to simulate.
#' @param n_sample Number of samples to generate.
#' @param svd_inertias A vector of singular value decomposition (SVD) inertias to use for generating similarity matrices.
#' @param chr_list A list of chromosome data.
#' @param chr_vcf_list A list of VCF matrices for different chromosomes.
#' @param phenot A dataframe containing phenotype data.
#' @param position_list A list of positions to use for cutting SNP matrices.
#' @param n_chr Number of chromosomes to sample for QTL generation.
#'
#' @return A list containing the simulation results, including predictions for different methods and signal sizes, and the corresponding QTL effects.
third_simulation_function = function(methods,sim_method,matrix_interval,signal_sizes,n_qtl,n_sample,svd_inertias,chr_list,chr_vcf_list,phenot,position_list,n_chr) {
  X_list = lmrbay::cut_snp_matrices(chr_list,position_list,interval=matrix_interval,keep_position=TRUE)
  matrix_position_list = X_list$position_list
  X_list = X_list$X_list
  svd_list = lapply(svd_inertias,function(svd_inertia){
    temp = lmrbay::get_snp_similarity_matrices(X_list,svd_inertia,method=sim_method)
    V_list = lapply(temp,function(chr){
      lapply(chr,function(svd)svd$V)
    })
    D_list = lapply(temp,function(chr){
      lapply(chr,function(svd)svd$D)
    })
    list(V_list = V_list,D_list = D_list,svd_inertia=svd_inertia)
  })
  samples = lapply(1:n_sample,function(i){
    print(paste0("sample",i))
    qtl_list = third_haplo_qtl_generation(n_qtl,phenot,chr_vcf_list,n_chr)
    total_effect = Reduce("+",lapply(qtl_list,function(qtl)qtl$y_effect))
    epsilon = rnorm(length(total_effect))
    epsilon = epsilon - mean(epsilon)
    epsilon = epsilon/sqrt(mean(epsilon^2))
    result = list()
    for(signal_size in signal_sizes) {
      noise_size = (1-signal_size)/(signal_size)
      y = total_effect+epsilon*sqrt(noise_size)*sqrt(mean(total_effect^2))
      for(svd in svd_list) {
        for(method in methods) {
          pred_name = paste(method,"signal",signal_size,"inertia",svd$svd_inertia,sep="_")
          result[[pred_name]] = get(paste0("rhm_",method))(y,svd$V_list,svd$D_list)
        }
      }
    }
    list(result=result,qtl_list=qtl_list,epsilon=epsilon,effect=total_effect)
  })
  chr_start = cumsum(c(0,sapply(position_list,function(chr)max(chr))))
  for(chr in 1:length(position_list)) {
    position_list[[chr]] = position_list[[chr]] + chr_start[chr]
    matrix_position_list[[chr]] = matrix_position_list[[chr]] + chr_start[chr]
  }
  return(list(samples=samples,position_list=position_list,matrix_position_list=matrix_position_list,chr_start=chr_start))
}
