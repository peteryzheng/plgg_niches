library(data.table)
library(parallel)

# pairwise categorical enrichment test
pairwise_categorical_enrichment_test = function(var1, var2, num_cores = 4, padj_method = 'BH'){
    if(!padj_method %in% c('holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr', 'none')){
        stop('padj_method must be one of the following: holm, hochberg, hommel, bonferroni, BH, BY, fdr, none')
    }
    # every single pair for pairwise comparison
    grid_df = expand.grid(var1 = sort(unique(var1)), var2 = sort(unique(var2))) 
    data.table(t(mcmapply(
        function(x,y){
            # individual fisher test for each pair
            fisher_test_var = fisher.test(table(var1 == x, var2 == y))
            return(c(
                var1 = as.character(x), var2 = as.character(y), 
                N_var1 = sum(var1 == x), N_var2 = sum(var2 == y),
                N_var1var2 = sum(var1 == x & var2 == y),
                p_value = fisher_test_var$p.value,
                odds_ratios = as.numeric(fisher_test_var$estimate)
            ))
        }, grid_df$var1, grid_df$var2, 
        mc.preschedule = TRUE, mc.set.seed = 55555,
        mc.cores = num_cores
    )))[
        , c(
            'N_var1','N_var2','N_var1var2','p_value','odds_ratios','q_value'
        ) := list(
            as.numeric(N_var1), as.numeric(N_var2), as.numeric(N_var1var2), 
            as.numeric(p_value), as.numeric(odds_ratios), 
            p.adjust(as.numeric(p_value), method = padj_method)
        )
    ]
}

# test multi-category categorical variable with multiple continuous variables
# an example of this can be -- test the enrichment of certain signature exposures in different cancer types
# the requirement for dcast df is that:
#   the first column is the individual id, 
#   the second column is the categorical variable that we are testing 
#   the rest of the columns are continuous variables
#       since the original goals was to test proportions, it's assumed that these continuous variables are related and affected by each other
#       and therefore we are using BY here to test for associations when we do multiple testing
pairwise_categorical_continuous_enrichment_test = function(dcast_df, num_cores = 4, padj_method = 'BY'){
    if(!padj_method %in% c('holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr', 'none')){
        stop('padj.method must be one of the following: holm, hochberg, hommel, bonferroni, BH, BY, fdr, none')
    }

    # test every single category with every single continuous variable
    individual_id_col = colnames(dcast_df)[1]
    categorical_col = colnames(dcast_df)[2]
    categories = sort(dcast_df[,.N,categorical_col][N >= 30][[categorical_col]])
    print(paste0(
        'The following categories are left out because there are less than 30 observations: ',
        paste0(sort(dcast_df[,.N,categorical_col][N < 30][[categorical_col]]), collapse = ', ')
    ))
    continuous_vars = colnames(dcast_df)[-(1:2)]
    grid_df = data.table(expand.grid(
        categorical_col = categories, 
        continuous_col = continuous_vars
    ))
    enrich_dt = data.table(t(mcmapply(
        function(x,y){
            t_test_var = t.test(
                dcast_df[[as.character(y)]][dcast_df[[categorical_col]] == x], 
                dcast_df[[as.character(y)]][dcast_df[[categorical_col]] != x]
            )
            return(c(
                categorical_col = as.character(x), continuous_cols = as.character(y),
                N_var1 = sum(dcast_df[[categorical_col]] == x), N_var2 = as.numeric(t_test_var$estimate[1]),
                N_var1var2 = sum(dcast_df[[categorical_col]] == x) * as.numeric(t_test_var$estimate[1]),
                p_value = t_test_var$p.value,
                delta = as.numeric(t_test_var$estimate[1] - t_test_var$estimate[2])
            ))
        }, grid_df$categorical_col, grid_df$continuous_col,
        mc.preschedule = TRUE, mc.set.seed = 55555,
        mc.cores = num_cores
    )))[
        , c(
            'N_var1','N_var2','N_var1var2','p_value','delta','q_value'
        ) := list(
            as.numeric(N_var1), as.numeric(N_var2), as.numeric(N_var1var2), 
            as.numeric(p_value), as.numeric(delta), 
            p.adjust(as.numeric(p_value), method = padj_method)
        )
    ]
    setnames(enrich_dt, 1, categorical_col)
    return(enrich_dt)
}

