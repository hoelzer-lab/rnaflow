/**************************************************
* TPM filter
***************************************************/
process tpm_filter {
    label 'python3'
    label 'smallTask'

    if (params.cloudProcess) { 
        publishDir "${params.output}/${params.tpm_filter_dir}", mode: 'copy', pattern: "**.counts.filtered.formated.tsv"
        publishDir "${params.output}/${params.tpm_filter_dir}", mode: 'copy', pattern: "**.counts.tpm.tsv"
    } else { 
        publishDir "${params.output}/${params.tpm_filter_dir}", pattern: "**.counts.filtered.formated.tsv" 
        publishDir "${params.output}/${params.tpm_filter_dir}", pattern: "**.counts.tpm.tsv"
    }

    input:
    val(sample)
    path(count)
    val(condition)

    output:
    val sample, emit: samples // [mock_rep1, mock_rep2, treat_rep1, treat_rep2]
    path "**.counts.filtered.formated.tsv", emit: filtered_counts // [mock_rep1.counts.filtered.formated.tsv, mock_rep2.counts.filtered.formated.tsv, treat_rep1.counts.filtered.formated.tsv, treat_rep2.counts.filtered.formated.tsv]
    path "**.counts.tpm.tsv", emit: counts_with_tpm
    path "tpm_stats.tsv", emit: stats

    script:
    // make lists of strings for the python script
    samples = sample.collect { "\"${it}\"" }
    counts = count.collect { "\"${it}\"" }
    conditions = condition.collect { "\"${it}\"" }
    
    shell:
    '''
    #!/usr/bin/env python3

    import pandas as pd
    import numpy as np

    num_unfiltered_features = None
    num_filtered_features = None

    df = pd.DataFrame()
    cols = []
    for i, sample in enumerate(!{samples}):
        col_counts = (!{conditions}[i], sample, 'counts')
        col_tpm = (!{conditions}[i], sample, 'tpm')
        cols.append(col_counts)
        cols.append(col_tpm)
        
        df_i = pd.read_csv(!{counts}[i], index_col=['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length'], sep='\\t', comment='#')
        df_i.columns = [col_counts]

        df_i_tpm = df_i.reset_index()
        sample_reads = df_i_tpm.loc[:, [col_counts]].copy()
        gene_len = df_i_tpm.loc[:, ['Length']]
        rate = sample_reads.values / gene_len.values
        tpm = rate / np.sum(rate).reshape(1, -1) * 1e6

        df_i[col_tpm] = tpm
        df = pd.concat([df, df_i], axis=1)

        # write table with count and tpm values
        df_i.columns = ['counts', 'tpm']
        df_i.to_csv(f"{sample}.counts.tpm.tsv", sep='\\t')

        if not num_unfiltered_features:
            num_unfiltered_features = df_i.shape[0]

    df.columns = pd.MultiIndex.from_tuples(df.columns)

    df_mean_tpm = pd.DataFrame()
    for cond in set(!{conditions}):
        mean = df.T.loc[cond, slice(None), 'tpm'].mean(axis=0)
        df_mean_tpm[cond] = mean

    df_filtered = df[np.sum(df_mean_tpm < !{params.tpm}, axis=1) != len(set(!{conditions}))]

    for cond, df_cond in df_filtered.groupby(level=[0,1,2], axis=1):
        if cond[-1] == 'counts':
            df_cond.columns = [cond[1]]
            df_cond.reset_index(inplace=True)
            df_cond.to_csv(f"{cond[1]}.counts.filtered.formated.tsv", sep='\\t', columns=['Geneid', cond[1]], header=False, index=False)
            
            if not num_filtered_features:
                num_filtered_features = df_cond.shape[0]
        
    with open('tpm_stats.tsv', 'w') as stats:
        stats.write(f'"Number of retained features"\\t{num_filtered_features}\\n')
        stats.write(f'"Number of filtered out features"\\t{num_unfiltered_features-num_filtered_features}\\n')        
    '''
}