/**************************************************
* TPM filter
***************************************************/
process tpm_filter {
    label 'python3'
    label 'smallTask'

    if (params.cloudProcess) { publishDir "${params.output}/${params.tpm_filter_dir}", mode: 'copy', pattern: "**.counts.filtered.formated" }
    else { publishDir "${params.output}/${params.tpm_filter_dir}", pattern: "**.counts.filtered.formated" }

    input:
    val(sample)
    path(count)
    val(condition)

    output:
    val sample, emit: samples // [mock_rep1, mock_rep2, treat_rep1, treat_rep2]
    path "**.counts.filtered.formated", emit: filtered_counts // [mock_rep1.counts.filtered.formated, mock_rep2.counts.filtered.formated, treat_rep1.counts.filtered.formated, treat_rep2.counts.filtered.formated]

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

    df.columns = pd.MultiIndex.from_tuples(df.columns)

    df_mean_tpm = pd.DataFrame()
    for cond in set(!{conditions}):
        mean = df.T.loc[cond, slice(None), 'tpm'].mean(axis=0)
        df_mean_tpm[cond] = mean

    df_filtered = df[np.sum(df_mean_tpm < !{params.tpm}, axis=1) != len(set(!{conditions}))]

    #out_sep_list = [None for i in !{samples}]

    for cond, df_cond in df_filtered.groupby(level=[0,1,2], axis=1):
        if cond[-1] == 'counts':
            df_cond.columns = [cond[1]]
            df_cond.reset_index(inplace=True)
            #out_sep_list[!{samples}.index(cond[1])] = df_cond
            df_cond.to_csv(f"{cond[1]}.counts.filtered.formated", sep='\\t', columns=['Geneid', cond[1]], header=False, index=False)
        
    #for i,df_i in enumerate(out_sep_list):
        #sample_name =!{samples}[i]
        #df_i.to_csv(f"{sample_name}.counts.filtered.formated", sep='\\t', columns=['Geneid', sample_name], header=False, index=False)
    '''
}