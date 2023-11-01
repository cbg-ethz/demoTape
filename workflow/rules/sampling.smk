rule sampling:
    input:
        config['general']['input-loom']
    output:
        os.path.join(OUT_DIR, 'samples', 'c{cells}.r{run}.filtered_variants.csv'),
        os.path.join(OUT_DIR, 'samples', 'c{cells}.r{run}.filtered_variants_gt.csv')
    threads: 1
    resources:
        mem_mb = 49152,
        runtime = 60,
    params:
        minGQ = config['COMPASS'].get('minGQ', 30),
        minDP = config['COMPASS'].get('minDP', 10),
        minVAF = config['COMPASS'].get('minaVAF', 0.2),
        minVarGeno = config['COMPASS'].get('minVarGeno', 0.5),
        minCellGeno = config['COMPASS'].get('minCellGeno', 0.5),
        minMutated = lambda w: \
            max(50, int(w.cells) * config['COMPASS'].get('minMutated', 0.01)),
        maxRefVAF = config['COMPASS'].get('maxRefVAF', 0.05),
        minHomVAF = config['COMPASS'].get('minHomVAF', 0.95),
        minHetVAF = config['COMPASS'].get('minHetVAF', 0.35),
        proximity = config['COMPASS'].get('proximity', '25 50 100 200')
    shell:
        """
        python {SCRIPT_DIR}/mosaic_preprocessing.py \
            -i {input} \
            -o {output[0]} \
            --minGQ {params.minGQ} \
            --minDP {params.minDP} \
            --minVAF {params.minVAF} \
            --minVarGeno {params.minVarGeno} \
            --minCellGeno {params.minCellGeno} \
            --minMutated {params.minMutated} \
            --max_ref_VAF {params.maxRefVAF} \
            --min_hom_VAF {params.minHomVAF} \
            --min_het_VAF {params.minHetVAF} \
            --proximity {params.proximity} \
            --cell_no {wildcards.cells}
        """


rule run_COMPASS:
    input:
        os.path.join(OUT_DIR, 'samples', 'c{cells}.r{run}.filtered_variants.csv')
    output:
        os.path.join(OUT_DIR, 'COMPASS', 'c{cells}.r{run}.d{dbt}_cellAssignments.tsv')
    threads: 4
    resources:
        mem_mb = 10240,
        runtime = 1440,
    params:
        exe = config['COMPASS']['exe'],
        chain_length = config['COMPASS'].get('chain_length', 20000),
        CNV = config['COMPASS'].get('CNV', 0),
        sex = config['COMPASS'].get('sex', 'female'),
    shell:
        """
        {params.exe} \
            -i {OUT_DIR}/samples/c{wildcards.cells}.r{wildcards.run}.filtered \
            -o {OUT_DIR}/COMPASS/c{wildcards.cells}.r{wildcards.run}.d{wildcards.dbt} \
            --CNV {params.CNV} \
            --nchains {threads} \
            --doubletrate {wildcards.dbt} \
            --chainlength {params.chain_length} \
            --sex {params.sex}
        """


rule run_BnpC:
    input:
        os.path.join(OUT_DIR, 'samples',
            'c{cells}.r{run}.filtered_variants_gt.csv')
    output:
        os.path.join(OUT_DIR, 'BnpC', 'c{cells}.r{run}', 'assignment.txt')
    threads: 4
    resources:
        mem_mb = 10240,
        runtime = 630,
    params:
        exe = config['BnpC']['exe'],
        chain_length = config['BnpC'].get('chain_length', 20000),
        runtime = f'-r {config["BnpC"]["runtime"]} ' \
            if config['BnpC'].get('runtime', False) else '',
        pp = ' '.join([str(i) for i in config['BnpC'].get('params_prior', [1, 1])])
    shell:
        """
        python {params.exe} \
            {input} \
            -o {OUT_DIR}/BnpC/c{wildcards.cells}.r{wildcards.run} \
            -s {params.chain_length} \
            {params.runtime} \
            -e posterior ML MAP \
            -n {threads} \
            --param_prior {params.pp} \
            -np \
            -v 0
        """


rule run_SCG:
    input:
        os.path.join(OUT_DIR, 'samples',
            'c{cells}.r{run}.filtered_variants_gt.csv')
    output:
        os.path.join(OUT_DIR, 'SCG', 'c{cells}.r{run}.d{dbt}', 'assignments.tsv')
    threads: 1
    resources:
        mem_mb = 10240,
        runtime = 1440,
    params:
        wrapper_exe = config['SCG']['wrapper_exe'],
        chain_length = config['SCG'].get('chain_length', 100000),
        max_cl = 20,
    shell:
        """
        python {params.wrapper_exe} \
            {input} \
            -o {OUT_DIR}/SCG/c{wildcards.cells}.r{wildcards.run}.d{wildcards.dbt}/ \
            -s {params.chain_length} \
            -k {params.max_cl} \
            --doublets {wildcards.dbt} \
            --silent
        """


def get_result_files(wildcards):
    res_files = []
    for cells in CELL_NO:
        for run in RUNS_NO:
            for dbt in DOUBLET_RATE:
                if config.get('COMPASS', {}).get('run', False):
                    res_files.append(os.path.join(OUT_DIR, 'COMPASS',
                        f'c{cells}.r{run}.d{dbt}_cellAssignments.tsv'))
                if config.get('SCG', {}).get('run', False):
                    res_files.append(os.path.join(OUT_DIR, 'SCG',
                        f'c{cells}.r{run}.d{dbt}', 'assignments.tsv'))
                if config.get('BnpC', {}).get('run', False):
                    res_files.append(os.path.join(OUT_DIR, 'BnpC',
                        f'c{cells}.r{run}', 'assignment.txt'))
    return res_files


rule compare_runs:
    input:
        get_result_files
    output:
        os.path.join(OUT_DIR, f'summary.r{RUNS_TOTAL}.c{CELL_STR}.tsv')
    threads: 1
    resources:
        mem_mb = 32768,
        runtime = 60,
    params:
        min_cells = [1, 25, 50]
    script:
        f'{SCRIPT_DIR}/evaluate_sampling.py'