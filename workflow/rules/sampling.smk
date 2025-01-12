rule sampling:
    input:
        loom = ancient(LOOM_FILE)
    output:
        var = os.path.join(OUT_DIR, 'samples',
            'c{cells}.r{run}.filtered_variants.csv'),
        gt = os.path.join(OUT_DIR, 'samples',
            'c{cells}.r{run}.filtered_variants_gt.csv')
    wildcard_constraints:
        cells = r'\d+',
        run = r'\d+'
    conda: os.path.join(ENV_DIR, 'analysis.yaml')
    threads: 1
    resources:
        mem_mb_per_cpu = 49152,
        runtime = 60,
    params:
        minGQ = config.get('mosaic', {}).get('minGQ', 30),
        minDP = config.get('mosaic', {}).get('minDP', 10),
        minVAF = config.get('mosaic', {}).get('minVAF', 0.2),
        minVarGeno = config.get('mosaic', {}).get('minVarGeno', 0.5),
        minCellGeno = config.get('mosaic', {}).get('minCellGeno', 0.5),
        minMutated = lambda w: \
            max(50, int(w.cells) * config.get('mosaic', {}).get('minMutated', 0.05)),
        maxRefVAF = config.get('mosaic', {}).get('maxRefVAF', 0.05),
        minHomVAF = config.get('mosaic', {}).get('minHomVAF', 0.95),
        minHetVAF = config.get('mosaic', {}).get('minHetVAF', 0.35),
        proximity = config.get('mosaic', {}).get('proximity', '25 50 100 200')
    shell:
        """
        python {SCRIPT_DIR}/mosaic_preprocessing.py \
            -i {input} \
            -o {output.var} \
            --full_output \
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

rule filter_nonRelevant_SNPs:
    input:
        var = os.path.join(OUT_DIR, 'samples', 
            'c{cells}.r{run}.filtered_variants.csv'),
        gt = os.path.join(OUT_DIR, 'samples',
            'c{cells}.r{run}.filtered_variants_gt.csv'),
        dbsnp = DBSNP
    output:
        var = os.path.join(OUT_DIR, 'samples', 
            'c{cells}.r{run}.relevant.filtered_variants.csv'),
        gt = os.path.join(OUT_DIR, 'samples', 
            'c{cells}.r{run}.relevant.filtered_variants_gt.csv'),
    wildcard_constraints:
        cells = r'\d+',
        run = r'\d+'
    conda: os.path.join(ENV_DIR, 'analysis.yaml')
    threads: 1
    resources:
        mem_mb_per_cpu = 8096,
        runtime = 60,
    shell:
        """
        python {SCRIPT_DIR}/filter_nonRelevant_SNPs.py \
            -v {input.var} \
            -g {input.gt} \
            -d {input.dbsnp}
        """


rule run_COMPASS:
    input:
        os.path.join(OUT_DIR, 'samples',
            'c{cells}.r{run}{rel}.filtered_variants.csv')
    output:
        os.path.join(OUT_DIR, 'COMPASS',
            'c{cells}.r{run}{rel}.d{dbt}_cellAssignments.tsv')
    wildcard_constraints:
        cells = r'\d+',
        run = r'\d+',
        rel = r'(\.relevant|)',
        dbt = r'[\d\.]+'
    threads: config['analysis']['algorithms']['COMPASS'].get('chains', 4)
    resources:
        cpus_per_task = config['analysis']['algorithms']['COMPASS'] \
            .get('chains', 4),
        mem_mb_per_cpu = 10240,
        runtime = 1440,
    params:
        exe = config['analysis']['algorithms']['COMPASS']['exe'],
        in_files = lambda w: os.path.join(OUT_DIR, 'samples',
            f'c{w.cells}.r{w.run}' + w.rel + f'.filtered'),
        out_files = lambda w: os.path.join(OUT_DIR, 'COMPASS',
            f'c{w.cells}.r{w.run}' + w.rel + f'.d{w.dbt}'),
        CNV = config['analysis']['algorithms']['COMPASS'].get('CNV', 0),
        chain_length = config['analysis']['algorithms']['COMPASS'] \
            .get('chain_length', 20000),
        sex = SEX
    shell:
        """
        {params.exe} \
            -i {params.in_files} \
            -o {params.out_files} \
            --CNV {params.CNV} \
            --nchains {resources.cpus_per_task} \
            --chainlength {params.chain_length} \
            --doubletrate {wildcards.dbt} \
            --sex {params.sex}
        """


rule run_BnpC:
    input:
        os.path.join(OUT_DIR, 'samples',
            'c{cells}.r{run}{rel}.filtered_variants_gt.csv')
    output:
        os.path.join(OUT_DIR, 'BnpC',
            'c{cells}.r{run}{rel}', 'assignment.txt')
    wildcard_constraints:
        cells = r'\d+',
        run = r'\d+',
        rel = r'(\.relevant|)'
    conda: os.path.join(ENV_DIR, 'BnpC.yaml')
    threads: config['analysis']['algorithms']['BnpC'].get('chains', 4)
    resources:
        cpus_per_task = config['analysis']['algorithms']['BnpC'] \
            .get('chains', 4),
        mem_mb_per_cpu = 10240,
        runtime = config['analysis']['algorithms']['BnpC'].get('runtime', 90) \
            + 120,
    params:
        exe = config['analysis']['algorithms']['BnpC']['exe'],
        out_dir = lambda w: os.path.join(OUT_DIR, 'BnpC', 
            f'c{w.cells}.r{w.run}' + w.rel),
        runtime = config['analysis']['algorithms']['BnpC'].get('runtime', 90),
        pp = ' '.join([str(i) for i in \
            config['analysis']['algorithms']['BnpC'].get('params_prior', [1, 1])]),
        ap = ' '.join([str(i) for i in \
            config['analysis']['algorithms']['BnpC'].get('dpa_prior', [1, 1])]),
        cup = config['analysis']['algorithms']['BnpC'].get('conc_update_prob', 0)
    shell:
        """
        python {params.exe} \
            {input} \
            -o {params.out_dir} \
            -r {params.runtime} \
            -e posterior MAP \
            -n {resources.cpus_per_task} \
            --param_prior {params.pp} \
            --DPa_prior {params.ap} \
            --conc_update_prob {params.cup} \
            -np
        """


rule run_SCG:
    input:
        os.path.join(OUT_DIR, 'samples',
            'c{cells}.r{run}{rel}.filtered_variants_gt.csv')
    output:
        os.path.join(OUT_DIR, 'SCG',
            'c{cells}.r{run}{rel}.d{dbt}', 'assignments.tsv')
    wildcard_constraints:
        cells = r'\d+',
        run = r'\d+',
        rel = r'(\.relevant|)',
        dbt = r'[\d\.]+'
    conda: os.path.join(ENV_DIR, 'scg.yaml')
    threads: 1
    resources:
        mem_mb_per_cpu = 10240,
        runtime = 1440,
    params:
        out_dir = lambda w: os.path.join(OUT_DIR, 'SCG', 
            f'c{w.cells}.r{w.run}' + w.rel + f'.d{w.dbt}'),
        wrapper_exe =  config['analysis']['algorithms']['SCG']['wrapper_exe'],
        chain_length =  config['analysis']['algorithms']['SCG'] \
            .get('chain_length', 100000),
        max_cl = config['analysis']['algorithms']['SCG'].get('clusters', 10)
    shell:
        """
        python {params.wrapper_exe} \
            {input} \
            -o {params.out_dir} \
            -s {params.chain_length} \
            -k {params.max_cl} \
            --doublets {wildcards.dbt} \
            --silent
        """


def get_algorithm_files(wildcards):
    res_files = []
    for cells in CELL_NO:
        for run in  range(RUNS):
            for rel in [r'.relevant', '']:
                out_base = f'c{cells:04d}.r{run:04d}' + f'{rel}'
                if 'BnpC' in ALGORITHMS:
                    res_files.append(os.path.join(OUT_DIR, 'BnpC',
                        out_base, 'assignment.txt'))
                for dbt in DBTS:
                    if 'COMPASS' in ALGORITHMS:
                        res_files.append(os.path.join(OUT_DIR, 'COMPASS',
                            f'{out_base}.d{dbt:0<4}_cellAssignments.tsv'))
                    if 'SCG' in ALGORITHMS:
                        res_files.append(os.path.join(OUT_DIR, 'SCG',
                            f'{out_base}.d{dbt:0<4}', 'assignments.tsv'))
    return res_files


rule compare_runs:
    input:
        get_algorithm_files
    output:
        os.path.join(OUT_DIR, 'summary_sampling.tsv')
    wildcard_constraints:
        cells = r'\d+',
        run = r'\d+',
        rel = r'(\.relevant|)',
        dbt = r'[\d\.]+'
    conda: os.path.join(ENV_DIR, 'analysis.yaml')
    threads: 1
    resources:
        mem_mb_per_cpu = 32768,
        runtime = 120,
    params:
        min_cells = [1, 25, 50]
    shell:
        """
        python {SCRIPT_DIR}/evaluate_sampling.py \
            -i {input} \
            -o {output} \
            -m {params.min_cells}
        """
