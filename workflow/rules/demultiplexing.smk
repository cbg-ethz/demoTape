rule demultiplexing_processing:
    input:
        IN_FILE_ABS
    output:
        os.path.join(OUT_DIR, f'{IN_FILE_NAME}_variants.csv')
    threads: 1
    resources:
        mem_mb = 16384,
        runtime = 240,
    params:
        minGQ = config.get('mosaic', {}).get('minGQ', 30),
        minDP = config.get('mosaic', {}).get('minDP', 10),
        minVAF = config.get('mosaic', {}).get('minaVAF', 0.2),
        minVarGeno = config.get('mosaic', {}).get('minVarGeno', 0.5),
        minCellGeno = config.get('mosaic', {}).get('minCellGeno', 0.5),
        minMutated = config.get('mosaic', {}).get('minMutated', 0.01),
        max_ref_VAF = config.get('mosaic', {}).get('max_ref_VAF', 0.05),
        min_hom_VAF = config.get('mosaic', {}).get('min_hom_VAF', 0.95),
        min_het_VAF = config.get('mosaic', {}).get('min_het_VAF', 0.35),
        proximity = ' '.join([str(i) for i in \
            config.get('mosaic', {}).get('proximity', [50, 100, 150, 200])])
    conda:
        '../envs/mosaic.yml'
    shell:
        """
        python {SCRIPT_DIR}/mosaic_preprocessing.py \
            -i {input} \
            -o {output} \
            --minGQ {params.minGQ} \
            --minDP {params.minDP} \
            --minVAF {params.minVAF} \
            --minVarGeno {params.minVarGeno} \
            --minCellGeno {params.minCellGeno} \
            --minMutated {params.minMutated} \
            --max_ref_VAF {params.max_ref_VAF} \
            --min_hom_VAF {params.min_hom_VAF} \
            --min_het_VAF {params.min_het_VAF} \
            --proximity {params.proximity} \
            --full_output
        """


rule demultiplexing_demoTape:
    input:
        os.path.join(OUT_DIR, f'{IN_FILE_NAME}_variants.csv')
    output:
        os.path.join(OUT_DIR, f'{IN_FILE_NAME}_variants.assignments.tsv')
    threads: 1
    resources:
        mem_mb = 4096,
        runtime = 90,
    conda:
        '../envs/sample_assignment.yml'
    shell:
        """
        python {SCRIPT_DIR}/demultiplex_distance.py \
            -i {input} \
            -n {SAMPLES} 
        """

rule demultiplexing_reanalyze_whitelist:
    input:
        loom = IN_FILE_ABS,
        assignment = os.path.join(OUT_DIR, f'{IN_FILE_NAME}_variants.assignments.tsv'),
        whitelist = os.path.join(SNP_DIR, f'{SAMPLE_NAME}.snps.stripped.tsv')
    output:
        [temp(os.path.join(OUT_DIR,'temp', f'{IN_FILE_NAME}.{i}_variants.csv')) \
            for i in range(SAMPLES)],
    threads: 1
    resources:
        mem_mb = 16384,
        runtime = 240,
    params:
        out_file_raw = os.path.join(OUT_DIR, 'temp', f'{IN_FILE_NAME}_variants.csv'),
        minGQ = config.get('mosaic', {}).get('minGQ', 30),
        minDP = config.get('mosaic', {}).get('minDP', 10),
        minVAF = config.get('mosaic', {}).get('minaVAF', 0.2),
        minVarGeno = config.get('mosaic', {}).get('minVarGeno', 0.5),
        minCellGeno = config.get('mosaic', {}).get('minCellGeno', 0.5),
        minMutated = config.get('mosaic', {}).get('minMutated', 50), # Set to 50, instead of 0.01, for smaller datasets!
        max_ref_VAF = config.get('mosaic', {}).get('max_ref_VAF', 0.05),
        min_hom_VAF = config.get('mosaic', {}).get('min_hom_VAF', 0.95),
        min_het_VAF = config.get('mosaic', {}).get('min_het_VAF', 0.35),
        proximity = ' '.join([str(i) for i in \
            config.get('mosaic', {}).get('proximity', [50, 100, 150, 200])])
    conda:
        '../envs/mosaic.yml'
    shell:
        """
        python {SCRIPT_DIR}/mosaic_preprocessing.py \
            -i {input.loom} \
            -a {input.assignment} \
            -o {params.out_file_raw} \
            --minGQ {params.minGQ} \
            --minDP {params.minDP} \
            --minVAF {params.minVAF} \
            --minVarGeno {params.minVarGeno} \
            --minCellGeno {params.minCellGeno} \
            --minMutated {params.minMutated} \
            --max_ref_VAF {params.max_ref_VAF} \
            --min_hom_VAF {params.min_hom_VAF} \
            --min_het_VAF {params.min_het_VAF} \
            --proximity {params.proximity} \
            --full_output \
            --prepare_cnv_file \
            --panel_annotation /cluster/work/bewi/members/jgawron/INTeRCePT/tapestri/panels/Myeloid/Myeloid_amplicon_annotation.csv \
            --whitelist {input.whitelist}
        """

rule demultiplex_sample_assignment:
    input:
        [os.path.join(OUT_DIR,'temp', f'{IN_FILE_NAME}.{i}_variants.csv') \
            for i in range(SAMPLES)],
    output:
        heatmap = os.path.join(OUT_DIR,f'{SAMPLE_NAME}_variants_distance_heatmap.png'),
        yaml = os.path.join(OUT_DIR,f'{SAMPLE_NAME}.assigment.yaml')
    threads: 1
    resources:
        mem_mb = 4096,
        runtime = 30,
    params:
       RNA_SNP_profile_dir = SNP_DIR,
       output_base = os.path.join(OUT_DIR,f'{SAMPLE_NAME}_variants_distance.png')
    conda:
        '../envs/sample_assignment.yml'
    shell:
        """
        python {SCRIPT_DIR}/assign_cells_via_scRNAseq.py \
        -i {input} \
        -p $(ls {params.RNA_SNP_profile_dir}/*.final.vcf) \
        -o {params.output_base}
        --assigment-output {output.yaml}
        """

rule demultiplexing_splitting:
    input:
        loom = IN_FILE_ABS,
        assignment = os.path.join(OUT_DIR, f'{IN_FILE_NAME}_variants.assignments.tsv')
    output:
        [os.path.join(OUT_DIR, f'{IN_FILE_NAME}.{i}_variants.csv') \
            for i in range(SAMPLES)],
        [os.path.join(OUT_DIR, f'{IN_FILE_NAME}.{i}_variants_gt.csv') \
            for i in range(SAMPLES)]
    threads: 1
    resources:
        mem_mb = 16384,
        runtime = 240,
    params:
        out_file_raw = os.path.join(OUT_DIR, f'{IN_FILE_NAME}_variants.csv'),
        minGQ = config.get('mosaic', {}).get('minGQ', 30),
        minDP = config.get('mosaic', {}).get('minDP', 10),
        minVAF = config.get('mosaic', {}).get('minaVAF', 0.2),
        minVarGeno = config.get('mosaic', {}).get('minVarGeno', 0.5),
        minCellGeno = config.get('mosaic', {}).get('minCellGeno', 0.5),
        minMutated = config.get('mosaic', {}).get('minMutated', 50), # Set to 50, instead of 0.01, for smaller datasets!
        max_ref_VAF = config.get('mosaic', {}).get('max_ref_VAF', 0.05),
        min_hom_VAF = config.get('mosaic', {}).get('min_hom_VAF', 0.95),
        min_het_VAF = config.get('mosaic', {}).get('min_het_VAF', 0.35),
        proximity = ' '.join([str(i) for i in \
            config.get('mosaic', {}).get('proximity', [50, 100, 150, 200])])
    conda:
        '../envs/mosaic.yml'
    shell:
        """
        python {SCRIPT_DIR}/mosaic_preprocessing.py \
            -i {input.loom} \
            -a {input.assignment} \
            -o {params.out_file_raw} \
            --minGQ {params.minGQ} \
            --minDP {params.minDP} \
            --minVAF {params.minVAF} \
            --minVarGeno {params.minVarGeno} \
            --minCellGeno {params.minCellGeno} \
            --minMutated {params.minMutated} \
            --max_ref_VAF {params.max_ref_VAF} \
            --min_hom_VAF {params.min_hom_VAF} \
            --min_het_VAF {params.min_het_VAF} \
            --proximity {params.proximity} \
            --full_output \
            --prepare_cnv_file \
            --panel_annotation /cluster/work/bewi/members/jgawron/INTeRCePT/tapestri/panels/Myeloid/Myeloid_amplicon_annotation.csv
        """


rule run_COMPASS:
    input:
        variants = os.path.join(OUT_DIR, f'{IN_FILE_NAME}.{{cl}}_variants.csv'),
        report = os.path.join(INPUT_DIRECTORY,f'{SAMPLE_NAME}.report.json')
    output:
        os.path.join(OUT_DIR, 'COMPASS', 'cl{cl}', 'r{run}.d{dbt}_cellAssignments.tsv')
    threads: 4
    resources:
        mem_mb = 10240,
        runtime = 1440,
    log: os.path.join(INPUT_DIRECTORY, 'logs', 'r{run}.cl{cl}.d{dbt}_COMPASS.log')
    params:
        in_file_raw = lambda w: os.path.join(OUT_DIR, f'{IN_FILE_NAME}.{w.cl}'),
        out_file_raw = lambda w: os.path.join(
            OUT_DIR, 'COMPASS', f'cl{w.cl}', f'r{w.run}.d{w.dbt}'),
        exe = config['COMPASS']['exe'],
        chain_length = config['COMPASS'].get('chain_length', 20000),
        CNV = config['COMPASS'].get('CNV', 0),
        sex = config['COMPASS'].get('sex', 'female'),
    shell:
        """
        ( ado=$(jq '.variant_calling.ado_rate' {input.report}); \
        {params.exe} \
            -i {params.in_file_raw} \
            -o {params.out_file_raw} \
            --CNV {params.CNV} \
            --nchains {threads} \
            --doubletrate {wildcards.dbt} \
            --chainlength {params.chain_length} \
            --sex {params.sex} \
            --dropoutrate $ado ) &> {log}
        """


rule run_BnpC:
    input:
        os.path.join(OUT_DIR, f'{IN_FILE_NAME}.{{cl}}_variants_gt.csv')
    output:
        os.path.join(OUT_DIR, 'BnpC', 'cl{cl}', 'r{run}', 'assignment.txt')
    threads: 4
    resources:
        mem_mb = 10240,
        runtime = 630,
    params:
        out_dir = lambda w: os.path.join(
            OUT_DIR, 'BnpC', f'cl{w.cl}', f'r{w.run}'),
        exe = config['BnpC']['exe'],
        chain_length = config['BnpC'].get('chain_length', 20000),
        runtime = f'-r {config["BnpC"]["runtime"]} ' \
            if config['BnpC'].get('runtime', False) else '',
        pp = ' '.join([str(i) for i in config['BnpC'].get('params_prior', [1, 1])])
    shell:
        """
        python {params.exe} \
            {input} \
            -o {params.out_dir} \
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
        os.path.join(OUT_DIR, f'{IN_FILE_NAME}.{{cl}}_variants_gt.csv')
    output:
        os.path.join(OUT_DIR, 'SCG', 'cl{cl}', 'r{run}.d{dbt}', 'assignments.tsv')
    threads: 1
    resources:
        mem_mb = 10240,
        runtime = 1440,
    params:
        out_dir = lambda w: os.path.join(
            OUT_DIR, 'SCG', f'cl{w.cl}', f'r{w.run}.d{w.dbt}'),
        wrapper_exe = config['SCG']['wrapper_exe'],
        chain_length = config['SCG'].get('chain_length', 100000),
        max_cl = 20,
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


def get_assignments(*args):
    all_files = []
    for cl in range(SAMPLES):
        for run in RUNS:
            all_files.append(os.path.join(OUT_DIR, 'BnpC', f'cl{cl}', f'r{run}', 'assignment.txt'))
            for dbt in DOUBLET_RATE:
                all_files.append(
                    os.path.join(OUT_DIR, 'SCG', f'cl{cl}' , f'r{run}.d{dbt}', 'assignments.tsv'))
                all_files.append(
                    os.path.join(OUT_DIR, 'COMPASS', f'cl{cl}', f'r{run}.d{dbt}_cellAssignments.tsv'))
    return all_files


rule generate_cooclusterin_matrix:
    input:
        get_assignments
    output:
        expand([os.path.join(OUT_DIR, 'cl{cl}_coclustering_COMPASS_d{dbt}.png'),
            os.path.join(OUT_DIR, 'cl{cl}_coclustering_SCG_d{dbt}.png')],
        cl=range(SAMPLES), dbt=DOUBLET_RATE),
        expand(os.path.join(OUT_DIR, 'cl{cl}_coclustering_BnpC_d-1.png'),
            cl=range(SAMPLES))
    threads: 1
    resources:
        mem_mb = 32768,
        runtime = 60,
    params:
        out_dir = OUT_DIR
    script:
        f'{SCRIPT_DIR}/plot_cooccurence_matrix.py'


