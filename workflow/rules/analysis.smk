rule mosaic_full:
    input:
        loom = ancient(LOOM_FILE)
    output:
        var = os.path.join(OUT_DIR, f'{DNA_PREFIX}.filtered_variants.csv'),
        gt = os.path.join(OUT_DIR, f'{DNA_PREFIX}.filtered_variants_gt.csv'),
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
        minMutated = config.get('mosaic', {}).get('minMutated', 0.05), # default: 0.01
        maxRefVAF = config.get('mosaic', {}).get('maxRefVAF', 0.05),
        minHomVAF = config.get('mosaic', {}).get('minHomVAF', 0.95),
        minHetVAF = config.get('mosaic', {}).get('minHetVAF', 0.35),
        proximity = config.get('mosaic', {}).get('proximity', '25 50 100 200')
    shell:
        """
        python {SCRIPT_DIR}/mosaic_preprocessing.py \
            -i {input.loom} \
            -o {output.var} \
            --minGQ {params.minGQ} \
            --minDP {params.minDP} \
            --minVAF {params.minVAF} \
            --minVarGeno {params.minVarGeno} \
            --minCellGeno {params.minCellGeno} \
            --minMutated {params.minMutated} \
            --max_ref_VAF {params.maxRefVAF} \
            --min_hom_VAF {params.minHomVAF} \
            --min_het_VAF {params.minHetVAF} \
            --proximity {params.proximity}
        """

# ------------------------- DEMULTIPLEXING -------------------------------------
if POOLED:
    rule demoTape_demultiplexing:
        input:
            os.path.join(OUT_DIR, f'{DNA_PREFIX}.filtered_variants.csv'),
        output:
            os.path.join(OUT_DIR, f'{DNA_PREFIX}.demoTape.assignments.tsv'),
            os.path.join(OUT_DIR, f'{DNA_PREFIX}.demoTape.heatmap.png')
        conda: os.path.join(ENV_DIR, 'analysis.yaml')
        threads: 1
        resources:
            mem_mb_per_cpu = 4096,
            runtime = 90,
        params:
            out_base = os.path.join(OUT_DIR, f'{DNA_PREFIX}.demoTape')
        shell:
            """
            python {SCRIPT_DIR}/demultiplex_distance.py \
                -i {input} \
                -o {params.out_base} \
                -n {POOLED_SAMPLES} \
                -op
            """


    rule mosaic_demultiplexed:
        input:
            loom = ancient(LOOM_FILE),
            var_all = os.path.join(OUT_DIR,
                f'{DNA_PREFIX}.filtered_variants.csv'),
            assignment = os.path.join(OUT_DIR,
                f'{DNA_PREFIX}.demoTape.assignments.tsv')
        output:
            expand([os.path.join(OUT_DIR, '{prefix}.filtered_variants.csv'),
                os.path.join(OUT_DIR, '{prefix}.filtered_variants_gt.csv')],
                prefix=PREFIXES[1:])
        wildcard_constraints:
            prefix = rf'{DNA_PREFIX}_\d(?=\.)'
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
            minMutated = 50, # Set to 50 for smaller datasets!
            maxRefVAF = config.get('mosaic', {}).get('maxRefVAF', 0.05),
            minHomVAF = config.get('mosaic', {}).get('minHomVAF', 0.95),
            minHetVAF = config.get('mosaic', {}).get('minHetVAF', 0.35),
            proximity = config.get('mosaic', {}).get('proximity', '25 50 100 200')
        shell:
            """
            python {SCRIPT_DIR}/mosaic_preprocessing.py \
                -i {input.loom} \
                -o {input.var_all} \
                -a {input.assignment} \
                --minGQ {params.minGQ} \
                --minDP {params.minDP} \
                --minVAF {params.minVAF} \
                --minVarGeno {params.minVarGeno} \
                --minCellGeno {params.minCellGeno} \
                --minMutated {params.minMutated} \
                --max_ref_VAF {params.maxRefVAF} \
                --min_hom_VAF {params.minHomVAF} \
                --min_het_VAF {params.minHetVAF} \
                --proximity {params.proximity}
            """


    if MATCHED_RNA:
        rule generate_loci_whitelist:
            input:
                os.path.join(OUT_DIR, f'{DNA_PREFIX}.filtered_variants.csv'),
                expand(os.path.join(OUT_DIR, '{prefix}.filtered_variants.csv'),
                    prefix=PREFIXES[1:])
            output:
                os.path.join(OUT_DIR, f'{DNA_PREFIX}.loci_whitelist.csv')
            threads: 1
            resources:
                mem_mb_per_cpu = 4096,
                runtime = 30
            run:
                import pandas as pd

                CHR_ORDER = {str(i): i for i in range(1, 23)}
                CHR_ORDER.update({'X':23, 'Y':24})

                loci_all = []
                for in_file in input:
                    df = pd.read_csv(in_file)
                    loci_all.extend(
                        df.apply(lambda x: f'{x["CHR"]}_{x["POS"]}', axis=1).to_list()
                    )
                loci_uniq = list(set(loci_all))
                loci_sorted = sorted(loci_uniq,
                    key=lambda x: CHR_ORDER[x.split('_')[0]])

                with open(output[0], 'w') as f:
                    f.write('CHR,POS\n')
                    for i in loci_sorted:
                        f.write(i.replace('_', ',') + '\n')


        rule mosaic_whitelist:
            input:
                loom = ancient(LOOM_FILE),
                var_all = os.path.join(OUT_DIR,
                    f'{DNA_PREFIX}.filtered_variants.csv'),
                assignment = os.path.join(OUT_DIR,
                    f'{DNA_PREFIX}.demoTape.assignments.tsv'),
                whitelist = os.path.join(OUT_DIR,
                    f'{DNA_PREFIX}.loci_whitelist.csv')
            output:
                expand(os.path.join(OUT_DIR, 
                        '{prefix}.whitelist.filtered_variants.csv'),
                    prefix=PREFIXES[1:])
            wildcard_constraints:
                prefix = rf'{DNA_PREFIX}_\d(?=\.)'
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
                minMutated = 50, # Set to 50 for smaller datasets!
                maxRefVAF = config.get('mosaic', {}).get('maxRefVAF', 0.05),
                minHomVAF = config.get('mosaic', {}).get('minHomVAF', 0.95),
                minHetVAF = config.get('mosaic', {}).get('minHetVAF', 0.35),
                proximity = config.get('mosaic', {}).get('proximity', '25 50 100 200')
            shell:
                """
                python {SCRIPT_DIR}/mosaic_preprocessing.py \
                    -i {input.loom} \
                    -o {input.var_all} \
                    -a {input.assignment} \
                    -wl {input.whitelist} \
                    --minGQ {params.minGQ} \
                    --minDP {params.minDP} \
                    --minVAF {params.minVAF} \
                    --minVarGeno {params.minVarGeno} \
                    --minCellGeno {params.minCellGeno} \
                    --minMutated {params.minMutated} \
                    --max_ref_VAF {params.maxRefVAF} \
                    --min_hom_VAF {params.minHomVAF} \
                    --min_het_VAF {params.minHetVAF} \
                    --proximity {params.proximity}
                """


        rule assign_sample_to_patient:
            input:
                var = expand(os.path.join(OUT_DIR,
                        '{prefix}.whitelist.filtered_variants.csv'),
                    prefix=PREFIXES[1:]),
                rna = RNA_SNPS
            output:
                assignment = os.path.join(OUT_DIR,
                    f'{DNA_PREFIX}.sample_patient_assigments.tsv'),
                profiles = os.path.join(OUT_DIR,
                    f'{DNA_PREFIX}.sample_patient_profiles.png'),
                distance = os.path.join(OUT_DIR,
                    f'{DNA_PREFIX}.sample_patient_distance.png'),
            conda: os.path.join(ENV_DIR, 'analysis.yaml')
            threads: 1
            resources:
                mem_mb_per_cpu = 4096,
                runtime = 30,
            shell:
                """
                python {SCRIPT_DIR}/assign_sample_to_patient.py \
                -i {input.var} \
                -r {input.rna} \
                -oa {output.assignment} \
                -op {output.profiles} \
                -od {output.distance}
                """

else:
    rule mock_demultiplexing_demoTape:
        input:
            os.path.join(OUT_DIR, f'{DNA_PREFIX}.filtered_variants.csv'),
        output:
            os.path.join(OUT_DIR, f'{DNA_PREFIX}.demoTape.assignments.tsv')
        threads: 1
        resources:
            mem_mb_per_cpu = 4096,
            runtime = 15,
        run:
            with open(input[0], 'r') as f:
                cells = f.readline().strip().split(',')[7:]

            cell_str = '\t'.join(cells)
            cl_str = '\t'.join(['0'] * len(cells))
            order_str = '\t'.join([str(i) for i in range(len(cells))])
            with open(output[0], 'w') as f:
                f.write(f'Barcode\t{cell_str}\nCluster\t{cl_str}\nOrder\t{order_str}')  


# ------------------------------------------------------------------------------


rule plot_VAF_hm:
    input:
        var = os.path.join(OUT_DIR, '{prefix}.filtered_variants.csv'),
    output:
        hm = os.path.join(OUT_DIR, '{prefix}.SNP_heatmap.png')
    wildcard_constraints:
        prefix = rf'{DNA_PREFIX}(_\d)*(?=\.)'
    conda: os.path.join(ENV_DIR, 'analysis.yaml')
    threads: 1
    resources:
        mem_mb_per_cpu = 16384,
        runtime = 30,
    params:
        cell_annotation = f'-ca {CELL_ANNOT}' if CELL_ANNOT else '',
    shell:
        """
        python {SCRIPT_DIR}/plot_VAF_heatmap.py \
            -i {input.var} \
            -o {output.hm} \
            {params.cell_annotation}
        """


rule plot_read_hm:
    input:
        reads = READS_FILE,
        panel_annotated = PANEL_INFO,
        sample_annotation = os.path.join(OUT_DIR, 
            f'{DNA_PREFIX}.demoTape.assignments.tsv')
    output:
        png = os.path.join(OUT_DIR, '{prefix}.reads_heatmap.png'),
        cnv = os.path.join(OUT_DIR, '{prefix}.reads_CNVclones.csv')
    wildcard_constraints:
        prefix = rf'{DNA_PREFIX}(_\d)*(?=\.)'
    conda: os.path.join(ENV_DIR, 'analysis.yaml')
    threads: 1
    resources:
        mem_mb_per_cpu = 16384,
        runtime = 30,
    params:
        cell_annotation = f'-ca {CELL_ANNOT}' if CELL_ANNOT else '',
    shell:
        """
        python {SCRIPT_DIR}/plot_read_heatmap.py \
            -i {input.reads} \
            -o {output.png} \
            -p {input.panel_annotated} \
            -sa {input.sample_annotation} \
            {params.cell_annotation}
        """


rule filter_nonRelevant_SNPs:
    input:
        var = os.path.join(OUT_DIR, '{prefix}.filtered_variants.csv'),
        gt = os.path.join(OUT_DIR, '{prefix}.filtered_variants_gt.csv'),
        dbsnp = DBSNP
    output:
        var = os.path.join(OUT_DIR, '{prefix}.relevant.filtered_variants.csv'),
        gt = os.path.join(OUT_DIR, '{prefix}.relevant.filtered_variants_gt.csv'),
        stats = os.path.join(OUT_DIR, '{prefix}.relevant.statistics.csv')
    wildcard_constraints:
        prefix = rf'{DNA_PREFIX}(_\d)*(?=\.)'
    conda: os.path.join(ENV_DIR, 'analysis.yaml')
    threads: 1
    resources:
        mem_mb_per_cpu = 8096,
        runtime = 60,
    params:
        cell_annotation = f'-ca {CELL_ANNOT}' if CELL_ANNOT else '',
        clinVar = f'-c {CLINVAR}' if CLINVAR else '',
    shell:
        """
        python {SCRIPT_DIR}/filter_nonRelevant_SNPs.py \
            -v {input.var} \
            -g {input.gt} \
            -d {input.dbsnp} \
            {params.cell_annotation} \
            {params.clinVar}
        """


rule plot_VAF_hm_relevant:
    input:
        var = os.path.join(OUT_DIR, '{prefix}.relevant.filtered_variants.csv'),
    output:
        hm = os.path.join(OUT_DIR, '{prefix}.relevant.SNP_heatmap.png')
    wildcard_constraints:
        prefix = rf'{DNA_PREFIX}(_\d)*(?=\.)'
    conda: os.path.join(ENV_DIR, 'analysis.yaml')
    threads: 1
    resources:
        mem_mb_per_cpu = 16384,
        runtime = 60,
    params:
        cell_annotation = f'-ca {CELL_ANNOT}' if CELL_ANNOT else '',
    shell:
        """
        python {SCRIPT_DIR}/plot_VAF_heatmap.py \
            -i {input.var} \
            -o {output.hm} \
            {params.cell_annotation}
        """


# # ------------------------------------------------------------------------------

def get_compass_in(wildcards):
    if wildcards.snp_set == 'relevant':
        var_file = os.path.join(OUT_DIR, 
            f'{wildcards.prefix}.relevant.filtered_variants.csv')
    else:
        var_file = os.path.join(OUT_DIR,
            f'{wildcards.prefix}.filtered_variants.csv')
    return {'var': var_file}


rule prepare_COMPASS_regions_files:
    input:
        unpack(get_compass_in),
        reads = READS_FILE,
        gene_annotation = GENE_INFO,
        panel_annotated = PANEL_INFO
    output:
        expand(os.path.join(OUT_DIR,
                '{{prefix}}.{{snp_set}}.{region}.filtered_variants.csv'),
            region=REGIONS),
        expand(os.path.join(OUT_DIR,
                '{{prefix}}.{{snp_set}}.{region}.filtered_regions.csv'),
            region=REGIONS)
    wildcard_constraints:
        prefix = rf'{DNA_PREFIX}(_\d)*(?=\.)'
    conda: os.path.join(ENV_DIR, 'analysis.yaml')
    threads: 1
    resources:
        mem_mb_per_cpu = 8096,
        runtime = 30,
    params:
        out_base = lambda w: os.path.join(OUT_DIR, f'{w.prefix}.{w.snp_set}'),
        cell_annotation = f'-ca {CELL_ANNOT}' if CELL_ANNOT else '',
    shell:
        """
        python {SCRIPT_DIR}/prepare_COMPASS_regions_file.py \
            -v {input.var} \
            -d {input.reads} \
            -ga {input.gene_annotation} \
            -p {input.panel_annotated} \
            -o {params.out_base} \
            {params.cell_annotation} \
        """


rule run_COMPASS:
    input:
        var = os.path.join(OUT_DIR,
            '{prefix}.{snp_set}.{region}.filtered_variants.csv'),
        reg = os.path.join(OUT_DIR,
            '{prefix}.{snp_set}.{region}.filtered_regions.csv')
    output:
        os.path.join(OUT_DIR, 'COMPASS',
            '{prefix}.{snp_set}.r{run}.d{dbt}.{region}_tree.gv')
    wildcard_constraints:
        prefix = rf'{DNA_PREFIX}(_\d)*(?=\.)',
        run = r'\d+',
        dbt = r'[\d\.]+',
    threads: config['analysis']['algorithms']['COMPASS'].get('chains', 4)
    resources:
        mem_mb_per_cpu = 10240,
        runtime = 900,
    params:
        exe = config['analysis']['algorithms']['COMPASS']['exe'],
        in_files = lambda w: os.path.join(OUT_DIR,
            f'{w.prefix}.{w.snp_set}.{w.region}.filtered'),
        out_files = lambda w: os.path.join(OUT_DIR, 'COMPASS',
            f'{w.prefix}.{w.snp_set}.r{w.run}.d{w.dbt}.{w.region}'),
        CNV = config['analysis']['algorithms']['COMPASS'].get('CNV', 1),
        chain_length = config['analysis']['algorithms']['COMPASS'] \
            .get('chain_length', 20000),
        sex = SEX
    shell:
        """
        {params.exe} \
            -i {params.in_files} \
            -o {params.out_files} \
            --CNV {params.CNV} \
            --nchains {threads} \
            --chainlength {params.chain_length} \
            --doubletrate {wildcards.dbt} \
            --sex {params.sex}
        """


rule plot_COMPASS_tree_simple:
    input:
        os.path.join(OUT_DIR, 'COMPASS',
            '{prefix}.{snp_set}.r{run}.d{dbt}.{region}_tree.gv')
    output:
        os.path.join(OUT_DIR, 'COMPASS',
            '{prefix}.{snp_set}.r{run}.d{dbt}.{region}_tree.png')
    wildcard_constraints:
        prefix = rf'{DNA_PREFIX}(_\d)*(?=\.)',
        run = r'\d+',
        dbt = r'[\d\.]+',
    threads: 1
#    retries: 3
    resources:
        mem_mb_per_cpu = 4096,
        runtime = 5
    shell:
        """
        python {SCRIPT_DIR}/../../sandbox/render_trees.py \
            -i {input}
        """


def get_clustering_gt_in(wildcards):
    if wildcards.snp_set == 'relevant':
        return os.path.join(OUT_DIR,
            f'{wildcards.prefix}.relevant.filtered_variants_gt.csv')
    return os.path.join(OUT_DIR,
        f'{wildcards.prefix}.filtered_variants_gt.csv')


rule run_BnpC:
    input:
        get_clustering_gt_in
    output:
        os.path.join(OUT_DIR, 'BnpC', '{prefix}.{snp_set}.r{run}',
            'assignment.txt')
    wildcard_constraints:
        prefix = rf'{DNA_PREFIX}(_\d)*(?=\.)',
        run = r'\d+',
    conda: os.path.join(ENV_DIR, 'BnpC.yaml')
    threads: config['analysis']['algorithms']['BnpC'].get('chains', 4)
    resources:
        mem_mb_per_cpu = 10240,
        runtime = 630,
    params:
        out_dir = lambda w: os.path.join(OUT_DIR, 'BnpC', 
            f'{w.prefix}.{w.snp_set}.r{w.run}'),
        exe = config['analysis']['algorithms']['BnpC']['exe'],
        chain_length = config['analysis']['algorithms']['BnpC'] \
            .get('chain_length', 20000),
        runtime = f'-r {config["analysis"]["algorithms"]["BnpC"]["runtime"]} ' \
            if config['analysis']['algorithms']['BnpC'].get('runtime', False) \
                else '',
        pp = ' '.join([str(i) for i in \
            config['analysis']['algorithms']['BnpC'].get('params_prior', [1, 1])])
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
        get_clustering_gt_in
    output:
        os.path.join(OUT_DIR, 'SCG', '{prefix}.{snp_set}.r{run}.d{dbt}',
            'assignments.tsv')
    wildcard_constraints:
        prefix = rf'{DNA_PREFIX}(_\d)*(?=\.)',
        run = r'\d+',
        dbt = r'[\d\.]+'
    conda: os.path.join(ENV_DIR, 'scg.yaml')
    threads: 1
    resources:
        mem_mb_per_cpu = 10240,
        runtime = 1440,
    params:
        out_dir = lambda w: os.path.join(OUT_DIR, 'SCG', 
            f'{w.prefix}.{w.snp_set}.r{w.run}.d{w.dbt}'),
        wrapper_exe =  config['analysis']['algorithms']['SCG']['wrapper_exe'],
        chain_length =  config['analysis']['algorithms']['SCG'] \
            .get('chain_length', 100000),
        max_cl = config['analysis']['algorithms']['SCG'].get('clusters', 10),
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


# ------------------------------------------------------------------------------
if RNA:
    rule prepare_cloneAlign:
        input:
            exp = RNA_FILE,
            cnv = os.path.join(OUT_DIR, f'{DNA_PREFIX}.reads_CNVclones.csv'),
            gene_annotation = GENE_INFO
        output:
            expand(os.path.join(OUT_DIR,
                    f'{DNA_PREFIX}.{{region}}.cloneAlign_expression.csv'),
                region=REGIONS),
            expand(os.path.join(OUT_DIR,
                    f'{DNA_PREFIX}.{{region}}.cloneAlign_CNVclones.csv'),
                region=REGIONS)
        conda: os.path.join(ENV_DIR, 'analysis.yaml')
        threads: 1
        resources:
            mem_mb_per_cpu = 32768,
            runtime = 60,
        shell:
            """
            python {SCRIPT_DIR}/prepare_cloneAlign_input.py \
                -exp {input.exp} \
                -cn {input.cnv} \
                -ga {input.gene_annotation} \
                -o {output[0]}
            """


    rule run_cloneAlign:
        input:
            exp = os.path.join(OUT_DIR,
                f'{DNA_PREFIX}.{{region}}.cloneAlign_expression.csv'),
            cnv = os.path.join(OUT_DIR,
                f'{DNA_PREFIX}.{{region}}.cloneAlign_CNVclones.csv')
        output:
            os.path.join(OUT_DIR, f'{DNA_PREFIX}.{{region}}.cloneAlign.tsv')
        log:
            os.path.join('logs', 'slurm', 'run_cloneAlign',
                f'{DNA_PREFIX}.{{region}}.cloneAlign.log')
        conda: os.path.join(CONDA_PREFIX, 'r_env')
        threads: 8
        resources:
            mem_mb_per_cpu = 32768,
            runtime = 180,
        params:
            min_counts_per_cell = lambda w: 1 if w.region == 'genes' else 100
        shell:
            """
            Rscript {SCRIPT_DIR}/run_cloneAlign.R \
                --exp {input.exp} \
                --cnv {input.cnv} \
                --mcpc {params.min_counts_per_cell} \
                --out_file {output} \
            >> {log} 2>&1
            """
