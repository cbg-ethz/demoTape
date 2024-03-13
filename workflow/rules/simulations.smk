if config.get('simulations', {}).get('run_early', False):
    rule multiplexing_simulate_early:
        input:
            ancient(os.path.join(SIM_RES_DIR, f'{PREFIX}.cells.loom'))
        output:
            os.path.join(OUT_DIR, '{frac_dir}', 'variants',
                'rep{n}.filtered_variants.csv')
        threads: 1
        resources:
            mem_mb = 32768,
            runtime = 90,
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
                config.get('mosaic', {}).get('proximity', [2, 3, 5, 10])])
        shell:
            """
            python {SCRIPT_DIR}/mosaic_preprocessing.py \
                -i {input} \
                -o {output} \
                --minGQ {params.minGQ} \
                --minDP {params.minDP} \
                --minVAF {params.minVAF} \
                --minVarGeno {params.minVarGeno},
                --minCellGeno {params.minCellGeno} \
                --minMutated {params.minMutated},
                --max_ref_VAF {params.max_ref_VAF} \
                --min_hom_VAF {params.min_hom_VAF} \
                --min_het_VAF {params.min_het_VAF} \
                --proximity {params.proximity} \
                --full_output
            """

else:
    rule multiplexing_simulate_late:
        input:
            config['general']['input-looms']
        output:
            var = os.path.join(OUT_DIR, '{frac_dir}', 'variants',
                'rep{n}.filtered_variants.csv'),
            RD = os.path.join(OUT_DIR, '{frac_dir}', 'variants',
                'rep{n}.filtered_variants_RD.csv'),
            AD = os.path.join(OUT_DIR, '{frac_dir}', 'variants',
                'rep{n}.filtered_variants_AD.csv'),
            RD_sm = os.path.join(OUT_DIR, '{frac_dir}', 'variants',
                'rep{n}.filtered_variants_RD.mtx'),
            AD_sm = os.path.join(OUT_DIR, '{frac_dir}', 'variants',
                'rep{n}.filtered_variants_AD.mtx'),
            barcodes = os.path.join(OUT_DIR, '{frac_dir}', 'variants',
                'rep{n}.filtered_variants_barcodes.csv'),
            vcf = os.path.join(OUT_DIR, '{frac_dir}', 'variants',
                'rep{n}.filtered_variants.vcf')
        threads: 1
        resources:
            mem_mb = 49152,
            runtime = 240,
        params:
            ratios = lambda w: re.search('([-0\.\d]+)-d0\.\d+', f'{w.frac_dir}') \
                .group(1).replace('-', ' '),
            doublets = lambda w: re.search('-d(0\.\d+)', f'{w.frac_dir}').group(1),
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
        shell:
            """
            python {SCRIPT_DIR}/mosaic_preprocessing.py \
                -i {input} \
                -o {output.var} \
                -r {params.ratios} \
                -d {params.doublets} \
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
        os.path.join(OUT_DIR, '{frac_dir}', 'variants',
            'rep{n}.filtered_variants.csv')
    output:
        os.path.join(OUT_DIR, '{frac_dir}', 'demoTape',
            'rep{n}.demoTape.assignments.tsv')
    threads: 1
    resources:
        mem_mb = 4096,
        runtime = 360,
    benchmark:
        os.path.join(OUT_DIR, 'benchmarks', '{frac_dir}.{n}.demoTape.txt')
    params:
        out_base = lambda w: os.path.join(OUT_DIR, f'{w.frac_dir}', 'demoTape',
            f'rep{w.n}.demoTape'),
    shell:
        """
        python {SCRIPT_DIR}/demultiplex_distance.py \
            -i {input} \
            -o {params.out_base} \
            -n {SAMPLES_N}
        """


rule demultiplexing_scSplit:
    input:
        RD = os.path.join(OUT_DIR, '{frac_dir}', 'variants',
            'rep{n}.filtered_variants_RD.csv'),
        AD = os.path.join(OUT_DIR, '{frac_dir}', 'variants',
            'rep{n}.filtered_variants_AD.csv'),
    output:
        os.path.join(OUT_DIR, '{frac_dir}', 'scSplit',
            'rep{n}', 'scSplit_result.csv')
    threads: 1
    resources:
        mem_mb = 16384,
        runtime = 180,
    benchmark:
        os.path.join(OUT_DIR, 'benchmarks', '{frac_dir}.{n}.scSplit.txt')
    params:
        exe = config['scSplit']['exe'],
        out_dir = lambda w: os.path.join(OUT_DIR, f'{w.frac_dir}', 'scSplit',
            f'rep{w.n}')
    shell:
        """
        python {params.exe} run  \
            -r {input.RD} \
            -a {input.AD} \
            -n {SAMPLES_N} \
            -o {params.out_dir}
        """


rule demultiplexing_souporcell:
    input:
        RD = os.path.join(OUT_DIR, '{frac_dir}', 'variants',
            'rep{n}.filtered_variants_RD.mtx'),
        AD = os.path.join(OUT_DIR, '{frac_dir}', 'variants',
            'rep{n}.filtered_variants_AD.mtx'),
        barcodes = os.path.join(OUT_DIR, '{frac_dir}', 'variants',
            'rep{n}.filtered_variants_barcodes.csv')
    output:
        os.path.join(OUT_DIR, '{frac_dir}', 'souporcell',
            'rep{n}', 'souporcell_result.tsv')
    threads: 4
    resources:
        mem_mb = 8192,
        runtime = 120,
    benchmark:
        os.path.join(OUT_DIR, 'benchmarks', '{frac_dir}.{n}.souporcell.txt')
    params:
        exe = config['souporcell']['exe'],
        exe_dbt = config['souporcell']['troublet'],
        min_alt = config['souporcell'].get('min_alt', 10),
        min_ref = config['souporcell'].get('min_ref', 10),
        restarts = config['souporcell'].get('restarts', 100)
    shell:
        """
        {params.exe} \
            --alt_matrix {input.AD} \
            --ref_matrix {input.RD} \
            --barcodes {input.barcodes} \
            --num_clusters {SAMPLES_N} \
            --threads {threads} \
            --min_alt {params.min_alt} \
            --min_ref {params.min_ref} \
            --restarts {params.restarts} \
        > {output}.tmp 2> /dev/null \
        && {params.exe_dbt} \
            --alts {input.AD} \
            --refs {input.RD} \
            --clusters {output}.tmp \
        > {output} 2> /dev/null
        """


rule demultiplexing_vireo:
    input:
        vcf = os.path.join(OUT_DIR, '{frac_dir}', 'variants',
            'rep{n}.filtered_variants.vcf'),
    output:
        os.path.join(OUT_DIR, '{frac_dir}', 'vireo',
            'rep{n}', 'vireo_result.tsv')
    threads: 4
    resources:
        mem_mb = 8192,
        runtime = 90,
    benchmark:
        os.path.join(OUT_DIR, 'benchmarks', '{frac_dir}.{n}.vireo.txt')
    params:
        out_dir = lambda w: os.path.join(OUT_DIR, f'{w.frac_dir}', 'vireo',
            f'rep{w.n}')
    shell:
        """
        vireo -c {input.vcf} \
            -N {SAMPLES_N} \
            -o {params.out_dir} \
            -t GT \
            --noPlot \
        && mv {params.out_dir}/donor_ids.tsv {params.out_dir}/vireo_result.tsv
        """


def get_demultiplexed_files(frac_dir):
    files = []
    for n in REPEATS:
        files.append(os.path.join(OUT_DIR, f'{frac_dir}', 'demoTape',
            f'rep{n}.demoTape.assignments.tsv'))
        if config.get('scSplit', {}).get('run', False):
            files.append(os.path.join(OUT_DIR, f'{frac_dir}', 'scSplit',
                f'rep{n}', 'scSplit_result.csv'))
        if config.get('souporcell', {}).get('run', False):
            files.append(os.path.join(OUT_DIR, f'{frac_dir}', 'souporcell',
                f'rep{n}', 'souporcell_result.tsv'))
        if config.get('vireo', {}).get('run', False):
            files.append(os.path.join(OUT_DIR, f'{frac_dir}', 'vireo',
                f'rep{n}', 'vireo_result.tsv'))
    return files


rule demultiplexing_summary:
    input:
        get_demultiplexed_files
    output:
        os.path.join(OUT_DIR, '{frac_dir}', 'summary.tsv')
    threads: 1
    resources:
        mem_mb = 2048,
        runtime = 30
    shell:
        """
        python {SCRIPT_DIR}/evaluate_simulations.py \
            -i {input} \
            -o {output}
        """


def get_summary_files(*args):
    files = []
    if not isinstance(config['simulations']['fractions'][0], list):
        config['simulations']['fractions'] = [config['simulations']['fractions']]
    if not isinstance(config['simulations']['doublets'], list):
        config['simulations']['doublets'] = [config['simulations']['doublets']]

    for fracs in config['simulations']['fractions']:
        for dbt in config['simulations']['doublets']:
            frac_str = '-'.join([str(i) for i in fracs]) + f'-d{dbt:.2f}'
            out_dir_frac = os.path.join(OUT_DIR, frac_str)
            files.append(os.path.join(out_dir_frac, 'summary.tsv'))

            if config.get('simulations', {}).get('run_early', False):
                files.append(os.path.join(out_dir_frac, 'multiplexed-' + frac_str,
                    'results', f'{PREFIX}.report.html'))

    return files


rule demultiplexing_merge_summaries:
    input:
        get_summary_files
    output:
        os.path.join(OUT_DIR, 'summary.tsv')
    threads: 1
    resources:
        mem_mb = 2048,
        runtime = 30
    shell:
        """
        python {SCRIPT_DIR}/merge_simulations.py \
            -i {input} \
            -o {output}
        """


rule generate_summary_plot:
    input:
        os.path.join(OUT_DIR, 'summary.tsv')
    output:
        os.path.join(OUT_DIR, 'summary.png')
    threads: 1
    resources:
        mem_mb = 4096,
        runtime = 30
    shell:
        """
        python {SCRIPT_DIR}/plot_simulations_summary.py -i {input}
        """


def get_benchmark_files(frac_dir):
    files = []
    if not isinstance(config['simulations']['fractions'][0], list):
        config['simulations']['fractions'] = [config['simulations']['fractions']]
    if not isinstance(config['simulations']['doublets'], list):
        config['simulations']['doublets'] = [config['simulations']['doublets']]
    algs = ['demoTape']
    for alg_other in ['scSplit', 'souporcell', 'vireo']:
        if config.get(alg_other, {}).get('run', False):
            algs.append(alg_other)  


    for fracs in config['simulations']['fractions']:
        for dbt in config['simulations']['doublets']:
            frac_str = '-'.join([str(i) for i in fracs]) + f'-d{dbt:.2f}'
            for n in REPEATS:
                for alg in algs:
                    files.append(os.path.join(OUT_DIR, 'benchmarks',
                        f'{frac_str}.{n}.{alg}.txt'))
    return files


rule merge_benchmarks:
    input:
        get_benchmark_files
    output:
        os.path.join(OUT_DIR, 'benchmark_summary.tsv')
    threads: 1
    resources:
        mem_mb = 4096,
        runtime = 30
    shell:
        """
        python {SCRIPT_DIR}/merge_simulations_benchmark.py -i {input} -o {output}
        """


rule generate_benchmark_plot:
    input:
        os.path.join(OUT_DIR, 'benchmark_summary.tsv')
    output:
        os.path.join(OUT_DIR, 'benchmark_summary.png')
    threads: 1
    resources:
        mem_mb = 4096,
        runtime = 30
    shell:
        """
        python {SCRIPT_DIR}/plot_simulations_benchmark.py -i {input} -o {output}
        """