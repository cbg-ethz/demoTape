rule simulate_multiplexing:
    input:
        LOOM_FILES
    output:
        var = os.path.join(OUT_DIR, '{frac}-d{dbt}{downsample}', 'variants',
            'rep{n}.filtered_variants.csv'),
        vcf = os.path.join(OUT_DIR, '{frac}-d{dbt}{downsample}', 'variants',
            'rep{n}.filtered_variants.vcf'),
        DP_T = os.path.join(OUT_DIR, '{frac}-d{dbt}{downsample}', 'variants',
            'rep{n}.filtered_variants_DP_T.csv'),
        RD = os.path.join(OUT_DIR, '{frac}-d{dbt}{downsample}', 'variants',
            'rep{n}.filtered_variants_RD.csv'),
        AD = os.path.join(OUT_DIR, '{frac}-d{dbt}{downsample}', 'variants',
            'rep{n}.filtered_variants_AD.csv'),
        AD_T = os.path.join(OUT_DIR, '{frac}-d{dbt}{downsample}', 'variants',
            'rep{n}.filtered_variants_AD_T.csv'),
        RD_mtx = os.path.join(OUT_DIR, '{frac}-d{dbt}{downsample}', 'variants',
            'rep{n}.filtered_variants_RD.mtx'),
        AD_mtx = os.path.join(OUT_DIR, '{frac}-d{dbt}{downsample}', 'variants', 
            'rep{n}.filtered_variants_AD.mtx'),
        barcodes = os.path.join(OUT_DIR, '{frac}-d{dbt}{downsample}', 'variants',
            'rep{n}.barcodes.csv')
    wildcard_constraints:
        frac = r'[\d\.-]+',
        dbt = r'[\d\.]+',
        downsample = r'(-dp[\d.]+)?'
    conda: os.path.join(ENV_DIR, 'analysis.yaml')
    threads: 1
    resources:
        mem_mb_per_cpu = 49152,
        runtime = 240,
    params:
        ratios = lambda w: w.frac.replace('-', ' '),
        downsampling_depth = lambda w: \
            f'-dp {w.downsample[3:]}' if w.downsample else '',
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
            -d {wildcards.dbt} \
            {params.downsampling_depth} \
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


def get_demotape_outdir(w):
    if w.downsample:
        subdir = f'{w.frac}-d{w.dbt}-{w.downsample[1:]}'
    else:
        subdir = f'{w.frac}-d{w.dbt}'

    return os.path.join(OUT_DIR, subdir, 'demoTape', f'rep{w.n}.demoTape')


rule demultiplexing_demoTape:
    input:
        os.path.join(OUT_DIR, '{frac}-d{dbt}{downsample}', 'variants',
            'rep{n}.filtered_variants.csv')
    output:
        os.path.join(OUT_DIR, '{frac}-d{dbt}{downsample}', 'demoTape',
            'rep{n}.demoTape.assignments.tsv')
    conda: os.path.join(ENV_DIR, 'analysis.yaml')
    wildcard_constraints:
        frac = r'[\d\.-]+',
        dbt = r'[\d\.]+',
        downsample = r'(-dp[\d\.]+)?'
    threads: 1
    resources:
        mem_mb_per_cpu = 4096,
        runtime = 360,
    benchmark:
        os.path.join(OUT_DIR, 'benchmarks',
            '{frac}-d{dbt}{downsample}.{n}.demoTape.txt')
    params:
        out_base = get_demotape_outdir,
        sample_no = len(LOOM_FILES)
    shell:
        """
        python {SCRIPT_DIR}/demultiplex_distance.py \
            -i {input} \
            -o {params.out_base} \
            -n {params.sample_no}
        """


def get_scSplit_outdir(w):
    if w.downsample:
        subdir = f'{w.frac}-d{w.dbt}-{w.downsample[1:]}'
    else:
        subdir = f'{w.frac}-d{w.dbt}'
    return os.path.join(OUT_DIR, subdir, 'scSplit', f'rep{w.n}')


rule demultiplexing_scSplit:
    input:
        RD = os.path.join(OUT_DIR, '{frac}-d{dbt}{downsample}', 'variants',
            'rep{n}.filtered_variants_RD.csv'),
        AD = os.path.join(OUT_DIR, '{frac}-d{dbt}{downsample}', 'variants',
            'rep{n}.filtered_variants_AD.csv'),
    output:
        os.path.join(OUT_DIR, '{frac}-d{dbt}{downsample}', 'scSplit',
            'rep{n}', 'scSplit_result.csv')
    wildcard_constraints:
        frac = r'[\d\.-]+',
        dbt = r'[\d\.]+',
        downsample = r'(-dp[\d\.]+)?'
    conda: os.path.join(ENV_DIR, 'scSplit.yaml') # requires to manually set exe file in config
    threads: 1
    resources:
        mem_mb_per_cpu = 16384,
        runtime = 180,
    benchmark:
        os.path.join(OUT_DIR, 'benchmarks',
            '{frac}-d{dbt}{downsample}.{n}.scSplit.txt')
    params:
        exe = config['simulations']['algorithms']['scSplit']['exe'],
        out_dir = get_scSplit_outdir,
        sample_no = len(LOOM_FILES)
    shell:
        """
        python {params.exe} run  \
            -r {input.RD} \
            -a {input.AD} \
            -n {params.sample_no} \
            -o {params.out_dir}
        """


rule demultiplexing_souporcell:
    input:
        RD = os.path.join(OUT_DIR, '{frac}-d{dbt}{downsample}', 'variants',
            'rep{n}.filtered_variants_RD.mtx'),
        AD = os.path.join(OUT_DIR, '{frac}-d{dbt}{downsample}', 'variants',
            'rep{n}.filtered_variants_AD.mtx'),
        barcodes = os.path.join(OUT_DIR, '{frac}-d{dbt}{downsample}', 'variants',
            'rep{n}.barcodes.csv')
    output:
        os.path.join(OUT_DIR, '{frac}-d{dbt}{downsample}', 'souporcell',
            'rep{n}', 'souporcell_result.tsv')
    wildcard_constraints:
        frac = r'[\d\.-]+',
        dbt = r'[\d\.]+',
        downsample = r'(-dp[\d\.]+)?'
    conda: os.path.join(ENV_DIR, 'souporcell.yaml') # Env requires manual action!
    threads: 4
    resources:
        mem_mb_per_cpu = 8192,
        runtime = 120,
    benchmark:
        os.path.join(OUT_DIR, 'benchmarks',
            '{frac}-d{dbt}{downsample}.{n}.souporcell.txt')
    params:
        exe = config['simulations']['algorithms']['souporcell']['exe'],
        exe_dbt = config['simulations']['algorithms']['souporcell']['troublet'],
        min_alt = config['simulations']['algorithms']['souporcell'] \
            .get('min_alt', 10),
        min_ref = config['simulations']['algorithms']['souporcell'] \
            .get('min_ref', 10),
        restarts = config['simulations']['algorithms']['souporcell'] \
            .get('restarts', 100),
        sample_no = len(LOOM_FILES)
    shell:
        """
        {params.exe} \
            --alt_matrix {input.AD} \
            --ref_matrix {input.RD} \
            --barcodes {input.barcodes} \
            --num_clusters {params.sample_no} \
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


def get_vireo_outdir(w):
    if w.downsample:
        subdir = f'{w.frac}-d{w.dbt}-{w.downsample[1:]}'
    else:
        subdir = f'{w.frac}-d{w.dbt}'
    return os.path.join(OUT_DIR, subdir, 'vireo', f'rep{w.n}')


rule demultiplexing_vireo:
    input:
        vcf = os.path.join(OUT_DIR, '{frac}-d{dbt}{downsample}', 'variants',
            'rep{n}.filtered_variants.vcf'),
    output:
        os.path.join(OUT_DIR, '{frac}-d{dbt}{downsample}', 'vireo',
            'rep{n}', 'vireo_result.tsv')
    wildcard_constraints:
        frac = r'[\d\.-]+',
        dbt = r'[\d\.]+',
        downsample = r'(-dp[\d.]+)?'
    conda: os.path.join(ENV_DIR, 'vireo.yaml')
    threads: 4
    resources:
        mem_mb_per_cpu = 8192,
        runtime = 90,
    benchmark:
        os.path.join(OUT_DIR, 'benchmarks',
            '{frac}-d{dbt}{downsample}.{n}.vireo.txt')
    params:
        out_dir = get_vireo_outdir,
        sample_no = len(LOOM_FILES)
    shell:
        """
        vireo -c {input.vcf} \
            -N {params.sample_no} \
            -o {params.out_dir} \
            -t GT \
            --noPlot \
        && mv {params.out_dir}/donor_ids.tsv {params.out_dir}/vireo_result.tsv
        """


rule run_doubletD:
    input:
        DP = os.path.join(OUT_DIR, '{frac}-d{dbt}{downsample}', 'variants',
            'rep{n}.filtered_variants_DP_T.csv'),
        AD = os.path.join(OUT_DIR, '{frac}-d{dbt}{downsample}', 'variants',
            'rep{n}.filtered_variants_AD_T.csv'),
    output:
        os.path.join(OUT_DIR, '{frac}-d{dbt}{downsample}', 'doubletD',
            'rep{n}.doubletD_result.tsv')
    wildcard_constraints:
        frac = r'[\d\.-]+',
        dbt = r'[\d\.]+',
        downsample = r'(-dp[\d\.]+)?'
    conda: os.path.join(ENV_DIR, 'doubletD.yaml')
    threads: 1
    resources:
        mem_mb_per_cpu = 8192,
        runtime = 30,
    benchmark:
        os.path.join(OUT_DIR, 'benchmarks',
            '{frac}-d{dbt}{downsample}.{n}.doubletD.txt')
    params:
        ADO_rate = 0.05
    shell:
        """
        doubletd --inputAlternate {input.AD} \
            --inputTotal {input.DP} \
            --delta {wildcards.dbt} \
            --beta {params.ADO_rate} \
            -o {output}
        """


def get_demultiplexed_files(w):
    files = []
    ident_str = f'{w.frac}-d{w.dbt}{w.downsample}'
    for n in REPEATS:
        for alg in ALGORITHMS:
            if alg == 'demoTape':
                files.append(os.path.join(OUT_DIR, f'{ident_str}', 'demoTape',
                    f'rep{n}.demoTape.assignments.tsv'))
            elif alg == 'scSplit':
                files.append(os.path.join(OUT_DIR, f'{ident_str}', 'scSplit',
                    f'rep{n}', 'scSplit_result.csv'))
            elif alg == 'souporcell':
                files.append(os.path.join(OUT_DIR, f'{ident_str}', 'souporcell',
                    f'rep{n}', 'souporcell_result.tsv'))
            elif alg == 'vireo':
                files.append(os.path.join(OUT_DIR, f'{ident_str}', 'vireo',
                    f'rep{n}', 'vireo_result.tsv'))
            elif alg == 'doubletD':
                files.append(os.path.join(OUT_DIR, f'{ident_str}', 'doubletD',
                    f'rep{n}.doubletD_result.tsv'))
    return files


rule demultiplexing_summary:
    input:
        get_demultiplexed_files
    output:
        os.path.join(OUT_DIR, '{frac}-d{dbt}{downsample}', 'summary.tsv')
    wildcard_constraints:
        frac = r'[\d\.-]+',
        dbt = r'[\d\.]+',
        downsample = r'(-dp[\d\.]+)?'
    conda: os.path.join(ENV_DIR, 'analysis.yaml')
    threads: 1
    resources:
        mem_mb_per_cpu = 2048,
        runtime = 30
    shell:
        """
        python {SCRIPT_DIR}/evaluate_simulations.py \
            -i {input} \
            -o {output}
        """


def get_summary_files(*args):
    files = []
    for downsample in DOWNSAMPLES:
        for fracs in FRACTIONS:
            for dbt in DBTS:
                ident_str = '-'.join([str(i) for i in fracs])
                ident_str += f'-d{dbt:.2f}'
                if downsample < 1:
                    ident_str += f'-dp{downsample:.2f}'
                out_dir_frac = os.path.join(OUT_DIR, ident_str)
                files.append(os.path.join(out_dir_frac, 'summary.tsv'))
    return files


rule demultiplexing_merge_summaries:
    input:
        get_summary_files
    output:
        os.path.join(OUT_DIR, 'summary.tsv')
    conda: os.path.join(ENV_DIR, 'analysis.yaml')
    threads: 1
    resources:
        mem_mb_per_cpu = 2048,
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
    conda: os.path.join(ENV_DIR, 'analysis.yaml')
    threads: 1
    resources:
        mem_mb_per_cpu = 4096,
        runtime = 30
    shell:
        """
        python {SCRIPT_DIR}/plot_simulations_summary.py -i {input}
        """


def get_benchmark_files(*args):
    files = []
    for downsample in DOWNSAMPLES:
        for fracs in FRACTIONS:
            for dbt in DBTS:
                ident_str = '-'.join([str(i) for i in fracs])
                ident_str += f'-d{dbt:.2f}'
                if downsample < 1:
                    ident_str += f'-dp{downsample:.2f}'

                for n in REPEATS:
                    for alg in ALGORITHMS:
                        files.append(os.path.join(OUT_DIR, 'benchmarks',
                            f'{ident_str}.{n}.{alg}.txt'))
    return files


rule merge_benchmarks:
    input:
        get_benchmark_files
    output:
        os.path.join(OUT_DIR, 'benchmark_summary.tsv')
    conda: os.path.join(ENV_DIR, 'analysis.yaml')
    threads: 1
    resources:
        mem_mb_per_cpu = 4096,
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
    conda: os.path.join(ENV_DIR, 'analysis.yaml')
    threads: 1
    resources:
        mem_mb_per_cpu = 4096,
        runtime = 30
    shell:
        """
        python {SCRIPT_DIR}/plot_simulations_benchmark.py -i {input} -o {output}
        """
