import sys
print("Python executable:", sys.executable)
from snakemake import snakemake
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-c', '--config_file_path', help='path to the config file')
args = parser.parse_args()

success = snakemake(
    snakefile="Snakefile",
    cores=16,                             # --cores
    keepgoing=True,                       # --keep-going
    configfiles=[args.config_file_path],  # --configfile
    nodes=500,                            # -j
    # unlock=True,                        # --unlock
    use_conda=True,                       # --use-conda
    force_incomplete=True,                # --force-incomplete
    latency_wait=60,                      # --latency-wait
    printshellcmds=True,                  # --printshellcmds
    cluster='bsub -q gsla-cpu -n {threads} '
            '-R "span[hosts=1]" '
            '-R "rusage[mem={params.memory}000]" '
            '-J {params.job_name} '
            '-o {params.job_out_dir}/{params.job_out_file}.o.txt '
            '-e {params.job_out_dir}/{params.job_out_file}.e.txt '
            '-W {params.run_time}'
)

if not success:
    raise RuntimeError("Snakemake workflow failed")