set -e
true
true
/opt/conda/envs/wf1/bin/spades-hammer /workspace/gitpod/nf-training/work/aa/b15235d6ffbdb044f5835075e795cd/ANN0831_R_assembled/corrected/configs/config.info
/opt/conda/envs/wf1/bin/python /opt/conda/envs/wf1/share/spades/spades_pipeline/scripts/compress_all.py --input_file /workspace/gitpod/nf-training/work/aa/b15235d6ffbdb044f5835075e795cd/ANN0831_R_assembled/corrected/corrected.yaml --ext_python_modules_home /opt/conda/envs/wf1/share/spades --max_threads 16 --output_dir /workspace/gitpod/nf-training/work/aa/b15235d6ffbdb044f5835075e795cd/ANN0831_R_assembled/corrected --gzip_output
true
true
/opt/conda/envs/wf1/bin/spades-core /workspace/gitpod/nf-training/work/aa/b15235d6ffbdb044f5835075e795cd/ANN0831_R_assembled/K21/configs/config.info
/opt/conda/envs/wf1/bin/spades-core /workspace/gitpod/nf-training/work/aa/b15235d6ffbdb044f5835075e795cd/ANN0831_R_assembled/K33/configs/config.info
/opt/conda/envs/wf1/bin/spades-core /workspace/gitpod/nf-training/work/aa/b15235d6ffbdb044f5835075e795cd/ANN0831_R_assembled/K55/configs/config.info
/opt/conda/envs/wf1/bin/spades-core /workspace/gitpod/nf-training/work/aa/b15235d6ffbdb044f5835075e795cd/ANN0831_R_assembled/K77/configs/config.info
/opt/conda/envs/wf1/bin/python /opt/conda/envs/wf1/share/spades/spades_pipeline/scripts/copy_files.py /workspace/gitpod/nf-training/work/aa/b15235d6ffbdb044f5835075e795cd/ANN0831_R_assembled/K77/before_rr.fasta /workspace/gitpod/nf-training/work/aa/b15235d6ffbdb044f5835075e795cd/ANN0831_R_assembled/before_rr.fasta /workspace/gitpod/nf-training/work/aa/b15235d6ffbdb044f5835075e795cd/ANN0831_R_assembled/K77/assembly_graph_after_simplification.gfa /workspace/gitpod/nf-training/work/aa/b15235d6ffbdb044f5835075e795cd/ANN0831_R_assembled/assembly_graph_after_simplification.gfa /workspace/gitpod/nf-training/work/aa/b15235d6ffbdb044f5835075e795cd/ANN0831_R_assembled/K77/final_contigs.fasta /workspace/gitpod/nf-training/work/aa/b15235d6ffbdb044f5835075e795cd/ANN0831_R_assembled/contigs.fasta /workspace/gitpod/nf-training/work/aa/b15235d6ffbdb044f5835075e795cd/ANN0831_R_assembled/K77/first_pe_contigs.fasta /workspace/gitpod/nf-training/work/aa/b15235d6ffbdb044f5835075e795cd/ANN0831_R_assembled/first_pe_contigs.fasta /workspace/gitpod/nf-training/work/aa/b15235d6ffbdb044f5835075e795cd/ANN0831_R_assembled/K77/strain_graph.gfa /workspace/gitpod/nf-training/work/aa/b15235d6ffbdb044f5835075e795cd/ANN0831_R_assembled/strain_graph.gfa /workspace/gitpod/nf-training/work/aa/b15235d6ffbdb044f5835075e795cd/ANN0831_R_assembled/K77/scaffolds.fasta /workspace/gitpod/nf-training/work/aa/b15235d6ffbdb044f5835075e795cd/ANN0831_R_assembled/scaffolds.fasta /workspace/gitpod/nf-training/work/aa/b15235d6ffbdb044f5835075e795cd/ANN0831_R_assembled/K77/scaffolds.paths /workspace/gitpod/nf-training/work/aa/b15235d6ffbdb044f5835075e795cd/ANN0831_R_assembled/scaffolds.paths /workspace/gitpod/nf-training/work/aa/b15235d6ffbdb044f5835075e795cd/ANN0831_R_assembled/K77/assembly_graph_with_scaffolds.gfa /workspace/gitpod/nf-training/work/aa/b15235d6ffbdb044f5835075e795cd/ANN0831_R_assembled/assembly_graph_with_scaffolds.gfa /workspace/gitpod/nf-training/work/aa/b15235d6ffbdb044f5835075e795cd/ANN0831_R_assembled/K77/assembly_graph.fastg /workspace/gitpod/nf-training/work/aa/b15235d6ffbdb044f5835075e795cd/ANN0831_R_assembled/assembly_graph.fastg /workspace/gitpod/nf-training/work/aa/b15235d6ffbdb044f5835075e795cd/ANN0831_R_assembled/K77/final_contigs.paths /workspace/gitpod/nf-training/work/aa/b15235d6ffbdb044f5835075e795cd/ANN0831_R_assembled/contigs.paths
true
/opt/conda/envs/wf1/bin/python /opt/conda/envs/wf1/share/spades/spades_pipeline/scripts/breaking_scaffolds_script.py --result_scaffolds_filename /workspace/gitpod/nf-training/work/aa/b15235d6ffbdb044f5835075e795cd/ANN0831_R_assembled/scaffolds.fasta --misc_dir /workspace/gitpod/nf-training/work/aa/b15235d6ffbdb044f5835075e795cd/ANN0831_R_assembled/misc --threshold_for_breaking_scaffolds 3
true
