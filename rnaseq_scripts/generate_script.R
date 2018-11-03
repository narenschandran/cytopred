# This is a hacky solution to generate the entire script, but it works.
# Make sure to edit the paths in the following lines as necessary
# At the very least, you would need to edit the salmon_index_dir

run_table_file_path="../SRP050272_SraRunTable.txt"
script_file="./script.sh"
srr_files_dir="../srr_files"
fastq_dump_dir="../fastq_dump"
fastqc_dir="../fastqc"
salmon_index_dir="/home/naren/Workspace/software/Salmon-0.8.2_linux_x86_64/reference_index/Homo_sapiens.GRCh38.cdna.all_index"
salmon_output_dir="../salmon_output"
num_threads=5

run_table <- read.table(run_table_file_path, sep="\t", header=T)
run_table <- run_table[order(run_table$Run_s),]
exp_ids <- as.character(unique(run_table$Experiment_s))

lapply(exp_ids, function(exp_id) {
    write(paste("#", exp_id, "commands"), file=script_file, append=T)
    write("#-----", file=script_file, append=T)
    srr_ids <- as.character(run_table$Run_s[run_table$Experiment_s==exp_id])
    srr_file_paths <- file.path(srr_files_dir, paste(srr_ids, ".sra", sep=""))
    sapply(srr_file_paths, function(srr_file_path) {
        fastq_dump_command <- paste("fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files", srr_file_path, "-O", fastq_dump_dir)
        write(fastq_dump_command, file=script_file, append=T)
    })
    fastq1_dump_file_paths <- file.path(fastq_dump_dir, paste(srr_ids, "_1.fastq", sep=""))
    exp_fastq1_path <- file.path(fastq_dump_dir, paste(exp_id, "_1.fastq", sep=""))
    fastq1_cat_command <- paste("cat", paste(fastq1_dump_file_paths, collapse=" "), ">", exp_fastq1_path)
    write(fastq1_cat_command, file=script_file, append=T)
    rm_fastq1_dump_command <- paste("rm", paste(fastq1_dump_file_paths, collapse=" "))
    write(rm_fastq1_dump_command, file=script_file, append=T)

    fastq2_dump_file_paths <- file.path(fastq_dump_dir, paste(srr_ids, "_2.fastq", sep=""))
    exp_fastq2_path <- file.path(fastq_dump_dir, paste(exp_id, "_2.fastq", sep=""))
    fastq2_cat_command <- paste("cat", paste(fastq2_dump_file_paths, collapse=" "), ">", exp_fastq2_path)
    write(fastq2_cat_command, file=script_file, append=T)
    rm_fastq2_dump_command <- paste("rm", paste(fastq2_dump_file_paths, collapse=" "))
    write(rm_fastq2_dump_command, file=script_file, append=T)

    #fastqc_command <- paste("fastqc", exp_fastq_path, "-o", fastqc_dir, "-f fastq", "-t", num_threads)
    #write(fastqc_command, file=script_file, append=T)
    salmon_command <- paste("salmon quant -i", salmon_index_dir, "-l A -1", exp_fastq1_path, "-2", exp_fastq2_path, "-o", file.path(salmon_output_dir, exp_id), "-p", num_threads)
    write(salmon_command,file=script_file, append=T)
    rm_exp_fastq_command <- paste("rm", exp_fastq1_path, exp_fastq2_path)
    write(rm_exp_fastq_command, file=script_file, append=T)
    write("\n", file=script_file, append=T)
})
