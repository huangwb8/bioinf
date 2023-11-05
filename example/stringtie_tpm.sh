


ws=$1

# Parameters
path_res=$ws/output/stringtie/STAR_Ensembl_104_unmasked
path_gtf=/data/reference/Ensembl/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.gtf
path_log=$ws/log/stringtie/STAR_Ensembl_104_unmasked
prepDE=/home/huangwb8/Downloads/stringtie/stringtie-2.1.6/prepDE.py

# Dir management
mkdir -p $path_res $path_log 

# Programe: 不需要预测新型转录本

##================Get TPM/FPKM of genes and transcripts=================##

# --fr	Assumes a stranded library fr-secondstrand.
# -e Limits the processing of read alignments to only estimate and output the assembled transcripts matching the reference transcripts given with the -G option (requires -G, recommended for -B/-b). 

ls $ws/output/align/STAR_ensembl_104_unmasked/*bam | while read id
do
    #  |head -n 1; grep -v D1
    name=$(basename ${id} .bam)
    mkdir -p ${path_res}/${name}
    
    if [ -f ${path_res}/${name}/gene_abund.tab ]
    then
    echo `date`" ${path_res}/${sample}: stringtie data exists. Ignore!"
    else
    echo `date`" ${path_res}/${sample}: run stringtie..."
    nohup stringtie ${id} -e -v \
        -A ${path_res}/${name}/gene_abund.tab \
        -C ${path_res}/${name}/cov_refs.gtf \
        -B -p 10 \
        -G ${path_gtf} \
        -o ${path_res}/${name}/merged.gtf --fr \
        > ${path_log}/${name}.log 2>&1 &
    fi
done

echo `date`" All done!"
