#configfile: "config/config.yaml"

PROTEINSETS = ['model', 'uniprot', 'orthodb']

localrules: preprocess_fold, preprocess_proteins, embl

# Define conditions for producing different output types 
is_rna_seq = ("rna_seq_r1" in config and "rna_seq_r2" in config and config["rna_seq_r1"] and config["rna_seq_r2"])
is_rna_seq_u = ("rna_seq_u" in config and config["rna_seq_u"])
is_iso_seq =  ("iso_seq" in config and config["iso_seq"])
is_on_seq = ("on_seq" in config and config["on_seq"])
is_locus_tag = ("locus_tag"  in config and config["locus_tag"])

# Create output paths if conditions met
rna_seq = ["stringtie/hisat2.sort.bam"] if is_rna_seq  else []
rna_seq_u = ["stringtie/hisat2_u.sort.bam"] if is_rna_seq_u else []
iso_seq = ["stringtie/minimap_iso.sort.bam"] if is_iso_seq  else []
on_seq = ["stringtie/minimap_on.sort.bam"] if is_on_seq else []
merged_iso_on = if is_iso_seq and is_on_seq ["stringtie/minimap.sort.bam"] else []

embl_list = if is_locus_tag [expand("{prefix}.embl.gz", prefix=config["prefix"])] else []

trans_seq = list()
if is_iso_seq or is_rna_seq or is_rna_seq_u or is_on_seq:
    trans_seq.append("stringtie/transcripts.fasta.transdecoder.genome.gff3")
    trans_seq.append(expand("stringtie/busco_stringtie_{lineage}", lineage=config["busco_lineage"]))


rule all:
    input:
        rna_seq,
        rna_seq_u,
        iso_seq,
        on_seq,
        trans_seq,
        embl_list,
        merged_iso_on,
        expand("functional/{prefix}.gff", prefix=config["prefix"]),
        expand("functional/{prefix}.proteins.fa", prefix=config["prefix"]),
        expand("evm/busco_evm_{lineage}", lineage=config["busco_lineage"]),
        expand("galba/busco_galba_{lineage}", lineage=config["busco_lineage"]),
        expand("functional/busco_{prefix}_{lineage}", prefix=config["prefix"], lineage=config["busco_lineage"]),
        expand("miniprot/busco_model_mp_{lineage}", lineage=config["busco_lineage"]),
        expand("{prefix}_EarlGrey/{prefix}_summaryFiles/{prefix}-families.fa.strained", prefix=config["prefix"]),

rule no_te:
    input:
        rna_seq,
        rna_seq_u,
        iso_seq,
        on_seq,
        trans_seq,
        embl_list,
        merged_iso_on,
        expand("functional/{prefix}.gff", prefix=config["prefix"]),
        expand("functional/{prefix}.proteins.fa", prefix=config["prefix"]),
        expand("evm/busco_evm_{lineage}", lineage=config["busco_lineage"]),
        expand("galba/busco_galba_{lineage}", lineage=config["busco_lineage"]),
        expand("functional/busco_{prefix}_{lineage}", prefix=config["prefix"], lineage=config["busco_lineage"]),
        expand("miniprot/busco_model_mp_{lineage}", lineage=config["busco_lineage"]),

rule all_mp:
    input:
        expand("miniprot/{protein_set}_mp4evm.gff3",  protein_set=PROTEINSETS)

# Do some processing, mainly folding the assembly since single line fasta does not work well.
# Complex headers can create problems. Best if user fixes that outside this, because handling
# all different cases would be complex
rule preprocess_fold:
    output:
        processed = expand("{prefix}.fold.fa", prefix=config["prefix"]),
    input: 
        assembly = config["assembly"],
    shell:
        r"""
        cat {input.assembly} | fold > {output.processed}
        """

rule preprocess_proteins:
    output:
        model_proteins = "proteins/model.fasta",
        uniprot_proteins ="proteins/uniprot.fasta",
        orthodb_proteins ="proteins/orthodb.fasta",
    input:
        model = config["model_proteins"],
        uniprot = config["uniprot_proteins"],
        orthodb = config["orthodb_proteins"],
    shell:
        r"""
        mkdir -p proteins
        cp -rf {input.model} {output.model_proteins}
        ln -s {input.uniprot} {output.uniprot_proteins}
        ln -s {input.orthodb} {output.orthodb_proteins}
        """


#1a: create softmasked assembly with RepeatModeler from TETools (https://github.com/Dfam-consortium/TETools). RepeatModeler takes a lot of time. Maybe do RED instead?
#earlgrey instead

rule earlgrey:
    output:
       repeats = "{prefix}_EarlGrey/{prefix}_summaryFiles/{prefix}-families.fa.strained"
    input:
       assembly = "{prefix}.fold.fa"
    params:
        threads = 64                        
    resources:
        time = "128:0:0",
        mem_per_cpu = "5G",
        partition = "bigmem",
        ntasks = 64   
    shell:
        r"""
        earlGrey -t {params.threads} -g {input.assembly} -s {wildcards.prefix} -o .
        """

rule softmask:
    output:
        masked = "{prefix}.softmasked.fa"
    input:
        assembly = "{prefix}.fold.fa"
    resources:
        time = "20:0:0",
        mem_per_cpu = "6G",
    shell:
        r"""
        python /cluster/projects/nn8013k/opt/redmask/redmask.py -i {input.assembly} -o {wildcards.prefix} 
        """
#1b: map proteins from a model species against the genome assembly, zebrafish for fish, human for mammals, etc., using miniprot, convert to EVM input
#1c: map SwissProt/UniProtKB against the genome assembly using miniprot, convert to EVM input
#1d: map OrthoDB v11 clade of proteins against the genome assembly using miniprot (can take some time and memory), convert to EVM input 

rule miniprot:
    output:
        gff = "miniprot/{protein_set}_mp_aln.gff",
    input:
        assembly = expand("{prefix}.fold.fa", prefix=config["prefix"]),
        proteins = "proteins/{protein_set}.fasta", 
    params:
        threads = 10
    resources:
        mem_per_cpu = "5G",
        ntasks = 10,
        time = "72:0:0"
    shell:
        r"""
        miniprot -It{params.threads} --gff {input.assembly} {input.proteins} > miniprot/{wildcards.protein_set}_mp_aln.gff 2> miniprot/miniprot_{wildcards.protein_set}_"`date +\%y\%m\%d_\%H\%M\%S`".err
        """

def get_partition(wildcards):
    if wildcards.protein_set == "orthodb":
        return "bigmem"
    else: 
        return "normal"

def get_mem(wildcards):
    if wildcards.protein_set == "orthodb":
        return "750G"
    else:
        return "100G"

def get_time(wildcards):
    if wildcards.protein_set == "orthodb":
        return "120:0:0"
    else:
        return "48:0:0"

rule process_mp:
    output:
        "miniprot/{protein_set}_mp.proteins.fa",
        gff = "miniprot/{protein_set}_mp4evm.gff3",
    input:
        gff = "miniprot/{protein_set}_mp_aln.gff",
        assembly = expand("{prefix}.fold.fa", prefix=config["prefix"]),
    params:
        threads = 1
    resources:
        #time = config["fcs_gx_time"],
        mem_per_cpu = get_mem,
        partition = get_partition,
        ntasks = 1,
        time = get_time,
    shell:
        r"""        
        #need a bit of temporary dir, standard on Saga might be too small for some purposes
        mkdir -p $USERWORK/tmp
        export TMPDIR=$USERWORK/tmp
        
        if [ {wildcards.protein_set} != orthodb ]
        then
            agat_sp_extract_sequences.pl --gff miniprot/{wildcards.protein_set}_mp_aln.gff -f {input.assembly} -p -o miniprot/{wildcards.protein_set}_mp.proteins.fa #might not work for vertebrata
            agat_sp_extract_sequences.pl --gff miniprot/{wildcards.protein_set}_mp_aln.gff -f {input.assembly} -t exon --merge -o miniprot/{wildcards.protein_set}_mp.mrna.fa #might not work for odb11_vertebrata
        else
	    touch miniprot/{wildcards.protein_set}_mp.proteins.fa
            touch miniprot/{wildcards.protein_set}_mp.mrna.fa
        fi
        python $EVM_HOME/EvmUtils/misc/miniprot_GFF_2_EVM_GFF3.py {input.gff}  > {output.gff}
        """    

#1e1: if RNA-Seq data, map against genome with HiSat2
rule hisat:
    output:
        "stringtie/hisat2.sort.bam",
    input:
        assembly = expand("{prefix}.fold.fa", prefix=config["prefix"]),
        r1 = config["rna_seq_r1"],
        r2 = config["rna_seq_r2"],
    params:
        threads = 20,
        prefix = config["prefix"]
    resources:
        time = "96:0:0",
        mem_per_cpu = "6G",
        ntasks = 20
    shell:
        r"""
        mkdir -p stringtie
        cd stringtie
        hisat2-build ../{input.assembly} {params.prefix} 2> hisat-build.err
        
        hisat2 -p {params.threads} --dta -x {params.prefix} \
        -1 {input.r1} \
        -2 {input.r2} \
        --rna-strandness RF \
        2> hisat.err | samtools view -buS - | \
        samtools sort -m 4G -@ {params.threads} -T tmp -O bam - > hisat2.sort.bam 2> samtools_hisat.err
        """

rule hisat_u:
    output:
        "stringtie/hisat2_u.sort.bam",
    input:
        assembly = expand("{prefix}.fold.fa", prefix=config["prefix"]),
        u = config["rna_seq_r1"],
    params:
        threads = 20,
        prefix = config["prefix"]
    resources:
        #time = config["fcs_gx_time"],
        time = "96:0:0",
        mem_per_cpu = "6G",
        #partition = config["fcs_gx_partition"],
        ntasks = 20
    shell:
        r"""
        mkdir -p stringtie
        cd stringtie
        hisat2-build ../{input.assembly} {params.prefix} 2> hisat-build.err

        hisat2 -p {params.threads} --dta -x {params.prefix} \
        -u {input.u} \
        --rna-strandness R \
        2> hisat.err | samtools view -buS - | \
        samtools sort -m 4G -@ {params.threads} -T tmp -O bam - > hisat2.sort.bam 2> samtools_hisat.err
        """

def get_mode(wildcards):
    if is_rna_seq and is_iso_seq and is_on_seq:
        return "onisorna"
    elif is_rna_seq and is_iso_seq:
        return "isorna"
    elif is_iso_seq and is_on_seq:
        return "oniso"
    elif is_rna_seq and is_on_seq:
        return "onrna"
    elif is_rna_seq:
        return "rna"
    elif is_rna_seq_u:
        return "rna_u"
    elif is_iso_seq:
        return "iso"
    elif is_on_seq:
        return "on"
    else:
        print ("something went wrong, should not get here")

#if iso data, run against genome with minimap2
rule isoseq:
    output:
        "stringtie/minimap_iso.sort.bam",
    input:
        assembly = expand("{prefix}.fold.fa", prefix=config["prefix"]),
        iso = config["iso_seq"],
    params:
        threads = 20,
        prefix = config["prefix"]
        mode = get_mode
    resources:
        time = "96:0:0",
        mem_per_cpu = "4G",
        ntasks = 20
    shell:
        """
        mkdir -p stringtie
        cd stringtie

        minimap2 -t {params.threads} -ax splice:hq -uf ../{input.assembly} {input.iso} | samtools view -buS - | \
        samtools sort -m 4G -@ {params.threads} -T tmp -O bam - > minimap_iso.sort.bam 2> samtools_minimap_iso.err        
        """

#if on data present, run against genome with minimap2
rule onseq:
    output:
        "stringtie/minimap_on.sort.bam",
    input:
        assembly = expand("{prefix}.fold.fa", prefix=config["prefix"]),
        on = config["on_seq"],
    params:
        threads = 20,
        prefix = config["prefix"]
        mode = get_mode
    resources:
        time = "96:0:0",
        mem_per_cpu = "4G",
        ntasks = 20
    shell:
        """
        mkdir -p stringtie
        cd stringtie

        minimap2 -t {params.threads} -ax splice -uf -k14 ../{input.assembly} {input.on} | samtools view -buS - | \
        samtools sort -m 4G -@ {params.threads} -T tmp -O bam - > minimap_on.sort.bam 2> samtools_minimap_on.err
        """


#TODO: what is the point of this?
def get_hisat(wildcards):
    if rna_seq:
        return "stringtie/hisat2.sort.bam"
    elif rna_seq_u:
        return "stringtie/hisat2.sort.bam"
    else:
        return list()

def get_minimap(wildcards):
    if iso_seq and on_seq:
        return "stringtie/minimap.sort.bam"
    elif iso_seq:
        return "stringtie/minimap_iso.sort.bam"
    elif on_seq:
        return "stringtie/minimap_on.sort.bam"
    else:
        return list()

#1f merge iso and on reads if both are present
rule merge_iso_on:
    output:
        "stringtie/minimap.sort.bam",
    input:
        on_minimap = "stringtie/minimap_on.sort.bam"
        iso_minimap = "stringtie/minimap_iso.sort.bam"     
    params:
        threads = 20,
        prefix = config["prefix"],
        mode = get_mode,
    resources:
        time = "96:0:0",
        mem_per_cpu = "4G",
        ntasks = 20
    shell:
        r"""
            samtools merge -o stringtie/minimap.sort.bam {input.on_minimap} {input.iso_minimap}
        """

#1e2: run StringTie to generate a GTF, convert to EVM input via TransDecoder
rule stringtie:
    output:
        "stringtie/transcripts.fasta.transdecoder.genome.gff3",
        "stringtie/stringtie.proteins.fa",
    input:
        assembly = expand("{prefix}.fold.fa", prefix=config["prefix"]),
        hisat = get_hisat,
        minimap = get_minimap,        
        pfam = config["pfam"],
        uniprot =  "proteins/uniprot.fasta",
    params:
        threads = 20,
        prefix = config["prefix"],
        mode = get_mode,
    resources:
        time = "96:0:0",
        mem_per_cpu = "4G",
        ntasks = 20
    shell:
        r"""
        cd stringtie

        if [ {params.mode} == onisorna]; then
            stringtie --mix --rf hisat2.sort.bam {input.minimap} > stringtie.gtf 2> stringtie.err
        elif [ {params.mode} == isorna]; then
            stringtie --mix --rf hisat2.sort.bam {input.minimap} > stringtie.gtf 2> stringtie.err
        elif [ {params.mode} == rna ];	then
            stringtie --rf hisat2.sort.bam > stringtie.gtf 2> stringtie.err
        elif [ {params.mode} == rna_u ];  then
            stringtie hisat2.sort.bam > stringtie.gtf 2> stringtie.err
        elif [[ {params.mode} == iso || {params.mode} == on ]];	then
            stringtie -L {input.minimap} > stringtie.gtf 2> stringtie.err
        else
            echo "Unknown setting"
        fi

        gffread -E stringtie.gtf -o- > stringtie.gff 
        gffread -w stringtie.fa -g ../{input.assembly} stringtie.gtf
        #agat does not work here for some reason
        
        gtf_genome_to_cdna_fasta.pl stringtie.gtf ../{input.assembly} > transcripts.fasta 2> gtf_genome_to_cdna_fasta.err

        gtf_to_alignment_gff3.pl stringtie.gtf > transcripts.gff3 2> transcripts.err

        TransDecoder.LongOrfs -t transcripts.fasta --output_dir . 1> transdecoder_longorfs.out 2> transdecoder_longorfs.err

        hmmsearch --cpu {params.threads} --domtblout pfam.domtblout {input.pfam} transcripts.fasta.transdecoder_dir/longest_orfs.pep 1> hmmscan.out 2> hmmscan.err

        diamond blastp \
        --query transcripts.fasta.transdecoder_dir/longest_orfs.pep  \
        --db ../{input.uniprot} \
	--outfmt 6 \
        --sensitive \
        --max-target-seqs 1 \
        --evalue 1e-5 \
        --threads {params.threads} \
        > blastp.outfmt6 2> diamond_blastp.err

        TransDecoder.Predict -t transcripts.fasta --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6 --output_dir . 1> transdecoder_predict.out 2> transdecoder_predict.err

        cdna_alignment_orf_to_genome_orf.pl \
        transcripts.fasta.transdecoder.gff3 \
        transcripts.gff3 \
        transcripts.fasta > transcripts.fasta.transdecoder.genome.gff3

        ln -s transcripts.fasta.transdecoder.pep stringtie.proteins.fa

        """



#2: run GALBA on the softmasked genome assembly using suitable proteins (zebrafish, human, etc).

rule galba:
    output: 
        "galba/galba.proteins.fa",
        annotation = "galba/galba.gtf",
    input:
        proteins =  "proteins/model.fasta",
        assembly = expand("{prefix}.softmasked.fa", prefix=config["prefix"])
    params:
        threads = 20,
        prefix = config["prefix"]
    resources:
        #time = config["fcs_gx_time"],
        time = "96:0:0",
        mem_per_cpu = "4G",
        #partition = config["fcs_gx_partition"],
        ntasks = 20
    shell:
        r"""
        rm -r galba
        mkdir -p galba 
        singularity exec -B $PWD:/data /cluster/projects/nn8013k/opt/galba/galba.sif cp -rf /usr/share/augustus/config /data/galba
        singularity exec -B $PWD:/data /cluster/projects/nn8013k/opt/galba/galba.sif  galba.pl --version > galba.version
        singularity exec -B $PWD:/data /cluster/projects/nn8013k/opt/galba/galba.sif  galba.pl --species={params.prefix} --threads={params.threads} --genome=/data/{input.assembly} \
         --verbosity=4 --prot_seq=/data/{input.proteins}  --workingdir=/data/galba \
         --augustus_args="--stopCodonExcludedFromCDS=False" --gff3 \
         --AUGUSTUS_CONFIG_PATH=/data/galba/config \
         1> galba/galba_"`date +\%y\%m\%d_\%H\%M\%S`".out 2> galba/galba_"`date +\%y\%m\%d_\%H\%M\%S`".err
        
        agat_sp_extract_sequences.pl --gff galba/galba.gtf -f {input.assembly}  -t cds -p -o galba/galba.proteins.fa
        """    

def get_stringtie(wildcards):
    if trans_seq:
        print("trans_seq")
        return "stringtie/transcripts.fasta.transdecoder.genome.gff3"
    else:
        print("nothing")
        return list()

#3: run EVM on everything
rule evm:
    output:
        "evm/evm.gff3",
        "evm/evm.proteins.fa",
    input:
        orthodb = "miniprot/orthodb_mp4evm.gff3",
        model =  "miniprot/model_mp4evm.gff3",
        uniprot =  "miniprot/uniprot_mp4evm.gff3",
        masked = expand("{prefix}.softmasked.fa", prefix=config["prefix"]),
        galba = "galba/galba.gtf",
        #expand('stringtie/transcripts.fasta.transdecoder.genome.gff3', proxy=[] if rna_seq else [None]),
        stringtie = get_stringtie,
    params:
        threads = 20,
        prefix = config["prefix"]
    resources:
        #time = config["fcs_gx_time"],
        mem_per_cpu = "8G",
        #partition = config["fcs_gx_partition"],
        ntasks = 20
    shell:
        r"""
        mkdir -p evm 
        $EVM_HOME/EvmUtils/misc/augustus_GTF_to_EVM_GFF3.pl {input.galba} > evm/galba.evm.gff3 
        #/cluster/projects/nn8013k/programs/miniconda3/envs/anno_pipeline/opt/evidencemodeler-1.1.1/EvmUtils/misc/augustus_GTF_to_EVM_GFF3.pl {input.galba} > evm/galba.evm.gff3 
        cat evm/galba.evm.gff3 > evm/gene_predictions.gff3 
        cat {input.orthodb} | awk '/gene/ {{printf "\n"}}1' >> evm/gene_predictions.gff3 
        cat {input.model} | awk '/gene/ {{printf "\n"}}1' >> evm/gene_predictions.gff3 
        cat {input.uniprot} | awk '/gene/ {{printf "\n"}}1' >> evm/gene_predictions.gff3
        cat {input.stringtie} | awk '/gene/ {{printf "\n"}}1' >> evm/gene_predictions.gff3
        #the awk thingy adds an empty line before each line matching gene to make the script below work

        printf "ABINITIO_PREDICTION\tAugustus\t1\n" > evm/weights.evm.txt 
        printf "OTHER_PREDICTION\tminiprot\t4\n" >> evm/weights.evm.txt
        printf "OTHER_PREDICTION\ttransdecoder\t6\n" >> evm/weights.evm.txt
        
        cd evm
        
        #/cluster/projects/nn8013k/programs/miniconda3/envs/anno_pipeline/lib/python3.8/site-packages/funannotate/aux_scripts/funannotate-runEVM.py \
        /cluster/projects/nn8013k/programs/miniconda3/envs/anno_pipeline3/lib/python3.8/site-packages/funannotate/aux_scripts/funannotate-runEVM.py \
        -w weights.evm.txt \
        -c {params.threads} \
        -d . \
        -g gene_predictions.gff3 \
        -f ../{input.masked} \
        -l evm.log \
        -m 10 \
        -i 1500 \
        -o evm.gff3 \
        --EVM_HOME $EVM_HOME \
        1> evm_"`date +\%y\%m\%d_\%H\%M\%S`".out 2> evm_"`date +\%y\%m\%d_\%H\%M\%S`".err 
               
        #--EVM_HOME /cluster/projects/nn8013k/programs/miniconda3/envs/anno_pipeline/opt/evidencemodeler-1.1.1 \

        #for some stupid reason, funannotate-runEVM.py does not create the GFF file...
	#"scripts/concat_gff.py"        

        agat_sp_extract_sequences.pl --gff evm.gff3 -f ../{input.masked}  -t cds -p -o evm.proteins.fa 1> agat_proteins_"`date +\%y\%m\%d_\%H\%M\%S`".out 2> agat_proteins_"`date +\%y\%m\%d_\%H\%M\%S`".err
        agat_sp_extract_sequences.pl --gff evm.gff3 -f ../{input.masked}  -t exon --merge -o evm.mrna.fa 1> agat_exons_"`date +\%y\%m\%d_\%H\%M\%S`".out 2> agat_exons_"`date +\%y\%m\%d_\%H\%M\%S`".err
        """
#4: filter annotated proteins based on matches to repeat proteins from Funannotate. And remove all proteins with XX in then (spanning gaps)
rule filter:
    output:
        "filter/filtered.proteins.fa",
        expand("filter/filtered_genes_sup{min_aa_size}.gff", min_aa_size=config["min_aa_size"]),
    input:
        evm_proteins = "evm/evm.proteins.fa",
        evm_gff = "evm/evm.gff3",
        masked = expand("{prefix}.softmasked.fa", prefix=config["prefix"]),
    params:
       min_aa_size = config["min_aa_size"],
    resources:
        #time = config["fcs_gx_time"],
        mem_per_cpu = "10G",
        #partition = config["fcs_gx_partition"],
        ntasks = 10
    shell:
        r"""
        mkdir -p filter
        diamond blastp --sensitive --query {input.evm_proteins} --threads 10 --out filter/repeats.tsv --db /cluster/projects/nn8013k/opt/funannotate/repeats.dmnd --evalue 1e-10 --max-target-seqs 1 --outfmt 6
        
        set +e
        cut -f 1 filter/repeats.tsv |sort -u > filter/kill.list
        #the line below adds any and all proteins with X to the kill list. X is because the nucleotide sequence contains Ns. For some reason it sometimes returns 1. Have to set +e above to avoid aborting the pipeline
        cat {input.evm_proteins} |awk '/^>/ {{printf("%s%s\t",(N>0?"\n":""),$0);N++;next;}} {{printf("%s",$0);}} END {{printf("\n");}}' |grep X |cut -f 1 -d " " | tr -d ">" >> filter/kill.list
        #[ -f filter/removed_repeats.gff ] && rm -rf filter/removed_repeats*
	echo "remove some old files 1"
        rm -rf filter/removed_repeats*
        agat_sp_filter_feature_from_kill_list.pl --gff {input.evm_gff} --kill_list filter/kill.list -o filter/removed_repeats.gff

        #[ -f filter/filtered_genes.gff ] && rm -rf filter/filtered_genes* 
        echo "remove some old files 2"
        rm -rf filter/filtered_genes*
        agat_sp_filter_by_ORF_size.pl --gff filter/removed_repeats.gff -s {params.min_aa_size} -o filter/filtered_genes.gff
        agat_sp_extract_sequences.pl --gff filter/filtered_genes_sup{params.min_aa_size}.gff -f {input.masked} -t cds -p -o filter/filtered.proteins.fa 1> filter/agat_proteins_"`date +\%y\%m\%d_\%H\%M\%S`".out 2> filter/agat_proteins_"`date +\%y\%m\%d_\%H\%M\%S`".err
        """

#5a: run InterProScan on the filtered genes
rule ipr:
    output:
       "ipr/ipr.tsv"
    input:
       filtered = "filter/filtered.proteins.fa"
    resources:
        #time = config["fcs_gx_time"],
        mem_per_cpu = "3G",
        #partition = config["fcs_gx_partition"],
        ntasks = 20
    shell:
       r"""
       mkdir -p ipr
       cat {input.filtered} |sed 's/*//g' > ipr/input.proteins.fa
       cd ipr
       module load StdEnv InterProScan/5.62-94.0-foss-2022a

       timepoint="`date +\%y\%m\%d_\%H\%M\%S`"
       
       mkdir -p $USERWORK/$timepoint       

       interproscan.sh --cpu 20 -dp -i input.proteins.fa  \
       -f XML,GFF3,TSV -goterms -iprlookup -b ipr -T $USERWORK/timepoint/ 1> interproscan.out 2> interproscan.err 

       """

#5b: BLAST filtered genes against uniprot
rule uniprot:
    output:
        "uniprot/diamond.blastp.out"
    input:
       filtered = "filter/filtered.proteins.fa",
       uniprot =  "proteins/uniprot.fasta",
    resources:
        #time = config["fcs_gx_time"],
        mem_per_cpu = "5G",
        #partition = config["fcs_gx_partition"],
        ntasks = 20
    shell:
        r"""
        diamond blastp \
        --query {input.filtered} \
        --db {input.uniprot} \
        --outfmt 6 \
        --sensitive \
        --max-target-seqs 1 \
        --evalue 1e-25 \
        --threads 20 \
        > uniprot/diamond.blastp.out
        """

#6: merge functional annotations together to produce a final set.
rule functional:
    output:
        expand("functional/{prefix}.gff", prefix=config["prefix"]),
        expand("functional/{prefix}.proteins.fa", prefix=config["prefix"]),
    input:
        uniprot_blast = "uniprot/diamond.blastp.out",
        ipr = "ipr/ipr.tsv",
        filtered = expand("filter/filtered_genes_sup{min_aa_size}.gff", min_aa_size=config["min_aa_size"]),
        uniprot =  "proteins/uniprot.fasta",
        masked = expand("{prefix}.softmasked.fa", prefix=config["prefix"]),
    resources:
        #time = config["fcs_gx_time"],
        mem_per_cpu = "40G",
        #partition = config["fcs_gx_partition"],
        ntasks = 1,
    params:
        min_aa_size = config["min_aa_size"],
        prefix = config["prefix"],
    shell:
        r"""
        mkdir -p functional
        cd functional
        agat_sp_manage_functional_annotation.pl --gff ../{input.filtered} -i ../{input.ipr} -b ../{input.uniprot_blast} \
        --ID FUNC \
        -o fun \
        --clean_name \
        -db ../{input.uniprot} \
        1> manage_functional_annotation.out 2> manage_functional_annotation.err

        agat_sp_statistics.pl --gff fun/filtered_genes_sup{params.min_aa_size}.gff -o gff_stats 1> gff_stats.out 2> gff_stats.err

        cp fun/filtered_genes_sup{params.min_aa_size}.gff {params.prefix}.gff 

        agat_sp_extract_sequences.pl --gff {params.prefix}.gff -f ../{input.masked} -t cds -p -o {params.prefix}.proteins.fa 1> agat_proteins_"`date +\%y\%m\%d_\%H\%M\%S`".out 2> agat_proteins_"`date +\%y\%m\%d_\%H\%M\%S`".err
        agat_sp_extract_sequences.pl --gff {params.prefix}.gff -f ../{input.masked} -t exon --merge -o {params.prefix}.mrna.fa 1> agat_mrna_"`date +\%y\%m\%d_\%H\%M\%S`".out 2> agat_mrna_"`date +\%y\%m\%d_\%H\%M\%S`".err
        """

#7: emblmygff. This program contacts NCBI to get more information about the species. Need possibly internet access
rule embl:
    input:
        annotation = expand("functional/{prefix}.gff", prefix=config["prefix"]),
        assembly = expand("{prefix}.fold.fa", prefix=config["prefix"]),
    output:
        expand("{prefix}.embl.gz", prefix=config["prefix"]),
    params:
        locus_tag = config["locus_tag"],
        species = config["species"],
        project_id = config["project_id"],
        rg = "EBP-Nor",
        prefix = config["prefix"],
        mito_code = config["mito_code"],
        plastid_code = config["plastid_code"],
        main_code = config["main_code"]
    shell:
        r"""

        #Need to create separate EMBL files for mitochondria and chloroplast because they have different translation tables

        mkdir -p embl
        cd embl
        samtools faidx ../{input.assembly}
        samtools faidx -c ../{input.assembly} MT > MT.fasta 
        samtools faidx -c ../{input.assembly} PL > PL.fasta 

        seqtk seq ../{input.assembly} | paste - - |grep -v MT |grep -v PL | tr -s "\t" "\n" |fold > wo_MT_PL.fasta
        
        grep ^MT ../{input.annotation} > MT.gff || true
        grep ^PL ../{input.annotation} > PL.gff || true

        grep -v ^MT ../{input.annotation} |grep -v ^PL > wo_MT_PL.gff

        EMBLmyGFF3 \
        wo_MT_PL.gff \
        wo_MT_PL.fasta \
        --topology linear \
        --molecule_type 'genomic DNA' \
        --transl_table {params.main_code}  \
        --rg {params.rg} \
        --species '{params.species}' \
        --locus_tag {params.locus_tag} \
        --project_id {params.project_id} \
        -o wo_MT_PL.embl 

        start_locus=$(grep locus_tag wo_MT_PL.embl |tail -n 1 |cut -d '_' -f 3 |tr -d "LOCUS" |tr -d '"')
        ((start_locus+=1))


        EMBLmyGFF3 \
        MT.gff \
        MT.fasta \
        --topology linear \
        --molecule_type 'genomic DNA' \
        --transl_table {params.mito_code} \
        --organelle mitochondrion \
        --rg {params.rg} \
        --species '{params.species}' \
        --locus_tag {params.locus_tag} \
        --project_id {params.project_id} \
        --locus_numbering_start $start_locus \
        -o MT.embl || true

        if [ -s MT.embl ]
        then 
            continue_locus=$(grep locus_tag MT.embl |tail -n 1 |cut -d '_' -f 3 |tr -d "LOCUS" |tr -d '"')
        else
            continue_locus=$start_locus
        fi
        ((continue_locus+=1))

        EMBLmyGFF3 \
        PL.gff \
        PL.fasta \
        --topology linear \
        --molecule_type 'genomic DNA' \
        --transl_table {params.plastid_code}  \
        --organelle plastid \
        --rg {params.rg} \
        --species '{params.species}' \
        --locus_tag {params.locus_tag} \
        --project_id {params.project_id} \
        --locus_numbering_start $continue_locus \
        -o PL.embl || true

        cd ..

        #EMBLmyGFF3 \
        #{input.annotation} \
        #{input.assembly} \
        #--topology linear \
        #--molecule_type 'genomic DNA' \
        #--transl_table 1  \
        #--rg {params.rg} \
        #--species '{params.species}' \
        #--locus_tag {params.locus_tag} \
        #--project_id {params.project_id} \
        #-o {params.prefix}.embl > embl_gff.out 2> embl_gff.err
        
        #EMBLmyGFF has a bug which puts OG on the same line as XX, meaning that downstream analyses don't see that the following is an organelle
        cat embl/*embl |sed "s/XXOG/XX\\nOG/g" > {params.prefix}.embl

        samtools faidx {input.assembly}
        grep SUPER {input.assembly}.fai |awk '{{f=$1 ; sub(/SUPER_/,"",f) ; printf ("%s\t%s\tLinear-Chromosome\n", $1, f)}}' > {params.prefix}.chrs
        grep MT {input.assembly}.fai |awk '{{printf ("%s\tMT\tLinear-Chromosome\tMitochondrion\n", $1)}}' >> {params.prefix}.chrs || true
        grep PL {input.assembly}.fai |awk '{{printf ("%s\tPltd\tLinear-Chromosome\tChloroplast\n", $1)}}' >> {params.prefix}.chrs || true
        # || true forces the command to return 0, instead of 1 if it doesn't find MT or PL in the fasta file. OK to not be successful if failing to find SUPER        

        cat {params.prefix}.chrs | gzip > {params.prefix}.chrs.gz
        cat {params.prefix}.embl | gzip > {params.prefix}.embl.gz
        """

#8: busco
rule busco:
    output:
        #directory(expand("{{dir}}/busco_{{analysis}}_{lineage}", lineage=config["busco_lineage"])),
        directory("{dir}/busco_{analysis}_{lineage}"),
        #"{dir}/busco_{analysis}_{lineage}/short_summary.specific.{lineage}_odb10.busco_{analysis}_{lineage}.txt"
        #functional/busco_inChrCarn3_arthropoda/short_summary.specific.arthropoda_odb10.busco_inChrCarn3_arthropoda.txt
    input:
        proteins = "{dir}/{analysis}.proteins.fa",
        #busco = str(expand("{busco_lib}", busco_lib=config["busco_lib"])) + "/{lineage}_odb10"
        #busco = expand("{config['busco_lib']}/{lineage}_odb10"
        #busco = expand("{busco_lib}/{{lineage}}_odb10",busco_lib=config["busco_lib"]),
        #busco = expand("{busco_lib}/{lineage}_odb10", busco_lib=config["busco_lib"], lineage=config["busco_lineage"])
    resources:
        mem_per_cpu = "7G",
        ntasks = 10
    params:
       busco_lib = config['busco_lib'] 
    shell:
        r"""
        busco --tar -c 10 -i {input.proteins} -l {params.busco_lib}/{wildcards.lineage}_odb10 -o busco_{wildcards.analysis}_{wildcards.lineage} --out_path {wildcards.dir} \
        -m proteins --offline 
        """
