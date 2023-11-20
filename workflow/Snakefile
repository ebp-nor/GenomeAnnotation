#configfile: "config/config.yaml"

PROTEINSETS = ['model', 'uniprot', 'orthodb']

localrules: preprocess, embl

rna_seq = list()
if config["rna_seq_r1"] and config["rna_seq_r2"]:
    rna_seq.append("stringtie/transcripts.fasta.transdecoder.genome.gff3")
    rna_seq.append("stringtie/stringtie.proteins.fa")
    lineage = config["busco_lineage"]
    
    rna_seq.append("stringtie/busco_stringtie_%s" % lineage)

rule all:
    input:
        rna_seq,
        expand("miniprot/{protein_set}_mp4evm.gff3",  protein_set=PROTEINSETS),
        expand("{prefix}.softmasked.fa", prefix=config["prefix"]),
        expand("{prefix}.fold.fa", prefix=config["prefix"]),
        "galba/augustus.hints.gff",
        "galba/galba.proteins.fa",
        "evm/evm.gff3",
        "filter/filtered.proteins.fa",
        "uniprot/diamond.blastp.out",
        "ipr/ipr.tsv",
        expand("functional/{prefix}.gff", prefix=config["prefix"]),
        expand("functional/{prefix}.proteins.fa", prefix=config["prefix"]),
        expand("evm/busco_evm_{lineage}", lineage=config["busco_lineage"]),
        expand("galba/busco_galba_{lineage}", lineage=config["busco_lineage"]),
        expand("functional/busco_{prefix}_{lineage}", prefix=config["prefix"], lineage=config["busco_lineage"]),
        expand("miniprot/busco_model_mp_{lineage}", lineage=config["busco_lineage"]),
        expand("{prefix}.embl", prefix=config["prefix"]),

rule all_mp:
    input:
        expand("miniprot/{protein_set}_mp4evm.gff3",  protein_set=PROTEINSETS)

#do some processing, mainly folding the assembly since single line fasta does not work well
#complex headers can create problems. Best if user fixes that outside this, because handling
#all different cases would be complex
rule preprocess:
    output:
        processed = expand("{prefix}.fold.fa", prefix=config["prefix"]),
        model_proteins = "proteins/model.fasta",
        uniprot_proteins ="proteins/uniprot.fasta",
        orthodb_proteins ="proteins/orthodb.fasta",
    input: 
        assembly = config["assembly"],
        model = config["model_proteins"],
        uniprot = config["uniprot_proteins"],
        orthodb = config["orthodb_proteins"],
    shell:
        r"""
        cat {input.assembly} | fold > {output.processed}
        mkdir -p proteins
        cp -rf {input.model} {output.model_proteins}
        ln -s {input.uniprot} {output.uniprot_proteins}
        ln -s {input.orthodb} {output.orthodb_proteins}
        """

#1a: create softmasked assembly with RepeatModeler from TETools (https://github.com/Dfam-consortium/TETools). RepeatModeler takes a lot of time. Maybe do RED instead?

rule softmask_rm:
    output:
       database = "database"
    input:
       assembly = "{prefix}.fold.fa"
#    threads: 20
    params:
        threads = 60                        
    resources:
        time = "128:0:0",
#        mem_per_cpu = "60G",
#        partition = "bigmem",
#        ntasks = 30   
    shell:
       r"""
       singularity exec --bind $PWD:/data --pwd /data /cluster/projects/nn8013k/opt/tetools/dfam-tetools-1.85.sif BuildDatabase -name {output.database} {input.assembly}
       singularity exec --bind $PWD:/data --pwd /data /cluster/projects/nn8013k/opt/tetools/dfam-tetools-1.85.sif RepeatModeler -database {output.database} -LTRStruct -threads {threads}
       """

rule softmask:
    output:
#        masked = expand("{prefix}.softmasked.fa", prefix=config["prefix"])
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
        #time = config["fcs_gx_time"],
        #mem_per_cpu = config["fcs_gx_mem"],
        #partition = config["fcs_gx_partition"],
        mem_per_cpu = "3G",
        ntasks = 10
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
        return "650G"
    else:
        return "100G"

def get_time(wildcards):
    if wildcards.protein_set == "orthodb":
        return "96:0:0"
    else:
        return "48:0:0"


rule process_mp:
    output:
        "miniprot/{protein_set}_mp.proteins.fa",
        gff = "miniprot/{protein_set}_mp4evm.gff3",
#        gff = expand("{protein_set}_mp4evm.gff3", protein_set=PROTEINSETS),
    input:
        gff = "miniprot/{protein_set}_mp_aln.gff",
        assembly = expand("{prefix}.fold.fa", prefix=config["prefix"]),
#        proteins = "proteins/{protein_set}.fasta",
        #proteins = "proteins/{protein_set}.fasta",
#       proteins = expand("proteins/{protein_set}.fasta",, protein_set=PROTEINSETS),    
#    threads: 10
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
        grep -v "^##" miniprot/{wildcards.protein_set}_mp_aln.gff |sort -k1,1 -k4,4n -k5,5n > miniprot/{wildcards.protein_set}_temp.mp.gff
        rm -f miniprot/{wildcards.protein_set}_temp.agat.gff
        agat_convert_sp_gxf2gxf.pl -g miniprot/{wildcards.protein_set}_temp.mp.gff -o miniprot/{wildcards.protein_set}_temp.agat.gff > miniprot/{wildcards.protein_set}_agat_convert.out 2> miniprot/{wildcards.protein_set}_agat_convert.err
        gffread -T -V -g {input.assembly} miniprot/{wildcards.protein_set}_temp.agat.gff > miniprot/{wildcards.protein_set}_temp.agat.gtf
        cat miniprot/{wildcards.protein_set}_temp.agat.gtf |gtf2gff.pl --gff3 --printExon --out=miniprot/{wildcards.protein_set}_temp.convert.gff
        rm -f miniprot/{wildcards.protein_set}_temp.agat2.gff
        agat_convert_sp_gxf2gxf.pl -g miniprot/{wildcards.protein_set}_temp.convert.gff -o miniprot/{wildcards.protein_set}_temp.agat2.gff > miniprot/{wildcards.protein_set}_agat_convert_2.out 2> miniprot/{wildcards.protein_set}_agat_convert_2.err
        python /cluster/projects/nn8013k/scripts/annotation/convert_to_evm_output.py -i miniprot/{wildcards.protein_set}_temp.agat2.gff -o {output.gff}
        """    

#1e: if RNA-Seq data, map against genome with HiSat2, and run StringTie to generate a GTF, convert to EVM input via TransDecoder

rule hisat_stringtie:
    output:
        "stringtie/transcripts.fasta.transdecoder.genome.gff3",
        "stringtie/stringtie.proteins.fa",
    input:
        assembly = expand("{prefix}.fold.fa", prefix=config["prefix"]),
        r1 = config["rna_seq_r1"],
        r2 = config["rna_seq_r2"],
        pfam = config["pfam"],
        uniprot =  "proteins/uniprot.fasta",
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
        mkdir -p stringtie
        cd stringtie
        hisat2-build ../{input.assembly} {params.prefix} 2> hisat-build.err
        
        hisat2 -p {params.threads} --dta -x {params.prefix} \
        -1 {input.r1} \
        -2 {input.r2} \
        --rna-strandness RF \
        2> hisat.err | samtools view -buS - | \
        samtools sort -m 4G -@ {params.threads} -T tmp -O bam - > hisat2.sort.bam 2> samtools_hisat.err

        stringtie --rf hisat2.sort.bam > stringtie.gtf 2> stringtie.err

        gffread -E stringtie.gtf -o- > stringtie.gff 
        gffread -w stringtie.fa -g ../{input.assembly} stringtie.gtf
        #agat does not work here for some reason
        #agat_sp_extract_sequences.pl --gff stringtie.gff -f ../{input.assembly} -p -o stringtie.proteins.fa 
        
        gtf_genome_to_cdna_fasta.pl stringtie.gtf ../{input.assembly} > transcripts.fasta 2> gtf_genome_to_cdna_fasta.err

        gtf_to_alignment_gff3.pl stringtie.gtf > transcripts.gff3 2> transcripts.err

        TransDecoder.LongOrfs -t transcripts.fasta --output_dir transdecoder_dir 1> transdecoder_longorfs.out 2> transdecoder_longorfs.err

        hmmsearch --cpu {params.threads} --domtblout pfam.domtblout {input.pfam} transdecoder_dir/longest_orfs.pep 1> hmmscan.out 2> hmmscan.err

        diamond blastp \
        --query transdecoder_dir/longest_orfs.pep  \
        --db ../{input.uniprot} \
	--outfmt 6 \
        --sensitive \
        --max-target-seqs 1 \
        --evalue 1e-5 \
        --threads 20 \
        > blastp.outfmt6 2> diamond_blastp.err

        TransDecoder.Predict -t transcripts.fasta --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6 --output_dir transdecoder_dir 1> transdecoder_predict.out 2> transdecoder_predict.err

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
        annotation = "galba/augustus.hints.gff",
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
        mkdir -p galba 
        singularity exec -B $PWD:/data /cluster/projects/nn8013k/opt/galba/galba.sif cp -rf /usr/share/augustus/config /data/galba
        singularity exec -B $PWD:/data /cluster/projects/nn8013k/opt/galba/galba.sif  galba.pl --species={params.prefix} --threads={params.threads} --genome=/data/{input.assembly} \
         --verbosity=4 --prot_seq=/data/{input.proteins}  --workingdir=/data/galba \
         --augustus_args="--stopCodonExcludedFromCDS=False" --gff3 \
         --AUGUSTUS_CONFIG_PATH=/data/galba/config \
         1> galba/galba_"`date +\%y\%m\%d_\%H\%M\%S`".out 2> galba/galba_"`date +\%y\%m\%d_\%H\%M\%S`".err
        
        agat_sp_extract_sequences.pl --gff galba/augustus.hints.gff -f {input.assembly}  -t cds -p -o galba/galba.proteins.fa
        """    

def get_stringtie(wildcards):
    if rna_seq:
        return "stringtie/transcripts.fasta.transdecoder.genome.gff3"
    else:
        return ""


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
        galba = "galba/augustus.hints.gff",
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
        /cluster/projects/nn8013k/programs/miniconda3/envs/anno_pipeline/opt/evidencemodeler-1.1.1/EvmUtils/misc/augustus_GTF_to_EVM_GFF3.pl {input.galba} > evm/galba.evm.gff3 
        cat evm/galba.evm.gff3 > evm/gene_predictions.gff3 
        cat {input.orthodb} >> evm/gene_predictions.gff3 
        cat {input.model} >> evm/gene_predictions.gff3 
        cat {input.uniprot} >> evm/gene_predictions.gff3
        cat {input.stringtie} >> evm/gene_predictions.gff3

        printf "ABINITIO_PREDICTION\tAugustus\t1\n" > evm/weights.evm.txt 
        printf "OTHER_PREDICTION\tminiprot\t4\n" >> evm/weights.evm.txt
        printf "OTHER_PREDICTION\ttransdecoder\t6\n" >> weights.evm.txt
        
        cd evm
        
        /cluster/projects/nn8013k/programs/miniconda3/envs/anno_pipeline/lib/python3.8/site-packages/funannotate/aux_scripts/funannotate-runEVM.py \
        -w weights.evm.txt \
        -c {params.threads} \
        -d . \
        -g gene_predictions.gff3 \
        -f ../{input.masked} \
        -l evm.log \
        -m 10 \
        -i 1500 \
        -o evm.gff3 \
        --EVM_HOME /cluster/projects/nn8013k/programs/miniconda3/envs/anno_pipeline/opt/evidencemodeler-1.1.1 \
        1> evm_"`date +\%y\%m\%d_\%H\%M\%S`".out 2> evm_"`date +\%y\%m\%d_\%H\%M\%S`".err 
               
        #for some stupid reason, funannotate-runEVM.py does not create the GFF file...
	#"scripts/concat_gff.py"        

        agat_sp_extract_sequences.pl --gff evm.gff3 -f ../{input.masked}  -t cds -p -o evm.proteins.fa 1> agat_proteins_"`date +\%y\%m\%d_\%H\%M\%S`".out 2> agat_proteins_"`date +\%y\%m\%d_\%H\%M\%S`".err
        agat_sp_extract_sequences.pl --gff evm.gff3 -f ../{input.masked}  -t exon --merge -o evm.mrna.fa 1> agat_exons_"`date +\%y\%m\%d_\%H\%M\%S`".out 2> agat_exons_"`date +\%y\%m\%d_\%H\%M\%S`".err
        """
#4: filter annotated proteins based on matches to repeat proteins from Funannotate
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

        cut -f 1 filter/repeats.tsv |sort -u > filter/kill.list

        rm -rf filter/removed_repeats*
        agat_sp_filter_feature_from_kill_list.pl --gff {input.evm_gff} --kill_list filter/kill.list -o filter/removed_repeats.gff

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
        expand("{prefix}.embl", prefix=config["prefix"]),
    params:
        locus_tag = config["locus_tag"],
        species = config["species"],
        project_id = config["project_id"],
        rg = "EBP-Nor",
        prefix = config["prefix"],
    shell:
        r"""
        EMBLmyGFF3 \
        {input.annotation} \
        {input.assembly} \
        --topology linear \
        --molecule_type 'genomic DNA' \
        --transl_table 1  \
        --rg {params.rg} \
        --species {params.species} \
        --locus_tag {params.locus_tag} \
        --project_id {params.project_id} \
        -o {params.prefix}.embl > create_hap1.out 2> create_hap1.err
        """

#8: busco
rule busco:
    output:
        directory(expand("{{dir}}/busco_{{analysis}}_{lineage}", lineage=config["busco_lineage"])),
    input:
        proteins = "{dir}/{analysis}.proteins.fa",
        busco = expand("{busco_lib}/{lineage}_odb10", busco_lib=config["busco_lib"], lineage=config["busco_lineage"])
    resources:
        #time = config["fcs_gx_time"],
        mem_per_cpu = "7G",
        #partition = config["fcs_gx_partition"],
        ntasks = 10
    params:
        lineage = config["busco_lineage"],
        busco_lib = config["busco_lib"],
    shell:
        r"""
        busco --tar -c 10 -i {input.proteins} -l {input.busco} -o busco_{wildcards.analysis}_{params.lineage} --out_path {wildcards.dir} \
        -m proteins --offline 
        """
