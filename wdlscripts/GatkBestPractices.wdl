workflow GatkBestPractices {

    # Common options for every GATK task
    File gatkJar                # GATK .jar file
    File refFasta               # Reference genome fasta file
    File refIndex               # Index for reference genome fasta
    File refSeqDict             # Dict file for reference genome
    File targetedIntervals      # BED file of targeted genomic intervals
    String xmx                  # JVM xmx option for GATK jobs
    String xms                  # JVM xms option for GATK jobs
    String xmn                  # JVM xmn option for GATK jobs

    # Other options
    File sampleTable            # Table of sample, bam file
    Array[Array[String]] sampleToBam = read_tsv(sampleTable)    # Sample, bam file
    File thousandGenomesSNPs    # VCF file of known SNPs from 1000 Genomes project
    File hapmap                 # VCF file of known variants from HapMap
    File dbSNP                  # VCF file of known variants from dbSNP
    File millsIndels            # VCF file of known indels from Mills
    File omni                   # VCF file of known variants from Omni

    # Per-sample data processing
    scatter (sampleBam in sampleToBam) {

        call realignerTargetCreator {
            input:
                gatk = gatkJar,
                ref = refFasta,
                refIdx = refIndex,
                refDict = refSeqDict,
                intervals = targetedIntervals,
                xmx = xmx,
                xms = xms,
                xmn = xmn,
                sample = sampleBam[0],
                bam = sampleBam[1],
                bamIdx = sampleBam[2],
                indels = millsIndels
        }

        call indelRealigner {
            input:
                gatk = gatkJar,
                ref = refFasta,
                refIdx = refIndex,
                refDict = refSeqDict,
                intervals = targetedIntervals,
                xmx = xmx,
                xms = xms,
                xmn = xmn,
                sample = sampleBam[0],
                bam = sampleBam[1],
                bamIdx = sampleBam[2],
                indels = millsIndels,
                realignerTargets = realignerTargetCreator.realignerTargets
        }



    }

}


task realignerTargetCreator {
# https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_indels_RealignerTargetCreator.php

    # Common options
    File gatk
    File ref
    File refIdx
    File refDict
    File intervals
    String xmx
    String xms
    String xmn

    # Other options
    String sample
    File bam
    File bamIdx
    File indels

    command {
        java \
        -Xmx${xmx} -Xms${xms} -Xmn${xmn} -jar \
        ${gatk} \
        -R ${ref} \
        -L ${intervals} \
        -T RealignerTargetCreator \
        -I ${bam} \
        --known ${indels} \
        -o ${sample}.interval_list
    }

    output {
        File realignerTargets = "${sample}.interval_list"
    }
}

task indelRealigner {
# https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_indels_IndelRealigner.php

    # Common options
    File gatk
    File ref
    File refIdx
    File refDict
    File intervals
    String xmx
    String xms
    String xmn

    # Other options
    File realignerTargets # Output from RealignerTargetCreator
    String sample
    File bam
    File bamIdx
    File indels

    command {
        java \
        -Xmx${xmx} -Xms${xms} -Xmn${xmn} -jar \
        ${gatk} \
        -R ${ref} \
        -L ${intervals} \
        -T IndelRealigner \
        -I ${bam} \
        -known ${indels} \
        -targetIntervals ${realignerTargets} \
        -o ${sample}.realign.bam
    }

    output {
        File realignedBam = "${sample}.realign.bam"
    }

}


}

