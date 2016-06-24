workflow GatkBestPractices {

    # Common options for every GATK task
    File gatk                   # GATK .jar file
    File ref                    # Reference genome fasta file
    File refIdx                 # Index for reference genome fasta
    File refDict                # Dict file for reference genome
    File intervals              # BED file of targeted genomic intervals
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
    # http://gatkforums.broadinstitute.org/wdl/discussion/3441/queue-how-to-connect-gatk-walkers
    scatter (sampleBam in sampleToBam) {

        call realignerTargetCreator {
            input:
                gatk = gatk,
                ref = ref,
                refIdx = refIdx,
                refDict = refDict,
                intervals = intervals,
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
                gatk = gatk,
                ref = ref,
                refIdx = refIdx,
                refDict = refDict,
                intervals = intervals,
                xmx = xmx,
                xms = xms,
                xmn = xmn,
                sample = sampleBam[0],
                bam = sampleBam[1],
                bamIdx = sampleBam[2],
                indels = millsIndels,
                realignerTargets = realignerTargetCreator.out
        }

        call baseRecalibratorBefore {
            input:
                gatk = gatk,
                ref = ref,
                refIdx = refIdx,
                refDict = refDict,
                intervals = intervals,
                xmx = xmx,
                xms = xms,
                xmn = xmn,
                realignedBam = indelRealigner.out,
                realignedBamIdx = indelRealigner.outIdx,
                knownSites = dbSNP,
                sample = sampleBam[0]
        }

        call baseRecalibratorAfter {
            input:
                gatk = gatk,
                ref = ref,
                refIdx = refIdx,
                refDict = refDict,
                intervals = intervals,
                xmx = xmx,
                xms = xms,
                xmn = xmn,
                realignedBam = indelRealigner.out,
                realignedBamIdx = indelRealigner.outIdx,
                knownSites = dbSNP,
                sample = sampleBam[0],
                bqsr = baseRecalibratorBefore.out
        }

        call analyzeCovariates {
            input:
                gatk = gatk,
                ref = ref,
                refIdx = refIdx,
                refDict = refDict,
                intervals = intervals,
                xmx = xmx,
                xms = xms,
                xmn = xmn,
                before = baseRecalibratorBefore.out,
                after = baseRecalibratorAfter.out,
                sample = sampleBam[0]
        }

        call printReads {
            input:
                gatk = gatk,
                ref = ref,
                refIdx = refIdx,
                refDict = refDict,
                intervals = intervals,
                xmx = xmx,
                xms = xms,
                xmn = xmn,
                sample = sampleBam[0],
                bam = sampleBam[1],
                bamIdx = sampleBam[2],
                bqsr = baseRecalibratorAfter.out
        }

    }

}

# **********************************************************************
# TODO map to reference
# https://www.broadinstitute.org/gatk/guide/bp_step.php?p=1
# **********************************************************************


# **********************************************************************
# TODO mark duplicates
# https://www.broadinstitute.org/gatk/guide/bp_step.php?p=1
# **********************************************************************


# **********************************************************************
# LOCAL REALIGNMENT AROUND INDELS
# https://www.broadinstitute.org/gatk/guide/article?id=38
# https://www.broadinstitute.org/gatk/guide/article?id=2800
# **********************************************************************


task realignerTargetCreator {
# RealignerTargetCreator: Define intervals to target for local realignment
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
        File out = "${sample}.interval_list"
    }
}

task indelRealigner {
# IndelRealigner: Perform local realignment of reads around indels
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
        File out = "${sample}.realign.bam"
        File outIdx = "${sample}.realign.bai"
    }

}

# **********************************************************************
# BASE QUALITY SCORE RECALIBRATION
# https://www.broadinstitute.org/gatk/guide/article?id=44
# https://www.broadinstitute.org/gatk/guide/article?id=2801
# **********************************************************************


task baseRecalibratorBefore {
# BaseRecalibrator: Generate base recalibration table to compensate for systematic errors in basecalling confidences
# Generate the first pass recalibration table file
# https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_bqsr_BaseRecalibrator.php

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
    File realignedBam # Output from IndelRealigner
    File realignedBamIdx # Output from IndelRealigner
    File knownSites
    String sample

    command {
        java \
        -Xmx${xmx} -Xms${xms} -Xmn${xmn} -jar \
        ${gatk} \
        -R ${ref} \
        -L ${intervals} \
        -T BaseRecalibrator \
        -I ${realignedBam} \
        -knownSites ${knownSites} \
        -o ${sample}.base_recalibrator_first_pass.out
    }

    output {
        File out = "${sample}.base_recalibrator_first_pass.out"
    }

}

task baseRecalibratorAfter {
# BaseRecalibrator: Generate base recalibration table to compensate for systematic errors in basecalling confidences
# Generate the second pass recalibration table file
# https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_bqsr_BaseRecalibrator.php

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
    File realignedBam # Output from IndelRealigner
    File realignedBamIdx # Output from IndelRealigner
    File knownSites
    String sample
    File bqsr # Output from BaseRecalibrator first pass

    command {
        java \
        -Xmx${xmx} -Xms${xms} -Xmn${xmn} -jar \
        ${gatk} \
        -R ${ref} \
        -L ${intervals} \
        -T BaseRecalibrator \
        -I ${realignedBam} \
        -knownSites ${knownSites} \
        -BQSR ${bqsr} \
        -o ${sample}.base_recalibrator_second_pass.out
    }

    output {
        File out = "${sample}.base_recalibrator_second_pass.out"
    }

}

task analyzeCovariates {
# AnalyzeCovariates: Create plots to visualize base recalibration results
# https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_bqsr_AnalyzeCovariates.php

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
    File before # Output from BaseRecalibrator first pass
    File after # Output from BaseRecalibrator second pass
    String sample

    command {
        java \
        -Xmx${xmx} -Xms${xms} -Xmn${xmn} -jar \
        ${gatk} \
        -R ${ref} \
        -L ${intervals} \
        -T AnalyzeCovariates \
        -before ${before} \
        -after ${after} \
        -plots analyze_covariates_${sample}.realign.BQSR.pdf \
        -csv analyze_covariates_${sample}.realign.BQSR.csv
    }

    output {
        File plots = "analyze_covariates_${sample}.realign.BQSR.pdf"
        File csv = "analyze_covariates_${sample}.realign.BQSR.csv"
    }

}

task printReads {
# PrintReads: Write out sequence read data (for filtering, merging, subsetting etc)
# https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_readutils_PrintReads.php

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
    File bam # Original bam
    File bamIdx
    File bqsr # Output from BaseRecalibrator second pass
    String sample

    command {
        java \
        -Xmx${xmx} -Xms${xms} -Xmn${xmn} -jar \
        ${gatk} \
        -R ${ref} \
        -L ${intervals} \
        -T PrintReads \
        -I ${bam} \
        -BQSR ${bqsr} \
        -O ${sample}.recalibrated.bam
    }

    output {
        File out = "${sample}.recalibrated.bam"
    }

}
