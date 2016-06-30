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

        call haplotypeCaller {
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
                bam = printReads.out,
                bamIdx = printReads.outIdx,
                emitRefConfidence = "GVCF",
                dbSNP = dbSNP
        }

    }

    call combineGVCFs {
        input:
            gatk = gatk,
            ref = ref,
            refIdx = refIdx,
            refDict = refDict,
            intervals = intervals,
            xmx = xmx,
            xms = xms,
            xmn = xmn,
            gVCFs = haplotypeCaller.out,
            dbSNP = dbSNP
    }

    call genotypeGVCFs {
        input:
            gatk = gatk,
            ref = ref,
            refIdx = refIdx,
            refDict = refDict,
            intervals = intervals,
            xmx = xmx,
            xms = xms,
            xmn = xmn,
            variant = combineGVCFs.out,
            dbSNP = dbSNP
    }

    call variantRecalibratorSNPs {
        input:
            gatk = gatk,
            ref = ref,
            refIdx = refIdx,
            refDict = refDict,
            intervals = intervals,
            xmx = xmx,
            xms = xms,
            xmn = xmn,
            vcf = genotypeGVCFs.out,
            nt = 4,
            hapmap = hapmap,
            omni = omni,
            thousandGenomes = thousandGenomesSNPs,
            dbSNP = dbSNP,
            hapmapParam = "known=false,training=true,truth=true,prior=15.0",
            omniParam = "known=false,training=true,truth=true,prior=12.0",
            thousandGenomesParam = "known=false,training=true,truth=false,prior=10.0",
            dbSNPParam = "known=true,training=false,truth=false,prior=2.0",
            an = ["QD", "MQ", "MQRankSum", "ReadPosRankSum", "FS", "SQR", "InbreedingCoeff"],
            mode = "SNP"
    }

    call variantRecalibratorIndels {
        input:
            gatk = gatk,
            ref = ref,
            refIdx = refIdx,
            refDict = refDict,
            intervals = intervals,
            xmx = xmx,
            xms = xms,
            xmn = xmn,
            vcf = genotypeGVCFs.out,
            nt = 4,
            maxGaussians = 4,
            dbSNP = dbSNP,
            millsIndels = millsIndels,
            dbSNPParam = "known=true,training=false,truth=false,prior=2.0",
            millsParam = "known=false,training=true,truth=true,prior=12.0",
            an = ["QD", "FS", "SQR", "ReadPosRankSum", "MQRankSum", "InbreedingCoeff"],
            mode = "INDEL"
    }

    call applyRecalibrationSNPs {
         input:
             gatk = gatk,
             ref = ref,
             refIdx = refIdx,
             refDict = refDict,
             intervals = intervals,
             xmx = xmx,
             xms = xms,
             xmn = xmn,
             vcf = genotypeGVCFs.out,
             tranchesFile = variantRecalibratorSNPs.tranches,
             recalFile = variantRecalibratorSNPs.recal,
             tsFilterLevel = 99.5,
             mode = "SNP"

    }

    call applyRecalibrationIndels {
         input:
             gatk = gatk,
             ref = ref,
             refIdx = refIdx,
             refDict = refDict,
             intervals = intervals,
             xmx = xmx,
             xms = xms,
             xmn = xmn,
             vcf = applyRecalibrationSNPs.out,
             tranchesFile = variantRecalibratorIndels.tranches,
             recalFile = variantRecalibratorIndels.recal,
             tsFilterLevel = 99.0,
             mode = "INDEL"

    }

    call calculateGenotypePosteriors {
         input:
             gatk = gatk,
             ref = ref,
             refIdx = refIdx,
             refDict = refDict,
             intervals = intervals,
             xmx = xmx,
             xms = xms,
             xmn = xmn,
             vcf = applyRecalibrationIndels.out
    }

    call variantFiltration {
         input:
             gatk = gatk,
             ref = ref,
             refIdx = refIdx,
             refDict = refDict,
             intervals = intervals,
             xmx = xmx,
             xms = xms,
             xmn = xmn,
             vcf = calculateGenotypePosteriors.out,
             filterName = "lowGQ",
             filterExpression = "\"GQ < 20.0\""
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
        File outIdx = "${sample}.recalibrated.bai"
    }

}


# **********************************************************************
# VARIANT DISCOVERY
# https://www.broadinstitute.org/gatk/guide/bp_step.php?p=2
# **********************************************************************

task haplotypeCaller {
# haplotypeCaller: Call germline SNPs and indels via local re-assembly of haplotypes
# For cohort mode, call variants per sample then combine with CombineGVCFs
# https://www.broadinstitute.org/gatk/guide/article?id=3893
# https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php

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
    File bam # processed bam
    File bamIdx
    String emitRefConfidence # mode for emitting reference confidence scores e.g. GVCF
    File dbSNP
    String sample

    command {
        java \
        -Xmx${xmx} -Xms${xms} -Xmn${xmn} -jar \
        ${gatk} \
        -R ${ref} \
        -L ${intervals} \
        -T HaplotypeCaller \
        -I ${bam} \
        --emitRefConfidence ${emitRefConfidence} \
        --dbsnp ${dbSNP} \
        -o ${sample}.recalibrated.raw.snps.indels.g.vcf
    }

    output {
        File out = "${sample}.recalibrated.raw.snps.indels.g.vcf"
    }


}

task combineGVCFs {
# CombineGVCFs: Combine per-sample gVCF files produced by HaplotypeCaller into a multi-sample gVCF file
# https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_variantutils_CombineGVCFs.php

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
    Array[File] gVCFs # Output from haplotypeCaller
    File dbSNP

    command {
        java \
        -Xmx${xmx} -Xms${xms} -Xmn${xmn} -jar \
        ${gatk} \
        -R ${ref} \
        -L ${intervals} \
        -T CombineGVCFs \
        -V ${sep = " -V " gVCFs} \
        -o "multisample.g.vcf"
    }

    output {
        File out = "multisample.g.vcf"
    }


}

task genotypeGVCFs {
# GenotypeGVCFs: Perform joint genotyping on gVCF files produced by HaplotypeCaller
# https://www.broadinstitute.org/gatk/guide/article?id=3893
# https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_variantutils_GenotypeGVCFs.php

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
    File variant # Output from combineGVCFs
    File dbSNP

    command {
        java \
        -Xmx${xmx} -Xms${xms} -Xmn${xmn} -jar \
        ${gatk} \
        -R ${ref} \
        -L ${intervals} \
        -T GenotypeGVCFs \
        -V ${variant} \
        -o "multisample.genotyped.vcf"
    }

    output {
        File out = "multisample.genotyped.vcf"
    }

}

task variantRecalibratorSNPs {
# VariantRecalibrator for SNPs: Build a recalibration model to score variant quality for filtering purposes (SNPs)
# https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_variantrecalibration_VariantRecalibrator.php
# https://www.broadinstitute.org/gatk/guide/article?id=39
# https://www.broadinstitute.org/gatk/guide/article?id=2805
# https://www.broadinstitute.org/gatk/guide/article?id=1259 for parameter recommendations

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
    File vcf # output from GenotypeGVCFs
    Int nt # number of threads
    File hapmap
    File omni
    File thousandGenomes
    File dbSNP
    String hapmapParam
    String omniParam
    String thousandGenomesParam
    String dbSNPParam
    Array[String] an # names of annotations which should be used for calculations
    String mode

    command {
        java \
        -Xmx${xmx} -Xms${xms} -Xmn${xmn} -jar \
        ${gatk} \
        -R ${ref} \
        -L ${intervals} \
        -T VariantRecalibrator \
        -input ${vcf} \
        -nt ${nt} \
        -resource:hapmap,${hapmapParam} ${hapmap} \
        -resource:omni,${omniParam} ${omni} \
        -resource:1000G,${thousandGenomesParam} ${thousandGenomes} \
        -resource:dbsnp,${dbSNPParam} ${dbSNP} \
        -an ${sep = " -an " an} \
        -mode ${mode}
    }

    output {
        File recal = "multisample.genotyped.SNP.recal"
        File tranches = "multisample.genotyped.SNP.tranches"
    }

}

task variantRecalibratorIndels {
# VariantRecalibrator for SNPs: Build a recalibration model to score variant quality for filtering purposes (SNPs)
# https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_variantrecalibration_VariantRecalibrator.php
# https://www.broadinstitute.org/gatk/guide/article?id=39
# https://www.broadinstitute.org/gatk/guide/article?id=2805
# https://www.broadinstitute.org/gatk/guide/article?id=1259 for parameter recommendations

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
    File vcf # output from GenotypeGVCFs
    Int nt # number of threads
    Int maxGaussians
    File dbSNP
    File millsIndels
    String dbSNPParam
    String millsParam
    Array[String] an # names of annotations which should be used for calculations
    String mode

    command {
        java \
        -Xmx${xmx} -Xms${xms} -Xmn${xmn} -jar \
        ${gatk} \
        -R ${ref} \
        -L ${intervals} \
        -T VariantRecalibrator \
        -input ${vcf} \
        -nt ${nt} \
        -resource:dbsnp,${dbSNPParam} ${dbSNP} \
        -resource:millsIndels,${millsParam} ${millsIndels} \
        -an ${sep = " -an " an} \
        -mode ${mode}
    }

    output {
        File recal = "multisample.genotyped.indel.recal"
        File tranches = "multisample.genotyped.indel.tranches"
    }

}

task applyRecalibrationSNPs {
# ApplyRecalibration for SNPs: Apply a score cutoff to filter variants based on a recalibration table (SNPs)
# https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_variantrecalibration_ApplyRecalibration.php
# https://www.broadinstitute.org/gatk/guide/article?id=1259 for parameter recommendations

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
    File vcf # output from GenotypeGVCFs
    File tranchesFile # output from VariantRecalibrator
    File recalFile # output from VariantRecalibrator
    Float tsFilterLevel
    String mode

    command {
        java \
        -Xmx${xmx} -Xms${xms} -Xmn${xmn} -jar \
        ${gatk} \
        -R ${ref} \
        -L ${intervals} \
        -T ApplyRecalibration \
        -input ${vcf} \
        --ts_filter_level ${tsFilterLevel} \
        -tranchesFile ${tranchesFile} \
        -recalFile ${recalFile} \
        -mode ${mode} \
        -o "multisample.genotyped.recalibrated.filtered.SNPs_only.vcf"
    }

    output {
        File out = "multisample.genotyped.recalibrated.filtered.SNPs_only.vcf"
    }

}

task applyRecalibrationIndels {
# ApplyRecalibration for indels: Apply a score cutoff to filter variants based on a recalibration table (SNPs)
# https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_variantrecalibration_ApplyRecalibration.php
# https://www.broadinstitute.org/gatk/guide/article?id=1259 for parameter recommendations

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
    File vcf # output from ApplyRecalibration for SNPs
    File tranchesFile # output from VariantRecalibrator
    File recalFile # output from VariantRecalibrator
    Float tsFilterLevel
    String mode

    command {
        java \
        -Xmx${xmx} -Xms${xms} -Xmn${xmn} -jar \
        ${gatk} \
        -R ${ref} \
        -L ${intervals} \
        -T ApplyRecalibration \
        -input ${vcf} \
        --ts_filter_level ${tsFilterLevel} \
        -tranchesFile ${tranchesFile} \
        -recalFile ${recalFile} \
        -mode ${mode} \
        -o "multisample.genotyped.recalibrated.filtered.vcf"
    }

    output {
        File out = "multisample.genotyped.recalibrated.filtered.vcf"
    }

}


# *********************************************************************
#                        CALLSET REFINEMENT
# https://www.broadinstitute.org/gatk/guide/bp_step.php?p=3
# https://www.broadinstitute.org/gatk/guide/article?id=4723
# *********************************************************************

task calculateGenotypePosteriors {
# CalculateGenotypePosteriors: Calculate genotype posterior likelihoods given panel data
# https://www.broadinstitute.org/gatk/guide/article?id=4727
# https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_variantutils_CalculateGenotypePosteriors.php
# Here, we follow the method to "refine the genotypes of a large panel based on the discovered allele frequency"

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
    File vcf # output from ApplyRecalibration for indels

    command {
        java \
        -Xmx${xmx} -Xms${xms} -Xmn${xmn} -jar \
        ${gatk} \
        -R ${ref} \
        -L ${intervals} \
        -T ApplyRecalibration \
        -V ${vcf} \
        -o "multisample.genotyped.recalibrated.filtered.withPosteriors.vcf"
    }

    output {
        File out = "multisample.genotyped.recalibrated.filtered.withPosteriors.vcf"
    }

}

task variantFiltration {
# VariantFiltration: Filter variant calls based on INFO and FORMAT annotations
# https://www.broadinstitute.org/gatk/guide/article?id=4727
# https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_filters_VariantFiltration.php
# Use recommendations here: https://www.broadinstitute.org/gatk/guide/article?id=4727

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
    File vcf # output from calculateGenotypePosteriors
    String filterName
    String filterExpression # in double quotes

    command {
        java \
        -Xmx${xmx} -Xms${xms} -Xmn${xmn} -jar \
        ${gatk} \
        -R ${ref} \
        -L ${intervals} \
        -T VariantFiltration \
        -V ${vcf} \
        -G_filter ${filterExpression} \
        -G_filterName ${filterName} \
        -o "multisample.genotyped.recalibrated.filtered.withPosteriors.Gfiltered.vcf"
    }

    output {
        File out = "multisample.genotyped.recalibrated.filtered.withPosteriors.Gfiltered.vcf"
    }


}