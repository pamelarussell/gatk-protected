package qscripts

import org.broadinstitute.gatk.queue.QScript
//import org.broadinstitute.gatk.queue.extensions.picard.MarkDuplicates
import org.broadinstitute.gatk.queue.extensions.gatk._

/**
  * Created by prussell on 3/3/16.
  * Implements GATK best practices
  */
class GatkBestPractices extends QScript {

  /**
    * *********************************************************************
    *                               INPUTS
    * *********************************************************************
    */

  /**
    * Reference genome fasta
    */
  @Input(doc="Reference genome for the bam files", shortName="R", fullName="REF_FASTA", required=true)
  var referenceFile: File = null

  /**
    * Bam files
    */
  @Input(doc="One or more bam files", shortName="I", fullName="INPUT_BAM", required=true)
  var bamFiles: List[File] = Nil

  /**
    * VCF files
    */
  @Input(doc="VCF file of known SNPs from 1000 Genomes project", fullName="THOUSAND_GENOMES_SNPS", required=true)
  var thousandGenomesSNPs: File = null
  @Input(doc="VCF file of known variants from HapMap", fullName="HAPMAP", required=true)
  var hapmap: File = null
  @Input(doc="VCF file of known variants from dbSNP", fullName="DBSNP", required=true)
  var dbSNP: File = null
  @Input(doc="VCF file of known indels from Mills", fullName="MILLS_INDELS", required=true)
  var millsIndels: File = null
  @Input(doc="VCF file of known variants from Omni", fullName="OMNI", required=true)
  var omni: File = null

  /**
    * Targeted intervals
    */
  @Input(doc="File of targeted intervals", shortName="ti", fullName="TARGETED_INTERVALS", required=true)
  var intervalsFile : File = null

  /**
    * Output directory
    */
  @Input(doc="Output directory", shortName="od", fullName="OUT_DIR", required=true)
  var outputDirectory : File = null

  /**
    * Output file prefix not including directory
    */
  @Input(doc="Output prefix not including directory", shortName="op", fullName="OUT_PREFIX", required=true)
  var outputPrefix : File = null

  /**
    * Common arguments
    */
  trait CommonArguments extends CommandLineGATK {
    this.reference_sequence = referenceFile
    this.intervals = Seq(intervalsFile)
  }

  def script() = {

    /**
      * *********************************************************************
      *                         MAP TO REFERENCE
      * https://www.broadinstitute.org/gatk/guide/bp_step.php?p=1
      * *********************************************************************
      */

    // TODO

    /**
      * *********************************************************************
      *                         MARK DUPLICATES
      * https://www.broadinstitute.org/gatk/guide/bp_step.php?p=1
      * *********************************************************************
      */

    // TODO


    /**
      * Container for the processed bam files
      */
    var processedFiles = Seq.empty[File]

    /**
      * **********************************************************************************************
      *                              PER-SAMPLE DATA PROCESSING
      * http://gatkforums.broadinstitute.org/wdl/discussion/3441/queue-how-to-connect-gatk-walkers
      * **********************************************************************************************
      */
    for(bam <- bamFiles) {

      /**
        * *********************************************************************
        * LOCAL REALIGNMENT AROUND INDELS
        * https://www.broadinstitute.org/gatk/guide/article?id=38
        * https://www.broadinstitute.org/gatk/guide/article?id=2800
        * *********************************************************************
        */

      /**
        * Local realignment around indels
        * RealignerTargetCreator: Define intervals to target for local realignment
        * https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_indels_RealignerTargetCreator.php
        */
      val realignerTargetCreator = new RealignerTargetCreator with CommonArguments
      realignerTargetCreator.input_file +:= bam
      realignerTargetCreator.known = Seq(millsIndels)
      realignerTargetCreator.out = swapExt(bam, "bam", "interval_list")
      add(realignerTargetCreator)

      /**
        * Local realignment around indels
        * IndelRealigner: Perform local realignment of reads around indels
        * https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_indels_IndelRealigner.php
        */
      val indelRealigner = new IndelRealigner with CommonArguments
      indelRealigner.targetIntervals = realignerTargetCreator.out
      indelRealigner.input_file +:= bam
      indelRealigner.knownAlleles = Seq(millsIndels)
      indelRealigner.out = swapExt(bam, "bam", "realign.bam")
      add(indelRealigner)

      /**
        * *********************************************************************
        * BASE QUALITY SCORE RECALIBRATION
        * https://www.broadinstitute.org/gatk/guide/article?id=44
        * https://www.broadinstitute.org/gatk/guide/article?id=2801
        * *********************************************************************
        */

      /**
        * Base quality score recalibration
        * BaseRecalibrator: Generate base recalibration table to compensate for systematic errors in basecalling confidences
        * https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_bqsr_BaseRecalibrator.php
        */

      // Generate the first pass recalibration table file
      val baseRecalibratorBefore = new BaseRecalibrator with CommonArguments
      baseRecalibratorBefore.input_file +:= indelRealigner.out
      baseRecalibratorBefore.out = swapExt(indelRealigner.out, "realign.bam", "base_recalibrator_first_pass.out")
      baseRecalibratorBefore.knownSites = Seq(dbSNP)
      add(baseRecalibratorBefore)

      // Generate the second pass recalibration table file
      val baseRecalibratorAfter = new BaseRecalibrator with CommonArguments
      baseRecalibratorAfter.BQSR = baseRecalibratorBefore.out
      baseRecalibratorAfter.input_file +:= indelRealigner.out
      baseRecalibratorAfter.out = swapExt(indelRealigner.out, "realign.bam", "base_recalibrator_second_pass.out")
      baseRecalibratorAfter.knownSites = Seq(dbSNP)
      add(baseRecalibratorAfter)

      /**
        * Base quality score recalibration
        * AnalyzeCovariates: Create plots to visualize base recalibration results
        * https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_bqsr_AnalyzeCovariates.php
        */
      val analyzeCovariates = new AnalyzeCovariates with CommonArguments
      analyzeCovariates.beforeReportFile = baseRecalibratorBefore.out
      analyzeCovariates.afterReportFile = baseRecalibratorAfter.out
      analyzeCovariates.plotsReportFile = new File(outputDirectory.getAbsolutePath + "/analyzeCovariates_" + swapExt(indelRealigner.out, "bam", "BQSR.pdf"))
      analyzeCovariates.intermediateCsvFile = new File(outputDirectory.getAbsolutePath + "/analyzeCovariates_" + swapExt(indelRealigner.out, "bam", "BQSR.csv"))
      analyzeCovariates.ignoreLMT = false
      add(analyzeCovariates)

      /**
        * Base quality score recalibration
        * PrintReads: Write out sequence read data (for filtering, merging, subsetting etc)
        * https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_readutils_PrintReads.php
        */
      val printReads = new PrintReads with CommonArguments
      printReads.input_file +:= bam
      printReads.BQSR = baseRecalibratorAfter.out
      printReads.out = swapExt(bam, "bam", "recalibrated.bam")
      add(printReads)

      processedFiles +:= printReads.out

    }

    /**
      * *********************************************************************
      *                        VARIANT DISCOVERY
      * https://www.broadinstitute.org/gatk/guide/bp_step.php?p=2
      * *********************************************************************
      */

    /**
      * Variant discovery
      * HaplotypeCaller: Call germline SNPs and indels via local re-assembly of haplotypes
      * For cohort mode, call variants per sample then combine with CombineGVCFs
      * https://www.broadinstitute.org/gatk/guide/article?id=3893
      * https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php
      */

    /**
      * Container for the processed bam files
      */
    var sampleGVCFs = Seq.empty[File]

    for(processedBam <- processedFiles) {

      val haplotypeCaller = new HaplotypeCaller with CommonArguments
      haplotypeCaller.input_file = Seq(processedBam)
      haplotypeCaller.emitRefConfidence = org.broadinstitute.gatk.tools.walkers.haplotypecaller.ReferenceConfidenceMode.GVCF
      haplotypeCaller.dbsnp = dbSNP
      haplotypeCaller.out = swapExt(processedBam, "bam", "raw.snps.indels.g.vcf")
      add(haplotypeCaller)
      sampleGVCFs +:= haplotypeCaller.out
     }

    /**
      * Variant discovery
      * CombineGVCFs: Combine per-sample gVCF files produced by HaplotypeCaller into a multi-sample gVCF file
      * https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_variantutils_CombineGVCFs.php
      */
    val combineGVCFs = new CombineGVCFs with CommonArguments
    combineGVCFs.variant = sampleGVCFs
    combineGVCFs.dbsnp = dbSNP
    combineGVCFs.out = "multisample.g.vcf"
    combineGVCFs.breakBandsAtMultiplesOf = 0 // Default 0
    combineGVCFs.convertToBasePairResolution = false // Default false
    add(combineGVCFs)

    /**
      * Variant discovery
      * GenotypeGVCFs: Perform joint genotyping on gVCF files produced by HaplotypeCaller
      * https://www.broadinstitute.org/gatk/guide/article?id=3893
      * https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_variantutils_GenotypeGVCFs.php
      */
    val genotypeGVCFs = new GenotypeGVCFs with CommonArguments
    genotypeGVCFs.variant = Seq(combineGVCFs.out)
    genotypeGVCFs.out = swapExt(combineGVCFs.out, "g.vcf", "genotyped.vcf")
    genotypeGVCFs.dbsnp = dbSNP
    add(combineGVCFs)


    /**
      * Variant discovery
      * VariantFiltration: Filter variant calls based on INFO and FORMAT annotations
      * https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_filters_VariantFiltration.php
      */
    // Filter with VQSR instead

    /**
      * Variant discovery
      * VariantRecalibrator for SNPs: Build a recalibration model to score variant quality for filtering purposes (SNPs)
      * https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_variantrecalibration_VariantRecalibrator.php
      * https://www.broadinstitute.org/gatk/guide/article?id=39
      * https://www.broadinstitute.org/gatk/guide/article?id=2805
      * https://www.broadinstitute.org/gatk/guide/article?id=1259 for parameter recommendations
      */
    val variantRecalibratorSNPs = new VariantRecalibrator with CommonArguments
    variantRecalibratorSNPs.input = Seq(genotypeGVCFs.out)
    variantRecalibratorSNPs.recalFile = swapExt(genotypeGVCFs.out, "vcf", "SNP.recal")
    variantRecalibratorSNPs.tranchesFile = swapExt(variantRecalibratorSNPs.recalFile, "recal", "tranches")
    variantRecalibratorSNPs.nt = 4
    variantRecalibratorSNPs.resource = Seq(
      new org.broadinstitute.gatk.queue.extensions.gatk.TaggedFile(hapmap.getAbsolutePath, "known=false,training=true,truth=true,prior=15.0"),
      new org.broadinstitute.gatk.queue.extensions.gatk.TaggedFile(omni.getAbsolutePath, "known=false,training=true,truth=true,prior=12.0"),
      new org.broadinstitute.gatk.queue.extensions.gatk.TaggedFile(thousandGenomesSNPs.getAbsolutePath, "known=false,training=true,truth=false,prior=10.0"),
      new org.broadinstitute.gatk.queue.extensions.gatk.TaggedFile(dbSNP.getAbsolutePath, "known=true,training=false,truth=false,prior=2.0")
    )
    variantRecalibratorSNPs.an = Seq("QD", "MQ", "MQRankSum", "ReadPosRankSum", "FS", "SQR", "InbreedingCoeff")
    variantRecalibratorSNPs.mode = org.broadinstitute.gatk.tools.walkers.variantrecalibration.VariantRecalibratorArgumentCollection.Mode.SNP
    add(variantRecalibratorSNPs)

    /**
      * Variant discovery
      * VariantRecalibrator for indels: Build a recalibration model to score variant quality for filtering purposes (indels)
      * https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_variantrecalibration_VariantRecalibrator.php
      * https://www.broadinstitute.org/gatk/guide/article?id=39
      * https://www.broadinstitute.org/gatk/guide/article?id=2805
      * https://www.broadinstitute.org/gatk/guide/article?id=1259 for parameter recommendations
      */
    val variantRecalibratorIndels = new VariantRecalibrator with CommonArguments
    variantRecalibratorIndels.input = Seq(genotypeGVCFs.out)
    variantRecalibratorIndels.recalFile = swapExt(genotypeGVCFs.out, "vcf", "indel.recal")
    variantRecalibratorIndels.tranchesFile = swapExt(variantRecalibratorIndels.recalFile, "recal", "tranches")
    variantRecalibratorIndels.nt = 4
    variantRecalibratorIndels.maxGaussians = 4
    variantRecalibratorIndels.resource = Seq(
      new org.broadinstitute.gatk.queue.extensions.gatk.TaggedFile(millsIndels.getAbsolutePath, "known=false,training=true,truth=true,prior=12.0"),
      new org.broadinstitute.gatk.queue.extensions.gatk.TaggedFile(dbSNP.getAbsolutePath, "known=true,training=false,truth=false,prior=2.0")
    )
    variantRecalibratorIndels.an = Seq("QD", "FS", "SQR", "ReadPosRankSum", "MQRankSum", "InbreedingCoeff")
    variantRecalibratorIndels.mode = org.broadinstitute.gatk.tools.walkers.variantrecalibration.VariantRecalibratorArgumentCollection.Mode.INDEL
    add(variantRecalibratorIndels)

    /**
      * Variant discovery
      * ApplyRecalibration for SNPs: Apply a score cutoff to filter variants based on a recalibration table (SNPs)
      * https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_variantrecalibration_ApplyRecalibration.php
      * https://www.broadinstitute.org/gatk/guide/article?id=1259 for parameter recommendations
      */
    val applyRecalibrationSNPs = new ApplyRecalibration with CommonArguments
    applyRecalibrationSNPs.input = Seq(genotypeGVCFs.out)
    applyRecalibrationSNPs.tranchesFile = variantRecalibratorSNPs.tranchesFile
    applyRecalibrationSNPs.recalFile = variantRecalibratorSNPs.recalFile
    applyRecalibrationSNPs.out = swapExt(applyRecalibrationSNPs.input.head, "vcf", "recalibrated.filtered.SNPs_only.vcf")
    applyRecalibrationSNPs.ts_filter_level = 99.5
    applyRecalibrationSNPs.mode = org.broadinstitute.gatk.tools.walkers.variantrecalibration.VariantRecalibratorArgumentCollection.Mode.SNP
    add(applyRecalibrationSNPs)

    /**
      * Variant discovery
      * ApplyRecalibration for indels: Apply a score cutoff to filter variants based on a recalibration table (indels)
      * https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_variantrecalibration_ApplyRecalibration.php
      * https://www.broadinstitute.org/gatk/guide/article?id=1259 for parameter recommendations
      */
    val applyRecalibrationIndels = new ApplyRecalibration with CommonArguments
    applyRecalibrationIndels.input = Seq(applyRecalibrationSNPs.out) // Use the VCF file generated by ApplyRecalibration for SNPs
    applyRecalibrationIndels.tranchesFile = variantRecalibratorIndels.tranchesFile
    applyRecalibrationIndels.recalFile = variantRecalibratorIndels.recalFile
    applyRecalibrationIndels.out = swapExt(variantRecalibratorIndels.input.head, "recalibrated.filtered.SNPs_only.vcf", "recalibrated.filtered.vcf")
    applyRecalibrationIndels.ts_filter_level = 99.0
    applyRecalibrationIndels.mode = org.broadinstitute.gatk.tools.walkers.variantrecalibration.VariantRecalibratorArgumentCollection.Mode.INDEL
    add(applyRecalibrationIndels)

    /**
      * *********************************************************************
      *                        CALLSET REFINEMENT
      * https://www.broadinstitute.org/gatk/guide/bp_step.php?p=3
      * https://www.broadinstitute.org/gatk/guide/article?id=4723
      * *********************************************************************
      */

    /**
      * Callset refinement
      * CalculateGenotypePosteriors: Calculate genotype posterior likelihoods given panel data
      * https://www.broadinstitute.org/gatk/guide/article?id=4727
      * https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_variantutils_CalculateGenotypePosteriors.php
      */
    // Here, we follow the method to "refine the genotypes of a large panel based on the discovered allele frequency"
    val calculateGenotypePosteriors = new CalculateGenotypePosteriors with CommonArguments
    calculateGenotypePosteriors.variant = applyRecalibrationIndels.out
    calculateGenotypePosteriors.out = swapExt(calculateGenotypePosteriors.variant, "vcf", "withPosteriors.vcf")
    add(calculateGenotypePosteriors)

    /**
      * Callset refinement
      * VariantFiltration: Filter variant calls based on INFO and FORMAT annotations
      * https://www.broadinstitute.org/gatk/guide/article?id=4727
      * https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_filters_VariantFiltration.php
      */
    // Use recommendations here: https://www.broadinstitute.org/gatk/guide/article?id=4727
    val variantFiltration = new VariantFiltration with CommonArguments
    variantFiltration.variant = calculateGenotypePosteriors.out
    variantFiltration.out = swapExt(variantFiltration.variant, "vcf", "Gfiltered.vcf")
    variantFiltration.filterExpression = Seq("GQ < 20.0")
    variantFiltration.filterName = Seq("lowGQ")
    add(variantFiltration)

    /**
      * Other callset refinement steps that we probably don't need
      * All listed in https://www.broadinstitute.org/gatk/guide/bp_step.php?p=3
      */

//    /**
//      * VariantAnnotator: Annotate variant calls with context information (annotate possible de novo mutations)
//      * https://www.broadinstitute.org/gatk/guide/article?id=4727
//      * https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_VariantAnnotator.php
//      */
//    // Can't call de novo mutations without pedigree information
//    // See https://www.broadinstitute.org/gatk/guide/article?id=4723 and https://www.broadinstitute.org/gatk/guide/article?id=4727
//
//    /**
//      * SelectVariants: Select a subset of variants from a larger callset
//      * https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_variantutils_SelectVariants.php
//      */
//
//    /**
//      * CombineVariants: Combine variant records from different sources
//      * https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_variantutils_CombineVariants.php
//      */
//
//    /**
//      * VariantEval: General-purpose tool for variant evaluation (% in dbSNP, genotype concordance, Ti/Tv ratios, and a lot more)
//      * https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_varianteval_VariantEval.php
//      */
//
//    /**
//      * VariantsToTable: Extract specific fields from a VCF file to a tab-delimited table
//      * https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_variantutils_VariantsToTable.php
//      */
//
//    /**
//      * GenotypeConcordance (Picard version): Genotype concordance
//      * https://broadinstitute.github.io/picard/command-line-overview.html#GenotypeConcordance
//      */



  }

}
