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
  @Input(doc="VCF file(s) with known indels", fullName="KNOWN_INDELS", required=true)
  var knownIndels: File = null
  @Input(doc="Database of known variants e.g. dbSNP", fullName="KNOWN_VARIANTS", required=true)
  var knownPolymorphicSites: File = null

  /**
    * Downsampling fraction
    */
  @Argument(doc="Downsample to fraction", fullName="DOWNSAMPLE_TO_FRACTION", required=false)
  var downsampleToFraction: Double = 1.0

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
    this.downsample_to_fraction = Some(downsampleToFraction)
    // TODO other common arguments?
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

    // TODO fix dependency issue to get picard package
    //val markDuplicates = new MarkDuplicates
    // TODO replace with actual output from mark duplicates
    val bamFilesDuplicatesMarked : List[File] = Nil
    // TODO
    // TODO use duplicates marked bam files as input to pipeline


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
        * Step 1 of 2: RealignerTargetCreator: Define intervals to target for local realignment
        * https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_indels_RealignerTargetCreator.php
        */
      val realignerTargetCreator = new RealignerTargetCreator with CommonArguments
      realignerTargetCreator.input_file +:= bam
      realignerTargetCreator.known = Seq(knownIndels)
      realignerTargetCreator.out = swapExt(bam, "bam", "interval_list")
      realignerTargetCreator.maxIntervalSize = int2intOption(500) // Default 500
      realignerTargetCreator.minReadsAtLocus = int2intOption(4) // Default 4
      realignerTargetCreator.mismatchFraction = double2doubleOption(0.0) // Default 0.0
      realignerTargetCreator.windowSize = int2intOption(10) // Default 10
      add(realignerTargetCreator)
      println("LOG:\tAdded RealignerTargetCreator. " +
        "Input: " + realignerTargetCreator.input_file + ". " +
        "Output: " + realignerTargetCreator.out + ".")

      /**
        * Local realignment around indels
        * Step 2 of 2: IndelRealigner: Perform local realignment of reads around indels
        * https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_indels_IndelRealigner.php
        */
      val indelRealigner = new IndelRealigner with CommonArguments
      indelRealigner.targetIntervals = realignerTargetCreator.out
      indelRealigner.input_file +:= bam
      indelRealigner.knownAlleles = Seq(knownIndels)
      indelRealigner.out = swapExt(bam, "bam", "realign.bam")
      indelRealigner.consensusDeterminationModel = null
      indelRealigner.LODThresholdForCleaning = double2doubleOption(5.0) // Default 5.0
      indelRealigner.nWayOut = null
      indelRealigner.entropyThreshold = double2doubleOption(0.15) // Default 0.15
      indelRealigner.maxConsensuses = int2intOption(30) // Default 30
      indelRealigner.maxIsizeForMovement = int2intOption(3000) // Default 3000
      indelRealigner.maxPositionalMoveAllowed = int2intOption(200) // Default 200
      indelRealigner.maxReadsForConsensuses = int2intOption(120) // Default 120
      indelRealigner.maxReadsForRealignment = int2intOption(20000) // Default 20000
      indelRealigner.maxReadsInMemory = int2intOption(150000) // Default 150000
      indelRealigner.noOriginalAlignmentTags = false
      add(indelRealigner)
      println("LOG:\tAdded IndelRealigner. " +
        "Input: " + indelRealigner.input_file + ". " +
        "Output: " + indelRealigner.out + ".")

      /**
        * *********************************************************************
        * BASE QUALITY SCORE RECALIBRATION
        * https://www.broadinstitute.org/gatk/guide/article?id=44
        * https://www.broadinstitute.org/gatk/guide/article?id=2801
        * *********************************************************************
        */

      /**
        * Base quality score recalibration
        * Step 1 of 3: BaseRecalibrator: Generate base recalibration table to compensate for systematic errors in basecalling confidences
        * https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_bqsr_BaseRecalibrator.php
        */

      // Generate the first pass recalibration table file
      val baseRecalibratorBefore = new BaseRecalibrator with CommonArguments
      baseRecalibratorBefore.input_file +:= indelRealigner.out
      baseRecalibratorBefore.out = swapExt(indelRealigner.out, "realign.bam", "base_recalibrator_first_pass.out")
      baseRecalibratorBefore.knownSites = Seq(knownPolymorphicSites)
      baseRecalibratorBefore.indels_context_size = int2intOption(3) // Default 3
      baseRecalibratorBefore.maximum_cycle_value = int2intOption(500) // Default 500
      baseRecalibratorBefore.mismatches_context_size = int2intOption(2) // Default 2
      baseRecalibratorBefore.solid_nocall_strategy = null
      baseRecalibratorBefore.solid_recal_mode = null
      baseRecalibratorBefore.list = false
      baseRecalibratorBefore.lowMemoryMode = false
      baseRecalibratorBefore.no_standard_covs = false
      baseRecalibratorBefore.sort_by_all_columns = false
      baseRecalibratorBefore.binary_tag_name = null
      baseRecalibratorBefore.bqsrBAQGapOpenPenalty = double2doubleOption(40.0) // Default 40.0
      baseRecalibratorBefore.deletions_default_quality = int2byteOption(45) // Default 45
      baseRecalibratorBefore.insertions_default_quality = int2byteOption(45) // Default 45
      baseRecalibratorBefore.low_quality_tail = int2byteOption(2) // Default 2
      baseRecalibratorBefore.mismatches_default_quality = int2byteOption(-1) // Default -1
      baseRecalibratorBefore.quantizing_levels = int2intOption(16) // Default 16
      baseRecalibratorBefore.run_without_dbsnp_potentially_ruining_quality = false
      add(baseRecalibratorBefore)
      println("LOG:\tAdded BaseRecalibrator. " +
        "Input: " + baseRecalibratorBefore.input_file + ". " +
        "Output: " + baseRecalibratorBefore.out + ".")

      // Generate the second pass recalibration table file
      val baseRecalibratorAfter = new BaseRecalibrator with CommonArguments
      baseRecalibratorAfter.BQSR = baseRecalibratorBefore.out
      baseRecalibratorAfter.input_file +:= indelRealigner.out
      baseRecalibratorAfter.out = swapExt(indelRealigner.out, "realign.bam", "base_recalibrator_second_pass.out")
      baseRecalibratorAfter.knownSites = Seq(knownPolymorphicSites)
      baseRecalibratorAfter.indels_context_size = int2intOption(3) // Default 3
      baseRecalibratorAfter.maximum_cycle_value = int2intOption(500) // Default 500
      baseRecalibratorAfter.mismatches_context_size = int2intOption(2) // Default 2
      baseRecalibratorAfter.solid_nocall_strategy = null
      baseRecalibratorAfter.solid_recal_mode = null
      baseRecalibratorAfter.list = false
      baseRecalibratorAfter.lowMemoryMode = false
      baseRecalibratorAfter.no_standard_covs = false
      baseRecalibratorAfter.sort_by_all_columns = false
      baseRecalibratorAfter.binary_tag_name = null
      baseRecalibratorAfter.bqsrBAQGapOpenPenalty = double2doubleOption(40.0) // Default 40.0
      baseRecalibratorAfter.deletions_default_quality = int2byteOption(45) // Default 45
      baseRecalibratorAfter.insertions_default_quality = int2byteOption(45) // Default 45
      baseRecalibratorAfter.low_quality_tail = int2byteOption(2) // Default 2
      baseRecalibratorAfter.mismatches_default_quality = int2byteOption(-1) // Default -1
      baseRecalibratorAfter.quantizing_levels = int2intOption(16) // Default 16
      baseRecalibratorAfter.run_without_dbsnp_potentially_ruining_quality = false
      add(baseRecalibratorAfter)
      println("LOG:\tAdded BaseRecalibrator. " +
        "Input: " + baseRecalibratorAfter.input_file + ". " +
        "Output: " + baseRecalibratorAfter.out + ".")

      /**
        * Base quality score recalibration
        * Step 2 of 3: AnalyzeCovariates: Create plots to visualize base recalibration results
        * https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_bqsr_AnalyzeCovariates.php
        */
      val analyzeCovariates = new AnalyzeCovariates with CommonArguments
      analyzeCovariates.beforeReportFile = baseRecalibratorBefore.out
      analyzeCovariates.afterReportFile = baseRecalibratorAfter.out
      analyzeCovariates.plotsReportFile = new File(outputDirectory.getAbsolutePath + "/analyzeCovariates_" + swapExt(indelRealigner.out, "bam", "BQSR.pdf"))
      analyzeCovariates.intermediateCsvFile = new File(outputDirectory.getAbsolutePath + "/analyzeCovariates_" + swapExt(indelRealigner.out, "bam", "BQSR.csv"))
      analyzeCovariates.ignoreLMT = false
      add(analyzeCovariates)
      println("LOG:\tAdded AnalyzeCovariates. " +
        "Input: " + analyzeCovariates.input_file + ". " +
        "Plots report file: " + analyzeCovariates.plotsReportFile + ". " +
        "Intermediate CSV file: " + analyzeCovariates.intermediateCsvFile + ".")

      /**
        * Base quality score recalibration
        * Step 3 of 3: PrintReads: Write out sequence read data (for filtering, merging, subsetting etc)
        * https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_readutils_PrintReads.php
        */
      val printReads = new PrintReads with CommonArguments
      printReads.input_file +:= bam
      printReads.BQSR = baseRecalibratorAfter.out
      printReads.out = swapExt(bam, "bam", "recalibrated.bam")
      add(printReads)
      println("LOG:\tAdded PrintReads. " +
        "Input: " + printReads.input_file + ". " +
        "Output: " + printReads.out + ".")

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
      * Step 1 of 6: HaplotypeCaller: Call germline SNPs and indels via local re-assembly of haplotypes
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
      haplotypeCaller.dbsnp = knownPolymorphicSites
      //haplotypeCaller.intervals = ??? // TODO do we want this?
      haplotypeCaller.out = swapExt(processedBam, "bam", "raw.snps.indels.g.vcf")
      add(haplotypeCaller)
      sampleGVCFs +:= haplotypeCaller.out
      println("LOG:\tAdded HaplotypeCaller. " +
        "Input: " + haplotypeCaller.input_file + ". " +
        "Output: " + haplotypeCaller.out + ".")
    }

    /**
      * Variant discovery
      * Step 2 of 6: CombineGVCFs: Combine per-sample gVCF files produced by HaplotypeCaller into a multi-sample gVCF file
      * https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_variantutils_CombineGVCFs.php
      */
    // TODO

    /**
      * Variant discovery
      * Step 3 of 6: GenotypeGVCFs: Perform joint genotyping on gVCF files produced by HaplotypeCaller
      * https://www.broadinstitute.org/gatk/guide/article?id=3893
      * https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_variantutils_GenotypeGVCFs.php
      */
    // TODO

    /**
      * Variant discovery
      * Step 4 of 6: VariantFiltration: Filter variant calls based on INFO and FORMAT annotations
      * https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_filters_VariantFiltration.php
      */
    // TODO

    /**
      * Variant discovery
      * Step 5 of 6: VariantRecalibrator: Build a recalibration model to score variant quality for filtering purposes
      * https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_variantrecalibration_VariantRecalibrator.php
      * https://www.broadinstitute.org/gatk/guide/article?id=39
      * https://www.broadinstitute.org/gatk/guide/article?id=2805
      *
      */
    // TODO

    /**
      * Variant discovery
      * Step 6 of 6: ApplyRecalibration: Apply a score cutoff to filter variants based on a recalibration table
      * https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_variantrecalibration_ApplyRecalibration.php
      * https://www.broadinstitute.org/gatk/guide/article?id=2806
      */
    // TODO

    /**
      * *********************************************************************
      *                        CALLSET REFINEMENT
      * https://www.broadinstitute.org/gatk/guide/bp_step.php?p=3
      * *********************************************************************
      */

    /**
      * Callset refinement
      * Step 1 of 8: CalculateGenotypePosteriors: Calculate genotype posterior likelihoods given panel data
      * https://www.broadinstitute.org/gatk/guide/article?id=4727
      * https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_variantutils_CalculateGenotypePosteriors.php
      */
    // TODO

    /**
      * Callset refinement
      * Step 2 of 8: VariantFiltration: Filter variant calls based on INFO and FORMAT annotations
      * https://www.broadinstitute.org/gatk/guide/article?id=4727
      * https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_filters_VariantFiltration.php
      */
    // TODO

    /**
      * Callset refinement
      * Step 3 of 8: VariantAnnotator: Annotate variant calls with context information
      * https://www.broadinstitute.org/gatk/guide/article?id=4727
      * https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_VariantAnnotator.php
      */
    // TODO

    /**
      * Callset refinement
      * Step 4 of 8: SelectVariants: Select a subset of variants from a larger callset
      * https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_variantutils_SelectVariants.php
      */
    // TODO

    /**
      * Callset refinement
      * Step 5 of 8: CombineVariants: Combine variant records from different sources
      * https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_variantutils_CombineVariants.php
      */
    // TODO

    /**
      * Callset refinement
      * Step 6 of 8: VariantEval: General-purpose tool for variant evaluation (% in dbSNP, genotype concordance, Ti/Tv ratios, and a lot more)
      * https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_varianteval_VariantEval.php
      */
    // TODO

    /**
      * Callset refinement
      * Step 7 of 8: VariantsToTable: Extract specific fields from a VCF file to a tab-delimited table
      * https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_variantutils_VariantsToTable.php
      */
    // TODO

    /**
      * Callset refinement
      * Step 8 of 8: GenotypeConcordance (Picard version): Genotype concordance
      * https://broadinstitute.github.io/picard/command-line-overview.html#GenotypeConcordance
      */
    // TODO



  }

}
