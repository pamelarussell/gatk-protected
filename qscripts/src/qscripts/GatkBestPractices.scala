package qscripts

import org.broadinstitute.gatk.queue.QScript
//import org.broadinstitute.gatk.queue.extensions.picard.MarkDuplicates
import org.broadinstitute.gatk.queue.extensions.gatk._

/**
  * Created by prussell on 3/3/16.
  * Implements GATK best practices
  */
class GatkBestPractices extends QScript{

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
  var knownIndels: List[File] = null
  @Input(doc="Database of known variants e.g. dbSNP", fullName="KNOWN_VARIANTS", required=true)
  var knownPolymorphicSites: List[File] = null

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
      * *********************************************************************
      *              LOCAL REALIGNMENT AROUND INDELS
      * https://www.broadinstitute.org/gatk/guide/article?id=38
      * https://www.broadinstitute.org/gatk/guide/article?id=2800
      * *********************************************************************
      */

    /**
      * Local realignment around indels
      * Step 1 of 2: RealignerTargetCreator: Define intervals to target for local realignment
      * https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_indels_RealignerTargetCreator.php
      */
    val realignerTargetCreator = new RealignerTargetCreator
    val realignerTargetCreatorOutput : File = new File(outputDirectory.getAbsolutePath + "/" +
      outputPrefix + "_realignment_targets.list")
    realignerTargetCreator.reference_sequence = referenceFile
    realignerTargetCreator.known = knownIndels
    realignerTargetCreator.out = realignerTargetCreatorOutput
    realignerTargetCreator.maxIntervalSize = int2intOption(500) // Default 500
    realignerTargetCreator.minReadsAtLocus = int2intOption(4) // Default 4
    realignerTargetCreator.mismatchFraction = double2doubleOption(0.0) // Default 0.0
    realignerTargetCreator.windowSize = int2intOption(10) // Default 10
    add(realignerTargetCreator)

    /**
      * Local realignment around indels
      * Step 2 of 2: IndelRealigner: Perform local realignment of reads around indels
      * https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_indels_IndelRealigner.php
      */
    val indelRealigner = new IndelRealigner
    val indelRealignerOutput : File = new File(outputDirectory.getAbsolutePath + "/" +
      outputPrefix + "_realigned_reads.bam")
    indelRealigner.reference_sequence = referenceFile
    indelRealigner.targetIntervals = realignerTargetCreatorOutput
    indelRealigner.input_file = bamFiles
    indelRealigner.knownAlleles = knownIndels
    indelRealigner.out = indelRealignerOutput
    val indelRealignerOutputAsList : List[File] = List(indelRealignerOutput)
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

    /**
      * *********************************************************************
      *              BASE QUALITY SCORE RECALIBRATION
      * https://www.broadinstitute.org/gatk/guide/article?id=44
      * https://www.broadinstitute.org/gatk/guide/article?id=2801
      * *********************************************************************
      */

    /**
      * Base quality score recalibration
      * Step 1 of 3: BaseRecalibrator: Generate base recalibration table to compensate for systematic errors in basecalling confidences
      * https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_bqsr_BaseRecalibrator.php
      */
    val baseRecalibrator = new BaseRecalibrator
    val baseRecalibratorOutput : File = new File(outputDirectory.getAbsolutePath + "/" +
      outputPrefix + "_base_recalibrator.out")
    baseRecalibrator.input_file = indelRealignerOutputAsList
    indelRealigner.reference_sequence = referenceFile
    baseRecalibrator.out = baseRecalibratorOutput
    baseRecalibrator.knownSites = knownPolymorphicSites
    baseRecalibrator.indels_context_size = int2intOption(3) // Default 3
    baseRecalibrator.maximum_cycle_value = int2intOption(500) // Default 500
    baseRecalibrator.mismatches_context_size = int2intOption(2) // Default 2
    baseRecalibrator.solid_nocall_strategy = null
    baseRecalibrator.solid_recal_mode = null
    baseRecalibrator.list = false
    baseRecalibrator.lowMemoryMode = false
    baseRecalibrator.no_standard_covs = false
    baseRecalibrator.sort_by_all_columns = false
    baseRecalibrator.binary_tag_name = null
    baseRecalibrator.bqsrBAQGapOpenPenalty = double2doubleOption(40.0) // Default 40.0
    baseRecalibrator.deletions_default_quality = int2byteOption(45) // Default 45
    baseRecalibrator.insertions_default_quality = int2byteOption(45) // Default 45
    baseRecalibrator.low_quality_tail = int2byteOption(2) // Default 2
    baseRecalibrator.mismatches_default_quality = int2byteOption(-1) // Default -1
    baseRecalibrator.quantizing_levels = int2intOption(16) // Default 16
    baseRecalibrator.run_without_dbsnp_potentially_ruining_quality = false
    add(baseRecalibrator)

    /**
      * Base quality score recalibration
      * Step 2 of 3: AnalyzeCovariates: Create plots to visualize base recalibration results
      * https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_bqsr_AnalyzeCovariates.php
      */
    // TODO

    /**
      * Base quality score recalibration
      * Step 3 of 3: PrintReads: Write out sequence read data (for filtering, merging, subsetting etc)
      * https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_readutils_PrintReads.php
      */
    // TODO

    /**
      * *********************************************************************
      *                        VARIANT DISCOVERY
      * https://www.broadinstitute.org/gatk/guide/bp_step.php?p=2
      * *********************************************************************
      */

    /**
      * Variant discovery
      * Step 1 of 6: HaplotypeCaller: Call germline SNPs and indels via local re-assembly of haplotypes
      * https://www.broadinstitute.org/gatk/guide/article?id=2803
      * https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php
      * https://www.broadinstitute.org/gatk/guide/article?id=3893
      */
    // TODO

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
