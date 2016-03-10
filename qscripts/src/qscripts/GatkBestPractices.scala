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

  @Input(doc="Reference genome for the bam files", shortName="R", fullName="REF", required=true)
  var referenceFile: File = null

  @Input(doc="One or more bam files", shortName="I", fullName="INPUT", required=true)
  var bamFiles: List[File] = Nil

  @Input(doc="VCF file(s) with known indels", fullName="KNOWN_INDELS", required=false)
  var knownIndels: List[File] = null

  @Input(doc="Output directory", shortName="od", fullName="OUTDIR", required=true)
  var outputDirectory : File = null

  @Input(doc="Output prefix not including directory", shortName="op", fullName="OUT_PREFIX", required=true)
  var outputPrefix : String = null

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

    /**
      * *********************************************************************
      *              LOCAL REALIGNMENT AROUND INDELS
      * https://www.broadinstitute.org/gatk/guide/article?id=38
      * https://www.broadinstitute.org/gatk/guide/article?id=2800
      * *********************************************************************
      */

    /**
      * Local realignment around indels
      * Step 1 of 2: Define intervals to target for local realignment
      * https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_indels_RealignerTargetCreator.php
      */
    val realignerTargetCreator = new RealignerTargetCreator
    val realignerTargetCreatorOutput : File = new File(outputDirectory.getAbsolutePath + "/" +
      outputPrefix + "_realignment_targets.list")
    realignerTargetCreator.reference_sequence = referenceFile
    realignerTargetCreator.input_file = bamFilesDuplicatesMarked
    realignerTargetCreator.known = knownIndels
    realignerTargetCreator.out = realignerTargetCreatorOutput
    add(realignerTargetCreator)

    /**
      * Local realignment around indels
      * Step 2 of 2: Perform local realignment of reads around indels
      * https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_indels_IndelRealigner.php
      */
    val indelRealigner = new IndelRealigner
    val indelRealignerOutput : File = new File(outputDirectory.getAbsolutePath + "/" +
      outputPrefix + "_realigned_reads.bam")
    indelRealigner.reference_sequence = referenceFile
    indelRealigner.targetIntervals = realignerTargetCreatorOutput
    indelRealigner.known = knownIndels
    indelRealigner.out = indelRealignerOutput
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
      * Step 1 of 3: Generate base recalibration table to compensate for systematic errors in basecalling confidences
      * https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_bqsr_BaseRecalibrator.php
      */
    // TODO

    /**
      * Base quality score recalibration
      * Step 2 of 3: Create plots to visualize base recalibration results
      * https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_bqsr_AnalyzeCovariates.php
      */
    // TODO

    /**
      * Base quality score recalibration
      * Step 3 of 3: Write out sequence read data (for filtering, merging, subsetting etc)
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
      * Step 1 of 6: Call germline SNPs and indels via local re-assembly of haplotypes
      * https://www.broadinstitute.org/gatk/guide/article?id=2803
      * https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php
      * https://www.broadinstitute.org/gatk/guide/article?id=3893
      */
    // TODO

    /**
      * Variant discovery
      * Step 2 of 6: Combine per-sample gVCF files produced by HaplotypeCaller into a multi-sample gVCF file
      * https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_variantutils_CombineGVCFs.php
      */
    // TODO

    /**
      * Variant discovery
      * Step 3 of 6: Perform joint genotyping on gVCF files produced by HaplotypeCaller
      * https://www.broadinstitute.org/gatk/guide/article?id=3893
      * https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_variantutils_GenotypeGVCFs.php
      */
    // TODO

    /**
      * Variant discovery
      * Step 4 of 6: Filter variant calls based on INFO and FORMAT annotations
      * https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_filters_VariantFiltration.php
      */
    // TODO

    /**
      * Variant discovery
      * Step 5 of 6: Build a recalibration model to score variant quality for filtering purposes
      * https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_variantrecalibration_VariantRecalibrator.php
      * https://www.broadinstitute.org/gatk/guide/article?id=39
      * https://www.broadinstitute.org/gatk/guide/article?id=2805
      *
      */
    // TODO

    /**
      * Variant discovery
      * Step 6 of 6: Apply a score cutoff to filter variants based on a recalibration table
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
      * Step 1 of : Calculate genotype posterior likelihoods given panel data
      * https://www.broadinstitute.org/gatk/guide/article?id=4727
      * https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_variantutils_CalculateGenotypePosteriors.php
      */
    // TODO

    /**
      * Callset refinement
      * Step 2 of : Filter variant calls based on INFO and FORMAT annotations
      * https://www.broadinstitute.org/gatk/guide/article?id=4727
      * https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_filters_VariantFiltration.php
      */
    // TODO

    /**
      * Callset refinement
      * Step 3 of : Annotate variant calls with context information
      * https://www.broadinstitute.org/gatk/guide/article?id=4727
      * https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_VariantAnnotator.php
      */
    // TODO

    /**
      * Callset refinement
      * Step 4 of : Annotate variant calls with context information
      * https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_VariantAnnotator.php
      */
    // TODO



  }

}
