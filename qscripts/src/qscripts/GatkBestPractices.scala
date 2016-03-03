package qscripts

import org.broadinstitute.gatk.queue.QScript
//import org.broadinstitute.gatk.queue.extensions.picard.MarkDuplicates
import org.broadinstitute.gatk.queue.extensions.gatk._

/**
  * Created by prussell on 3/3/16.
  * Implements GATK best practices
  */
class GatkBestPractices extends QScript{

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
      * Map to reference
      * https://www.broadinstitute.org/gatk/guide/bp_step.php?p=1
      */
    // TODO

    /**
      * Mark duplicates
      * https://www.broadinstitute.org/gatk/guide/bp_step.php?p=1
      */
    // TODO fix dependency issue to get picard package
    //val markDuplicates = new MarkDuplicates
    // TODO replace with actual output from mark duplicates
    val bamFilesDuplicatesMarked : List[File] = Nil

    /**
      * Create target list of intervals to be realigned
      * https://www.broadinstitute.org/gatk/guide/article?id=38
      * https://www.broadinstitute.org/gatk/guide/article?id=2800
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
      * Perform realignment of target intervals
      * https://www.broadinstitute.org/gatk/guide/article?id=38
      * https://www.broadinstitute.org/gatk/guide/article?id=2800
      */
    val indelRealigner = new IndelRealigner
    val indelRealignerOutput : File = new File(outputDirectory.getAbsolutePath + "/" +
      outputPrefix + "_realigned_reads.bam")
    indelRealigner.reference_sequence = referenceFile
    indelRealigner.targetIntervals = realignerTargetCreatorOutput
    indelRealigner.known = knownIndels
    indelRealigner.out = indelRealignerOutput
    add(indelRealigner)

  }

}
