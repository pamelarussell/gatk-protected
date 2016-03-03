package org.broadinstitute.gatk.queue.qscripts

import org.broadinstitute.gatk.queue.QScript
import org.broadinstitute.gatk.queue.extensions.picard
import org.broadinstitute.gatk.queue.extensions.gatk._

/**
  * Created by prussell on 3/3/16.
  * Implements GATK best practices
  */
class GatkBestPractices extends QScript {

  @Input(doc="Reference genome for the bam files", shortName="R", fullName="REF", required=true)
  var referenceFile: File = null

  @Input(doc="One or more bam files", shortName="I", fullName="INPUT", required=true)
  var bamFiles: List[File] = Nil

  def script() {

    val m = new MarkDuplicates
    val r = new RealignerTargetCreator
    val u = new UnifiedGenotyper

  }

}
