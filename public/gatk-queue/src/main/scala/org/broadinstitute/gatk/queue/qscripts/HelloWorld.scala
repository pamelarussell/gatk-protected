package org.broadinstitute.gatk.queue.qscripts

import org.broadinstitute.gatk.queue.QScript

/**
  * Created by prussell on 2/26/16.
  */
class HelloWorld extends QScript {
    def script = {
      add(new CommandLineFunction {
        def commandLine = "echo hello world"
      })
    }
}
