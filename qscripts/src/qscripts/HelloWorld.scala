package qscripts

import org.broadinstitute.gatk.queue.QScript

/**
  * Created by prussell on 3/3/16.
  */
class HelloWorld extends QScript {

  def script() = {
    add(new CommandLineFunction {
      def commandLine = "echo Hello World!"
    })
  }

}
