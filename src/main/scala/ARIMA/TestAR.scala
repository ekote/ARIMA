package ARIMA

import breeze.linalg._
import breeze.stats._
import breeze.numerics._

object TestAR {
  

	def main(args: Array[String]) {
      val x = Readcsv.read()
	  val y = Array[Double](1,2,3,4,5)
	  val model = AutoRegression.train(x, 2)
    }
}