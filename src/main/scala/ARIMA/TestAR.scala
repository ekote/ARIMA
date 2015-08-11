package ARIMA

import breeze.linalg._
import breeze.stats._
import breeze.numerics._
import IO.Readcsv._

object TestAR {
  

	def main(args: Array[String]) {
      val x = csvReader("D://Data/sharefolder/testdata.csv")
	  val y = Array[Double](1,2,3,4,5)
	  val model = AutoRegression.train(x, 2)
	  println(model.phi)
	  println(model.predict(10) toList)
    }
}