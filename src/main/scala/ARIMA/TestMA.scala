package ARIMA

import breeze.linalg._
import IO.Readcsv._

object TestMA {
 
	def main(args: Array[String]) {
	  
      val x = csvReader("D://Data/sharefolder/testdata.csv")
//	  val x = Array[Double](8.0,10.0,7.0,6.0,9,8,6,5,7,4)
	  val model = MovingAverage.train(x,q = 10, recursion_level = 17)
	  println(model.theta)
      println(model.predict(10) toList)
	  
    }
}