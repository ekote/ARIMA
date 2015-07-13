package ARIMA

import breeze.linalg._

object TestMA {
 
	def main(args: Array[String]) {
	  
      val x = Readcsv.read()
//	  val x = Array[Double](8.0,10.0,7.0,6.0,9,8,6,5,7,4)
	  val model = MovingAverage.train(x,q = 2, recursion_level = 17)
//      println(model.predict(20) toList)
	  
    }
}