package ARIMA

import breeze.linalg._

object TestARMA {
  def main(args: Array[String]) {
    val x = Readcsv.read()
    //	  val x = Array[Double](8.0,10.0,7.0,6.0,9,8,6,5,7,4)
    val y = Array.range(-10, 21) map { v => v.toDouble }
//    val model = ARIMA.train(x, order = (3,2,3),xv = List("diff","1","trend","1"), method = "MLE", demean = true)
    val model = AutoFit.autofit(x)
    //  val model = ARIMA.train(x, order = (1,0,1),xv = List(), method = "hannan", demean = true)
    //      println(model.predict(20) toList)
    model.print_ARIMA_Model
    model.print_Predictions(12)
    model.print_LjungBox
  }
}
