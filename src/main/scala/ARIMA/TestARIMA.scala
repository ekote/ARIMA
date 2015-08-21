/**
 * @author  Jian Wang
 * @version 1.0
 * @date    10/06/2015
 */

package ARIMA

import breeze.linalg._
import IO.Readcsv._
import IO.DateFrequency._

// test ARIMA model
object TestARIMA {
  def main(args: Array[String]) {
    val x = csvReader("D://Data/sharefolder/testdata.csv")
    val x_weekly = FromDaily(x.slice(0,280),7)
    val x_weekly2 = FromDaily(x.slice(280,397),7)
    //	  val x = Array[Double](8.0,10.0,7.0,6.0,9,8,6,5,7,4)
    val y = Array.range(-10, 21) map { v => v.toDouble }
//accepted test    val model = ARIMA.train(x, order = (3,1,2),xv = List("diff","1"), method = "MLE", demean = true)
//accepted test    val model = ARIMA.train(y, order = (1,1,0),xv = List("diff","1"), MaxLag = 20, method = "MLE", demean = true)
//accepted test   val model = ARIMA.train(y, order = (1,0,0),method = "MLE", demean = true, start_params = DenseVector[Double](0.9899))
//accepted test      val model = ARIMA.train(x, order = (1,1,1),xv = List("diff","1"), method = "MLE", demean = true)
    val model = ARIMA.train(x,order=(1,0,0),xv = List("season","13"), method = "MLE",MeanMLEQ = true, demean = true)
    println(x_weekly toList)
    println(x_weekly2 toList)
    //a = arma(e,1,1) Error in arima(x, c(p, 0, q)) : non-stationary AR part from CSS
    // When using CSS (conditional sum of squares), it is possible for the autoregressive coefficients to be non-stationary
    // so in R, MLE will fail sometimes. 
    // in the case of non-stationary (even if R do not fail), we will have very different coefficient estimation from R because of LBGFS algorithm
    //    val model = ARIMA.train(x, order = (1,0,1),xv = List(), method = "MLE", demean = true, start_params = DenseVector(0.1,0.2))
      //      println(model.predict(20) toList)
    model.print_ARIMA_Model
    model.print_Predictions(16)
    model.print_train_data
//    model.print_LjungBox
//    model.print_res_and_fits

  }
}
