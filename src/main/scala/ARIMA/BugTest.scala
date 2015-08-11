package ARIMA

import breeze.linalg._
import breeze.stats._
import breeze.numerics._
import IO.Readcsv._
import IO.DateFrequency._

object BugTest {


  def main(args: Array[String]) {
    val x = csvReader("D://Data/sharefolder/testdata.csv")
    val y = Array[Double](1, 2, 3, 4, 5)
    val A = DenseMatrix.zeros[Double](5, 5)
    val a = DenseVector[Double](y)
    val zeta = DenseVector[Double](Array(0.4, 0.3, 0.2, 1, 1, 1))
    val phi = DenseVector[Double](0.4, 0.3, 1.3)

    val theta = zeta(5 to 3 by -1)
    val c = DenseVector[Double](1, 2, 3,4,5,6,7)
    //    println(phi dot zeta(3 to 5))
    println(FromDaily(y,2) toList)
    println(DenseVector[Double](0.0))
  }
}